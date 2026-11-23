# 06-13_gblup_prs.R - Cross-validated genome-wide polygenic prediction (GBLUP) of CSLB scores. See README.md.
library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_keep <- as.integer(Sys.getenv("TEST_KEEP", "0"))
tag <- if (test_keep > 0) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("06-13", tag))
pic_prefix <- file.path(pic_dir, paste0("06-13", tag))
work_dir <- paste0(out_prefix, "_work")
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
plink2 <- "/data/00_tools/plink2_linux_avx2_20231123/plink2"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
null_perm <- as.integer(Sys.getenv("NULL_PERM", if (test_keep > 0) "1" else "20"))
K <- if (test_keep > 0) 3L else 10L
R <- if (test_keep > 0) 1L else 5L
rel_main <- 0.0884; rel_sens <- 0.0442
king_f <- file.path(base_dir, "GeneticData/v2023/DogAgingProject_GeneticData_CuratedRelease_2023.kinship.tsv")
grm_full <- file.path(in_dir, "04-1_grm", "grm_all")
bfile <- file.path(in_dir, "03-3_maf0.01_hwe")
prune <- file.path(in_dir, "03-4_pca", "prune.prune.in")
pheno_f <- file.path(in_dir, "04-2_pheno.txt"); qcov_f <- file.path(in_dir, "03-5_qcovar_PC5.txt"); cov_f <- file.path(in_dir, "03-5_covar.txt")
for (f in c(gcta, plink2, paste0(grm_full, c(".grm.bin", ".grm.id")), paste0(bfile, c(".bed", ".bim", ".fam")), prune, pheno_f, qcov_f, cov_f, king_f))
  if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261113)
sh <- function(cmd, args, log) invisible(system2(cmd, args, stdout = log, stderr = log))
wr <- function(d, f) { fwrite(d, f, sep = " ", col.names = FALSE, quote = FALSE); f }
qn <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5)); cn <- c("Sex", "Intact", "Vision", "Cohort_BHS")

ids <- fread(paste0(grm_full, ".grm.id"), header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
grm <- grm_full
if (test_keep > 0) {
  kf <- wr(ids[seq(1, .N, by = test_keep)], file.path(work_dir, "test_keep.txt")); grm <- file.path(work_dir, "grm_test")
  if (!file.exists(paste0(grm, ".grm.bin"))) sh(gcta, c("--grm", grm_full, "--keep", kf, "--make-grm", "--out", grm, "--threads", n_cores), file.path(work_dir, "grm_test.stdout"))
  ids <- fread(paste0(grm, ".grm.id"), header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
  cat(sprintf("TEST_KEEP=%d: %d마리 (결과 해석 금지)\n", test_keep, nrow(ids)))
}
ph <- fread(pheno_f, header = FALSE, colClasses = list(character = 1:2)); setnames(ph, c("FID", "IID", "y", "score"))
qc <- fread(qcov_f, header = FALSE, colClasses = list(character = 1:2)); setnames(qc, c("FID", "IID", qn))
cv <- fread(cov_f, header = FALSE, colClasses = "character"); setnames(cv, c("FID", "IID", cn))
d <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")), list(ids, ph[, .(FID, IID, y)], qc, cv))
d[, old := age_cslb >= 8]
idf <- wr(ids[, .(IID)], file.path(work_dir, "ids_for_king.txt"))
pairs <- fread(cmd = sprintf("awk -F'\\t' 'NR==FNR{k[$1]=1;next} FNR>1 && ($2 in k) && ($4 in k) && $8>%g{print $2\"\\t\"$4\"\\t\"$8}' %s %s",
                             rel_sens, shQuote(idf), shQuote(king_f)), header = FALSE, col.names = c("a", "b", "r"), colClasses = c("character", "character", "numeric"))
pairs <- rbind(pairs, pairs[, .(a = b, b = a, r)])
cat(sprintf("개 %d마리 (8세 이상 %d), 친척 쌍 KING > %.4f: %s · > %.4f: %s\n", nrow(d), sum(d$old), rel_sens, format(nrow(pairs) / 2, big.mark = ","),
            rel_main, format(pairs[r > rel_main, .N] / 2, big.mark = ",")))

pbfile <- file.path(work_dir, "prune_bfile")
if (!file.exists(paste0(pbfile, ".bed")))
  sh(plink2, c("--bfile", bfile, "--dog", "--extract", prune, "--keep", wr(ids[, .(FID, IID)], file.path(work_dir, "all_ids.txt")),
               "--make-bed", "--out", pbfile, "--threads", n_cores), paste0(pbfile, ".stdout"))
if (!file.exists(paste0(pbfile, ".bed"))) stop("작은 bfile 만들기 실패: ", pbfile, ".stdout")

fit_predict <- function(name, train, test, yvec = NULL) {
  out <- file.path(work_dir, name); prs_f <- paste0(out, ".sscore")
  if (!file.exists(prs_f)) {
    tr <- d[IID %in% train]
    if (!is.null(yvec)) tr[, y := yvec[match(IID, d$IID)]]
    keep_c <- cn[vapply(cn, function(j) uniqueN(tr[[j]]) > 1, TRUE)]
    keep_q <- qn[vapply(qn, function(j) isTRUE(sd(tr[[j]]) > 0), TRUE)]
    pf <- wr(tr[, .(FID, IID, y)], paste0(out, ".pheno")); qf <- wr(tr[, c("FID", "IID", keep_q), with = FALSE], paste0(out, ".qcovar"))
    cf <- wr(tr[, c("FID", "IID", keep_c), with = FALSE], paste0(out, ".covar")); kf <- wr(tr[, .(FID, IID)], paste0(out, ".keep"))
    sh(gcta, c("--reml", "--grm", grm, "--keep", kf, "--pheno", pf, "--qcovar", qf, "--covar", cf, "--reml-pred-rand", "--reml-no-lrt",
               "--reml-maxit", "200", "--out", out, "--threads", n_cores), paste0(out, ".gcta.stdout"))
    if (!file.exists(paste0(out, ".indi.blp"))) { warning(name, ": GREML 실패 (", out, ".gcta.stdout)"); return(NULL) }
    sh(gcta, c("--bfile", pbfile, "--keep", kf, "--blup-snp", paste0(out, ".indi.blp"), "--autosome-num", "38",
               "--out", out, "--threads", n_cores), paste0(out, ".blupsnp.stdout"))
    if (!file.exists(paste0(out, ".snp.blp"))) { warning(name, ": blup-snp 실패"); return(NULL) }
    tf <- wr(d[IID %in% test, .(FID, IID)], paste0(out, ".test"))
    sh(plink2, c("--bfile", pbfile, "--dog", "--keep", tf, "--score", paste0(out, ".snp.blp"), "1", "2", "3", "cols=+scoresums",
                 "--out", out, "--threads", n_cores), paste0(out, ".score.stdout"))
    if (!file.exists(prs_f)) { warning(name, ": plink2 --score 실패"); return(NULL) }
    unlink(paste0(out, c(".pheno", ".qcovar", ".covar")))
  }
  s <- fread(prs_f); setnames(s, 1:2, c("FID", "IID"))
  hsq <- if (file.exists(paste0(out, ".hsq"))) { h <- fread(paste0(out, ".hsq"), fill = TRUE); as.numeric(h[Source == "V(G)/Vp", Variance]) } else NA_real_
  data.table(IID = as.character(s$IID), prs = s[[grep("SCORE1_SUM", names(s))]], h2_train = hsq)
}
excl_rel <- function(test, thr) unique(pairs[r > thr & a %in% test, b])

make_folds <- function(seed) { set.seed(seed); f <- integer(nrow(d))
  for (g in c(TRUE, FALSE)) { w <- which(d$old == g); f[w] <- sample(rep_len(seq_len(K), length(w))) }; f }
folds <- sapply(seq_len(R), function(r) make_folds(20261113 + r))
fwrite(cbind(d[, .(IID, old)], setNames(as.data.table(folds), paste0("rep", seq_len(R)))), paste0(out_prefix, "_folds.csv"))
run_cv <- function(prefix, foldvec, thr, yvec = NULL) {
  rbindlist(lapply(seq_len(K), function(kk) {
    test <- d$IID[foldvec == kk]; drop <- excl_rel(test, thr)
    train <- setdiff(d$IID[foldvec != kk], drop)
    cat(sprintf("  %s fold %d: 학습 %d (KING > %.4f 친척 %d마리 뺌), 검증 %d\n", prefix, kk, length(train), thr, length(intersect(drop, d$IID[foldvec != kk])), length(test)))
    r <- fit_predict(sprintf("%s_f%02d", prefix, kk), train, test, yvec)
    if (is.null(r)) NULL else r[, fold := kk]
  }))
}
eval_r2 <- function(oof, yvec = NULL, label) {
  e <- merge(d, oof, by = "IID"); if (!is.null(yvec)) e[, y := yvec[match(IID, d$IID)]]
  e[, prs_z := (prs - mean(prs)) / sd(prs), by = fold]
  rbindlist(lapply(c("all", "age8", "young8"), function(g) {
    s <- if (g == "all") e else if (g == "age8") e[old == TRUE] else e[old == FALSE]
    if (nrow(s) < 20) return(NULL)
    keep_c <- cn[vapply(cn, function(j) uniqueN(s[[j]]) > 1, TRUE)]
    keep_q <- qn[vapply(qn, function(j) isTRUE(sd(s[[j]]) > 0), TRUE)]
    fx <- c(keep_q, keep_c, if (uniqueN(s$fold) > 1) "factor(fold)")
    f0 <- as.formula(paste("y ~", paste(fx, collapse = " + "))); f1 <- update(f0, . ~ . + prs_z)
    m0 <- lm(f0, data = s); m1 <- lm(f1, data = s)
    ff <- if (uniqueN(s$fold) > 1) "factor(fold)" else "1"
    r2_prs_only <- summary(lm(as.formula(paste("y ~ prs_z +", ff)), data = s))$r.squared - (if (ff == "1") 0 else summary(lm(y ~ factor(fold), data = s))$r.squared)
    r0 <- lm(as.formula(paste("y ~", paste(c(keep_q, keep_c), collapse = " + "))), data = s); r1 <- update(r0, . ~ . + prs)
    data.table(analysis = label, group = g, n = nrow(s), r2_prs_only = r2_prs_only, r2_base = summary(m0)$r.squared, r2_full = summary(m1)$r.squared,
               incr_r2 = summary(m1)$r.squared - summary(m0)$r.squared, beta_prs = coef(m1)[["prs_z"]], p = anova(m0, m1)$`Pr(>F)`[2],
               cor_prs_resid = cor(s$prs_z, resid(m0)), h2_train_mean = mean(s$h2_train, na.rm = TRUE),
               incr_r2_raw = summary(r1)$r.squared - summary(r0)$r.squared, p_raw = anova(r0, r1)$`Pr(>F)`[2])
  }))
}
t0 <- Sys.time()
cat("(가) 반복 교차검증\n")
ev <- rbindlist(lapply(seq_len(R), function(rr) {
  oof <- run_cv(sprintf("cv_r%d", rr), folds[, rr], rel_main)
  fwrite(oof[, .(IID, fold, prs, h2_train)], paste0(out_prefix, sprintf("_oof_prs_rep%d.csv", rr)))
  eval_r2(oof, NULL, sprintf("rep%d (KING > %.4f 제외)", rr, rel_main))[, rep := rr]
}))
cat(sprintf("(가) 민감도: KING > %.4f 제외 (반복 1)\n", rel_sens))
oof_s <- run_cv("cvsens_r1", folds[, 1], rel_sens)
ev <- rbind(ev, eval_r2(oof_s, NULL, sprintf("민감도 rep1 (KING > %.4f 제외)", rel_sens))[, rep := 1L], fill = TRUE)
fwrite(ev, paste0(out_prefix, "_eval.csv"))

cat("(나) 젊은 개로 학습 → 노령견 예측\n")
old_ids <- d[old == TRUE, IID]; drop <- excl_rel(old_ids, rel_main)
train_y <- setdiff(d[old == FALSE, IID], drop)
cat(sprintf("  학습 young8 %d (노령견과 KING > %.4f인 %d마리 뺌), 검증 age8 %d\n", length(train_y), rel_main, length(intersect(drop, d[old == FALSE, IID])), length(old_ids)))
y2o <- fit_predict("young2old", train_y, old_ids); if (!is.null(y2o)) y2o[, fold := 1L]
y2o_ev <- if (!is.null(y2o)) eval_r2(y2o, NULL, "young8 → age8")[group == "age8"] else data.table()
fwrite(y2o_ev, paste0(out_prefix, "_young2old.csv"))

cat(sprintf("null: y를 나이 집단 안에서 섞어 교차검증 %d번\n", null_perm))
nul <- rbindlist(lapply(seq_len(null_perm), function(pp) {
  set.seed(9000 + pp); yp <- d$y; for (g in c(TRUE, FALSE)) { w <- which(d$old == g); yp[w] <- d$y[w][sample.int(length(w))] }
  oof <- run_cv(sprintf("null_p%d", pp), folds[, 1], rel_main, yp)
  eval_r2(oof, yp, sprintf("null perm %d", pp))[, perm := pp]
}))
fwrite(nul, paste0(out_prefix, "_null.csv"))
cat(sprintf("GBLUP 끝 (%.1f분)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

main_ev <- ev[grepl("^rep", analysis)]
verdict <- rbindlist(lapply(c("all", "age8", "young8"), function(g) {
  m <- main_ev[group == g]; nm <- nul[group == g]
  data.table(group = g, reps = nrow(m), reps_p05 = sum(m$p < 0.05), mean_incr_r2 = mean(m$incr_r2), mean_cor = mean(m$cor_prs_resid),
             null_n = nrow(nm), null_p05 = sum(nm$p < 0.05), null_max_incr_r2 = if (nrow(nm)) max(nm$incr_r2) else NA_real_,
             p_emp = if (nrow(nm)) (sum(nm$incr_r2 >= mean(m$incr_r2)) + 1) / (nrow(nm) + 1) else NA_real_,
             predicted = nrow(m) == R && all(m$p < 0.05) && nrow(nm) > 0 && mean(m$incr_r2) > max(nm$incr_r2))
}))
pd <- rbind(main_ev[, .(type = "cross-validation (5 repeats)", group, incr_r2)], nul[, .(type = "null (permuted y)", group, incr_r2)],
            if (nrow(y2o_ev)) y2o_ev[, .(type = "young8 -> age8", group, incr_r2)])
g1 <- ggplot(pd, aes(group, incr_r2, colour = type)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6), size = 2) +
  scale_colour_manual(values = c(`cross-validation (5 repeats)` = "#1f4e9c", `null (permuted y)` = "grey60", `young8 -> age8` = "#b2182b")) +
  labs(x = NULL, y = "incremental R2 of PRS over covariates (out-of-fold)", colour = NULL,
       title = "Genome-wide GBLUP polygenic score for CSLB", subtitle = "relatives (KING > 0.0884, 2nd degree) of test dogs removed from training; covariates incl. PC1-5") +
  theme_minimal(base_size = 10)
ggsave(paste0(pic_prefix, "_r2.png"), g1, width = 7, height = 4.5, dpi = 300)
f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("06-13 유전체 전체 GBLUP 다유전 점수 (PRS)로 CSLB 예측 — 교차검증 (탐색; 설정·판정은 실행 전 헤더에 고정)\n")
if (test_keep > 0) cat(sprintf("*** TEST_KEEP=%d 시험 결과. 해석하지 말 것 ***\n", test_keep))
cat(sprintf("SNP %s개 (GRM과 같은 LD pruning 세트), %d-fold × %d번 반복, 친척 (KING > %.4f) 학습 제외, 공변량 = 06-6과 같음 (PC1 ~ 5 포함, 품종 이름 없음)\n\n",
            format(length(readLines(prune)), big.mark = ","), K, R, rel_main))
cat("■ (가) 교차검증: 반복 × 집단별 증분 R² (fold 안 표준화 PRS를 공변량 + fold 모형에 더했을 때; _raw = 첫 실행 방식, 참고)\n")
print(f3(ev[, .(analysis, group, n, r2_prs_only, r2_base, r2_full, incr_r2, beta_prs, p, cor_prs_resid, h2_train_mean, incr_r2_raw, p_raw)]))
cat("  r2_prs_only = CSLB ~ PRS (공변량 없이) · r2_base = CSLB ~ 공변량 · r2_full = CSLB ~ 공변량 + PRS · incr_r2 = r2_full − r2_base (PRS 고유 예측력)\n")
cat("\n■ (나) 젊은 개 → 노령견\n"); if (nrow(y2o_ev)) print(f3(y2o_ev[, .(analysis, n, incr_r2, beta_prs, p, cor_prs_resid, h2_train_mean)])) else cat("  결과 없음\n")
cat("\n■ null (y를 나이 집단 안에서 섞음)\n"); print(f3(nul[, .(perm, group, n, incr_r2, p, incr_r2_raw, p_raw)]))
cat(sprintf("  null p < 0.05 비율: 수정 방식 %.2f · 첫 실행 방식 %.2f (기대 0.05)\n", mean(nul$p < 0.05), mean(nul$p_raw < 0.05)))
cat("\n■ 판정 (예측됨 = 5번 반복 모두 p < 0.05 & 평균 증분 R² > null 최대값; p_emp = null 중 관측 평균 이상 비율)\n"); print(f3(verdict))
cat(sprintf("  (나) young8 → age8: %s\n", if (nrow(y2o_ev)) sprintf("증분 R² %.2g, p %.2g → %s", y2o_ev$incr_r2, y2o_ev$p, if (y2o_ev$p < 0.05) "예측됨 (탐색적)" else "예측 안 됨") else "결과 없음"))
cat("  사전 기대: 증분 R² 약 1e-4 (h² 0.018, 학습 약 6,400마리), 검출 한계 약 1e-3. 예측 안 됨이면 06-6 작은 h²의 예측 쪽 확인\n")
sink()
cat(sprintf("저장: %s_folds.csv · _oof_prs_rep*.csv · _eval.csv · _young2old.csv · _null.csv · _summary.txt / 그림 %s_r2.png\n", out_prefix, pic_prefix))

