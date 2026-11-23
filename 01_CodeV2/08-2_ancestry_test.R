# 08-2_ancestry_test.R - Variance in CSLB explained by breed ancestry (REML, LRT). See README.md.
library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_mode <- Sys.getenv("TEST", "0") == "1"
tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("08-2", tag)); pic_prefix <- file.path(pic_dir, paste0("08-2", tag))
work_dir <- paste0(out_prefix, "_work")
p81 <- if (test_mode) file.path(qc_dir, "08-1_test") else file.path(in_dir, "08-1")
kanc <- paste0(p81, "_kanc"); apc_f <- paste0(p81, "_ancestry_pcs.txt"); wide_f <- paste0(p81, "_ancestry_wide.csv")
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
grm <- file.path(in_dir, "04-1_grm", "grm_all")
a8 <- file.path(in_dir, "05-1_inputs", paste0("age8_", c("pheno", "qcovar", "covar"), ".txt")); names(a8) <- c("pheno", "qcovar", "covar")
main_f <- c(pheno = file.path(in_dir, "04-2_pheno.txt"), qcovar = file.path(in_dir, "03-5_qcovar_PC5.txt"), covar = file.path(in_dir, "03-5_covar.txt"))
for (f in c(paste0(kanc, ".grm.bin"), apc_f, wide_f, gcta, paste0(grm, ".grm.bin"), a8, main_f)) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
qn <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5)); cn <- c("Sex", "Intact", "Vision", "Cohort_BHS")
q_nopc <- qn[1:5]
wr <- function(d, f) { fwrite(d, f, sep = " ", col.names = FALSE, quote = FALSE); f }
rd <- function(f, nm) { d <- fread(f, header = FALSE, colClasses = list(character = 1:2)); setnames(d, c("FID", "IID", nm)); d }
build <- function(fl) Reduce(function(a, b) merge(a, b, by = c("FID", "IID")),
  list(rd(fl["pheno"], if (grepl("04-2", fl["pheno"])) c("y", "score") else "y"), rd(fl["qcovar"], qn), fread(fl["covar"], header = FALSE, colClasses = "character", col.names = c("FID", "IID", cn))))
kid <- fread(paste0(kanc, ".grm.id"), header = FALSE, colClasses = "character")$V2
A <- build(a8)[IID %in% kid]; M <- build(main_f)[IID %in% kid]
apc <- rd(apc_f, paste0("aPC", 1:20)); A <- merge(A, apc, by = c("FID", "IID")); M <- merge(M, apc, by = c("FID", "IID"))
cat(sprintf("8세 이상 %d마리 · 전체 %d마리\n", nrow(A), nrow(M)))

reml <- function(name, d, grms, qcols, lrt = TRUE) {
  out <- file.path(work_dir, name)
  if (!file.exists(paste0(out, ".hsq"))) {
    kq <- qcols[vapply(qcols, function(j) isTRUE(sd(d[[j]]) > 0), TRUE)]; kc <- cn[vapply(cn, function(j) uniqueN(d[[j]]) > 1, TRUE)]
    wr(d[, .(FID, IID, y)], paste0(out, ".pheno")); wr(d[, c("FID", "IID", kq), with = FALSE], paste0(out, ".qcovar"))
    wr(d[, c("FID", "IID", kc), with = FALSE], paste0(out, ".covar")); wr(d[, .(FID, IID)], paste0(out, ".keep"))
    g <- if (length(grms) == 1) c("--grm", grms) else { mf <- paste0(out, ".mgrm"); writeLines(grms, mf); c("--mgrm", mf, if (!lrt) "--reml-no-lrt", "--reml-no-constrain") }
    if (length(grms) == 1 && !lrt) g <- c(g, "--reml-no-lrt")
    system2(gcta, c("--reml", g, "--keep", paste0(out, ".keep"), "--pheno", paste0(out, ".pheno"), "--qcovar", paste0(out, ".qcovar"), "--covar", paste0(out, ".covar"),
                    "--reml-maxit", "200", "--out", out, "--threads", n_cores), stdout = paste0(out, ".stdout"), stderr = paste0(out, ".stdout"))
  }
  if (!file.exists(paste0(out, ".hsq"))) { warning(name, ": REML 실패 (", out, ".stdout)"); return(NULL) }
  h <- fread(paste0(out, ".hsq"), fill = TRUE, header = TRUE)
  g <- function(t) suppressWarnings(as.numeric(h[Source == t, ]))
  if (length(grms) == 1) {
    data.table(model = name, n = nrow(d), component = "anc", prop = g("V(G)/Vp")[2], se = g("V(G)/Vp")[3], p = g("Pval")[2], test = if (lrt) "LRT" else "Wald")
  } else {
    rbind(data.table(model = name, n = nrow(d), component = "anc", prop = g("V(G1)/Vp")[2], se = g("V(G1)/Vp")[3]),
          data.table(model = name, n = nrow(d), component = "SNP_GRM", prop = g("V(G2)/Vp")[2], se = g("V(G2)/Vp")[3]),
          data.table(model = name, n = nrow(d), component = "anc + SNP (sum)", prop = g("Sum of V(G)/Vp")[2], se = g("Sum of V(G)/Vp")[3]))[, `:=`(p = pnorm(prop / se, lower.tail = FALSE), test = "Wald")]
  }
}
res <- rbindlist(list(
  reml("P_age8_anc", A, kanc, q_nopc),
  reml("S1_age8_anc_plus_GRM", A, c(kanc, grm), q_nopc, lrt = FALSE),
  reml("S4_age8_anc_noWeight", A, kanc, setdiff(q_nopc, "wt_cov")),
  reml("S5_all_anc", M, kanc, q_nopc),
  reml("S6_age8_anc_withPC", A, kanc, qn)), fill = TRUE)

kc <- cn[vapply(cn, function(j) uniqueN(A[[j]]) > 1, TRUE)]
m0 <- lm(as.formula(paste("y ~", paste(c(q_nopc, kc), collapse = " + "))), data = A)
m1 <- update(m0, as.formula(paste(". ~ . +", paste(paste0("aPC", 1:20), collapse = " + "))))
s2 <- data.table(model = "S2_age8_aPC20_block", n = nrow(A), component = "aPC1-20 (fixed)", prop = summary(m1)$r.squared - summary(m0)$r.squared,
                 se = NA_real_, p = anova(m0, m1)$`Pr(>F)`[2], test = "F")
wide <- fread(wide_f, colClasses = c(IID = "character", FID = "character")); bcols <- setdiff(names(wide), c("FID", "IID"))
wm <- as.matrix(wide[, ..bcols]); top <- data.table(IID = wide$IID, breed = bcols[max.col(wm, ties.method = "first")], pmax = apply(wm, 1, max))
A[, resid := resid(m0)]
pb <- merge(A[, .(IID, resid)], top[pmax >= 0.85], by = "IID")
keepb <- pb[, .N, by = breed][N >= 10, breed]; pb <- pb[breed %in% keepb]
s3 <- if (length(keepb) >= 2) {
  av <- anova(lm(resid ~ breed, data = pb)); msb <- av$`Mean Sq`[1]; msw <- av$`Mean Sq`[2]; n0 <- (nrow(pb) - sum(table(pb$breed)^2) / nrow(pb)) / (length(keepb) - 1)
  data.table(model = "S3_age8_purebred_breed_means", n = nrow(pb), component = sprintf("%d breeds (ICC)", length(keepb)), prop = (msb - msw) / (msb + (n0 - 1) * msw),
             se = NA_real_, p = av$`Pr(>F)`[1], test = "F (one-way)")
} else data.table(model = "S3_age8_purebred_breed_means", n = nrow(pb), component = "품종 2개 미만", prop = NA_real_, se = NA_real_, p = NA_real_, test = "F")
res <- rbind(res, s2, s3, fill = TRUE)
fwrite(res, paste0(out_prefix, "_results.csv"))
bm <- pb[, .(n = .N, mean_resid = mean(resid), se = sd(resid) / sqrt(.N)), by = breed][order(mean_resid)]
fwrite(bm, paste0(out_prefix, "_breed_means.csv"))
if (nrow(bm)) {
  g1 <- ggplot(bm, aes(mean_resid, reorder(breed, mean_resid))) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_pointrange(aes(xmin = mean_resid - 1.96 * se, xmax = mean_resid + 1.96 * se), colour = "#1f4e9c", size = 0.3) +
    labs(x = "CSLB (RINT) residual after covariates (no genotype PCs), mean ± 95% CI", y = NULL, title = "Dogs aged 8+: genomic purebreds (ancestry ≥ 0.85), breeds with ≥ 10 dogs") +
    theme_minimal(base_size = 9)
  ggsave(paste0(pic_prefix, "_breed_means.png"), g1, width = 7, height = max(3, 0.18 * nrow(bm) + 1), dpi = 300)
}
P <- res[model == "P_age8_anc"]; H4 <- isTRUE(P$p < 0.05)
f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("08-2 가설 H4: 품종 계통과 노령기 CSLB (2026-11-13 추가, 판정 규칙은 실행 전 헤더에 고정)\n"); if (test_mode) cat("*** TEST=1 시험. 해석하지 말 것 ***\n")
cat("\n■ 결과 (prop = 계통 (또는 SNP GRM) 분산 비율 · 블록 partial R² · ICC)\n"); print(f3(res))
cat(sprintf("\n■ 판정 (주 검정 P_age8_anc): V(anc)/Vp %.3g ± %.3g, LRT p %.3g → %s\n", P$prop, P$se, P$p,
            if (H4) "H4 지지 (Holm 보정 전; '품종 계통에 따른 차이'까지만 해석)" else "H4 지지 안 됨"))
cat("  해석 경계: V(anc)/Vp는 유전력이 아님 (품종 유전 · 체격 · 수명 · 생존 선택 · 보호자 응답 습관이 함께 들어 있음)\n",
    "  S6 (PC 포함)과 주 검정의 차이 = 지금까지의 PC 보정이 지운 품종 (+ 지역) 몫. S1은 두 성분 분리가 불안정해 합 (anc + SNP)을 함께 봄\n",
    "  S3 (품종 평균 분산분석)은 품종 안 친척 때문에 관대한 쪽 → 기술용\n", sep = "")
sink()
cat(sprintf("저장: %s_results.csv · _breed_means.csv · _summary.txt / 그림 %s_breed_means.png\n", out_prefix, pic_prefix))

