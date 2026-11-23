# 11-2_purified_h2.R - SNP heritability after excluding the flagged dogs. See README.md.
library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_keep <- as.integer(Sys.getenv("TEST_KEEP", "0"))
n_rand <- c(age8 = as.integer(Sys.getenv("RAND_AGE8", "50")), all = as.integer(Sys.getenv("RAND_ALL", "10")))
tag <- if (test_keep > 0) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("11-2", tag)); pic_prefix <- file.path(pic_dir, paste0("11-2", tag))
work_dir <- paste0(out_prefix, "_work")
flags_f <- c(file.path(in_dir, "11-1_flags.csv"), file.path(qc_dir, "11-1_flags.csv")); flags_f <- flags_f[file.exists(flags_f)][1]
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
grm <- file.path(in_dir, "04-1_grm", "grm_all")
pheno_f <- file.path(in_dir, "04-2_pheno.txt"); qcov_f <- file.path(in_dir, "03-5_qcovar_PC5.txt"); cov_f <- file.path(in_dir, "03-5_covar.txt")
for (f in c(flags_f, gcta, paste0(grm, ".grm.bin"), pheno_f, qcov_f, cov_f)) if (is.na(f) || !file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
wr <- function(d, f) { fwrite(d, f, sep = " ", col.names = FALSE, quote = FALSE); f }

fl <- fread(flags_f, colClasses = c(FID = "character", IID = "character"))
need <- c("excl_primary", "excl_primary_any", "excl_sensory", "excl_severe", "excl_broad")
if (!all(need %in% names(fl))) stop("11-1 출력이 개정 전 형식: 11-1을 다시 실행 (", flags_f, ")")
ph <- fread(pheno_f, header = FALSE, colClasses = list(character = 1:2)); setnames(ph, c("FID", "IID", "rint", "score"))
qcv <- fread(qcov_f, header = FALSE, colClasses = list(character = 1:2)); cvv <- fread(cov_f, header = FALSE, colClasses = "character")
d <- merge(ph[, .(FID, IID, score)], fl[, c("IID", "age8", need), with = FALSE], by = "IID")
d <- merge(d, qcv[, .(IID = V2, age = V3)], by = "IID")
d[, age_bin := pmin(floor(age), 15)]
if (test_keep > 0) d <- d[seq(1, .N, by = test_keep)]
set.seed(20261113)
samples <- list()
for (g in c("all", "age8")) {
  base <- if (g == "all") d else d[age8 == TRUE]
  samples[[paste0(g, "_unpurified")]] <- base$IID
  samples[[paste0(g, "_primary")]] <- base[excl_primary == FALSE, IID]
  samples[[paste0(g, "_primary_anytime")]] <- base[excl_primary_any == FALSE, IID]
  samples[[paste0(g, "_sensory_only")]] <- base[excl_sensory == FALSE, IID]
  samples[[paste0(g, "_severe_sensory")]] <- base[excl_severe == FALSE, IID]
  samples[[paste0(g, "_broad")]] <- base[excl_broad == FALSE, IID]
  drop_n <- base[excl_primary == TRUE, .N, by = age_bin]
  for (r in seq_len(n_rand[[g]])) {
    out <- unlist(lapply(seq_len(nrow(drop_n)), function(i) { pool <- base[age_bin == drop_n$age_bin[i], IID]; sample(pool, min(length(pool), drop_n$N[i])) }))
    samples[[sprintf("%s_random_drop_%02d", g, r)]] <- setdiff(base$IID, out)
  }
}
reml <- function(name, keep_ids) {
  out <- file.path(work_dir, name)
  if (!file.exists(paste0(out, ".hsq"))) {
    s <- d[IID %in% keep_ids][, y := rint(score)]
    wr(s[, .(FID, IID, y)], paste0(out, ".pheno")); wr(s[, .(FID, IID)], paste0(out, ".keep"))
    q <- qcv[V2 %in% s$IID]; q <- q[, c(TRUE, TRUE, vapply(q[, -(1:2)], function(x) isTRUE(sd(x) > 0), TRUE)), with = FALSE]
    k <- cvv[V2 %in% s$IID]; k <- k[, c(TRUE, TRUE, vapply(k[, -(1:2)], function(x) uniqueN(x) > 1, TRUE)), with = FALSE]
    wr(q, paste0(out, ".qcovar")); wr(k, paste0(out, ".covar"))
    system2(gcta, c("--reml", "--grm", grm, "--keep", paste0(out, ".keep"), "--pheno", paste0(out, ".pheno"), "--qcovar", paste0(out, ".qcovar"), "--covar", paste0(out, ".covar"),
                    "--reml-maxit", "200", "--out", out, "--threads", n_cores), stdout = paste0(out, ".stdout"), stderr = paste0(out, ".stdout"))
  }
  if (!file.exists(paste0(out, ".hsq"))) { warning(name, ": GREML 실패"); return(NULL) }
  h <- fread(paste0(out, ".hsq"), fill = TRUE); g <- function(t) suppressWarnings(as.numeric(h[Source == t, ]))
  data.table(sample = name, n = g("n")[2], h2 = g("V(G)/Vp")[2], se = g("V(G)/Vp")[3], p_LRT = g("Pval")[2])
}
cat(sprintf("GREML %d번 (8세 이상 무작위 %d · 전체 무작위 %d)\n", length(samples), n_rand[["age8"]], n_rand[["all"]]))
res <- rbindlist(Map(reml, names(samples), samples))
res[, `:=`(group = sub("_.*", "", sample), version = sub("_[0-9]+$", "", sub("^(all|age8)_", "", sample)))]
fwrite(res, paste0(out_prefix, "_h2.csv"))

rd <- res[version == "random_drop", .(n_rand = .N, h2_q50 = median(h2), h2_q95 = quantile(h2, 0.95), h2_mean = mean(h2)), by = group]
judge <- function(g) {
  pr <- res[sample == paste0(g, "_primary")]; un <- res[sample == paste0(g, "_unpurified")]; rr <- res[group == g & version == "random_drop", h2]
  q95 <- if (length(rr)) quantile(rr, 0.95) else NA_real_; emp <- if (length(rr)) (1 + sum(rr >= pr$h2)) / (length(rr) + 1) else NA_real_
  data.table(group = g, h2_unpurified = un$h2, p_unpurified = un$p_LRT, h2_primary = pr$h2, se_primary = pr$se, p_primary = pr$p_LRT,
             rand_q95 = q95, emp_p_vs_random = emp, reveal = isTRUE(pr$h2 > un$h2 && pr$p_LRT < 0.05 && pr$h2 > q95))
}
jd <- rbindlist(lapply(c("age8", "all"), judge))
fwrite(jd, paste0(out_prefix, "_judge.csv"))

fix <- res[version != "random_drop"]
fix[, version := factor(version, levels = c("unpurified", "primary", "primary_anytime", "sensory_only", "severe_sensory", "broad"))]
g1 <- ggplot() + geom_hline(yintercept = 0, colour = "grey60") +
  geom_boxplot(data = res[version == "random_drop"], aes(x = "random_drop\n(age-matched)", y = h2), colour = "grey45", outlier.size = 0.6, width = 0.5) +
  geom_pointrange(data = fix, aes(version, h2, ymin = pmax(0, h2 - 1.96 * se), ymax = h2 + 1.96 * se), colour = "#1f4e9c") + facet_wrap(~ group) +
  scale_x_discrete(limits = c(levels(fix$version), "random_drop\n(age-matched)")) +
  labs(x = NULL, y = "SNP h² (95% CI, lower bound truncated at 0)", title = "CSLB SNP heritability after removing dogs with sensory / neurologic / pain / systemic conditions",
       subtitle = "Box: same number of dogs removed at random within 1-year age bins") +
  theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(paste0(pic_prefix, "_h2.png"), g1, width = 9, height = 4.5, dpi = 300)
f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("11-2 가설 H8: 정제 CSLB의 SNP 유전력 (측정 가설 · 추정; 해석 규칙은 실행 전 헤더에 고정, 2026-11-14 개정)\n"); if (test_keep > 0) cat("*** 시험. 해석하지 말 것 ***\n")
cat("\n■ 표본별 h² (표본 안에서 RINT 다시 계산; 무작위 제외는 아래 요약)\n"); print(f3(fix[order(group, version), .(group, version, n, h2, se, p_LRT)]))
cat("\n■ 무작위 제외 (주 정제와 같은 수 · 같은 나이 분포) h² 분포\n"); print(f3(rd))
cat("\n■ 판정 (정제 후 > 정제 전 · 정제 후 p < 0.05 · 정제 후 > 무작위 95 백분위 — 셋 다)\n"); print(f3(jd))
a8j <- jd[group == "age8"]
cat(sprintf("■ 8세 이상 주 판정: %s%s\n", if (a8j$reveal) "정제가 유전 신호를 드러냄 (탐색적)" else "정제로 유전 신호가 드러나지 않음",
            if (isTRUE(a8j$p_unpurified < 0.05)) { if (a8j$reveal) " (정제 전 h²도 이미 p < 0.05 — '드러냄'이 아니라 '키움'으로 읽음)" else " (참고: 정제 전 h²는 p < 0.05)" } else ""))
sink()
cat(sprintf("저장: %s_h2.csv · _judge.csv · _summary.txt / 그림 %s_h2.png\n", out_prefix, pic_prefix))

