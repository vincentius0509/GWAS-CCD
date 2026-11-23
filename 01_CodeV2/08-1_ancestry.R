# 08-1_ancestry.R - Builds a breed-ancestry similarity matrix from DAP ancestry proportions. See README.md.
library(data.table)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_keep <- as.integer(Sys.getenv("TEST_KEEP", "0"))
tag <- if (test_keep > 0) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("08-1", tag))
anc_f <- file.path(base_dir, "GeneticData/v2023/DogAgingProject_GeneticData_CuratedRelease_2023.ancestry.tsv")
id_f <- file.path(in_dir, "04-1_grm", "grm_all.grm.id"); qcov_f <- file.path(in_dir, "03-5_qcovar_PC5.txt")
for (f in c(anc_f, id_f, qcov_f)) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
n_apc <- 20

ids <- fread(id_f, header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
if (test_keep > 0) ids <- ids[seq(1, .N, by = test_keep)]
a <- fread(anc_f, colClasses = c(dog = "character"))[, pct := as.numeric(pct)][dog %in% ids$IID]
P <- dcast(a, dog ~ pop, value.var = "pct", fill = 0)
P <- P[match(ids$IID, dog)]
if (anyNA(P$dog)) stop("품종 비율이 없는 개: ", sum(is.na(P$dog)))
pm <- as.matrix(P[, -1]); rownames(pm) <- P$dog
pm <- pm[, colSums(pm > 0) >= 1, drop = FALSE]
cat(sprintf("개 %d마리 × 품종 %d개, 비율 합 %.2f ~ %.2f (나머지는 배정 안 됨)\n", nrow(pm), ncol(pm), min(rowSums(pm)), max(rowSums(pm))))
fwrite(cbind(ids, as.data.table(pm)), paste0(out_prefix, "_ancestry_wide.csv"))

Z <- scale(pm); Z[!is.finite(Z)] <- 0
K <- tcrossprod(Z) / ncol(pm)
v <- K[upper.tri(K, diag = TRUE)]
con <- file(paste0(out_prefix, "_kanc.grm.bin"), "wb"); writeBin(as.numeric(v), con, size = 4); close(con)
con <- file(paste0(out_prefix, "_kanc.grm.N.bin"), "wb"); writeBin(rep(as.numeric(ncol(pm)), length(v)), con, size = 4); close(con)
fwrite(ids, paste0(out_prefix, "_kanc.grm.id"), sep = "\t", col.names = FALSE)

pc <- prcomp(pm, center = TRUE, scale. = FALSE)
ve <- pc$sdev^2 / sum(pc$sdev^2)
apc <- data.table(ids, pc$x[, 1:n_apc]); setnames(apc, c("FID", "IID", paste0("aPC", 1:n_apc)))
fwrite(apc, paste0(out_prefix, "_ancestry_pcs.txt"), sep = " ", col.names = FALSE)
gpc <- fread(qcov_f, header = FALSE, colClasses = list(character = 1:2)); setnames(gpc, 1:2, c("FID", "IID"))
gpc <- gpc[match(ids$IID, IID), .SD, .SDcols = (ncol(gpc) - 4):ncol(gpc)]; setnames(gpc, paste0("PC", 1:5))
cc <- cor(as.matrix(gpc), pc$x[, 1:10])

sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("08-1 품종 계통 행렬 (가설 H4 준비)\n"); if (test_keep > 0) cat(sprintf("*** TEST_KEEP=%d 시험. 해석하지 말 것 ***\n", test_keep))
cat(sprintf("개 %d마리 × 품종 %d개. 유전체 순종형 (가장 큰 품종 비율 ≥ 0.85) %.1f%%\n", nrow(pm), ncol(pm), 100 * mean(apply(pm, 1, max) >= 0.85)))
cat(sprintf("K_anc: 대각 평균 %.3f, 비대각 평균 %.4f · SD %.4f\n", mean(diag(K)), mean(K[upper.tri(K)]), sd(K[upper.tri(K)])))
cat(sprintf("계통 PC 설명 분산: PC1 ~ 5 %s, PC1 ~ %d 합 %.1f%%\n", paste(sprintf("%.1f%%", 100 * ve[1:5]), collapse = " · "), n_apc, 100 * sum(ve[1:n_apc])))
cat("유전형 PC1 ~ 5 (03-5)와 계통 PC1 ~ 10의 상관 (절댓값 최대; 유전형 PC가 계통을 얼마나 담는지)\n")
print(round(apply(abs(cc), 1, max), 3))
sink()
cat(sprintf("저장: %s_ancestry_wide.csv · _ancestry_pcs.txt · _kanc.grm.* · _summary.txt\n", out_prefix))

