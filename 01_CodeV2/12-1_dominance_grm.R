# 12-1_dominance_grm.R - Additive and dominance GRMs. See README.md.
library(data.table)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_mode <- Sys.getenv("TEST", "0") == "1"
tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("12-1", tag))
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
bfile <- file.path(in_dir, "03-3_maf0.01_hwe"); prune <- file.path(in_dir, "03-4_pca", "prune.prune.in"); grm_a <- file.path(in_dir, "04-1_grm", "grm_all")
for (f in c(gcta, paste0(bfile, ".bed"), prune, paste0(grm_a, c(".grm.bin", ".grm.id")))) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

ids <- fread(paste0(grm_a, ".grm.id"), header = FALSE, colClasses = "character")
if (test_mode) ids <- ids[seq(1, .N, by = 8)]
kf <- paste0(out_prefix, "_keep.txt"); fwrite(ids, kf, sep = " ", col.names = FALSE)
mk <- function(out, what, extract = NULL) {
  done <- paste0(out, if (what == "--make-grm-d") ".d" else "", ".grm.bin")
  if (!file.exists(done)) {
    cat(sprintf("[%s] %s → %s\n", format(Sys.time(), "%F %T"), what, basename(out)))
    system2(gcta, c("--bfile", bfile, if (!is.null(extract)) c("--extract", extract), "--keep", kf, "--autosome-num", "38", "--autosome",
                    if (test_mode) c("--chr", "38"), what, "--out", out, "--thread-num", n_cores), stdout = paste0(out, ".stdout"), stderr = paste0(out, ".stdout"))
  }
  if (!file.exists(done)) stop("GRM 실패: ", out, ".stdout")
  sub("\\.grm\\.bin$", "", done)
}
g_dA <- mk(paste0(out_prefix, "_dense_A"), "--make-grm")
g_dD <- mk(paste0(out_prefix, "_dense_D"), "--make-grm-d")
g_pD <- mk(paste0(out_prefix, "_pruned_D"), "--make-grm-d", extract = prune)
g_pA <- if (test_mode) mk(paste0(out_prefix, "_pruned_A"), "--make-grm", extract = prune) else grm_a

rd_grm <- function(p) { n <- nrow(fread(paste0(p, ".grm.id"), header = FALSE)); readBin(paste0(p, ".grm.bin"), "numeric", n = n * (n + 1) / 2, size = 4) }
n_snp <- function(p) { l <- readLines(paste0(p, ".stdout"), warn = FALSE); x <- grep("SNPs have been processed|SNPs remain|SNPs to be included", l, value = TRUE); tail(x, 1) }
n <- nrow(ids); dg <- cumsum(seq_len(n)); off <- setdiff(seq_len(n * (n + 1) / 2), dg)
summ <- function(set, pa, pd) {
  xa <- rd_grm(pa); xd <- rd_grm(pd)
  data.table(set = set, diag_A_mean = mean(xa[dg]), diag_D_mean = mean(xd[dg]), diag_D_sd = sd(xd[dg]), off_D_sd = sd(xd[off]),
             cor_offdiag_A_D = cor(xa[off], xd[off]), cor_diag_A_D = cor(xa[dg], xd[dg]))
}
S <- rbind(summ("dense (전체 SNP, 주)", g_dA, g_dD), summ("pruned (민감도)", g_pA, g_pD))
f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("12-1 상가 · dominance GRM (GCTA --make-grm / --make-grm-d; 가설 H9 준비)\n"); if (test_mode) cat("*** TEST=1 시험 (chr38 · 일부 개). 해석하지 말 것 ***\n")
cat(sprintf("개 %d마리\n", n)); for (p in c(g_dA, g_dD, g_pD)) cat(sprintf("  %s: %s\n", basename(p), n_snp(sub("\\.d$", "", p))))
print(f3(S))
cat("비대각 상관이 높을수록 (예: > 0.85) 두 성분을 나누기 어려움 — 12-2가 다시 계산해 경고\n")
sink()
fwrite(S, paste0(out_prefix, "_grm_summary.csv"))
cat(sprintf("저장: %s_dense_A.grm.* · _dense_D.d.grm.* · _pruned_D.d.grm.* · _summary.txt\n", out_prefix))

