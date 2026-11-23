# 07-1_roh.R - Runs of homozygosity and F_ROH per dog. See README.md.
07
#   --homozyg-window-snp 50 · --homozyg-window-het 2 · --homozyg-window-missing 5 · --homozyg-window-threshold 0.05

library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_mode <- Sys.getenv("TEST", "0") == "1"
tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("07-1", tag)); pic_prefix <- file.path(pic_dir, paste0("07-1", tag))
work_dir <- paste0(out_prefix, "_work")
plink <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
bfile <- file.path(in_dir, "03-3_maf0.01_hwe"); prune <- file.path(in_dir, "03-4_pca", "prune.prune.in")
dog_file <- file.path(base_dir, "DogOverview/v2023/DAP_2023_DogOverview_v1.1.csv")
for (f in c(plink, gcta, paste0(bfile, c(".bed", ".bim", ".fam")), prune, dog_file)) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
sh <- function(cmd, args, log) { st <- system2(cmd, args, stdout = log, stderr = log); invisible(st) }
roh_args <- c("--homozyg", "--homozyg-kb", "500", "--homozyg-snp", "50", "--homozyg-density", "50", "--homozyg-gap", "1000",
              "--homozyg-window-snp", "50", "--homozyg-window-het", "2", "--homozyg-window-missing", "5", "--homozyg-window-threshold", "0.05")
cls <- list(short = c(0.5, 2), medium = c(2, 8), long = c(8, Inf))           # Mb

# 1. ROH (plink 1.9)
sub_args <- character(0)
if (test_mode) {
  fam <- fread(paste0(bfile, ".fam"), header = FALSE)
  kf <- file.path(work_dir, "test_keep.txt"); fwrite(fam[seq(1, .N, by = 8), .(V1, V2)], kf, sep = " ", col.names = FALSE)
  sub_args <- c("--chr", "38", "--keep", kf)
}
roh_out <- file.path(work_dir, "roh")
if (!file.exists(paste0(roh_out, ".hom"))) {
  t0 <- Sys.time(); cat("plink --homozyg 실행\n")
  sh(plink, c("--bfile", bfile, "--dog", if (test_mode) sub_args else "--autosome", roh_args, "--threads", n_cores, "--out", roh_out), paste0(roh_out, ".stdout"))
  cat(sprintf("  끝 (%.1f분)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (!file.exists(paste0(roh_out, ".hom"))) stop("ROH 실패: ", roh_out, ".stdout")
hom <- fread(paste0(roh_out, ".hom")); indiv <- fread(paste0(roh_out, ".hom.indiv"))
bim <- fread(paste0(bfile, ".bim"), header = FALSE, select = c(1, 4), col.names = c("chr", "bp"))
bim <- bim[chr %in% 1:38]; if (test_mode) bim <- bim[chr == 38]
auto_len <- bim[, .(len = (max(bp) - min(bp)) / 1e3), by = chr][, sum(len)]   # kb
hom[, mb := KB / 1e3]
fr <- indiv[, .(FID = as.character(FID), IID = as.character(IID), n_roh = NSEG, roh_kb = KB)]
fr[, F_ROH := roh_kb / auto_len]
for (k in names(cls)) {
  s <- hom[mb >= cls[[k]][1] & mb < cls[[k]][2], .(kb = sum(KB), n = .N), by = .(IID = as.character(IID))]
  fr <- merge(fr, s[, .(IID, tmp = kb / auto_len, tmpn = n)], by = "IID", all.x = TRUE)
  setnames(fr, c("tmp", "tmpn"), c(paste0("F_ROH_", k), paste0("n_roh_", k)))
}
for (j in grep("^(F_ROH_|n_roh_)", names(fr), value = TRUE)) set(fr, which(is.na(fr[[j]])), j, 0)

ibc_out <- file.path(work_dir, "ibc")
if (!file.exists(paste0(ibc_out, ".ibc")))
  sh(gcta, c("--bfile", bfile, "--extract", prune, "--autosome-num", "38", "--autosome", if (test_mode) c("--keep", kf), "--ibc", "--out", ibc_out,
             "--thread-num", n_cores), paste0(ibc_out, ".stdout"))
ibc <- fread(paste0(ibc_out, ".ibc"), colClasses = list(character = 1:2))
fr <- merge(fr, ibc[, .(IID, Fhat1, Fhat2, Fhat3)], by = "IID", all.x = TRUE)
dog <- fread(dog_file, select = c("dog_id", "Breed_Status", "Breed"), colClasses = c(dog_id = "character"))
fr <- merge(fr, dog, by.x = "IID", by.y = "dog_id", all.x = TRUE)
fwrite(fr, paste0(out_prefix, "_froh.csv"))

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
desc <- fr[, .(n = .N, F_ROH_mean = mean(F_ROH), F_ROH_median = median(F_ROH), F_ROH_q90 = quantile(F_ROH, 0.9),
               short = mean(F_ROH_short), medium = mean(F_ROH_medium), long = mean(F_ROH_long), Fhat3_mean = mean(Fhat3, na.rm = TRUE)), by = Breed_Status]
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("07-1 게놈 근친도 F_ROH (plink 1.9 --homozyg; 가설 H3 준비)\n")
if (test_mode) cat("*** TEST=1: chr38 · 8마리마다 1마리. 해석하지 말 것 ***\n")
cat(sprintf("개 %d마리, 상염색체 길이 (SNP가 덮은 구간) %.0f Mb, ROH %s개\n", nrow(fr), auto_len / 1e3, format(nrow(hom), big.mark = ",")))
cat("ROH 설정: ", paste(roh_args[-1], collapse = " "), "\n\n", sep = "")
cat("■ 순종 / 잡종별 F_ROH (전체 · 길이별 평균)과 GRM 기반 F̂3\n"); print(f3(desc))
cat(sprintf("\n■ F_ROH와 F̂3의 상관: Pearson %.3f · Spearman %.3f\n", cor(fr$F_ROH, fr$Fhat3, use = "complete.obs"),
            cor(fr$F_ROH, fr$Fhat3, use = "complete.obs", method = "spearman")))
cat("  (F̂3는 표본 전체 allele 빈도 기준이라 혼합 표본에서는 음수도 나옴. 상관이 높으면 ROH 설정이 근친도를 잘 반영)\n")
sink()
pd <- copy(fr)[!is.na(Breed_Status)]
g1 <- ggplot(pd, aes(F_ROH, fill = Breed_Status)) + geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  labs(x = "F_ROH (total ROH length / autosome length)", y = "dogs", fill = NULL, title = "Genomic inbreeding (F_ROH)") + theme_minimal(base_size = 10)
ggsave(paste0(pic_prefix, "_froh.png"), g1, width = 7, height = 4.5, dpi = 300)
cat(sprintf("저장: %s_froh.csv · _summary.txt / 그림 %s_froh.png\n", out_prefix, pic_prefix))

