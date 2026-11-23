# 01-1_missing_threshold.R - Examines variant missingness distributions to choose the --geno threshold. See README.md.
library(data.table)
library(ggplot2)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")    
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
dir.create(pic_dir, recursive = TRUE, showWarnings = FALSE)
in_prefix <- file.path(qc_dir, "01-1_raw")
out_prefix <- file.path(qc_dir, "01-1")
pic_prefix <- file.path(pic_dir, "01-1")
f_kin <- file.path(base_dir, "GeneticData/v2023/DogAgingProject_GeneticData_CuratedRelease_2023.kinship.tsv")
f_info <- file.path(qc_dir, "00_all_total.csv")

geno_cand <- c(0.02, 0.05, 0.1, 0.2)
maf_ref <- 0.05
kin_dup <- 0.354

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

vmiss <- fread(paste0(in_prefix, ".vmiss"), select = c("#CHROM", "ID", "F_MISS"), col.names = c("CHR", "SNP", "F_MISS"))
frq <- fread(paste0(in_prefix, ".afreq"), select = c("ID", "ALT_FREQS"), col.names = c("SNP", "ALT_FREQS"))
pvar <- fread(paste0(in_prefix, ".pvar"), select = c("ID", "POS"), col.names = c("SNP", "BP"))
stopifnot(identical(vmiss$SNP, frq$SNP), identical(vmiss$SNP, pvar$SNP))
vmiss[, `:=`(MAF = pmin(frq$ALT_FREQS, 1 - frq$ALT_FREQS), BP = pvar$BP)]
rm(frq, pvar); gc()

cat("SNP 수:", nrow(vmiss), " 염색체 수:", uniqueN(vmiss$CHR), " MAF 계산 불가(전원 결측):", vmiss[, sum(is.na(MAF))], "\n")

cat(sprintf("전체 genotyping rate (분석 대상, 필터 전): %.3f\n", 1 - mean(vmiss$F_MISS)))

print(summary(vmiss$F_MISS))
print(quantile(vmiss$F_MISS, c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)))
print(vmiss[, .(all_missing = sum(F_MISS == 1), no_missing = sum(F_MISS == 0))])

vmiss_bin <- vmiss[, .(n = .N), by = .(bin = pmin(floor(F_MISS / 0.01) * 0.01, 0.99))]
p_hist <- ggplot(vmiss_bin, aes(x = bin + 0.005, y = n / 1e6)) +
  geom_col(width = 0.01, fill = "steelblue", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = geno_cand, color = "red", linetype = "dashed") +
  annotate("text", x = geno_cand, y = Inf, label = geno_cand, vjust = c(1.5, 3.5, 1.5, 3.5), hjust = -0.1, color = "red") +
  labs(title = "SNP missing rate (autosomes, before any filter)",
       x = "Missing rate (F_MISS)", y = "Number of SNPs (million)") +
  theme_qc
ggsave(paste0(pic_prefix, "_vmiss_hist.png"), p_hist, width = 9, height = 5, dpi = 300)

vmiss_cum <- data.table(thr = seq(0, 1, by = 0.01))
vmiss_cum[, kept := findInterval(thr, sort(vmiss$F_MISS))]
p_cum <- ggplot(vmiss_cum, aes(x = thr, y = kept / 1e6)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_vline(xintercept = geno_cand, color = "red", linetype = "dashed") +
  labs(title = "SNPs retained by --geno threshold",
       x = "--geno threshold (remove SNPs with missing rate > threshold)", y = "SNPs retained (million)") +
  theme_qc
ggsave(paste0(pic_prefix, "_vmiss_cumulative.png"), p_cum, width = 9, height = 5, dpi = 300)

geno_tab <- rbindlist(lapply(geno_cand, function(t) {
  k <- vmiss[F_MISS <= t]
  by_chr <- k[, .N, by = CHR]
  data.table(geno = t, n_kept = nrow(k), pct_kept = round(100 * nrow(k) / nrow(vmiss), 1),
             n_kept_maf = k[MAF >= maf_ref, .N], min_chr_n = min(by_chr$N), min_chr = by_chr[which.min(N), CHR])
}))
print(geno_tab)
#    geno  n_kept pct_kept n_kept_maf min_chr_n min_chr
#    <num>   <int>    <num>      <int>     <int>   <int>
# 1:  0.02 2607181     11.6    2401021     21530      36
# 2:  0.05 4901806     21.9    4166515     46848      36
# 3:  0.10 6883442     30.7    5541606     74273      36
# 4:  0.20 8897853     39.7    6748178    106859      36

chr_tab <- dcast(rbindlist(lapply(geno_cand, function(t) vmiss[F_MISS <= t, .N, by = CHR][, geno := t])),
                 CHR ~ geno, value.var = "N", fill = 0)
setnames(chr_tab, -1, paste0("geno_", names(chr_tab)[-1]))
chr_tab <- merge(vmiss[, .(n_total = .N, median_F_MISS = round(median(F_MISS), 3)), by = CHR], chr_tab, by = "CHR")
print(chr_tab)

win_size <- 1e6
all_win <- unique(vmiss[, .(CHR, win = BP %/% win_size)])
dens_tab <- rbindlist(lapply(c(1, geno_cand), function(t) {
  k <- vmiss[F_MISS <= t, .(CHR, BP)]
  kept_win <- k[, .N, by = .(CHR, win = BP %/% win_size)]
  dens <- merge(all_win, kept_win, by = c("CHR", "win"), all.x = TRUE)[is.na(N), N := 0]
  gap <- k[, .(gap = diff(BP)), by = CHR]
  data.table(geno = t, n_win = nrow(dens), win_empty = sum(dens$N == 0), win_lt50 = sum(dens$N < 50),
             median_snp_per_Mb = median(dens$N), max_gap_kb = round(max(gap$gap) / 1e3), max_gap_chr = gap[which.max(gap), CHR])
}))
print(dens_tab)

vmiss[, maf_bin := cut(MAF, c(-Inf, 0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
                       labels = c("0 (monomorphic)", "(0,0.01]", "(0.01,0.05]", "(0.05,0.1]",
                                  "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", "(0.4,0.5]"))]
maf_tab <- vmiss[, c(list(n = .N, median_F_MISS = round(median(F_MISS), 3)),
                     setNames(lapply(geno_cand, function(t) round(100 * mean(F_MISS <= t), 1)), paste0("pct_le_", geno_cand))),
                 by = maf_bin][order(maf_bin)]
print(maf_tab)

smiss <- fread(paste0(in_prefix, ".smiss"), select = c("IID", "F_MISS"))
info <- fread(f_info, select = c("IID", "Breed_Status", "admix_class", "genetic_type", "coi", "top_pct"))
smiss <- merge(smiss, info, by = "IID", all.x = TRUE)
stopifnot(nrow(smiss) == 7166L, !anyNA(smiss$admix_class))
print(summary(smiss$F_MISS))
print(smiss[, .(n = .N, gt_0.05 = sum(F_MISS > 0.05), gt_0.1 = sum(F_MISS > 0.1), gt_0.2 = sum(F_MISS > 0.2))])
smiss_grp <- smiss[, .(n = .N, median_F_MISS = round(median(F_MISS), 4), q90_F_MISS = round(quantile(F_MISS, 0.9), 4)),
                   by = .(admix_class)][order(admix_class)]
print(smiss_grp)
cat(sprintf("개체 결측률과 COI 상관(Spearman): %.3f, top_pct: %.3f\n",
            cor(smiss$F_MISS, smiss$coi, method = "spearman"), cor(smiss$F_MISS, smiss$top_pct, method = "spearman")))
p_imiss <- ggplot(smiss, aes(x = F_MISS, fill = admix_class)) +
  geom_histogram(binwidth = 0.005, boundary = 0, color = "black", linewidth = 0.1) +
  facet_wrap(~ admix_class, ncol = 1, scales = "free_y") +
  labs(title = "Sample missing rate by genetic group (before --geno, reference only)",
       x = "Missing rate (F_MISS)", y = "Number of dogs") +
  theme_qc + theme(legend.position = "none")
ggsave(paste0(pic_prefix, "_smiss_raw_by_group.png"), p_imiss, width = 9, height = 7, dpi = 300)

kin <- fread(f_kin, select = c("IID1", "IID2", "KINSHIP"))
ids <- smiss$IID
kin <- kin[IID1 %in% ids & IID2 %in% ids]
kin_tab <- data.table(relation = c("duplicate/MZ (>= 0.354)", "1st degree [0.177, 0.354)", "2nd degree [0.0884, 0.177)", "3rd degree [0.0442, 0.0884)"),
                      n_pairs = c(kin[KINSHIP >= kin_dup, .N], kin[KINSHIP >= 0.177 & KINSHIP < kin_dup, .N],
                                  kin[KINSHIP >= 0.0884 & KINSHIP < 0.177, .N], kin[KINSHIP >= 0.0442 & KINSHIP < 0.0884, .N]))
print(kin_tab)
dup_pairs <- kin[KINSHIP >= kin_dup]
if (nrow(dup_pairs)) print(dup_pairs[order(-KINSHIP)])

fwrite(geno_tab, paste0(out_prefix, "_geno_threshold_summary.csv"))
fwrite(chr_tab, paste0(out_prefix, "_geno_threshold_by_chr.csv"))
fwrite(dens_tab, paste0(out_prefix, "_geno_threshold_density.csv"))
fwrite(maf_tab, paste0(out_prefix, "_missing_by_maf.csv"))
fwrite(smiss_grp, paste0(out_prefix, "_smiss_raw_by_group.csv"))
fwrite(kin_tab, paste0(out_prefix, "_kinship_summary.csv"))
fwrite(dup_pairs, paste0(out_prefix, "_duplicate_pairs.csv"))

