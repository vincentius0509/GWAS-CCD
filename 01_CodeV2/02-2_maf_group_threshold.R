# 02-2_maf_group_threshold.R - Compares the MAF sets on QC metrics and removes SNPs whose minor alleles are confined to one ancestry group. See README.md.
#   B : top_share > 0.8

library(data.table)
library(ggplot2)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
work_dir <- file.path(qc_dir, "02-1_group_freq")
out_prefix <- file.path(qc_dir, "02-2")
pic_prefix <- file.path(pic_dir, "02-2")

maf_list <- c(0.01, 0.05, 0.1)
poly_cut <- 0.05

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

frq <- fread(file.path(qc_dir, "02-1_geno0.05.acount"), select = c("ID", "ALT_CTS", "OBS_CT"), col.names = c("SNP", "ALT", "OBS"))
miss <- fread(file.path(qc_dir, "01-4_mind0.05_missing.vmiss"), select = c("ID", "F_MISS"), col.names = c("SNP", "F_MISS"))
frq[miss, F_MISS := i.F_MISS, on = "SNP"]
stopifnot(!anyDuplicated(frq$SNP), !anyNA(frq$F_MISS))
frq[, MAF := pmin(ALT, OBS - ALT) / OBS]
rm(miss)

set_ids <- lapply(setNames(maf_list, maf_list), function(m) fread(file.path(qc_dir, sprintf("02-1_maf%s.pvar", m)), select = "ID")$ID)
for (m in maf_list) frq[, (paste0("in_", m)) := SNP %chin% set_ids[[as.character(m)]]]
stopifnot(all(sapply(maf_list, function(m) sum(frq[[paste0("in_", m)]]) == length(set_ids[[as.character(m)]]))))
cat("SNP 수: --geno 재적용 전", nrow(frq), "/ 세트별", sapply(maf_list, function(m) sum(frq[[paste0("in_", m)]])), "\n")

base <- frq[match(set_ids[["0.01"]], SNP)]
cat(sprintf("MAF 0.01 세트 최소 MAF %.4f, 최소 소수 대립유전자 수 %d (= MAF 0.01 x 2 x 7,115 ≈ 142)\n",
            min(base$MAF), base[, min(pmin(ALT, OBS - ALT))]))
maf_tab <- base[, .(n = .N, median_F_MISS = round(median(F_MISS), 4), mean_F_MISS = round(mean(F_MISS), 4)),
                keyby = .(maf_bin = cut(MAF, c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), include.lowest = TRUE))]
maf_tab[, pct := round(100 * n / sum(n), 1)]
print(maf_tab)

glist <- fread(file.path(work_dir, "group_list.tsv"))
cat("집단 수:", nrow(glist), " 개체 수:", sum(glist$N), "\n")   # 67, 6942
n_snp <- nrow(base)
alt_all <- base$ALT; obs_all <- base$OBS
alt_minor <- alt_all <= obs_all - alt_all
minor_all <- pmin(alt_all, obs_all - alt_all)

top_share <- numeric(n_snp); top_grp <- rep(NA_integer_, n_snp); top_ratio <- rep(NA_real_, n_snp); n_poly <- integer(n_snp)
for (j in seq_len(nrow(glist))) {
  cnt <- fread(cmd = paste("zcat", shQuote(file.path(work_dir, paste0("G.", glist$GRP[j], ".counts.gz")))), header = FALSE, col.names = c("ALT", "OBS"))
  stopifnot(nrow(cnt) == n_snp)
  minor_g <- fifelse(alt_minor, cnt$ALT, cnt$OBS - cnt$ALT)
  share <- minor_g / minor_all
  up <- share > top_share
  top_share[up] <- share[up]; top_grp[up] <- j
  top_ratio[up] <- fifelse(cnt$OBS[up] > 0, share[up] / (cnt$OBS[up] / obs_all[up]), NA_real_)
  maf_g <- pmin(cnt$ALT, cnt$OBS - cnt$ALT) / cnt$OBS
  n_poly <- n_poly + (!is.na(maf_g) & maf_g >= poly_cut)
}
snp_stat <- data.table(base[, .(SNP, MAF, in_0.01, in_0.05, in_0.1)],
                       top_share = round(top_share, 4), top_grp = glist$GRP[top_grp], top_ratio = round(top_ratio, 2), n_poly)
snp_stat[, `:=`(CHR = as.integer(sub(":.*", "", SNP)), BP = as.integer(sub("^[^:]+:([0-9]+):.*", "\\1", SNP)))]
snp_stat[, `:=`(flag_A = top_share > 0.9, flag_B = top_share > 0.8, flag_C = n_poly == 0)]
stopifnot(nrow(snp_stat) == n_snp, !anyNA(snp_stat$top_grp), !anyNA(snp_stat$CHR), !anyNA(snp_stat$BP))
if (snp_stat[, any(!is.finite(top_ratio))]) cat("top_ratio 계산 불가(NA) SNP:", snp_stat[!is.finite(top_ratio), .N], "\n")

share_maf_tab <- snp_stat[, .(n = .N, pct_A = round(100 * mean(flag_A), 3), pct_B = round(100 * mean(flag_B), 3), pct_C = round(100 * mean(flag_C), 3),
                              median_top_share = round(median(top_share), 3)),
                          keyby = .(maf_bin = cut(MAF, c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5), include.lowest = TRUE))]
print(share_maf_tab)

chosen <- "A"
snp_stat[, excl := switch(chosen, A = flag_A, B = flag_B, C = flag_C)]

all_win <- unique(snp_stat[, .(CHR, win = BP %/% 1e6)])
set_tab <- rbindlist(lapply(maf_list, function(m) {
  s <- snp_stat[get(paste0("in_", m)) == TRUE]
  k <- s[excl == FALSE]
  dens <- merge(all_win, k[, .N, by = .(CHR, win = BP %/% 1e6)], by = c("CHR", "win"), all.x = TRUE)[is.na(N), N := 0]
  data.table(maf = m, n_snp = nrow(s), n_A = sum(s$flag_A), n_B = sum(s$flag_B), n_C = sum(s$flag_C),
             pct_A = round(100 * mean(s$flag_A), 3), n_keep = nrow(k), min_chr_n = k[, .N, by = CHR][, min(N)],
             win_empty = sum(dens$N == 0), median_snp_per_Mb = median(dens$N),
             median_F_MISS = round(median(frq[get(paste0("in_", m)) == TRUE, F_MISS]), 4))
}))
print(set_tab)

excl_grp <- rbindlist(lapply(maf_list, function(m) snp_stat[get(paste0("in_", m)) == TRUE & excl == TRUE, .N, by = top_grp][, maf := m]))
if (nrow(excl_grp)) print(dcast(excl_grp, top_grp ~ maf, value.var = "N", fill = 0)) else cat("기준", chosen, "에 걸린 SNP 없음 (세 세트 모두)\n")
print(snp_stat[order(-top_share)][1:10, .(SNP, MAF = round(MAF, 4), top_share, top_grp, top_ratio, n_poly)])
print(snp_stat[, .(q99 = quantile(top_share, 0.99), q999 = quantile(top_share, 0.999), max = max(top_share))])

share_bin <- rbindlist(lapply(maf_list, function(m)
  snp_stat[get(paste0("in_", m)) == TRUE, .N, by = .(bin = pmin(floor(top_share / 0.02) * 0.02, 0.98))][, maf := paste("MAF >=", m)]))
p_share <- ggplot(share_bin, aes(x = bin + 0.01, y = N / 1e3)) +
  geom_col(width = 0.02, fill = "steelblue", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = c(0.8, 0.9), color = "red", linetype = "dashed") +
  facet_wrap(~ maf, ncol = 1, scales = "free_y") +
  labs(title = "Top genetic-group share of minor allele copies",
       x = "Top group share of minor allele copies", y = "Number of SNPs (thousand)") +
  theme_qc
ggsave(paste0(pic_prefix, "_top_group_share_hist.png"), p_share, width = 9, height = 9, dpi = 300)

for (m in maf_list)
  fwrite(snp_stat[get(paste0("in_", m)) == TRUE & excl == TRUE, .(SNP)],
         paste0(out_prefix, "_maf", m, "_group_private_", chosen, "_exclude.txt"), col.names = FALSE)
fwrite(maf_tab, paste0(out_prefix, "_maf_bin_missing.csv"))
fwrite(share_maf_tab, paste0(out_prefix, "_group_share_by_maf.csv"))
fwrite(set_tab, paste0(out_prefix, "_maf_set_summary.csv"))
fwrite(snp_stat[flag_B | flag_C, .(CHR, SNP, BP, MAF = round(MAF, 4), top_share, top_grp, top_ratio, n_poly, flag_A, flag_B, flag_C)],
       paste0(out_prefix, "_group_private_flagged_snps.csv.gz"))

