# 03-2_hwe_threshold.R - Compares pooled vs. within-group HWE and defines the exclusion rule (p < 1e-6 in >= 2 groups). See README.md.
library(data.table)
library(ggplot2)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")   
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
work_dir <- file.path(qc_dir, "03-1_group_hwe")
base_prefix <- file.path(qc_dir, "02-1_maf0.01")
out_prefix <- file.path(qc_dir, "03-2")
pic_prefix <- file.path(pic_dir, "03-2")

maf_list <- c(0.01, 0.05, 0.1)
hwe_cut <- 1e-6

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

pvar <- fread(paste0(base_prefix, ".pvar"), select = c("#CHROM", "POS", "ID"), col.names = c("CHR", "BP", "SNP"))
n_snp <- nrow(pvar)
cat("SNP 수:", n_snp, "\n")   # 4851571

hwe_all <- fread(file.path(qc_dir, "03-1_hwe_all.hardy"), select = c("ID", "O(HET_A1)", "E(HET_A1)", "P"),
                 col.names = c("SNP", "O_HET", "E_HET", "P"))
hwe_all <- hwe_all[match(pvar$SNP, SNP)]
head(hwe_all)

stopifnot(identical(hwe_all$SNP, pvar$SNP))

glist <- fread(file.path(work_dir, "group_list.tsv"))
cat("집단 수:", nrow(glist), " 개체 수:", sum(glist$N), "\n")   # 45, 2624

n_fail <- n_fail_def <- n_fail_exc <- integer(n_snp)
min_p <- rep(2, n_snp); min_p_grp <- rep(NA_integer_, n_snp)
grp_tab <- copy(glist)
for (j in seq_len(nrow(glist))) {
  h <- fread(cmd = paste("zcat", shQuote(file.path(work_dir, paste0("G.", glist$GRP[j], ".hwe.gz")))),
             header = FALSE, col.names = c("P", "dHET"), na.strings = c("NA", "nan"))
  stopifnot(nrow(h) == n_snp)
  fail <- !is.na(h$P) & h$P < hwe_cut
  n_fail <- n_fail + fail
  n_fail_def <- n_fail_def + (fail & h$dHET < 0)
  n_fail_exc <- n_fail_exc + (fail & h$dHET > 0)
  p <- fifelse(is.na(h$P), 1, h$P)
  up <- p < min_p; min_p[up] <- p[up]; min_p_grp[up] <- j
  grp_tab[j, `:=`(n_fail = sum(fail), pct_deficit = round(100 * mean(h$dHET[fail] < 0), 1))]
}

snp_stat <- data.table(pvar, P_all = hwe_all$P, dHET_all = hwe_all$O_HET - hwe_all$E_HET,
                       n_fail, n_fail_def, n_fail_exc,
                       min_p = signif(min_p, 3), min_p_grp = glist$GRP[min_p_grp])
stopifnot(nrow(snp_stat) == n_snp, !anyNA(snp_stat$min_p_grp))
snp_stat[, `:=`(flag_P = P_all < hwe_cut, flag_W1 = n_fail >= 1, flag_W2 = n_fail >= 2)]

cat("우연히 걸리는 기대 수 (집단 수 x SNP 수 x 1e-6):", round(nrow(glist) * n_snp * hwe_cut), "\n")

crit_tab <- data.table(criterion = c("P: pooled p < 1e-6", "W1: >= 1 group p < 1e-6", "W2: >= 2 groups p < 1e-6"),
                       n_exclude = snp_stat[, c(sum(flag_P), sum(flag_W1), sum(flag_W2))])
crit_tab[, `:=`(pct_exclude = round(100 * n_exclude / n_snp, 2), n_keep = n_snp - n_exclude)]
print(crit_tab)

#                   criterion n_exclude pct_exclude  n_keep
#                      <char>     <int>       <num>   <int>
# 1:       P: pooled p < 1e-6   4806029       99.06   45542
# 2:  W1: >= 1 group p < 1e-6     42653        0.88 4808918
# 3: W2: >= 2 groups p < 1e-6      8835        0.18 4842736

print(snp_stat[flag_P == TRUE, .(n = .N, pct_deficit = round(100 * mean(dHET_all < 0), 1))])
cat(sprintf("전체 검정 F (= 1 - O/E) 중앙값: %.3f\n", snp_stat[, median(1 - (dHET_all + hwe_all$E_HET) / hwe_all$E_HET, na.rm = TRUE)]))
print(snp_stat[flag_P == TRUE, .(n = .N, also_W1 = sum(flag_W1), also_W2 = sum(flag_W2))])

print(snp_stat[n_fail > 0, .(n = .N, deficit_only = sum(n_fail_exc == 0), excess_only = sum(n_fail_def == 0)), keyby = .(n_fail = pmin(n_fail, 10))])

print(grp_tab[order(-n_fail)])

p_bin <- rbind(snp_stat[, .(src = "Pooled (N = 7,115)", bin = pmin(floor(-log10(pmax(P_all, 1e-300))), 30))],
               snp_stat[, .(src = "Min within-group", bin = pmin(floor(-log10(pmax(min_p, 1e-300))), 30))])[, .N, by = .(src, bin)]
p_hwe <- ggplot(p_bin, aes(x = bin + 0.5, y = N / 1e3)) +
  geom_col(width = 1, fill = "steelblue", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = -log10(hwe_cut), color = "red", linetype = "dashed") +
  facet_wrap(~ src, ncol = 1, scales = "free_y") +
  labs(title = "HWE p-value: pooled vs within genetic group", x = "-log10(p)  (30 = 30 or more)", y = "Number of SNPs (thousand)") +
  theme_qc
ggsave(paste0(pic_prefix, "_hwe_pooled_vs_group.png"), p_hwe, width = 9, height = 8, dpi = 300)

chosen <- "W2"
snp_stat[, excl := switch(chosen, P = flag_P, W1 = flag_W1, W2 = flag_W2)]

for (m in maf_list) {
  ids <- fread(file.path(qc_dir, sprintf("02-1_maf%s.pvar", m)), select = "ID")$ID
  stopifnot(all(ids %chin% snp_stat$SNP))
  snp_stat[, (paste0("in_", m)) := SNP %chin% ids]
}

all_win <- unique(snp_stat[, .(CHR, win = BP %/% 1e6)])
maf_compare <- rbindlist(lapply(maf_list, function(m) {
  s <- snp_stat[get(paste0("in_", m)) == TRUE]
  keep <- s[excl == FALSE]
  dens <- merge(all_win, keep[, .N, by = .(CHR, win = BP %/% 1e6)], by = c("CHR", "win"), all.x = TRUE)[is.na(N), N := 0]
  data.table(maf = m, n_snp = nrow(s), n_P = sum(s$flag_P), n_W1 = sum(s$flag_W1), n_W2 = sum(s$flag_W2),
             n_exclude = sum(s$excl), pct_exclude = round(100 * mean(s$excl), 2), n_final = nrow(keep),
             min_chr_n = min(keep[, .N, by = CHR]$N), median_snp_per_Mb = median(dens$N), win_empty = sum(dens$N == 0))
}))
print(maf_compare)

for (m in maf_list) {
  fwrite(snp_stat[get(paste0("in_", m)) == TRUE & excl == TRUE, .(SNP)],
         paste0(out_prefix, "_maf", m, "_hwe_", chosen, "_exclude.txt"), col.names = FALSE)
}
fwrite(crit_tab, paste0(out_prefix, "_hwe_criteria.csv"))
fwrite(grp_tab, paste0(out_prefix, "_hwe_by_group.csv"))
fwrite(maf_compare, paste0(out_prefix, "_hwe_maf_compare.csv"))
fwrite(snp_stat[flag_W1 == TRUE, .(CHR, SNP, BP, P_all = signif(P_all, 3), n_fail, n_fail_def, n_fail_exc, min_p, min_p_grp)],
       paste0(out_prefix, "_hwe_flagged_snps.csv.gz"))

