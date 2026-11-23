# 01-3_group_missing_threshold.R - Flags SNPs with missingness concentrated in one ancestry group (> 10% and binomial p < 1e-6). See README.md.
library(data.table)
library(ggplot2)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
work_dir <- file.path(qc_dir, "01-3_group_missing")
in_prefix <- file.path(qc_dir, "01-2_geno0.05")
out_prefix <- file.path(qc_dir, "01-3")
pic_prefix <- file.path(pic_dir, "01-3")

f_cut <- 0.1
p_cut <- 1e-6

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

glist <- fread(file.path(work_dir, "group_list.tsv"))
cat("집단 수:", nrow(glist), " 개체 수:", sum(glist$N), "\n")   # 67, 6993

pvar <- fread(paste0(in_prefix, ".pvar"), select = c("#CHROM", "POS", "ID"), col.names = c("CHR", "BP", "SNP"))
tot <- fread(paste0(in_prefix, "_missing.vmiss"), select = c("ID", "MISSING_CT", "OBS_CT"), col.names = c("SNP", "N_MISS", "N_OBS"))
stopifnot(identical(tot$SNP, pvar$SNP))
N_all <- tot$N_OBS[1]   # 7166

nmiss <- sapply(glist$GRP, function(g) fread(cmd = paste("zcat", shQuote(file.path(work_dir, paste0("G.", g, ".nmiss.gz")))), header = FALSE)$V1)
cat("행렬:", dim(nmiss), "\n")   # 4901806 x 67
stopifnot(nrow(nmiss) == nrow(pvar))

n_snp <- nrow(pvar)
max_f <- rep(-1, n_snp); max_f_grp <- integer(n_snp)
min_p <- rep(2, n_snp);  min_p_grp <- integer(n_snp)
flag_A <- flag_B <- flag_C <- logical(n_snp)
grp_tab <- copy(glist)

for (j in seq_len(nrow(glist))) {
  Ng <- glist$N[j]
  m <- nmiss[, j]
  f <- m / Ng
  p_other <- (tot$N_MISS - m) / (N_all - Ng)
  p <- pbinom(m - 1, Ng, p_other, lower.tail = FALSE)     # P(X >= m)

  a <- f > f_cut; b <- p < p_cut; c <- a & b
  flag_A <- flag_A | a; flag_B <- flag_B | b; flag_C <- flag_C | c

  up <- f > max_f; max_f[up] <- f[up]; max_f_grp[up] <- j
  up <- p < min_p; min_p[up] <- p[up]; min_p_grp[up] <- j

  k <- floor(f_cut * Ng)
  grp_tab[j, `:=`(mean_fmiss = round(mean(f), 4), n_A = sum(a), n_B = sum(b), n_C = sum(c),
                  exp_A_chance = round(sum(pbinom(k, Ng, p_other, lower.tail = FALSE))))]
}

crit_tab <- data.table(criterion = c("A: group F_MISS > 0.1", "B: binomial p < 1e-6", "C: A and B in same group"),
                       n_exclude = c(sum(flag_A), sum(flag_B), sum(flag_C)))
crit_tab[, `:=`(pct_exclude = round(100 * n_exclude / n_snp, 2), n_keep = n_snp - n_exclude)]
print(crit_tab)
#                   criterion n_exclude pct_exclude  n_keep
#                      <char>     <int>       <num>   <int>
# 1:    A: group F_MISS > 0.1   2755736       56.22 2146070
# 2:     B: binomial p < 1e-6     78534        1.60 4823272
# 3: C: A and B in same group     44595        0.91 4857211

print(grp_tab[order(-n_C)][1:20])
cat("잡종 집단 평균 결측률:\n"); print(grp_tab[GRP %in% c("intermediate", "highly_admixed")])
cat(sprintf("참조 품종 집단 평균 결측률 범위: %.4f ~ %.4f\n",
            grp_tab[!GRP %in% c("intermediate", "highly_admixed"), min(mean_fmiss)], grp_tab[!GRP %in% c("intermediate", "highly_admixed"), max(mean_fmiss)]))

max_bin <- data.table(bin = pmin(floor(max_f / 0.02) * 0.02, 0.98))[, .N, by = bin]
p_max <- ggplot(max_bin, aes(x = bin + 0.01, y = N / 1e3)) +
  geom_col(width = 0.02, fill = "steelblue", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = f_cut, color = "red", linetype = "dashed") +
  labs(title = sprintf("Max within-group SNP missing rate (%d genetic groups, N >= 10)", nrow(glist)),
       x = "Max group missing rate per SNP", y = "Number of SNPs (thousand)") +
  theme_qc
ggsave(paste0(pic_prefix, "_max_group_fmiss_hist.png"), p_max, width = 9, height = 5, dpi = 300)

p_bin <- data.table(bin = pmin(floor(-log10(pmax(min_p, 1e-300))), 30))[, .N, by = bin]
p_minp <- ggplot(p_bin, aes(x = bin + 0.5, y = N / 1e3)) +
  geom_col(width = 1, fill = "darkorange", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = -log10(p_cut), color = "red", linetype = "dashed") +
  labs(title = "Min binomial p per SNP (excess missingness in a genetic group)",
       x = "-log10(min p)  (30 = 30 or more)", y = "Number of SNPs (thousand)") +
  theme_qc
ggsave(paste0(pic_prefix, "_min_binom_p_hist.png"), p_minp, width = 9, height = 5, dpi = 300)

chosen <- "C"
flag <- switch(chosen, A = flag_A, B = flag_B, C = flag_C)

snp_stat <- data.table(pvar, max_f = round(max_f, 4), max_f_grp = glist$GRP[max_f_grp],
                       min_p = signif(min_p, 3), min_p_grp = glist$GRP[min_p_grp],
                       flag_A, flag_B, flag_C, excl = flag)
stopifnot(nrow(snp_stat) == n_snp, !anyNA(snp_stat$max_f_grp), !anyNA(snp_stat$min_p_grp))
cat("제외 SNP 수 (기준", chosen, "):", sum(snp_stat$excl), "\n")

print(snp_stat[, .(n_total = .N, n_keep = sum(!excl)), by = CHR][order(n_keep)][1:5])
dens <- merge(unique(snp_stat[, .(CHR, win = BP %/% 1e6)]),
              snp_stat[excl == FALSE, .N, by = .(CHR, win = BP %/% 1e6)], by = c("CHR", "win"), all.x = TRUE)[is.na(N), N := 0]
print(dens[, .(n_win = .N, win_empty = sum(N == 0), median_snp_per_Mb = median(N))])

print(snp_stat[excl == TRUE, .N, by = min_p_grp][order(-N)][1:10])
#                         min_p_grp     N
#                            <char> <int>
#  1:                    leonberger 11282
#  2: cavalier_king_charles_spaniel  3822
#  3:           miniature_schnauzer  2593
#  4:           australian_shepherd  1947
#  5:                           pug  1939
#  6:                  newfoundland  1642
#  7:                      shih_tzu  1398
#  8:          pembroke_welsh_corgi  1373
#  9:                 cairn_terrier  1336
# 10:                         boxer  1288

fwrite(snp_stat[excl == TRUE, .(SNP)], paste0(out_prefix, "_group_miss_exclude.txt"), col.names = FALSE)
fwrite(crit_tab, paste0(out_prefix, "_group_miss_criteria.csv"))
fwrite(grp_tab, paste0(out_prefix, "_group_miss_by_group.csv"))
fwrite(snp_stat[flag_A | flag_B], paste0(out_prefix, "_group_miss_flagged_snps.csv.gz"))

