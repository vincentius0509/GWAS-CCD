# 01-4_mind_threshold.R - Chooses the sample missingness (--mind) threshold. See README.md.
library(data.table)
library(ggplot2)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")    
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_prefix <- file.path(qc_dir, "01-4_groupmiss")
miss_prefix <- paste0(in_prefix, "_missing")
prev_prefix <- file.path(qc_dir, "01-2_geno0.05_missing")
out_prefix <- file.path(qc_dir, "01-4")
pic_prefix <- file.path(pic_dir, "01-4")

mind_cand <- c(0.02, 0.05, 0.1)

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

smiss <- fread(paste0(miss_prefix, ".smiss"), select = c("IID", "MISSING_CT", "OBS_CT", "F_MISS"))
cat("개체 수:", nrow(smiss), " SNP 수:", smiss$OBS_CT[1], "\n")   # 7166, 4857211

vmiss <- fread(paste0(miss_prefix, ".vmiss"), select = "F_MISS")
cat("SNP 결측률 최댓값:", max(vmiss$F_MISS), "(0.05 이하:", max(vmiss$F_MISS) <= 0.05, ")\n")
rm(vmiss); gc()

psam <- fread(paste0(in_prefix, ".psam"))
setnames(psam, 1, "FID")
psam[, Sex_lab := fifelse(SEX == 1, "Male", fifelse(SEX == 2, "Female", NA_character_))]
info <- fread(file.path(qc_dir, "00_all_total.csv"), select = c("IID", "Breed_Status", "admix_class", "top_pop", "coi"))
info[, GRP := fifelse(admix_class == "single_breed", top_pop, admix_class)]

smiss <- merge(smiss, psam[, .(FID, IID, Sex_lab)], by = "IID", all.x = TRUE)
smiss <- merge(smiss, info, by = "IID", all.x = TRUE)
print(colSums(is.na(smiss[, .(Sex_lab, admix_class, coi)])))

prev <- fread(paste0(prev_prefix, ".smiss"), select = c("IID", "F_MISS"))
setnames(prev, "F_MISS", "F_MISS_prev")
smiss <- merge(smiss, prev, by = "IID")
print(smiss[, .(median_prev = median(F_MISS_prev), median_now = median(F_MISS),
                max_prev = max(F_MISS_prev), max_now = max(F_MISS), cor = cor(F_MISS_prev, F_MISS))])

print(summary(smiss$F_MISS))
print(quantile(smiss$F_MISS, c(0.5, 0.9, 0.95, 0.99, 0.999)))

mind_tab <- data.table(mind = mind_cand,
                       n_remove = sapply(mind_cand, function(t) smiss[F_MISS > t, .N]),
                       n_keep = sapply(mind_cand, function(t) smiss[F_MISS <= t, .N]),
                       n_remove_prev = sapply(mind_cand, function(t) smiss[F_MISS_prev > t, .N]))
mind_tab[, pct_remove := round(100 * n_remove / nrow(smiss), 2)]
print(mind_tab)

outlier_ref <- median(smiss$F_MISS) + 3 * IQR(smiss$F_MISS)
cat(sprintf("참고 이상치 경계(중앙값 + 3 IQR): %.4f, 넘는 개체 %d마리\n", outlier_ref, smiss[F_MISS > outlier_ref, .N]))

p_hist <- ggplot(smiss, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.001, boundary = 0, fill = "steelblue", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = mind_cand[mind_cand <= max(smiss$F_MISS) * 1.2], color = "red", linetype = "dashed") +
  labs(title = "Sample missing rate (after all SNP missingness filters)", x = "Missing rate (F_MISS)", y = "Number of dogs") +
  theme_qc
ggsave(paste0(pic_prefix, "_smiss_hist.png"), p_hist, width = 9, height = 5, dpi = 300)

print(smiss[, .(n = .N, median = round(median(F_MISS), 4), q90 = round(quantile(F_MISS, 0.9), 4), max = round(max(F_MISS), 4)), by = Sex_lab])
print(wilcox.test(F_MISS ~ Sex_lab, data = smiss))

print(smiss[, .(n = .N, median = round(median(F_MISS), 4), max = round(max(F_MISS), 4)), by = admix_class])
grp_tab <- smiss[, .(n = .N, median = round(median(F_MISS), 4), max = round(max(F_MISS), 4),
                     n_gt_0.02 = sum(F_MISS > 0.02), n_gt_0.05 = sum(F_MISS > 0.05)), by = GRP][order(-median)]
print(grp_tab[n >= 10])
cat(sprintf("유전 집단이 설명하는 결측률 변이(R², 10마리 이상 집단): %.3f\n",
            summary(lm(F_MISS ~ GRP, data = smiss[GRP %in% grp_tab[n >= 10, GRP]]))$r.squared)) # 0.111
print(cor.test(smiss$F_MISS, smiss$coi, method = "spearman", exact = FALSE))
#        Spearman's rank correlation rho

# data:  smiss$F_MISS and smiss$coi
# S = 8.0488e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.3123632 

p_coi <- ggplot(smiss, aes(x = coi, y = F_MISS, color = admix_class)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(aes(group = 1), method = "loess", formula = y ~ x, se = FALSE, color = "black") +
  labs(title = "Sample missing rate vs COI", x = "COI (release)", y = "Missing rate (F_MISS)", color = NULL) +
  theme_qc + theme(legend.position = "bottom")
ggsave(paste0(pic_prefix, "_smiss_vs_coi.png"), p_coi, width = 7, height = 5.5, dpi = 300)

mind <- 0.05

remove_dogs <- smiss[F_MISS > mind][order(-F_MISS), .(FID, IID, F_MISS, Sex_lab, Breed_Status, admix_class, GRP, coi)]
cat("제외 개체 수 (--mind", mind, "):", nrow(remove_dogs), "\n")
print(remove_dogs)

if (nrow(remove_dogs)) {
  print(remove_dogs[, .N, by = Sex_lab])
  rm_grp <- merge(remove_dogs[, .(n_remove = .N), by = GRP], smiss[, .(n_grp = .N), by = GRP], by = "GRP")
  rm_grp[, pct_remove := round(100 * n_remove / n_grp, 1)]
  print(rm_grp[order(-n_remove)])
}

#    Sex_lab     N
#     <char> <int>
# 1:  Female    27
# 2:    Male    24
#                              GRP n_remove n_grp pct_remove
#                           <char>    <int> <int>      <num>
# 1:                highly_admixed       21  2282        0.9
# 2:                  intermediate       20  1809        1.1
# 3:              golden_retriever        3   388        0.8
# 4:                 cairn_terrier        2    16       12.5
# 5:                  basset_hound        1    21        4.8
# 6:                        beagle        1    57        1.8
# 7: cavalier_king_charles_spaniel        1    50        2.0
# 8:                     dalmatian        1    16        6.2
# 9:          pembroke_welsh_corgi        1    69        1.4

fwrite(mind_tab, paste0(out_prefix, "_mind_threshold_summary.csv"))
fwrite(grp_tab, paste0(out_prefix, "_smiss_by_group.csv"))
fwrite(remove_dogs, paste0(out_prefix, "_mind", mind, "_removed_dogs.csv"))

