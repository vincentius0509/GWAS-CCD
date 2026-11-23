# 03-5_pca_check.R - Compares the PC methods and writes the GCTA covariate files. See README.md.
library(data.table)
library(ggplot2)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")    
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
work_dir <- file.path(qc_dir, "03-4_pca")
in_prefix <- file.path(qc_dir, "03-3_maf0.05_hwe")
out_prefix <- file.path(qc_dir, "03-5")
pic_prefix <- file.path(pic_dir, "03-5")
f_anc <- file.path(base_dir, "GeneticData/v2023/DogAgingProject_GeneticData_CuratedRelease_2023.ancestry.tsv")

labels <- c("pcair", "all", "qmat")
ref_labels <- c("pcair", "all")
n_pc <- 30
pc_sets <- seq(5, 30, by = 5)
pc_names <- paste0("PC", 1:n_pc)

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

psam <- fread(paste0(in_prefix, ".psam")); setnames(psam, 1, "FID")
pcs <- lapply(setNames(labels, labels), function(l) {
  p <- fread(file.path(work_dir, l, "covar_pc.txt")); setnames(p, 1, "FID")
  p <- p[match(psam$IID, IID)]; stopifnot(identical(p$IID, psam$IID), !anyNA(p[, ..pc_names])); p
})
ref_ids <- lapply(setNames(ref_labels, ref_labels), function(l) fread(file.path(work_dir, l, "ref_select.txt"), header = FALSE)$V2)
eigval <- rbind(rbindlist(lapply(ref_labels, function(l) data.table(label = l, PC = 1:n_pc, eigval = scan(file.path(work_dir, l, "ref_pca.eigenval"), quiet = TRUE)))),
                fread(file.path(work_dir, "qmat", "eigenval.txt"))[PC <= n_pc, .(label = "qmat", PC, eigval = var)])
part <- fread(file.path(work_dir, "pcair_partition.tsv"))
cat("개체:", nrow(psam), " 기준 집단: pcair", length(ref_ids$pcair), "/ all", length(ref_ids$all), "\n")

info <- fread(file.path(qc_dir, "00_all_total.csv"), select = c("IID", "Breed_Status", "admix_class", "top_pop", "top_pct", "coi"))
smiss <- fread(file.path(qc_dir, "01-4_mind0.05_missing.smiss"), select = c("IID", "F_MISS"))
ind <- merge(psam[, .(FID, IID, SEX)], info, by = "IID", all.x = TRUE, sort = FALSE)
ind <- merge(ind, smiss, by = "IID", all.x = TRUE, sort = FALSE)[match(psam$IID, IID)]
ind[, GRP := fifelse(admix_class == "single_breed", top_pop, admix_class)]
stopifnot(identical(ind$IID, psam$IID), !anyNA(ind[, .(admix_class, coi, F_MISS, SEX)]))

anc <- fread(f_anc, select = c("dog", "pop", "pct"))[, `:=`(dog = as.character(dog), pct = as.numeric(pct))]

A <- dcast(anc[dog %in% as.character(psam$IID)], dog ~ pop, value.var = "pct", fill = 0)[match(as.character(psam$IID), dog)]
stopifnot(identical(A$dog, as.character(psam$IID)))
A <- as.matrix(A[, -1])
major_pops <- names(which(colMeans(A) >= 0.01))
cat("평균 비율 1% 이상 참조 집단:", length(major_pops), "개\n")

ind[, pcair_ref := IID %in% ref_ids$pcair]
part_tab <- rbind(ind[, .(n = .N, n_ref = sum(pcair_ref)), by = .(group = admix_class)],
                  ind[admix_class == "single_breed", .(n = .N, n_ref = sum(pcair_ref)), by = .(group = GRP)][order(-n)][1:15])
part_tab[, pct_ref := round(100 * n_ref / n, 1)]
print(part_tab)
print(part[, .(n = .N, median_n_rel = median(n_rel), median_n_div = median(n_div)), by = set])

lam <- fread(file.path(work_dir, "sim_lambda.tsv"))
lam_tab <- lam[, .(lambda_mean = round(mean(lambda), 3), lambda_sd = round(sd(lambda), 3)), by = .(label, n_pc)]
print(lam_tab[label == "none"])
print(dcast(lam_tab[label != "none"], label ~ n_pc, value.var = "lambda_mean"))
p_lam <- ggplot(lam_tab[label != "none"], aes(x = n_pc, y = lambda_mean, color = label)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = lambda_mean - lambda_sd, ymax = lambda_mean + lambda_sd), width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = sprintf("Simulated ancestry-structured phenotype: lambda (no PC = %.2f)", lam_tab[label == "none", lambda_mean]),
       x = "Number of PCs", y = "Genomic inflation (mean of 20 sims)", color = "PC method") +
  theme_qc

ggsave(paste0(pic_prefix, "_pca_lambda.png"), p_lam, width = 8, height = 5, dpi = 300)

r2_tab <- rbindlist(lapply(labels, function(l) rbindlist(lapply(pc_sets, function(k) {
  X <- as.matrix(pcs[[l]][, pc_names[1:k], with = FALSE])
  r2 <- sapply(major_pops, function(p) summary(lm(A[, p] ~ X))$r.squared)
  data.table(label = l, n_pc = k, mean_r2 = round(mean(r2), 3), min_r2 = round(min(r2), 3), min_pop = major_pops[which.min(r2)],
             r2_top_pct = round(summary(lm(ind$top_pct ~ X))$r.squared, 3), r2_coi = round(summary(lm(ind$coi ~ X))$r.squared, 3))
}))))
print(r2_tab)

shrink_tab <- rbindlist(lapply(ref_labels, function(l) {
  d <- data.table(ind[, .(IID, GRP, admix_class)], pcs[[l]][, pc_names[1:10], with = FALSE])
  d[, in_ref := IID %in% ref_ids[[l]]]
  ok <- d[admix_class == "single_breed", .(n_ref = sum(in_ref), n_proj = sum(!in_ref)), by = GRP][n_ref >= 10 & n_proj >= 10, GRP]
  if (!length(ok)) return(data.table(label = l, PC = pc_names[1:10], n_grp = 0L, median_var_ratio = NA_real_, median_shift_sd = NA_real_))
  rbindlist(lapply(pc_names[1:10], function(p) {
    s <- d[GRP %in% ok, .(var_ratio = var(get(p)[!in_ref]) / var(get(p)[in_ref]),
                          shift_sd = (mean(get(p)[!in_ref]) - mean(get(p)[in_ref])) / sd(get(p)[in_ref])), by = GRP]
    data.table(label = l, PC = p, n_grp = nrow(s), median_var_ratio = round(median(s$var_ratio), 3), median_shift_sd = round(median(s$shift_sd), 3))
  }))
}))
print(shrink_tab)
mix_shrink <- data.table(ind[, .(IID, admix_class)], pcs$pcair[, pc_names[1:10], with = FALSE])[admix_class != "single_breed"]
mix_shrink[, in_ref := IID %in% ref_ids$pcair]
cat("pcair 잡종형: 기준", mix_shrink[, sum(in_ref)], "/ 투영", mix_shrink[, sum(!in_ref)], " PC1~10 분산 비(투영/기준):",
    paste(sapply(pc_names[1:10], function(p) round(var(mix_shrink[in_ref == FALSE][[p]]) / var(mix_shrink[in_ref == TRUE][[p]]), 2)), collapse = " "), "\n")

cc_tab <- rbindlist(lapply(list(c("pcair", "all"), c("pcair", "qmat"), c("all", "qmat")), function(pr) rbindlist(lapply(c(5, 10, 20, 30), function(k) {
  cc <- cancor(as.matrix(pcs[[pr[1]]][, pc_names[1:k], with = FALSE]), as.matrix(pcs[[pr[2]]][, pc_names[1:k], with = FALSE]))$cor
  data.table(pair = paste(pr, collapse = " vs "), n_pc = k, cancor_mean = round(mean(cc), 3), cancor_min = round(min(cc), 3))
}))))

print(dcast(cc_tab, pair ~ n_pc, value.var = "cancor_mean"))

tech_tab <- rbindlist(lapply(labels, function(l) rbindlist(lapply(pc_names[1:10], function(p) {
  x <- pcs[[l]][[p]]
  data.table(label = l, PC = p, r_sex = round(cor(x, ind$SEX), 3), r_fmiss = round(cor(x, ind$F_MISS), 3),
             r_coi = round(cor(x, ind$coi), 3), r_top_pct = round(cor(x, ind$top_pct), 3))
}))))
print(tech_tab[, .(max_abs_r_sex = max(abs(r_sex)), max_abs_r_fmiss = max(abs(r_fmiss)), max_abs_r_coi = max(abs(r_coi))), by = label])

eigval[, pct := 100 * eigval / sum(eigval), by = label]
p_scree <- ggplot(eigval, aes(x = PC, y = pct, color = label)) + geom_line() + geom_point(size = 1.5) +
  labs(title = "Scree plot (% of top-30 eigenvalues)", x = "PC", y = "% of variance (top 30)", color = "PC method") + theme_qc
ggsave(paste0(pic_prefix, "_pca_scree.png"), p_scree, width = 8, height = 5, dpi = 300)

top_grps <- ind[admix_class == "single_breed", .N, by = GRP][order(-N)][1:10, GRP]
plot_dt <- rbindlist(lapply(labels, function(l) data.table(label = factor(l, levels = labels), ind[, .(GRP, admix_class)], PC1 = pcs[[l]]$PC1, PC2 = pcs[[l]]$PC2)))
plot_dt[, grp_lab := factor(fcase(GRP %in% top_grps, GRP, admix_class == "single_breed", "other single_breed", default = admix_class),
                            levels = c(top_grps, "other single_breed", "intermediate", "highly_admixed"))]
p_pc12 <- ggplot(plot_dt[order(grp_lab %in% c("intermediate", "highly_admixed"), decreasing = TRUE)], aes(x = PC1, y = PC2, color = grp_lab)) +
  geom_point(size = 0.5, alpha = 0.5) + facet_wrap(~ label, scales = "free") +
  scale_color_manual(values = c(scales::hue_pal()(10), "grey40", "grey70", "grey85")) +
  labs(title = "PC1 vs PC2 by genetic group", color = NULL) + theme_qc + theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 4, override.aes = list(size = 3, alpha = 1)))
ggsave(paste0(pic_prefix, "_pca_PC1_PC2.png"), p_pc12, width = 18, height = 7.5, dpi = 300)

chosen_ref <- "pcair"
n_pc_main <- 10
sens_refs <- setdiff(labels, chosen_ref)

qcov <- fread(file.path(qc_dir, "00_qcovar.txt"), header = FALSE,
              col.names = c("FID", "IID", "age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age"))
cov <- fread(file.path(qc_dir, "00_covar.txt"), header = FALSE, col.names = c("FID", "IID", "Sex", "Intact", "Vision", "Cohort_BHS"))
stopifnot(all(psam$IID %in% qcov$IID), all(psam$IID %in% cov$IID))
qcov <- qcov[match(psam$IID, IID)]; cov <- cov[match(psam$IID, IID)]
for (k in pc_sets) {
  fwrite(cbind(qcov, round(pcs[[chosen_ref]][, pc_names[1:k], with = FALSE], 6)), sprintf("%s_qcovar_PC%d.txt", out_prefix, k), sep = " ", col.names = FALSE)
}
for (l in sens_refs) {
  fwrite(cbind(qcov, round(pcs[[l]][, pc_names[1:n_pc_main], with = FALSE], 6)), sprintf("%s_qcovar_%s_PC%d.txt", out_prefix, l, n_pc_main), sep = " ", col.names = FALSE)
}
fwrite(cov, paste0(out_prefix, "_covar.txt"), sep = " ", col.names = FALSE)
cat("GCTA 공변량 저장:", nrow(psam), "마리 / 주 방법", chosen_ref, "PC", paste(pc_sets, collapse = "/"), "개 (주 분석", n_pc_main, "개) / 민감도",
    paste(sens_refs, collapse = ", "), "PC", n_pc_main, "개\n")

fwrite(part_tab, paste0(out_prefix, "_pcair_partition_by_group.csv"))
fwrite(lam_tab, paste0(out_prefix, "_pca_lambda.csv"))
fwrite(r2_tab, paste0(out_prefix, "_pca_structure_r2.csv"))
fwrite(shrink_tab, paste0(out_prefix, "_pca_shrinkage.csv"))
fwrite(cc_tab, paste0(out_prefix, "_pca_cancor.csv"))
fwrite(tech_tab, paste0(out_prefix, "_pca_tech.csv"))

