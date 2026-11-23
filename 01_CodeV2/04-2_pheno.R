# 04-2_pheno.R - Writes GCTA phenotype files (CSLB rank-based inverse-normal, raw score, simulated phenotypes). See README.md.
#   04-2_pheno.txt                : FID IID cslb_rint cslb_score
#   04-2_sim_pheno.txt            : FID IID SIM1..SIM20

library(data.table)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")    
out_prefix <- file.path(qc_dir, "04-2")
pc_sets <- seq(5, 30, by = 5)

rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))

fam <- fread(file.path(qc_dir, "03-3_maf0.01_hwe.fam"), header = FALSE, select = 1:2, col.names = c("FID", "IID"))
grm_id <- fread(file.path(qc_dir, "04-1_grm", "grm_all.grm.id"), header = FALSE, col.names = c("FID", "IID"))
stopifnot(identical(fam$IID, grm_id$IID))
cat("최종 개체:", nrow(fam), "\n")   # 7115

d <- fread(file.path(qc_dir, "00_all_total.csv"), select = c("IID", "cslb_score", "set_all"))
d <- d[match(fam$IID, IID)]
stopifnot(identical(d$IID, fam$IID), all(d$set_all), !anyNA(d$cslb_score))
d[, cslb_rint := rint(cslb_score)]
cat(sprintf("CSLB 원점수: 평균 %.2f, 중앙값 %.0f, 범위 %d~%d / RINT 평균 %.3f, SD %.3f\n",
            mean(d$cslb_score), median(d$cslb_score), min(d$cslb_score), max(d$cslb_score), mean(d$cslb_rint), sd(d$cslb_rint)))
fwrite(data.table(fam$FID, fam$IID, round(d$cslb_rint, 6), d$cslb_score), paste0(out_prefix, "_pheno.txt"), sep = " ", col.names = FALSE)

sim <- fread(file.path(qc_dir, "03-4_pca", "sim_pheno.txt"))
setnames(sim, 1, "FID")
sim <- sim[match(fam$IID, IID)]
stopifnot(identical(sim$IID, fam$IID), !anyNA(sim))
fwrite(sim, paste0(out_prefix, "_sim_pheno.txt"), sep = " ", col.names = FALSE)
cat("가상 표현형:", ncol(sim) - 2, "개\n")

pc <- fread(file.path(qc_dir, "03-4_pca", "pcair", "covar_pc.txt"))
setnames(pc, 1, "FID")
pc <- pc[match(fam$IID, IID)]
stopifnot(identical(pc$IID, fam$IID))
for (k in pc_sets) fwrite(pc[, c("FID", "IID", paste0("PC", 1:k)), with = FALSE], sprintf("%s_sim_qcovar_PC%d.txt", out_prefix, k), sep = " ", col.names = FALSE)
cat("PC 공변량 파일:", paste(pc_sets, collapse = "/"), "개\n")
