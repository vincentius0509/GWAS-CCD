# 09-1_size_pgs.R - Body-size polygenic score from published body-size markers (Plassais et al. 2019). See README.md.
library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
out_prefix <- file.path(qc_dir, "09-1"); pic_prefix <- file.path(pic_dir, "09-1")
work_dir <- paste0(out_prefix, "_work")
plink2 <- "/data/00_tools/plink2_linux_avx2_20231123/plink2"
raw_b <- file.path(base_dir, "GeneticData/v2023/DogAgingProject_GeneticData_CuratedRelease_2023_ed")
id_f <- file.path(in_dir, "04-1_grm", "grm_all.grm.id"); qcov_f <- file.path(in_dir, "03-5_qcovar_PC5.txt")
for (f in c(plink2, paste0(raw_b, c(".bed", ".bim", ".fam")), id_f, qcov_f)) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)

mk <- fread(text = "gene,chr,pos,af_modern,beta,af_small,af_medium,af_large
LCORL,3,91269525,0.17,11.15,0.00,0.01,0.68
IGF1,15,41221438,0.44,-7.72,0.85,0.35,0.18
HMGA2,10,8351907,0.21,-8.92,0.87,0.03,0.05
ACSL4,X,82919525,0.08,12.50,0.00,0.03,0.55
ZNF608,11,13945821,0.02,27.33,0.00,0.00,0.16
ADAMTS9-AS,20,26661051,0.02,21.46,0.00,0.00,0.21
R3HDM1,19,38303408,0.01,14.55,0.00,0.00,0.18
SMAD2,7,43782633,0.24,-5.95,0.57,0.16,0.05
GHR_1,4,67040898,0.23,-4.53,0.68,0.19,0.04
IGSF1,X,102292529,0.38,4.63,0.19,0.34,0.96
STC2,4,39182836,0.18,-5.29,0.47,0.10,0.11
HNF4G,29,23802662,0.01,9.97,0.00,0.00,0.20
IGF1R,3,41849479,0.04,-7.70,0.13,0.01,0.00
GHR_2,4,67040939,0.01,-10.17,0.06,0.00,0.00", colClasses = c(chr = "character"))
bim <- fread(paste0(raw_b, ".bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "a1", "a2"), colClasses = list(character = 1))
bim[, chrn := fifelse(chr == "39", "X", chr)]
mk <- merge(mk, bim[, .(chr = chrn, pos, id)], by = c("chr", "pos"), all.x = TRUE)
cat(sprintf("문헌 marker 14개 중 원본 유전형에 있는 것 %d개 (없음: %s)\n", sum(!is.na(mk$id)), paste(mk[is.na(id), gene], collapse = ", ")))
mk <- mk[!is.na(id)]

ids <- fread(id_f, header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
kf <- file.path(work_dir, "keep.txt"); fwrite(ids, kf, sep = " ", col.names = FALSE)
ef <- file.path(work_dir, "extract.txt"); writeLines(mk$id, ef)
rawp <- file.path(work_dir, "size_markers")
if (!file.exists(paste0(rawp, ".raw")))
  system2(plink2, c("--bfile", raw_b, "--dog", "--keep", kf, "--extract", ef, "--export", "A", "--out", rawp), stdout = paste0(rawp, ".stdout"), stderr = paste0(rawp, ".stdout"))
if (!file.exists(paste0(rawp, ".raw"))) stop("marker 추출 실패: ", rawp, ".stdout")
rw <- fread(paste0(rawp, ".raw"), colClasses = list(character = c("FID", "IID")))
qc <- fread(qcov_f, header = FALSE, colClasses = list(character = 1:2)); setnames(qc, c("FID", "IID", "age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5)))
rw <- merge(rw, qc[, .(IID, wt_cov)], by = "IID")
male <- rw$SEX == 1

res <- rbindlist(lapply(seq_len(nrow(mk)), function(i) {
  col <- grep(paste0("^", gsub("([.:])", "\\\\\\1", mk$id[i]), "_"), names(rw), value = TRUE)
  x <- rw[[col]]; counted <- sub(".*_", "", col)
  if (mk$chr[i] == "X") x[male] <- x[male] * ifelse(max(x[male], na.rm = TRUE) <= 1, 2, 1)
  f_small <- mean(x[rw$wt_cov < 10], na.rm = TRUE) / 2; f_large <- mean(x[rw$wt_cov > 41], na.rm = TRUE) / 2
  d_lit <- mk$af_small[i] - mk$af_large[i]; d_our <- f_small - f_large
  flip <- sign(d_our) != sign(d_lit)
  der <- if (flip) 2 - x else x
  keep <- abs(d_our) >= 0.05 && sign(d_lit) != 0
  data.table(gene = mk$gene[i], chr = mk$chr[i], pos = mk$pos[i], id = mk$id[i], counted = counted, flipped = flip, call_rate = mean(!is.na(x)),
             af_ours = mean(der, na.rm = TRUE) / 2, af_lit_modern = mk$af_modern[i],
             af_small_ours = if (flip) 1 - f_small else f_small, af_small_lit = mk$af_small[i],
             af_large_ours = if (flip) 1 - f_large else f_large, af_large_lit = mk$af_large[i], beta_kg = mk$beta[i], used = keep,
             col = col)[, dos := list(list(der))]
}))
D <- data.table(FID = rw$FID, IID = rw$IID)
for (i in seq_len(nrow(res))) { x <- res$dos[[i]]; x[is.na(x)] <- mean(x, na.rm = TRUE); set(D, j = res$gene[i], value = x) }
auto <- res[used == TRUE & chr != "X", gene]; allm <- res[used == TRUE, gene]
D[, PGS_auto := as.vector(as.matrix(.SD) %*% res[match(auto, gene), beta_kg]), .SDcols = auto]
D[, PGS_withX := as.vector(as.matrix(.SD) %*% res[match(allm, gene), beta_kg]), .SDcols = allm]
D <- merge(D, qc[, .(IID, wt_cov)], by = "IID")
fwrite(D, paste0(out_prefix, "_size_pgs.csv"))
fwrite(res[, !c("dos", "col")], paste0(out_prefix, "_markers.csv"))

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("09-1 체격 유전 점수 (Plassais 2019 Table 4 marker · 문헌 kg 가중치; 가설 H5 준비)\n\n■ marker (af = 파생 allele 빈도; ours = 이 자료, lit = 문헌)\n")
print(f3(res[, .(gene, chr, pos, flipped, call_rate, af_ours, af_lit_modern, af_small_ours, af_small_lit, af_large_ours, af_large_lit, beta_kg, used)]))
cat(sprintf("\n■ 점수에 쓴 marker: 상염색체 %d개 (%s) · X 포함 %d개\n", length(auto), paste(auto, collapse = ", "), length(allm)))
cat(sprintf("■ 확인: PGS와 보호자 보고 체중의 상관 — PGS_auto r %.3f (R² %.3f) · PGS_withX r %.3f (R² %.3f)\n",
            cor(D$PGS_auto, D$wt_cov), cor(D$PGS_auto, D$wt_cov)^2, cor(D$PGS_withX, D$wt_cov), cor(D$PGS_withX, D$wt_cov)^2))
cat("  (문헌은 품종 평균 체중 기준 95%; 개별 반려견 · 혼종 · 비만 등이 섞인 보호자 보고 체중이라 그보다 낮은 것이 정상. LCORL 등 3개 없음)\n")
sink()
g1 <- ggplot(D, aes(PGS_auto, wt_cov)) + geom_point(size = 0.4, alpha = 0.25) + geom_smooth(method = "lm", formula = y ~ x, colour = "#b2182b") +
  labs(x = "body-size PGS (autosomal Plassais 2019 markers, kg scale)", y = "owner-reported weight (kg)", title = "External body-size score vs weight") + theme_minimal(base_size = 10)
ggsave(paste0(pic_prefix, "_pgs_weight.png"), g1, width = 6.5, height = 4.5, dpi = 300)
cat(sprintf("저장: %s_size_pgs.csv · _markers.csv · _summary.txt / 그림 %s_pgs_weight.png\n", out_prefix, pic_prefix))

