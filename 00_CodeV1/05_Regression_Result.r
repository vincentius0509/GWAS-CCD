getwd()

library(qqman)
library(data.table)

# 04-1_DuplicateID_Filter.sh 실행 // 한번 실행해두면 다시 실행할 필요는 없음
# 04-2_Regression.sh 실행

# # 함수: λGC 계산
# calculate_lambda_gc <- function(file_path) {
#   if (!file.exists(file_path)) return(NA)
#   gwas <- fread(file_path)
#   gwas <- gwas[!is.na(P) & P > 0]
#   if (nrow(gwas) == 0) return(NA)
#   gwas[, chisq := qchisq(1 - P, df = 1)]
#   lambda_gc <- median(gwas$chisq) / qchisq(0.5, df = 1)
#   return(lambda_gc)
# }

calculate_lambda_gc <- function(file_path) {
  cat("Calculating:", file_path, "\n")
  if (!file.exists(file_path)) return(NA)

  gwas <- tryCatch(fread(file_path), error = function(e) {
    cat("❌ Error reading", file_path, "\n")
    return(NULL)
  })
  if (is.null(gwas) || !"P" %in% colnames(gwas)) {
    cat("⚠️ No valid 'P' column in", file_path, "\n")
    return(NA)
  }

  # 숫자로 변환 후 필터링
  suppressWarnings(gwas[, P := as.numeric(P)])
  gwas <- gwas[is.finite(P) & P > 0 & P < 1]
  if (nrow(gwas) == 0) {
    cat("⚠️ No valid P-values in", file_path, "\n")
    return(NA)
  }

  gwas[, chisq := qchisq(1 - P, df = 1)]
  lambda_gc <- median(gwas$chisq) / qchisq(0.5, df = 1)
  return(lambda_gc)
}

calculate_lambda_gc("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_5_linear.phe1.glm.linear")
# calculate_lambda_gc("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE2_001_1e6_5_logistic.phe2.glm.logistic.hybrid")

# 입력 파일 목록 생성
maf_vals <- c("001", "005", "01")
hwe_vals <- c("1e6")
pc_vals <- c(5, 10, 15, 20, 25, 30)

# 결과 저장용 테이블
results <- data.table()
base_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure"

# 루프: 모든 조합 반복
for (maf in maf_vals) {
  for (hwe in hwe_vals) {
    for (pc in pc_vals) {
      linear_path <- file.path(base_dir,sprintf("05_GWAS_PHE1_%s_%s_%d_linear.phe1.glm.linear", maf, hwe, pc))
      # logistic_path <- file.path(base_dir,sprintf("05_GWAS_PHE2_%s_%s_%d_logistic.phe2.glm.logistic.hybrid", maf, hwe, pc))

      lambda_linear <- calculate_lambda_gc(linear_path)
      # lambda_logistic <- calculate_lambda_gc(logistic_path)

      results <- rbind(
        results,
        data.table(
          MAF = maf,
          HWE = hwe,
          PC = pc,
          Lambda_Linear = lambda_linear
          # Lambda_Logistic = lambda_logistic
        ),
        fill = TRUE
      )
    }
  }
}

# 결과 저장
fwrite(results, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_Genomic_Inflation_Factor.csv")

head(results)

# Lambda_Linear, Lambda_Logistic 모두 0.95 이상 1.05 미만인 행 찾기
results[Lambda_Linear >= 0.95 & Lambda_Linear < 1.05 ]
# 모든 경우의 수에 대해 Lambda_Linear가 0.95 이상 1.05 미만

#       MAF    HWE    PC Lambda_Linear Lambda_Logistic
#    <char> <char> <num>         <num>           <num>
# 1:    001    1e4    15      1.011530       0.9597764
# 2:    001    1e6    15      1.010754       0.9644586

plink_001_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_5_linear.phe1.glm.linear")
plink_001_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear")
plink_001_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_15_linear.phe1.glm.linear")
plink_001_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_20_linear.phe1.glm.linear")
plink_001_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_25_linear.phe1.glm.linear")
plink_001_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_30_linear.phe1.glm.linear")

plink_005_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_5_linear.phe1.glm.linear")
plink_005_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_10_linear.phe1.glm.linear")
plink_005_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_15_linear.phe1.glm.linear")
plink_005_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_20_linear.phe1.glm.linear")
plink_005_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_25_linear.phe1.glm.linear")
plink_005_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_30_linear.phe1.glm.linear")

plink_01_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_5_linear.phe1.glm.linear")
plink_01_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_10_linear.phe1.glm.linear")
plink_01_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_15_linear.phe1.glm.linear") 
plink_01_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_20_linear.phe1.glm.linear")
plink_01_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_25_linear.phe1.glm.linear")
plink_01_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_30_linear.phe1.glm.linear")

genomewide_threshold <- 5e-8
suggestive_threshold <- 1e-5
low_threshold <- 1e-4

# λGC 값 수동 설정
head(results)
round(results$Lambda_Linear[1],7)
round(results$Lambda_Linear[2],7)
round(results$Lambda_Linear[3],7)
round(results$Lambda_Linear[18],7)

lambda_001_1e6_5 <- round(results$Lambda_Linear[1], 7)
lambda_001_1e6_10 <- round(results$Lambda_Linear[2], 7)
lambda_001_1e6_15 <- round(results$Lambda_Linear[3], 7)
lambda_001_1e6_20 <- round(results$Lambda_Linear[4], 7)
lambda_001_1e6_25 <- round(results$Lambda_Linear[5], 7)
lambda_001_1e6_30 <- round(results$Lambda_Linear[6], 7)
lambda_005_1e6_5 <- round(results$Lambda_Linear[7], 7)
lambda_005_1e6_10 <- round(results$Lambda_Linear[8], 7)
lambda_005_1e6_15 <- round(results$Lambda_Linear[9], 7)
lambda_005_1e6_20 <- round(results$Lambda_Linear[10], 7)
lambda_005_1e6_25 <- round(results$Lambda_Linear[11], 7)
lambda_005_1e6_30 <- round(results$Lambda_Linear[12], 7)
lambda_01_1e6_5 <- round(results$Lambda_Linear[13], 7)
lambda_01_1e6_10 <- round(results$Lambda_Linear[14], 7)
lambda_01_1e6_15 <- round(results$Lambda_Linear[15], 7)
lambda_01_1e6_20 <- round(results$Lambda_Linear[16], 7)
lambda_01_1e6_25 <- round(results$Lambda_Linear[17], 7)
lambda_01_1e6_30 <- round(results$Lambda_Linear[18], 7)

head(plink_001_1e6_5)
colnames(plink_001_1e6_5)[colnames(plink_001_1e6_5) == "#CHROM"] <- "CHR"
plink_001_1e6_5$CHR <- as.numeric(gsub("X", "39", plink_001_1e6_5$CHR))
table(plink_001_1e6_5$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_5_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_001_1e6_5, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=5]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_5_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_001_1e6_5$P, main = sprintf("QQ Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=5]\nλGC = %.5f", lambda_001_1e6_5), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_001_1e6_10)
colnames(plink_001_1e6_10)[colnames(plink_001_1e6_10) == "#CHROM"] <- "CHR"
plink_001_1e6_10$CHR <- as.numeric(gsub("X", "39", plink_001_1e6_10$CHR))
table(plink_001_1e6_10$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_10_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_001_1e6_10, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=10]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_10_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_001_1e6_10$P, main = sprintf("QQ Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=10]\nλGC = %.5f", lambda_001_1e6_10), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_001_1e6_15)
colnames(plink_001_1e6_15)[colnames(plink_001_1e6_15) == "#CHROM"] <- "CHR"
plink_001_1e6_15$CHR <- as.numeric(gsub("X", "39", plink_001_1e6_15$CHR))
table(plink_001_1e6_15$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_15_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_001_1e6_15, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=15]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)") 
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_15_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_001_1e6_15$P, main = sprintf("QQ Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=15]\nλGC = %.5f", lambda_001_1e6_15), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_001_1e6_20)
colnames(plink_001_1e6_20)[colnames(plink_001_1e6_20) == "#CHROM"] <- "CHR"
plink_001_1e6_20$CHR <- as.numeric(gsub("X", "39", plink_001_1e6_20$CHR))
table(plink_001_1e6_20$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_20_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_001_1e6_20, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=20]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_20_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_001_1e6_20$P, main = sprintf("QQ Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=20]\nλGC = %.5f", lambda_001_1e6_20), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_001_1e6_25)
colnames(plink_001_1e6_25)[colnames(plink_001_1e6_25) == "#CHROM"] <- "CHR"
plink_001_1e6_25$CHR <- as.numeric(gsub("X", "39", plink_001_1e6_25$CHR))
table(plink_001_1e6_25$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_25_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_001_1e6_25, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=25]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_25_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_001_1e6_25$P, main = sprintf("QQ Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=25]\nλGC = %.5f", lambda_001_1e6_25), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_001_1e6_30)
colnames(plink_001_1e6_30)[colnames(plink_001_1e6_30) == "#CHROM"] <- "CHR"
plink_001_1e6_30$CHR <- as.numeric(gsub("X", "39", plink_001_1e6_30$CHR))
table(plink_001_1e6_30$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_30_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_001_1e6_30, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=30]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_30_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_001_1e6_30$P, main = sprintf("QQ Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=30]\nλGC = %.5f", lambda_001_1e6_30), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_005_1e6_5)
colnames(plink_005_1e6_5)[colnames(plink_005_1e6_5) == "#CHROM"] <- "CHR"
plink_005_1e6_5$CHR <- as.numeric(gsub("X", "39", plink_005_1e6_5$CHR))
table(plink_005_1e6_5$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_5_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_005_1e6_5, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=5]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_5_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_005_1e6_5$P, main = sprintf("QQ Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=5]\nλGC = %.5f", lambda_005_1e6_5), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_005_1e6_10)
colnames(plink_005_1e6_10)[colnames(plink_005_1e6_10) == "#CHROM"] <- "CHR"
plink_005_1e6_10$CHR <- as.numeric(gsub("X", "39", plink_005_1e6_10$CHR))
table(plink_005_1e6_10$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_10_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_005_1e6_10, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=10]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_10_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_005_1e6_10$P, main = sprintf("QQ Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=10]\nλGC = %.5f", lambda_005_1e6_10), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_005_1e6_15)
colnames(plink_005_1e6_15)[colnames(plink_005_1e6_15) == "#CHROM"] <- "CHR"
plink_005_1e6_15$CHR <- as.numeric(gsub("X", "39", plink_005_1e6_15$CHR))
table(plink_005_1e6_15$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_15_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_005_1e6_15, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=15]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_15_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_005_1e6_15$P, main = sprintf("QQ Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=15]\nλGC = %.5f", lambda_005_1e6_15), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_005_1e6_20)
colnames(plink_005_1e6_20)[colnames(plink_005_1e6_20) == "#CHROM"] <- "CHR"
plink_005_1e6_20$CHR <- as.numeric(gsub("X", "39", plink_005_1e6_20$CHR))
table(plink_005_1e6_20$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_20_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_005_1e6_20, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=20]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_20_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_005_1e6_20$P, main = sprintf("QQ Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=20]\nλGC = %.5f", lambda_005_1e6_20), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_005_1e6_25)
colnames(plink_005_1e6_25)[colnames(plink_005_1e6_25) == "#CHROM"] <- "CHR"
plink_005_1e6_25$CHR <- as.numeric(gsub("X", "39", plink_005_1e6_25$CHR))
table(plink_005_1e6_25$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_25_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_005_1e6_25, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=25]", 
      col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_25_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_005_1e6_25$P, main = sprintf("QQ Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=25]\nλGC = %.5f", lambda_005_1e6_25), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_005_1e6_30)
colnames(plink_005_1e6_30)[colnames(plink_005_1e6_30) == "#CHROM"] <- "CHR"
plink_005_1e6_30$CHR <- as.numeric(gsub("X", "39", plink_005_1e6_30$CHR))
table(plink_005_1e6_30$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_30_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_005_1e6_30, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=30]", 
      col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_005_1e6_30_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_005_1e6_30$P, main = sprintf("QQ Plot for PHE1 [MAF=0.05, HWE=1e-6, PC=30]\nλGC = %.5f", lambda_005_1e6_30), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_01_1e6_5)
colnames(plink_01_1e6_5)[colnames(plink_01_1e6_5) == "#CHROM"] <- "CHR"
plink_01_1e6_5$CHR <- as.numeric(gsub("X", "39", plink_01_1e6_5$CHR))
table(plink_01_1e6_5$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_5_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_01_1e6_5, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=5]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_5_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_01_1e6_5$P, main = sprintf("QQ Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=5]\nλGC = %.5f", lambda_01_1e6_5), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_01_1e6_10)
colnames(plink_01_1e6_10)[colnames(plink_01_1e6_10) == "#CHROM"] <- "CHR"
plink_01_1e6_10$CHR <- as.numeric(gsub("X", "39", plink_01_1e6_10$CHR))
table(plink_01_1e6_10$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_10_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_01_1e6_10, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=10]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_10_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_01_1e6_10$P, main = sprintf("QQ Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=10]\nλGC = %.5f", lambda_01_1e6_10), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_01_1e6_15)
colnames(plink_01_1e6_15)[colnames(plink_01_1e6_15) == "#CHROM"] <- "CHR"
plink_01_1e6_15$CHR <- as.numeric(gsub("X", "39", plink_01_1e6_15$CHR))
table(plink_01_1e6_15$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_15_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_01_1e6_15, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=15]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_15_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_01_1e6_15$P, main = sprintf("QQ Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=15]\nλGC = %.5f", lambda_01_1e6_15), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_01_1e6_20)
colnames(plink_01_1e6_20)[colnames(plink_01_1e6_20) == "#CHROM"] <- "CHR"
plink_01_1e6_20$CHR <- as.numeric(gsub("X", "39", plink_01_1e6_20$CHR))
table(plink_01_1e6_20$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_20_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_01_1e6_20, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=20]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_20_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_01_1e6_20$P, main = sprintf("QQ Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=20]\nλGC = %.5f", lambda_01_1e6_20), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_01_1e6_25)
colnames(plink_01_1e6_25)[colnames(plink_01_1e6_25) == "#CHROM"] <- "CHR"
plink_01_1e6_25$CHR <- as.numeric(gsub("X", "39", plink_01_1e6_25$CHR))
table(plink_01_1e6_25$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_25_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_01_1e6_25, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=25]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_25_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_01_1e6_25$P, main = sprintf("QQ Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=25]\nλGC = %.5f", lambda_01_1e6_25), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

head(plink_01_1e6_30)
colnames(plink_01_1e6_30)[colnames(plink_01_1e6_30) == "#CHROM"] <- "CHR"
plink_01_1e6_30$CHR <- as.numeric(gsub("X", "39", plink_01_1e6_30$CHR))
table(plink_01_1e6_30$CHR)
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_30_linear.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
manhattan(plink_01_1e6_30, chr = "CHR", bp = "POS", p = "P", snp = "ID", main = "Manhattan Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=30]", 
          col = c("blue", "orange"), genomewideline = -log10(genomewide_threshold), suggestiveline = -log10(suggestive_threshold), ylim = c(0, 15), xlab = "Chromosome", ylab = "-log10(p-value)")
dev.off()
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_01_1e6_30_linear_QQ.png", width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)
qq(plink_01_1e6_30$P, main = sprintf("QQ Plot for PHE1 [MAF=0.1, HWE=1e-6, PC=30]\nλGC = %.5f", lambda_01_1e6_30), xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", col = "blue", cex = 0.5)
dev.off()

genomewide_threshold
head(plink_001_1e6_5)

cat("Genome-wide significance threshold: ", nrow(plink_001_1e6_5[plink_001_1e6_5$P < genomewide_threshold,]), "\n") # 4
cat("Genome-wide significant threshold: ", nrow(plink_001_1e6_10[plink_001_1e6_10$P < genomewide_threshold,]), "\n") # 6
cat("Genome-wide significant threshold: ", nrow(plink_001_1e6_15[plink_001_1e6_15$P < genomewide_threshold,]), "\n") # 6
cat("Genome-wide significant threshold: ", nrow(plink_001_1e6_20[plink_001_1e6_20$P < genomewide_threshold,]), "\n") # 4
cat("Genome-wide significant threshold: ", nrow(plink_001_1e6_25[plink_001_1e6_25$P < genomewide_threshold,]), "\n") # 3
cat("Genome-wide significant threshold: ", nrow(plink_001_1e6_30[plink_001_1e6_30$P < genomewide_threshold,]), "\n") # 3

cat("Genome-wide significant threshold: ", nrow(plink_005_1e6_5[plink_005_1e6_5$P < genomewide_threshold,]), "\n") # 0
cat("Genome-wide significant threshold: ", nrow(plink_005_1e6_10[plink_005_1e6_10$P < genomewide_threshold,]), "\n") # 0
cat("Genome-wide significant threshold: ", nrow(plink_005_1e6_15[plink_005_1e6_15$P < genomewide_threshold,]), "\n") # 0
cat("Genome-wide significant threshold: ", nrow(plink_005_1e6_20[plink_005_1e6_20$P < genomewide_threshold,]), "\n") # 0
cat("Genome-wide significant threshold: ", nrow(plink_005_1e6_25[plink_005_1e6_25$P < genomewide_threshold,]), "\n") # 0
cat("Genome-wide significant threshold: ", nrow(plink_005_1e6_30[plink_005_1e6_30$P < genomewide_threshold,]), "\n") # 0

cat("Genome-wide significant threshold: ", nrow(plink_01_1e6_5[plink_01_1e6_5$P < genomewide_threshold,]), "\n") # 0
cat("Genome-wide significant threshold: ", nrow(plink_01_1e6_10[plink_01_1e6_10$P < genomewide_threshold,]), "\n") # 0
cat("Genome-wide significant threshold: ", nrow(plink_01_1e6_15[plink_01_1e6_15$P < genomewide_threshold,]), "\n") # 0
cat("Genome-wide significant threshold: ", nrow(plink_01_1e6_20[plink_01_1e6_20$P < genomewide_threshold,]), "\n") # 0
cat("Genome-wide significant threshold: ", nrow(plink_01_1e6_25[plink_01_1e6_25$P < genomewide_threshold,]), "\n") # 0  
cat("Genome-wide significant threshold: ", nrow(plink_01_1e6_30[plink_01_1e6_30$P < genomewide_threshold,]), "\n") # 0

suggestive_threshold

cat("Suggestive significance threshold: ", nrow(plink_001_1e6_5[plink_001_1e6_5$P < suggestive_threshold,]), "\n") # 88
cat("Suggestive significance threshold: ", nrow(plink_001_1e6_10[plink_001_1e6_10$P < suggestive_threshold,]), "\n") # 89
cat("Suggestive significance threshold: ", nrow(plink_001_1e6_15[plink_001_1e6_15$P < suggestive_threshold,]), "\n") # 62
cat("Suggestive significance threshold: ", nrow(plink_001_1e6_20[plink_001_1e6_20$P < suggestive_threshold,]), "\n") # 49
cat("Suggestive significance threshold: ", nrow(plink_001_1e6_25[plink_001_1e6_25$P < suggestive_threshold,]), "\n") # 47
cat("Suggestive significance threshold: ", nrow(plink_001_1e6_30[plink_001_1e6_30$P < suggestive_threshold,]), "\n") # 25

cat("Suggestive significance threshold: ", nrow(plink_005_1e6_5[plink_005_1e6_5$P < suggestive_threshold,]), "\n") # 2
cat("Suggestive significance threshold: ", nrow(plink_005_1e6_10[plink_005_1e6_10$P < suggestive_threshold,]), "\n") # 2
cat("Suggestive significance threshold: ", nrow(plink_005_1e6_15[plink_005_1e6_15$P < suggestive_threshold,]), "\n") # 2
cat("Suggestive significance threshold: ", nrow(plink_005_1e6_20[plink_005_1e6_20$P < suggestive_threshold,]), "\n") # 2
cat("Suggestive significance threshold: ", nrow(plink_005_1e6_25[plink_005_1e6_25$P < suggestive_threshold,]), "\n") # 1
cat("Suggestive significance threshold: ", nrow(plink_005_1e6_30[plink_005_1e6_30$P < suggestive_threshold,]), "\n") # 1

cat("Suggestive significance threshold: ", nrow(plink_01_1e6_5[plink_01_1e6_5$P < suggestive_threshold,]), "\n") # 0
cat("Suggestive significance threshold: ", nrow(plink_01_1e6_10[plink_01_1e6_10$P < suggestive_threshold,]), "\n") # 0
cat("Suggestive significance threshold: ", nrow(plink_01_1e6_15[plink_01_1e6_15$P < suggestive_threshold,]), "\n") # 0
cat("Suggestive significance threshold: ", nrow(plink_01_1e6_20[plink_01_1e6_20$P < suggestive_threshold,]), "\n") # 0
cat("Suggestive significance threshold: ", nrow(plink_01_1e6_25[plink_01_1e6_25$P < suggestive_threshold,]), "\n") # 0
cat("Suggestive significance threshold: ", nrow(plink_01_1e6_30[plink_01_1e6_30$P < suggestive_threshold,]), "\n") # 0

low_threshold
cat("Low significance threshold: ", nrow(plink_001_1e6_5[plink_001_1e6_5$P < low_threshold,]), "\n") # 622
cat("Low significance threshold: ", nrow(plink_001_1e6_10[plink_001_1e6_10$P < low_threshold,]), "\n") # 563
cat("Low significance threshold: ", nrow(plink_001_1e6_15[plink_001_1e6_15$P < low_threshold,]), "\n") # 423
cat("Low significance threshold: ", nrow(plink_001_1e6_20[plink_001_1e6_20$P < low_threshold,]), "\n") # 348
cat("Low significance threshold: ", nrow(plink_001_1e6_25[plink_001_1e6_25$P < low_threshold,]), "\n") # 294
cat("Low significance threshold: ", nrow(plink_001_1e6_30[plink_001_1e6_30$P < low_threshold,]), "\n") # 175

cat("Low significance threshold: ", nrow(plink_005_1e6_5[plink_005_1e6_5$P < low_threshold,]), "\n") # 14
cat("Low significance threshold: ", nrow(plink_005_1e6_10[plink_005_1e6_10$P < low_threshold,]), "\n") # 12
cat("Low significance threshold: ", nrow(plink_005_1e6_15[plink_005_1e6_15$P < low_threshold,]), "\n") # 12
cat("Low significance threshold: ", nrow(plink_005_1e6_20[plink_005_1e6_20$P < low_threshold,]), "\n") # 12
cat("Low significance threshold: ", nrow(plink_005_1e6_25[plink_005_1e6_25$P < low_threshold,]), "\n") # 11
cat("Low significance threshold: ", nrow(plink_005_1e6_30[plink_005_1e6_30$P < low_threshold,]), "\n") # 13

cat("Low significance threshold: ", nrow(plink_01_1e6_5[plink_01_1e6_5$P < low_threshold,]), "\n") # 0
cat("Low significance threshold: ", nrow(plink_01_1e6_10[plink_01_1e6_10$P < low_threshold,]), "\n") # 0
cat("Low significance threshold: ", nrow(plink_01_1e6_15[plink_01_1e6_15$P < low_threshold,]), "\n") # 0 
cat("Low significance threshold: ", nrow(plink_01_1e6_20[plink_01_1e6_20$P < low_threshold,]), "\n") # 0
cat("Low significance threshold: ", nrow(plink_01_1e6_25[plink_01_1e6_25$P < low_threshold,]), "\n") # 1
cat("Low significance threshold: ", nrow(plink_01_1e6_30[plink_01_1e6_30$P < low_threshold,]), "\n") # 0


# MAF: 0.01, HWE: 1e-6, PC: 10
# 만족하는 SNP의 수가 더 많은 plink_001_1e6_10
suggestive_threshold # 6개
genomewide_threshold # 89개

head(plink_001_1e6_10)

colnames(plink_001_1e6_10)[colnames(plink_001_1e6_10) == "ID"] <- "SNP"
colnames(plink_001_1e6_10)

linear_genomewide_threshold <- plink_001_1e6_10[plink_001_1e6_10$P < genomewide_threshold, .(CHR, SNP, POS, A1, P, BETA)]
head(linear_genomewide_threshold)
dim(linear_genomewide_threshold) # 6  6

linear_suggestive_threshold <- plink_001_1e6_10[plink_001_1e6_10$P < suggestive_threshold, .(CHR, SNP, POS, A1, P, BETA)]
head(linear_suggestive_threshold)
dim(linear_suggestive_threshold) # 89 6

write.csv(linear_genomewide_threshold, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_Regression_Result_001_1e6_10_GW.csv", row.names = FALSE)
write.csv(linear_suggestive_threshold, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_Regression_Result_001_1e6_10_SG.csv", row.names = FALSE)

# low_threshold
low_threshold <- 1e-4
linear_low_threshold <- plink_001_1e6_10[plink_001_1e6_10$P < low_threshold, .(CHR, SNP, POS, A1, P, BETA)]
head(linear_low_threshold)
dim(linear_low_threshold) # 563 6

# Manhattan Plot
png("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/05_GWAS_PHE1_001_1e6_10_linear_new.png",
    width = 12, height = 8, units = "in", res = 300, bg = "white", pointsize = 12)

manhattan(plink_001_1e6_10,
          chr = "CHR", bp = "POS", p = "P", snp = "SNP",
          main = "Manhattan Plot for PHE1 [MAF=0.01, HWE=1e-6, PC=10]",
          col = c("blue", "orange"),
          genomewideline = -log10(genomewide_threshold),
          suggestiveline = -log10(suggestive_threshold),
          ylim = c(0, 15),
          xlab = "Chromosome", ylab = "-log10(p-value)")

abline(h = -log10(low_threshold), col = "darkgreen", lty = 2, lwd = 2)

dev.off()

write.csv(linear_low_threshold, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_Regression_Result_001_1e6_10_LT.csv", row.names = FALSE)


############################################################################################################################################################
# 논문용 figure
############################################################################################################################################################

# 패키지
library(qqman)        # manhattan(), qq()
# library(svglite)      # ggsave 없이도 svg()
library(grDevices)




pal <- c("#2A4F9E", "#FF9D00")     # 논문용 색
thr_gw <- 5e-8                     # genome-wide
thr_sg <- 1e-5                     # suggestive
thr_ex <- 1e-4                     # exploratory (추가)


prep_gwas_df <- function(df, chr_col="#CHROM", bp_col="POS", snp_col="ID", p_col="P") {
  # 표준 컬럼명으로 통일
  if (chr_col != "CHR" && "#CHROM" %in% names(df)) {
    names(df)[names(df) == "#CHROM"] <- "CHR"
  } else if (!"CHR" %in% names(df)) {
    names(df)[names(df) == chr_col] <- "CHR"
  }
  names(df)[names(df) == bp_col]  <- "BP"
  names(df)[names(df) == snp_col] <- "SNP"
  names(df)[names(df) == p_col]   <- "P"

  # X 염색체 → 39 (canfam 관례)
  df$CHR <- gsub("^X$", "39", df$CHR)
  df$CHR <- suppressWarnings(as.numeric(df$CHR))

  # p-value/좌표 유효값만
  df <- df[is.finite(df$P) & !is.na(df$P) & df$P > 0 & df$P <= 1, ]
  df <- df[is.finite(df$CHR) & !is.na(df$CHR) & is.finite(df$BP) & !is.na(df$BP), ]
  df
}

lambda_gc <- function(p){
  chisq <- qchisq(1-p, df=1)
  median(chisq,na.rm=TRUE)/0.456
}

plink_001_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear")

plink_001_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_5_linear.phe1.glm.linear")
plink_001_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_15_linear.phe1.glm.linear")
plink_001_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_20_linear.phe1.glm.linear")
plink_001_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_25_linear.phe1.glm.linear")
plink_001_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_30_linear.phe1.glm.linear")

plink_005_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_5_linear.phe1.glm.linear")
plink_005_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_10_linear.phe1.glm.linear")
plink_005_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_15_linear.phe1.glm.linear")
plink_005_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_20_linear.phe1.glm.linear")
plink_005_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_25_linear.phe1.glm.linear")
plink_005_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_005_1e6_30_linear.phe1.glm.linear")

plink_01_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_5_linear.phe1.glm.linear")
plink_01_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_10_linear.phe1.glm.linear")
plink_01_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_15_linear.phe1.glm.linear") 
plink_01_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_20_linear.phe1.glm.linear")
plink_01_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_25_linear.phe1.glm.linear")
plink_01_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_01_1e6_30_linear.phe1.glm.linear")

plot_manhattan_qq_base <- function(df,
                                   title_main,
                                   out_prefix,
                                   ymax = NULL,
                                   lambda_override = NULL,   # <- PPT 값 직접 입력
                                   write_pdf = TRUE,         # PDF 저장 여부 (문제 있으면 FALSE)
                                   annotate_top = 0) {

  df <- prep_gwas_df(df)

  if (is.null(ymax)) {
    ymax <- ceiling(max(-log10(df$P), na.rm = TRUE) * 1.05)
    ymax <- max(ymax, 10)
  }

  ## ---------- Manhattan (PNG) ----------
  png(paste0(out_prefix, "_manhattan.png"),
      width = 3000, height = 1800, res = 300, bg = "white", pointsize = 12)
  manhattan(df,
            chr="CHR", bp="BP", p="P", snp="SNP",
            col=c("#2A4F9E", "#FF9D00"),    # blue / orange
            genomewideline = -log10(thr_gw),
            suggestiveline = -log10(thr_sg),
            ylim = c(0, ymax),
            cex = 0.45, cex.axis = 0.9,
            main = title_main,
            ylab = expression(-log[10](italic(p))),
            xlab = "Chromosome")
  abline(h = -log10(thr_ex), lty = 2, col = "#555555")  # 1e-4 라인
  if (annotate_top > 0) {
    ord <- order(df$P)
    top <- df[ord[seq_len(min(annotate_top, nrow(df)))], ]
    with(top, text(x = CHR, y = -log10(P), labels = SNP, cex = 0.6, pos = 3))
  }
  dev.off()

  ## (선택) PDF도 저장
  if (isTRUE(write_pdf)) {
    pdf(paste0(out_prefix, "_manhattan.pdf"), width = 9, height = 5.5, useDingbats = FALSE)
    manhattan(df,
              chr="CHR", bp="BP", p="P", snp="SNP",
              col=c("#2A4F9E", "#FF9D00"),
              genomewideline = -log10(thr_gw),
              suggestiveline = -log10(thr_sg),
              ylim = c(0, ymax),
              cex = 0.45, cex.axis = 0.9,
              main = title_main,
              ylab = expression(-log[10](italic(p))),
              xlab = "Chromosome")
    abline(h = -log10(thr_ex), lty = 2, col = "#555555")
    dev.off()
  }

  ## ---------- QQ (PNG) ----------
  lambda_to_use <- if (!is.null(lambda_override)) lambda_override else lambda_gc(df$P)
  qq_title <- sprintf("QQ plot of the GWAS (MAF=0.1, HWE=1e-6, PC=30)\nlambdaGC = %.3f", lambda_to_use) ###########################

  png(paste0(out_prefix, "_qq.png"),
      width = 2000, height = 2000, res = 300, bg = "white", pointsize = 12)
  qq(df$P,
     main = qq_title,
     xlab = "Expected -log10(p)",
     ylab = "Observed -log10(p)",
     col = "#2A4F9E", cex = 0.6)
  abline(0, 1, col = "#666666", lwd = 1.2)
  dev.off()

  if (isTRUE(write_pdf)) {
    pdf(paste0(out_prefix, "_qq.pdf"), width = 5.5, height = 5.5, useDingbats = FALSE)
    qq(df$P,
       main = qq_title,
       xlab = "Expected -log10(p)",
       ylab = "Observed -log10(p)",
       col = "#2A4F9E", cex = 0.6)
    abline(0, 1, col = "#666666", lwd = 1.2)
    dev.off()
  }

  invisible(lambda_to_use)
}

########################
lam10 <- plot_manhattan_qq_base(
  plink_01_1e6_30,
  title_main = "Manhattan plot of the GWAS (MAF=0.1, HWE=1e-6, PC=30)", 
  out_prefix = "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/00_ArticleFigure/01_Supplementary1Q",
  lambda_override = 1.02914250399434,   # ← PPT값 직접 입력
  write_pdf = TRUE,                # PDF도 만들려면 TRUE
  annotate_top = 0                 # 상위 SNP 라벨(원하면 0)
)

plot_combined_manhattan_qq_horizontal <- function(df,
                                                  title_suffix,
                                                  out_prefix,
                                                  lambda_gc_value,
                                                  thr_gw = 5e-8, thr_sg = 1e-5, thr_ex = 1e-4,
                                                  ymax = NULL) {

  df <- prep_gwas_df(df)

  if (is.null(ymax)) {
    ymax <- ceiling(max(-log10(df$P), na.rm = TRUE) * 1.05)
    ymax <- max(ymax, 10)
  }

  ## ------------------ PDF (가로 1x2 패널) ------------------
  pdf(paste0(out_prefix, "_horizontal.pdf"), width = 12, height = 5, useDingbats = FALSE)
  par(mfrow = c(1,2), mar = c(4.2,4.6,3.0,1.0), mgp = c(2.2,0.8,0))

  ## Panel A - Manhattan
  manhattan(df,
            chr = "CHR", bp = "BP", p = "P", snp = "SNP",
            col = c("#2A4F9E", "#FF9D00"),
            suggestiveline = -log10(thr_sg),
            genomewideline = -log10(thr_gw),
            ylim = c(0, ymax),
            cex = 0.42, cex.axis = 0.85,
            main = paste("Manhattan plot of the GWAS", "\n", title_suffix),
            ylab = expression(-log[10](italic(p))),
            xlab = "Chromosome")
  abline(h = -log10(thr_ex), lty = 2, col = "#666666")

  ## Panel B - QQ
  qq(df$P,
     main = paste0("QQ plot of the GWAS ", title_suffix, "\n",
                   "lambdaGC = ", sprintf("%.3f", lambda_gc_value)),
     xlab = "Expected -log10(p)",
     ylab = "Observed -log10(p)",
     col = "#2A4F9E", cex = 0.6)
  abline(0, 1, col = "#666666", lwd = 1.2)

  dev.off()

  ## ------------------ PNG (고해상도) ------------------
  png(paste0(out_prefix, "_horizontal.png"),
      width = 3600, height = 1500, res = 300)
  par(mfrow = c(1,2), mar = c(4.2,4.6,3.0,1.0), mgp = c(2.2,0.8,0))

  manhattan(df,
            chr = "CHR", bp = "BP", p = "P", snp = "SNP",
            col = c("#2A4F9E", "#FF9D00"),
            suggestiveline = -log10(thr_sg),
            genomewideline = -log10(thr_gw),
            ylim = c(0, ymax),
            cex = 0.42, cex.axis = 0.85,
            main = paste("Manhattan plot of the GWAS", "\n", title_suffix),
            ylab = expression(-log[10](italic(p))),
            xlab = "Chromosome")
  abline(h = -log10(thr_ex), lty = 2, col = "#666666")

  qq(df$P,
     main = paste0("QQ plot of the GWAS ", title_suffix, "\n",
                   "lambdaGC = ", sprintf("%.3f", lambda_gc_value)),
     xlab = "Expected -log10(p)",
     ylab = "Observed -log10(p)",
     col = "#2A4F9E", cex = 0.6)
  abline(0, 1, col = "#666666", lwd = 1.2)

  dev.off()
}


## ===== 실행 예시 (네가 요청한 출력 경로) =====
plot_combined_manhattan_qq_horizontal(
  df = plink_01_1e6_30,
  title_suffix = "(MAF=0.1, HWE=1e-6, PC=30)",
  out_prefix = "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/00_ArticleFigure/01_Supplementary1Q_Combined",
  lambda_gc_value = 1.029   # PPT/표 값 그대로
)




plink_001_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear")
head(plink_001_1e6_10)

plink_001_1e6_10[plink_001_1e6_10$ID == "5:25529170:A:T",]
plink_001_1e6_10[plink_001_1e6_10$ID == "23:32288824:G:A",]
plink_001_1e6_10[plink_001_1e6_10$ID == "24:12787767:A:G",]
plink_001_1e6_10[plink_001_1e6_10$ID == "24:12787779:T:A",]
plink_001_1e6_10[plink_001_1e6_10$ID == "27:44543118:A:T",]
