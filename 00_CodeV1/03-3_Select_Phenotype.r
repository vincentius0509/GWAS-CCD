getwd()
load("../00_RDataFile/03_GWAS_Pure/03-2_Add_Phenotype.rdata")

library(data.table)
library(reshape2)

## data_001_1e6_5_cov_filtered2
eigenval_001_1e6_5 <- fread("../../03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_5.eigenval")[[1]]
head(eigenval_001_1e6_5)

variance_001_1e6_5 <- eigenval_001_1e6_5 / sum(eigenval_001_1e6_5)*100
cumvar_001_1e6_5 <- cumsum(variance_001_1e6_5)
pcs_5 <- seq_along(eigenval_001_1e6_5)

# output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
# file_name <- "03_PCA_001_1e6_5_EP.png"
# output_path <- file.path(output_dir, file_name)
# 
# png(filename = output_path, width = 800, height = 600, res = 60)
# plot(
#     pcs_5, eigenval_001_1e6_5,
#     type = "b", pch = 19, lwd = 2, col = "#1565c0",
#     xlab = "PC Number", ylab = "Eigenvalue",
#     main = "MAF 0.01, HWE 1E-6, PC 5"
# )
# dev.off()
# 
# output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
# file_name <- "03_PCA_001_1e6_5_CV.png"
# output_path <- file.path(output_dir, file_name)
# 
# png(filename = output_path, width = 800, height = 600, res = 60)
# plot(pcs_5, cumvar_001_1e6_5,
#      type = "b", pch = 19, lwd = 2, col = "#00695c",
#      xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
#      ylim = c(0, 100),
#      main = "MAF 0.01, HWE 1E-6, PC 5 [THRESHOLD: 85%]"
# )
# abline(h = 85, lty = 2, col = "red")
# dev.off()

dim(data_001_1e6_5_cov_filtered2)
head(data_001_1e6_5_cov_filtered2)
cov <- data_001_1e6_5_cov_filtered2
exclude_vars <- c("FID", "IID")
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2, ]  # 자기 자신 제외
cor_df$abs_corr <- abs(cor_df$value)

# data.table로 변환
setDT(cor_df)

# 이제 := 사용 가능
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

# unique 쌍만 남기면서 필터링
cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]

# 결과 출력
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
# Var1   Var2      value  abs_corr
# <fctr> <fctr>      <num>     <num>
#   1:    Golden_Retriever    PC3  0.8022424 0.8022424
# 2: German_Shepherd_Dog    PC4 -0.8506943 0.8506943


table(data_001_1e6_5_cov_filtered2$Golden_Retriever) # 3161 366
table(data_001_1e6_5_cov_filtered2$German_Shepherd_Dog) # 3360 167

# data_001_1e6_5_cov_filtered2에서 PC5, Golden_Retriever, German_Shepherd_Dog
dim(data_001_1e6_5_cov_filtered2) # 3527    184
table(data_001_1e6_5_cov_filtered2$Sex_cov)

exclude <- c("PC5","Golden_Retriever", "German_Shepherd_Dog")
data_001_1e6_5_cov_final <- data_001_1e6_5_cov_filtered2[, !..exclude, with=FALSE]
dim(data_001_1e6_5_cov_final) # 3527    181

exclude_vars <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_001_1e6_5_cov_final), exclude_vars)
length(exclude_vars) # 10
length(breed_vars) # 171

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_001_1e6_5_cov_final[[breed]], data_001_1e6_5_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 29
data_001_1e6_5_cov_chrX <- data_001_1e6_5_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_001_1e6_5_cov_chrX) # 3527 152

colnames(data_001_1e6_5_cov_chrX)
exclude_vars <- c("FID", "IID", paste0("PC", 1:4), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 10

breed_vars <- setdiff(names(data_001_1e6_5_cov_chrX), exclude_vars)
length(breed_vars) # 142

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_001_1e6_5_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds) # 41
data_001_1e6_5_cov_num <- data_001_1e6_5_cov_chrX[, !..low_count_breeds, with = FALSE]
dim(data_001_1e6_5_cov_num) # 3527  111

fwrite(data_001_1e6_5_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_001_1e6_5.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)



## data_001_1e6_10_cov_filtered2
eigenval_001_1e6_10 <- fread("../../03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_10.eigenval")[[1]]
head(eigenval_001_1e6_10)

variance_001_1e6_10 <- eigenval_001_1e6_10 / sum(eigenval_001_1e6_10)*100
cumvar_001_1e6_10 <- cumsum(variance_001_1e6_10)
pcs_10 <- seq_along(eigenval_001_1e6_10)

# output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
# file_name <- "03_PCA_001_1e6_10_EP.png"
# output_path <- file.path(output_dir, file_name)
# 
# png(filename = output_path, width = 800, height = 600, res = 60)
# plot(
#     pcs_10, eigenval_001_1e6_10,
#     type = "b", pch = 19, lwd = 2, col = "#1565c0",
#     xlab = "PC Number", ylab = "Eigenvalue",
#     main = "MAF 0.01, HWE 1E-6, PC 10"
# )
# dev.off()
# 
# file_name <- "03_PCA_001_1e6_10_CV.png"
# output_path <- file.path(output_dir, file_name)
# 
# png(filename = output_path, width = 800, height = 600, res = 60)
# plot(pcs_10, cumvar_001_1e6_10,
#      type = "b", pch = 19, lwd = 2, col = "#00695c",
#      xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
#      ylim = c(0, 100),
#      main = "MAF 0.01, HWE 1E-6, PC 10 [THRESHOLD: 85%]"
# )
# abline(h = 85, lty = 2, col = "red")
# dev.off()

# PC 8까지 87% 설명, PC9 & PC 10 제거
dim(data_001_1e6_10_cov_filtered2) # 3527   189
cov <- data_001_1e6_10_cov_filtered2
exclude_vars <- c("FID", "IID")
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)

cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]

print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])

#                   Var1   Var2      value  abs_corr
#                 <fctr> <fctr>      <num>     <num>
# 1:    Golden_Retriever    PC3  0.8022424 0.8022424
# 2: German_Shepherd_Dog    PC4 -0.8506943 0.8506943

# PC3와 상관관계가 높은 Golden_Retriever, German_Shepherd_Dog를 제거
dim(data_001_1e6_10_cov_filtered2) # 3527   189
table(data_001_1e6_10_cov_filtered2$Sex_cov)

exclude <- c("PC9", "PC10", "Golden_Retriever", "German_Shepherd_Dog")
data_001_1e6_10_cov_final <- data_001_1e6_10_cov_filtered2[, !..exclude, with=FALSE]
dim(data_001_1e6_10_cov_final) # 3527 185

exclude_vars <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",
    "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_001_1e6_10_cov_final), exclude_vars)
length(exclude_vars) # 14
length(breed_vars) # 171

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_001_1e6_10_cov_final[[breed]], data_001_1e6_10_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

# 29개 품종이 편향됨
length(biased_breed_vector)
data_001_1e6_10_cov_chrX <- data_001_1e6_10_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_001_1e6_10_cov_chrX) # 3527 156

colnames(data_001_1e6_10_cov_chrX)
exclude_vars <- c("FID", "IID", paste0("PC", 1:8), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 14

breed_vars <- setdiff(names(data_001_1e6_10_cov_chrX), exclude_vars)
length(breed_vars) # 142

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_001_1e6_10_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds) # 41
data_001_1e6_10_cov_num <- data_001_1e6_10_cov_chrX[, !..low_count_breeds, with = FALSE]
dim(data_001_1e6_10_cov_num) # 3527  115

fwrite(data_001_1e6_10_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_001_1e6_10.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

## data_001_1e6_15_cov_filtered2
eigenval_001_1e6_15 <- fread("../../03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_15.eigenval")[[1]]
length(eigenval_001_1e6_15)

variance_001_1e6_15 <- eigenval_001_1e6_15/sum(eigenval_001_1e6_15)*100
cumvar_001_1e6_15 <- cumsum(variance_001_1e6_15)
pcs_15 <- seq_along(eigenval_001_1e6_15)

# output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
# file_name <- "03_PCA_001_1e6_15_EP.png"
# output_path <- file.path(output_dir, file_name)
# 
# png(filename = output_path, width = 800, height = 600, res = 60)
# plot(
#     pcs_15, eigenval_001_1e6_15,
#     type = "b", pch = 19, lwd = 2, col = "#1565c0",
#     xlab = "PC Number", ylab = "Eigenvalue",
#     main = "MAF 0.01, HWE 1E-6, PC 15"
# )
# dev.off()
# 
# file_name <- "03_PCA_001_1e6_15_CV.png"
# output_path <- file.path(output_dir, file_name)
# 
# png(filename = output_path, width = 800, height = 600, res = 60)
# plot(pcs_15, cumvar_001_1e6_15,
#      type = "b", pch = 19, lwd = 2, col = "#00695c",
#      xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
#      ylim = c(0, 100),
#      main = "MAF 0.01, HWE 1E-6, PC 15 [THRESHOLD: 85%]"
# )
# abline(h = 85, lty = 2, col = "red")
# dev.off()

# PC 12까지 88% 정도 설명 (PC13, PC14, PC15 제외)
dim(data_001_1e6_15_cov_filtered2) # 3527   194
cov <- data_001_1e6_15_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 13:15))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)

cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]

print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#                   Var1   Var2      value  abs_corr
#                 <fctr> <fctr>      <num>     <num>
# 1:    Golden_Retriever    PC3  0.8022424 0.8022424
# 2: German_Shepherd_Dog    PC4 -0.8506943 0.8506943

# PC13~PC15 제거
dim(data_001_1e6_15_cov_filtered2) # 3527   194
exclude <- c(paste0("PC", 13:15))
data_001_1e6_15_cov_final <- data_001_1e6_15_cov_filtered2[, !..exclude, with=FALSE]
dim(data_001_1e6_15_cov_final) # 3527 191

exclude_vars <- c("FID", "IID", paste0("PC", 1:12), 
    "Sex_cov","Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_001_1e6_15_cov_final), exclude_vars)
length(exclude_vars) # 18
length(breed_vars) # 173

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_001_1e6_15_cov_final[[breed]], data_001_1e6_15_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 29

data_001_1e6_15_cov_chrX <- data_001_1e6_15_cov_final[, !..biased_breed_vector, with=FALSE]
dim(data_001_1e6_15_cov_chrX) # 3527    162

# VIF_TOO_HIGH => Golden Retriever, German_Shepherd_Dog
exclude <- c("Golden_Retriever", "German_Shepherd_Dog")
data_001_1e6_15_cov_VIF <- data_001_1e6_15_cov_chrX[, !..exclude, with = FALSE]
dim(data_001_1e6_15_cov_VIF) # 3527 160

colnames(data_001_1e6_15_cov_VIF)
exclude_vars <- c("FID", "IID", paste0("PC", 1:12), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 18

breed_vars <- setdiff(names(data_001_1e6_15_cov_VIF), exclude_vars)
length(breed_vars) # 142

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_001_1e6_15_cov_VIF[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}
length(low_count_breeds)
data_001_1e6_15_cov_num <- data_001_1e6_15_cov_VIF[, !..low_count_breeds, with = FALSE]
dim(data_001_1e6_15_cov_num) # 3527 119

fwrite(data_001_1e6_15_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_001_1e6_15.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)


## data_001_1e6_20_cov_filtered2
eigenval_001_1e6_20 <- fread("../../03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_20.eigenval")[[1]]
length(eigenval_001_1e6_20)

variance_001_1e6_20 <- eigenval_001_1e6_20 / sum(eigenval_001_1e6_20) *100
cumvar_001_1e6_20 <- cumsum(variance_001_1e6_20)
pcs_20 <- seq_along(eigenval_001_1e6_20)

# output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
# file_name <- "03_PCA_001_1e6_20_EP.png"
# output_path <- file.path(output_dir, file_name)
# 
# png(filename = output_path, width = 800, height = 600, res = 60)
# plot(
#     pcs_20, eigenval_001_1e6_20,
#     type = "b", pch = 19, lwd = 2, col = "#1565c0",
#     xlab = "PC Number", ylab = "Eigenvalue",
#     main = "MAF 0.01, HWE 1E-6, PC 20"
# )
# dev.off()
# 
# file_name <- "03_PCA_001_1e6_20_CV.png"
# output_path <- file.path(output_dir, file_name)
# 
# png(filename = output_path, width = 800, height = 600, res = 60)
# plot(pcs_20, cumvar_001_1e6_20,
#      type = "b", pch = 19, lwd = 2, col = "#00695c",
#      xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
#      ylim = c(0, 100),
#      main = "MAF 0.01, HWE 1E-6, PC 20 [THRESHOLD: 85%]"
# )
# abline(h = 85, lty = 2, col = "red")
# dev.off()

cumvar_001_1e6_20[16] # 87.9829
# PC16까지 누적 기여율이 88% (PC17, PC18, PC19, PC20 제외)
dim(data_001_1e6_20_cov_filtered2) # 3527   199
cov <- data_001_1e6_20_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 17:20))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)

cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]

print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#                   Var1   Var2      value  abs_corr
#                 <fctr> <fctr>      <num>     <num>
# 1:    Golden_Retriever    PC3  0.8022424 0.8022424
# 2: German_Shepherd_Dog    PC4 -0.8506943 0.8506943

dim(data_001_1e6_20_cov_filtered2) # 3527   199
exclude <- c(paste0("PC", 17:20))
data_001_1e6_20_cov_final <- data_001_1e6_20_cov_filtered2[, !..exclude, with = FALSE]
dim(data_001_1e6_20_cov_final) # 3527 195

exclude_vars <- c("FID", "IID", paste0("PC",1:16),
    "Sex_cov","Age_cov","Intact_cov","pa_score_cov")
breed_vars <- setdiff(names(data_001_1e6_20_cov_final), exclude_vars)
length(exclude_vars)

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_001_1e6_20_cov_final[[breed]], data_001_1e6_20_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 29

data_001_1e6_20_cov_chrX <- data_001_1e6_20_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_001_1e6_20_cov_chrX) # 3527    166

# Golden_Retriever, German_Shepherd_Dog (PC3와 상관관계)
exclude <- c("Golden_Retriever", "German_Shepherd_Dog")
data_001_1e6_20_cov_VIF <- data_001_1e6_20_cov_chrX[, !..exclude, with = FALSE]
dim(data_001_1e6_20_cov_VIF) # 3527 164

cumvar_001_1e6_20[14] # 81.63868
# VIF_TOO_HIGH in PC15 => PC15, PC16 제거
exclude <- c("PC15", "PC16")
data_001_1e6_20_cov_PC <- data_001_1e6_20_cov_VIF[, !..exclude, with = FALSE]
dim(data_001_1e6_20_cov_PC) # 3527  162

colnames(data_001_1e6_20_cov_PC) # PC1~PC14
exclude_vars <- c("FID", "IID", paste0("PC", 1:14), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 20

breed_vars <- setdiff(names(data_001_1e6_20_cov_PC), exclude_vars)
length(breed_vars) # 142

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_001_1e6_20_cov_PC[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds) # 41
data_001_1e6_20_cov_num <- data_001_1e6_20_cov_PC[, !..low_count_breeds, with = FALSE]
dim(data_001_1e6_20_cov_num) # 3527  121
fwrite(data_001_1e6_20_cov_num, "../../03_GWAS_Pure/14_regression_local/04_cov_001_1e6_20.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)


## data_001_1e6_25_cov_filtered2
eigenval_001_1e6_25  <- fread("../../03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_25.eigenval")[[1]]
length(eigenval_001_1e6_25)

variance_001_1e6_25 <- eigenval_001_1e6_25 / sum(eigenval_001_1e6_25) * 100
cumvar_001_1e6_25 <- cumsum(variance_001_1e6_25)
pcs_25 <- seq_along(eigenval_001_1e6_25)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_001_1e6_25_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
    pcs_25, eigenval_001_1e6_25,
    type = "b", pch = 19, lwd = 2, col = "#1565c0",
    xlab = "PC Number", ylab = "Eigenvalue",
    main = "MAF 0.01, HWE 1E-6, PC 25"
)
dev.off()

file_name <- "03_PCA_001_1e6_25_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_25, cumvar_001_1e6_25,
     type = "b", pch = 19, lwd = 2, col = "#00695c",
     xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
     ylim = c(0, 100),
     main = "MAF 0.01, HWE 1E-6, PC 25 [THRESHOLD: 85%]"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_001_1e6_25[19] # 85.14651
# PC19까지 누적기여율이 85%
dim(data_001_1e6_25_cov_filtered2) # 3527   204
cov <- data_001_1e6_25_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 20:25))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)

cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]

print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#                   Var1   Var2      value  abs_corr
#                 <fctr> <fctr>      <num>     <num>
# 1:    Golden_Retriever    PC3  0.8022424 0.8022424
# 2: German_Shepherd_Dog    PC4 -0.8506943 0.8506943
# 3:   Doberman_Pinscher   PC18  0.8308268 0.8308268

dim(data_001_1e6_25_cov_filtered2) # 3527   204
exclude <- c(paste0("PC", 20:25))
data_001_1e6_25_cov_final <- data_001_1e6_25_cov_filtered2[, !..exclude, with = FALSE]
dim(data_001_1e6_25_cov_final) # 3527   198

exclude_vars <- c("FID", "IID", paste0("PC",1:19),
    "Sex_cov","Age_cov","Intact_cov","pa_score_cov")
breed_vars <- setdiff(names(data_001_1e6_25_cov_final), exclude_vars)
length(exclude_vars) # 25

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_001_1e6_25_cov_final[[breed]], data_001_1e6_25_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 29

data_001_1e6_25_cov_chrX <- data_001_1e6_25_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_001_1e6_25_cov_chrX) # 3527 169

#                   Var1   Var2      value  abs_corr
#                 <fctr> <fctr>      <num>     <num>
# 1:    Golden_Retriever    PC3  0.8022424 0.8022424
# 2: German_Shepherd_Dog    PC4 -0.8506943 0.8506943
# 3:   Doberman_Pinscher   PC18  0.8308268 0.8308268
# Golden_Retriever, German_Shepherd_Dog 제거

# 품종 변수 추출: 이용 가능한 품종 컬럼만 추출
breed_vars_chrX <- setdiff(names(data_001_1e6_25_cov_chrX), 
                          c("FID", "IID", paste0("PC",1:19),
                            "Sex_cov","Age_cov","Intact_cov","pa_score_cov"))

# 품종별로 PC15와의 상관계수 계산
cor_results <- sapply(breed_vars_chrX, function(breed){
  cor(data_001_1e6_25_cov_chrX[[breed]], data_001_1e6_25_cov_chrX[["PC15"]], use="complete.obs")
})

# 상관계수 내림차순으로 정렬
cor_results_sorted <- sort(cor_results, decreasing=TRUE)

# PC15와 상관계수가 0.3(절댓값) 이상인 품종 추출 (threshold는 상황에 따라 조절)
high_corr_breeds <- names(cor_results_sorted[abs(cor_results_sorted) >= 0.3])
cor_results_sorted[high_corr_breeds]
#         Miniature_Schnauzer                  Great_Dane 
#                   0.3902890                   0.3537308 
# West_Highland_White_Terrier        Bernese_Mountain_Dog 
#                  -0.3382569                  -0.5656438 

exclude <- c("Golden_Retriever", "Bernese_Mountain_Dog", "German_Shepherd_Dog", "Doberman_Pinscher")
data_001_1e6_25_cov_VIF <- data_001_1e6_25_cov_chrX[, !..exclude, with = FALSE]
dim(data_001_1e6_25_cov_VIF) # 3527 165

colnames(data_001_1e6_25_cov_VIF) # PC1~PC19
exclude_vars <- c("FID", "IID", paste0("PC", 1:19), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 25

breed_vars <- setdiff(names(data_001_1e6_25_cov_VIF), exclude_vars)
length(breed_vars)

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_001_1e6_25_cov_VIF[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds) # 41
data_001_1e6_25_cov_num <- data_001_1e6_25_cov_VIF[, !..low_count_breeds, with = FALSE]
dim(data_001_1e6_25_cov_num) # 3527 124

fwrite(data_001_1e6_25_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_001_1e6_25.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

## data_001_1e6_30_cov_filtered2
eigenval_001_1e6_30  <- fread("../../03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_30.eigenval")[[1]]
length(eigenval_001_1e6_30)

variance_001_1e6_30 <- eigenval_001_1e6_30 / sum(eigenval_001_1e6_30)*100
cumvar_001_1e6_30 <- cumsum(variance_001_1e6_30)
pcs_30 <- seq_along(eigenval_001_1e6_30)


output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_001_1e6_30_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
    pcs_30, eigenval_001_1e6_30,
    type = "b", pch = 19, lwd = 2, col = "#1565c0",
    xlab = "PC Number", ylab = "Eigenvalue",
    main = "MAF 0.01, HWE 1E-6, PC 25"
)
dev.off()

file_name <- "03_PCA_001_1e6_30_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_30, cumvar_001_1e6_30,
     type = "b", pch = 19, lwd = 2, col = "#00695c",
     xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
     ylim = c(0, 100),
     main = "MAF 0.01, HWE 1E-6, PC 25 [THRESHOLD: 85%]"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_001_1e6_30[23]
# PC23까지 누적기여율이 85%

dim(data_001_1e6_30_cov_filtered2) # 3527   209
cov <- data_001_1e6_30_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 24:30))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)

cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]

print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#                   Var1   Var2      value  abs_corr
#                 <fctr> <fctr>      <num>     <num>
# 1:    Golden_Retriever    PC3  0.8022424 0.8022424
# 2: German_Shepherd_Dog    PC4 -0.8506943 0.8506943
# 3:   Doberman_Pinscher   PC18  0.8308268 0.8308268

dim(data_001_1e6_30_cov_filtered2) # 3527   209
exclude <- c(paste0("PC", 24:30))
data_001_1e6_30_cov_final <- data_001_1e6_30_cov_filtered2[, !..exclude, with = FALSE]
dim(data_001_1e6_30_cov_final) # 3527   202

exclude_vars <- c("FID", "IID", paste0("PC", 1:23),
    "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_001_1e6_30_cov_final), exclude_vars)
length(exclude_vars) # 29

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_001_1e6_25_cov_final[[breed]], data_001_1e6_25_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 29

data_001_1e6_30_cov_chrX <- data_001_1e6_30_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_001_1e6_30_cov_chrX) # 3527    173

breed_vars_chrX <- setdiff(names(data_001_1e6_30_cov_chrX), 
    c("FID", "IID", paste0("PC", 1:23), 
    "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov"))

cor_results <- sapply(breed_vars_chrX, function(breed) {
    cor(data_001_1e6_30_cov_chrX[[breed]], data_001_1e6_30_cov_chrX[["PC15"]], use = "complete.obs")
})

cor_results_sorted <- sort(cor_results, decreasing = TRUE)
high_corr_breeds <- names(cor_results_sorted[abs(cor_results_sorted) >= 0.3])
cor_results_sorted[high_corr_breeds]

# PC2
# Labrador_Retriever   Golden_Retriever 
#          0.6684225          0.5180913 

# PC15
#         Miniature_Schnauzer                  Great_Dane 
#                   0.3902890                   0.3537308 
# West_Highland_White_Terrier        Bernese_Mountain_Dog 
#                  -0.3382569                  -0.5656438 

exclude <- c("Labrador_Retriever", "Bernese_Mountain_Dog", "Golden_Retriever", "German_Shepherd_Dog", "Doberman_Pinscher")
data_001_1e6_30_cov_VIF <- data_001_1e6_30_cov_chrX[, !..exclude, with = FALSE]
dim(data_001_1e6_30_cov_VIF) # 3527 168

colnames(data_001_1e6_30_cov_VIF) # PC1~PC23
exclude_vars <- c("FID", "IID", paste0("PC", 1:23), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 29

breed_vars <- setdiff(names(data_001_1e6_30_cov_VIF),exclude_vars)
length(breed_vars) # 139

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_001_1e6_30_cov_VIF[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds) # 41
data_001_1e6_30_cov_num <- data_001_1e6_30_cov_VIF[, !..low_count_breeds, with = FALSE]
dim(data_001_1e6_30_cov_num) # 3527 127

fwrite(data_001_1e6_30_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_001_1e6_30.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)


## data_005_1e6_5_cov_filtered2
dim(data_005_1e6_5_cov_filtered2) # 3513    181
eigenval_005_1e6_5  <- fread("../../03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_5.eigenval")[[1]]
length(eigenval_005_1e6_5)

variance_005_1e6_5 <- eigenval_005_1e6_5 / sum(eigenval_005_1e6_5)*100
cumvar_005_1e6_5 <- cumsum(variance_005_1e6_5)
pcs_5 <- seq_along(eigenval_005_1e6_5)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_005_1e6_5_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_5, eigenval_005_1e6_5,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.05, HWE 1E-6, PC 5"
)
dev.off()

file_name <- "03_PCA_005_1e6_5_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_5, cumvar_005_1e6_5,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.05, HWE 1E-6, PC 5"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_005_1e6_5[4] # 88%

dim(data_005_1e6_5_cov_filtered2) # 3513    181
cov <- data_005_1e6_5_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 5))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])

dim(data_005_1e6_5_cov_filtered2) # 3513    181
exclude <- c("PC5")
data_005_1e6_5_cov_final <- data_005_1e6_5_cov_filtered2[, !..exclude, with = FALSE]
dim(data_005_1e6_5_cov_final) # 3513    180

exclude_vars <- c("FID", "IID", paste0("PC", 1:4),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_005_1e6_5_cov_final), exclude_vars)
length(exclude_vars) # 10

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_005_1e6_5_cov_final[[breed]], data_005_1e6_5_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
  biased_breeds[[breed]] <- breed_1_sex_count
  biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_005_1e6_5_cov_chrX <- data_005_1e6_5_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_005_1e6_5_cov_chrX) # 3513 153

colnames(data_005_1e6_5_cov_chrX) # PC1~PC4
exclude_vars <- c("FID", "IID", paste0("PC", 1:4), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars)

breed_vars <- setdiff(names(data_005_1e6_5_cov_chrX), exclude_vars)
length(breed_vars)

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_005_1e6_5_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds) # 42
data_005_1e6_5_cov_num <- data_005_1e6_5_cov_chrX[, !..low_count_breeds, with = FALSE]
dim(data_005_1e6_5_cov_num) # 3513  111

fwrite(data_005_1e6_5_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_005_1e6_5.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_005_1e6_10_cov_filtered2
dim(data_005_1e6_10_cov_filtered2) # 3513   186
eigenval_005_1e6_10  <- fread("../../03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_10.eigenval")[[1]]
length(eigenval_005_1e6_10)

variance_005_1e6_10 <- eigenval_005_1e6_10 / sum(eigenval_005_1e6_10)*100
cumvar_005_1e6_10 <- cumsum(variance_005_1e6_10)
pcs_10 <- seq_along(eigenval_005_1e6_10)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_005_1e6_10_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_10, eigenval_005_1e6_10,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.05, HWE 1E-6, PC 10"
)
dev.off()

file_name <- "03_PCA_005_1e6_10_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_10, cumvar_005_1e6_10,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.05, HWE 1E-6, PC 10"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_005_1e6_10[8] # 89%

dim(data_005_1e6_10_cov_filtered2) # 3513   186
cov <- data_005_1e6_10_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 9:10))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9367307 0.9367307

exclude <- c(paste0("PC", 9:10))
data_005_1e6_10_cov_final <- data_005_1e6_10_cov_filtered2[, !..exclude, with = FALSE]
dim(data_005_1e6_10_cov_final) # 3513   184

exclude_vars <- c("FID", "IID", paste0("PC", 1:8),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_005_1e6_10_cov_final), exclude_vars)
length(exclude_vars) # 14

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_005_1e6_10_cov_final[[breed]], data_005_1e6_10_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_005_1e6_10_cov_chrX <- data_005_1e6_10_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_005_1e6_10_cov_chrX) # 3513    157

colnames(data_005_1e6_10_cov_chrX)
exclude_vars <- c("FID", "IID", paste0("PC", 1:8), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 14

breed_vars <- setdiff(names(data_005_1e6_10_cov_chrX), exclude_vars)
length(breed_vars)

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_005_1e6_10_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds) # 42
data_005_1e6_10_cov_num <- data_005_1e6_10_cov_chrX[, !..low_count_breeds, with = FALSE]
dim(data_005_1e6_10_cov_num) # 3527  115

fwrite(data_005_1e6_10_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_005_1e6_10.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_005_1e6_15_cov_filtered2
dim(data_005_1e6_15_cov_filtered2) # 3513   191
eigenval_005_1e6_15  <- fread("../../03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_15.eigenval")[[1]]
length(eigenval_005_1e6_15)

variance_005_1e6_15 <- eigenval_005_1e6_15 / sum(eigenval_005_1e6_15)*100
cumvar_005_1e6_15 <- cumsum(variance_005_1e6_15)
pcs_15 <- seq_along(eigenval_005_1e6_15)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_005_1e6_15_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_15, eigenval_005_1e6_15,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.05, HWE 1E-6, PC 15"
)
dev.off()

file_name <- "03_PCA_005_1e6_15_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_15, cumvar_005_1e6_15,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.05, HWE 1E-6, PC 15"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_005_1e6_15[11] # 87%
cumvar_005_1e6_15[10] # 83%

dim(data_005_1e6_15_cov_filtered2) # 3513   191
cov <- data_005_1e6_15_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 12:15))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9367307 0.9367307

exclude <- c(paste0("PC", 12:15))
data_005_1e6_15_cov_final <- data_005_1e6_15_cov_filtered2[, !..exclude, with = FALSE]
dim(data_005_1e6_15_cov_final) # 3513   187

exclude_vars <- c("FID", "IID", paste0("PC", 1:11),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_005_1e6_15_cov_final), exclude_vars)
length(exclude_vars) # 17

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_005_1e6_15_cov_final[[breed]], data_005_1e6_15_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_005_1e6_15_cov_chrX <- data_005_1e6_15_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_005_1e6_15_cov_chrX) # 3513    160

colnames(data_005_1e6_15_cov_chrX)
exclude_vars <- c("FID", "IID", paste0("PC", 1:11), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 17

breed_vars <- setdiff(names(data_005_1e6_15_cov_chrX), exclude_vars)
length(breed_vars)

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_005_1e6_15_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds) # 42
data_005_1e6_15_cov_num <- data_005_1e6_15_cov_chrX[, !..low_count_breeds, with = FALSE]
dim(data_005_1e6_15_cov_num)

fwrite(data_005_1e6_15_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_005_1e6_15.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_005_1e6_20_cov_filtered2
eigenval_005_1e6_20  <- fread("../../03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_20.eigenval")[[1]]
length(eigenval_005_1e6_20)

variance_005_1e6_20 <- eigenval_005_1e6_20 / sum(eigenval_005_1e6_20)*100
cumvar_005_1e6_20 <- cumsum(variance_005_1e6_20)
pcs_20 <- seq_along(eigenval_005_1e6_20)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_005_1e6_20_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_20, eigenval_005_1e6_20,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.05, HWE 1E-6, PC 20"
)
dev.off()

file_name <- "03_PCA_005_1e6_20_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_20, cumvar_005_1e6_20,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.05, HWE 1E-6, PC 20"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_005_1e6_20[13] # 86%

dim(data_005_1e6_20_cov_filtered2) # 3513   196
cov <- data_005_1e6_20_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 14:20))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9367307 0.9367307

exclude <- c(paste0("PC", 14:20))
data_005_1e6_20_cov_final <- data_005_1e6_20_cov_filtered2[, !..exclude, with = FALSE]
dim(data_005_1e6_20_cov_final) # 3513   189

exclude_vars <- c("FID", "IID", paste0("PC", 1:13),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_005_1e6_20_cov_final), exclude_vars)
length(exclude_vars) # 19

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_005_1e6_20_cov_final[[breed]], data_005_1e6_20_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_005_1e6_20_cov_chrX <- data_005_1e6_20_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_005_1e6_20_cov_chrX) # 162

colnames(data_005_1e6_20_cov_chrX)
exclude_vars_20 <- c("FID", "IID", paste0("PC", 1:13), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars_20) # 19

breed_vars_20 <- setdiff(names(data_005_1e6_20_cov_chrX), exclude_vars_20)
length(breed_vars_20)

low_count_breeds_20 <- c()
for (breed in breed_vars_20) {
  count <- sum(data_005_1e6_20_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds_20 <- c(low_count_breeds_20, breed)
  }
}

length(low_count_breeds_20)
data_005_1e6_20_cov_num <- data_005_1e6_20_cov_chrX[, !..low_count_breeds_20, with = FALSE]
dim(data_005_1e6_20_cov_num)

fwrite(data_005_1e6_20_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_005_1e6_20.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_005_1e6_25_cov_filtered2
dim(data_005_1e6_25_cov_filtered2) # 3513   201
eigenval_005_1e6_25  <- fread("../../03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_25.eigenval")[[1]]
length(eigenval_005_1e6_25)

variance_005_1e6_25 <- eigenval_005_1e6_25 / sum(eigenval_005_1e6_25)*100
cumvar_005_1e6_25 <- cumsum(variance_005_1e6_25)
pcs_25 <- seq_along(eigenval_005_1e6_25)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_005_1e6_25_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_25, eigenval_005_1e6_25,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.05, HWE 1E-6, PC 25"
)
dev.off()

file_name <- "03_PCA_005_1e6_25_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_25, cumvar_005_1e6_25,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.05, HWE 1E-6, PC 25"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_005_1e6_25[16] # 86%
cumvar_005_1e6_25[15]

dim(data_005_1e6_25_cov_filtered2) # 3513   201
cov <- data_005_1e6_25_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 17:25))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9367307 0.9367307

exclude <- c(paste0("PC", 17:25))
data_005_1e6_25_cov_final <- data_005_1e6_25_cov_filtered2[, !..exclude, with = FALSE]
dim(data_005_1e6_25_cov_final) # 3513   192

exclude_vars <- c("FID", "IID", paste0("PC", 1:16),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_005_1e6_25_cov_final), exclude_vars)
length(exclude_vars) # 22

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_005_1e6_25_cov_final[[breed]], data_005_1e6_25_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_005_1e6_25_cov_chrX <- data_005_1e6_25_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_005_1e6_25_cov_chrX) # 3513    165

colnames(data_005_1e6_25_cov_chrX)
exclude_vars_25 <- c("FID", "IID", paste0("PC", 1:16), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars_25) # 22

breed_vars_25 <- setdiff(names(data_005_1e6_25_cov_chrX), exclude_vars_25)
length(breed_vars_25)

low_count_breeds_25 <- c()
for (breed in breed_vars_25) {
  count <- sum(data_005_1e6_25_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds_25 <- c(low_count_breeds_25, breed)
  }
}

length(low_count_breeds_25)
data_005_1e6_25_cov_num <- data_005_1e6_25_cov_chrX[, !..low_count_breeds_25, with = FALSE]
dim(data_005_1e6_25_cov_num)

fwrite(data_005_1e6_25_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_005_1e6_25.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_005_1e6_30_cov_filtered2 
dim(data_005_1e6_30_cov_filtered2) # 3513   206

eigenval_005_1e6_30  <- fread("../../03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_30.eigenval")[[1]]
length(eigenval_005_1e6_30)

variance_005_1e6_30 <- eigenval_005_1e6_30 / sum(eigenval_005_1e6_30)*100
cumvar_005_1e6_30 <- cumsum(variance_005_1e6_30)
pcs_30 <- seq_along(eigenval_005_1e6_30)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_005_1e6_30_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_30, eigenval_005_1e6_30,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.05, HWE 1E-6, PC 30"
)
dev.off()

file_name <- "03_PCA_005_1e6_30_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_30, cumvar_005_1e6_30,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.05, HWE 1E-6, PC 30"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_005_1e6_30[19] # 85%

dim(data_005_1e6_30_cov_filtered2) # 3513   206
cov <- data_005_1e6_30_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 20:30))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9367307 0.9367307

exclude <- c(paste0("PC", 20:30))
data_005_1e6_30_cov_final <- data_005_1e6_30_cov_filtered2[, !..exclude, with = FALSE]
dim(data_005_1e6_30_cov_final) # 3513   195

exclude_vars <- c("FID", "IID", paste0("PC", 1:30),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_005_1e6_30_cov_final), exclude_vars)
length(exclude_vars) # 36

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_005_1e6_30_cov_final[[breed]], data_005_1e6_30_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_005_1e6_30_cov_chrX <- data_005_1e6_30_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_005_1e6_30_cov_chrX) # 3527    168

colnames(data_005_1e6_30_cov_chrX)
exclude_vars_30 <- c("FID", "IID", paste0("PC", 1:19), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars_30) # 25

breed_vars_30 <- setdiff(names(data_005_1e6_30_cov_chrX), exclude_vars_30)
length(breed_vars_30)

low_count_breeds_30 <- c()
for (breed in breed_vars_30) {
  count <- sum(data_005_1e6_30_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds_30 <- c(low_count_breeds_30, breed)
  }
}

length(low_count_breeds_30)
data_005_1e6_30_cov_num <- data_005_1e6_30_cov_chrX[, !..low_count_breeds_30, with = FALSE]
dim(data_005_1e6_30_cov_num) # 3513 126

fwrite(data_005_1e6_30_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_005_1e6_30.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_01_1e6_5_cov_filtered2
eigenval_01_1e6_5  <- fread("../../03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_5.eigenval")[[1]]
length(eigenval_01_1e6_5)

variance_01_1e6_5 <- eigenval_01_1e6_5 / sum(eigenval_01_1e6_5)*100
cumvar_01_1e6_5 <- cumsum(variance_01_1e6_5)
pcs_5 <- seq_along(eigenval_01_1e6_5)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_01_1e6_5_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_5, eigenval_01_1e6_5,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.1, HWE 1E-6, PC 5"
)
dev.off()

file_name <- "03_PCA_01_1e6_5_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_5, cumvar_01_1e6_5,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.1, HWE 1E-6, PC 5"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_01_1e6_5[4] # 89%

dim(data_01_1e6_5_cov_filtered2) # 3512 181
cov <- data_01_1e6_5_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC",5))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9664447 0.9664447

exclude <- c(paste0("PC",5))
data_01_1e6_5_cov_final <- data_01_1e6_5_cov_filtered2[, !..exclude, with = FALSE]
dim(data_01_1e6_5_cov_final) # 3512 180

exclude_vars <- c("FID", "IID", paste0("PC", 1:4),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_01_1e6_5_cov_final), exclude_vars)
length(exclude_vars) # 10

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_01_1e6_5_cov_final[[breed]], data_01_1e6_5_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_01_1e6_5_cov_chrX <- data_01_1e6_5_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_01_1e6_5_cov_chrX) # 3512  153

colnames(data_01_1e6_5_cov_chrX)
exclude_vars <- c("FID", "IID", paste0("PC", 1:4), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 10

breed_vars <- setdiff(names(data_01_1e6_5_cov_chrX), exclude_vars)
length(breed_vars) # 143

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_01_1e6_5_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds)
data_01_1e6_5_cov_num <- data_01_1e6_5_cov_chrX[, !..low_count_breeds_30, with = FALSE]
dim(data_01_1e6_5_cov_num) # 3512 111

fwrite(data_01_1e6_5_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_01_1e6_5.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_01_1e6_10_cov_filtered2
eigenval_01_1e6_10  <- fread("../../03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_10.eigenval")[[1]]
length(eigenval_01_1e6_10)

variance_01_1e6_10 <- eigenval_01_1e6_10 / sum(eigenval_01_1e6_10)*100
cumvar_01_1e6_10 <- cumsum(variance_01_1e6_10)
pcs_10 <- seq_along(eigenval_01_1e6_10)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_01_1e6_10_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_10, eigenval_01_1e6_10,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.1, HWE 1E-6, PC 10"
)
dev.off()

file_name <- "03_PCA_01_1e6_10_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_10, cumvar_01_1e6_10,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.1, HWE 1E-6, PC 10"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_01_1e6_10[7] # 85

dim(data_01_1e6_10_cov_filtered2) # 3512    186
cov <- data_01_1e6_10_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 8:10))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9664447 0.9664447

exclude <- c(paste0("PC", 8:10))
data_01_1e6_10_cov_final <- data_01_1e6_10_cov_filtered2[, !..exclude, with = FALSE]
dim(data_01_1e6_10_cov_final) # 3512    183

exclude_vars <- c("FID", "IID", paste0("PC", 1:7),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_01_1e6_10_cov_final), exclude_vars)
length(exclude_vars) # 13
length(breed_vars)

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_01_1e6_10_cov_final[[breed]], data_01_1e6_10_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_01_1e6_10_cov_chrX <- data_01_1e6_10_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_01_1e6_10_cov_chrX) # 3512 156

colnames(data_01_1e6_10_cov_chrX)
exclude_vars <- c("FID", "IID", paste0("PC", 1:7), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 13

breed_vars <- setdiff(names(data_01_1e6_10_cov_chrX), exclude_vars)
length(breed_vars)

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_01_1e6_10_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds)
data_01_1e6_10_cov_num <- data_01_1e6_10_cov_chrX[, !..low_count_breeds, with = FALSE]
dim(data_01_1e6_10_cov_num) # 3512  114

fwrite(data_01_1e6_10_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_01_1e6_10.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_01_1e6_15_cov_filtered2
eigenval_01_1e6_15  <- fread("../../03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_15.eigenval")[[1]]
length(eigenval_01_1e6_15)

variance_01_1e6_15 <- eigenval_01_1e6_15 / sum(eigenval_01_1e6_15)*100
cumvar_01_1e6_15 <- cumsum(variance_01_1e6_15)
pcs_15 <- seq_along(eigenval_01_1e6_15)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_01_1e6_15_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_15, eigenval_01_1e6_15,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.1, HWE 1E-6, PC 15"
)
dev.off()

file_name <- "03_PCA_01_1e6_15_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_15, cumvar_01_1e6_15,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.1, HWE 1E-6, PC 15"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_01_1e6_15[10] # 85%

dim(data_01_1e6_15_cov_filtered2) # 3512    191
cov <- data_01_1e6_15_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 11:15))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9664447 0.9664447

exclude <- c(paste0("PC", 11:15))
data_01_1e6_15_cov_final <- data_01_1e6_15_cov_filtered2[, !..exclude, with = FALSE]
dim(data_01_1e6_15_cov_final) # 3512    186

exclude_vars <- c("FID", "IID", paste0("PC", 1:10),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_01_1e6_15_cov_final), exclude_vars)
length(exclude_vars)

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_01_1e6_15_cov_final[[breed]], data_01_1e6_15_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector)

data_01_1e6_15_cov_chrX <- data_01_1e6_15_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_01_1e6_15_cov_chrX) # 3512 159

colnames(data_01_1e6_15_cov_chrX) # ~PC10
exclude_vars <- c("FID", "IID", paste0("PC", 1:10), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 16

breed_vars <- setdiff(names(data_01_1e6_15_cov_chrX), exclude_vars)
length(breed_vars)

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_01_1e6_15_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds)
data_01_1e6_15_cov_num <- data_01_1e6_15_cov_chrX[, !..low_count_breeds, with = FALSE]
dim(data_01_1e6_15_cov_num) # 3512  117

fwrite(data_01_1e6_15_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_01_1e6_15.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_01_1e6_20_cov_filtered2
eigenval_01_1e6_20  <- fread("../../03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_20.eigenval")[[1]]
length(eigenval_01_1e6_20)

variance_01_1e6_20 <- eigenval_01_1e6_20 / sum(eigenval_01_1e6_20)*100
cumvar_01_1e6_20 <- cumsum(variance_01_1e6_20)
pcs_20 <- seq_along(eigenval_01_1e6_20)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_01_1e6_20_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_20, eigenval_01_1e6_20,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.1, HWE 1E-6, PC 20"
)
dev.off()

file_name <- "03_PCA_01_1e6_20_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_20, cumvar_01_1e6_20,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.1, HWE 1E-6, PC 20"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_01_1e6_20[13] # 85

dim(data_01_1e6_20_cov_filtered2) # 3512 196
cov <- data_01_1e6_20_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 14:20))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9664447 0.9664447

exclude <- c(paste0("PC", 14:20))
data_01_1e6_20_cov_final <- data_01_1e6_20_cov_filtered2[, !..exclude, with = FALSE]
dim(data_01_1e6_20_cov_final) # 3512    189

exclude_vars <- c("FID", "IID", paste0("PC", 1:13),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_01_1e6_20_cov_final), exclude_vars)
length(exclude_vars) # 19

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_01_1e6_20_cov_final[[breed]], data_01_1e6_20_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_01_1e6_20_cov_chrX <- data_01_1e6_20_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_01_1e6_20_cov_chrX) # 3512 162

colnames(data_01_1e6_20_cov_chrX) # ~PC13
exclude_vars <- c("FID", "IID", paste0("PC", 1:13), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 19

breed_vars <- setdiff(names(data_01_1e6_20_cov_chrX), exclude_vars)
length(breed_vars)

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_01_1e6_20_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds)
data_01_1e6_20_cov_num <- data_01_1e6_20_cov_chrX[, !..low_count_breeds, with = FALSE]
dim(data_01_1e6_20_cov_num) # 3512  120

fwrite(data_01_1e6_20_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_01_1e6_20.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_01_1e6_25_cov_filtered2
eigenval_01_1e6_25  <- fread("../../03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_25.eigenval")[[1]]
length(eigenval_01_1e6_25)

variance_01_1e6_25 <- eigenval_01_1e6_25 / sum(eigenval_01_1e6_25)*100
cumvar_01_1e6_25 <- cumsum(variance_01_1e6_25)
pcs_25 <- seq_along(eigenval_01_1e6_25)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_01_1e6_25_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_25, eigenval_01_1e6_25,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.1, HWE 1E-6, PC 25"
)
dev.off()

file_name <- "03_PCA_01_1e6_25_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_25, cumvar_01_1e6_25,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.1, HWE 1E-6, PC 25"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_01_1e6_25[16] # 85

dim(data_01_1e6_25_cov_filtered2) # 3512    201
cov <- data_01_1e6_25_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 17:25))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9664447 0.9664447


exclude <- c(paste0("PC", 17:25))
data_01_1e6_25_cov_final <- data_01_1e6_25_cov_filtered2[, !..exclude, with = FALSE]
dim(data_01_1e6_25_cov_final) # 3512    192

exclude_vars <- c("FID", "IID", paste0("PC", 1:16),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_01_1e6_25_cov_final), exclude_vars)
length(exclude_vars) # 22

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_01_1e6_25_cov_final[[breed]], data_01_1e6_25_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_01_1e6_25_cov_chrX <- data_01_1e6_25_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_01_1e6_25_cov_chrX) # 3512 165

colnames(data_01_1e6_25_cov_chrX)
exclude_vars <- c("FID", "IID", paste0("PC", 1:16), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 22

breed_vars <- setdiff(names(data_01_1e6_25_cov_chrX), exclude_vars)
length(breed_vars)

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_01_1e6_25_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds)
data_01_1e6_25_cov_num <- data_01_1e6_25_cov_chrX[, !..low_count_breeds, with = FALSE]
dim(data_01_1e6_25_cov_num) # 3512  123

fwrite(data_01_1e6_25_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_01_1e6_25.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

# data_01_1e6_30_cov_filtered2
eigenval_01_1e6_30  <- fread("../../03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_30.eigenval")[[1]]
length(eigenval_01_1e6_30)

variance_01_1e6_30 <- eigenval_01_1e6_30 / sum(eigenval_01_1e6_30)*100
cumvar_01_1e6_30 <- cumsum(variance_01_1e6_30)
pcs_30 <- seq_along(eigenval_01_1e6_30)

output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"
file_name <- "03_PCA_01_1e6_30_EP.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(
  pcs_30, eigenval_01_1e6_30,
  type = "b", pch = 19, lwd = 2, col = "#1565c0",
  xlab = "PC Number", ylab = "Eigenvalue",
  main = "MAF 0.1, HWE 1E-6, PC 30"
)
dev.off()

file_name <- "03_PCA_01_1e6_30_CV.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(pcs_30, cumvar_01_1e6_30,
   type = "b", pch = 19, lwd = 2, col = "#00695c",
   xlab = "PC Number", ylab = "Cumulative Variance Explained (%)",
   ylim = c(0, 100),
   main = "MAF 0.1, HWE 1E-6, PC 30"
)
abline(h = 85, lty = 2, col = "red")
dev.off()

cumvar_01_1e6_30[19] # 85%

dim(data_01_1e6_30_cov_filtered2) # 3512    206
cov <- data_01_1e6_30_cov_filtered2
exclude_vars <- c("FID","IID",paste0("PC", 20:30))
check_vars <- setdiff(names(cov), exclude_vars)
cor_mat <- cor(cov[, ..check_vars], use = "complete.obs")
threshold <- 0.8
cor_df <- melt(cor_mat)
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2,]
cor_df$abs_corr <- abs(cor_df$value)

setDT(cor_df)
cor_df[, Var_min := pmin(as.character(Var1), as.character(Var2))]
cor_df[, Var_max := pmax(as.character(Var1), as.character(Var2))]

cor_df_unique <- cor_df[abs_corr >= threshold]
cor_df_unique <- cor_df_unique[!duplicated(paste(Var_min, Var_max, sep="--"))]
print(cor_df_unique[, .(Var1, Var2, value, abs_corr)])
#       Var1   Var2      value  abs_corr
#     <fctr> <fctr>      <num>     <num>
# 1: Sex_cov    PC1 -0.9664447 0.9664447

exclude <- c(paste0("PC", 20:30))
data_01_1e6_30_cov_final <- data_01_1e6_30_cov_filtered2[, !..exclude, with = FALSE]
dim(data_01_1e6_30_cov_final) # 3512    195

exclude_vars <- c("FID", "IID", paste0("PC", 1:19),
  "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
breed_vars <- setdiff(names(data_01_1e6_30_cov_final), exclude_vars)
length(exclude_vars) # 25

biased_breeds <- list()
biased_breed_vector <- c()

for (breed in breed_vars) {
  tab <- table(data_01_1e6_30_cov_final[[breed]], data_01_1e6_30_cov_final$Sex_cov)
  breed_1_sex_count <- tab["1", , drop=FALSE]
  nonzero_sexes <- sum(breed_1_sex_count > 0)
  if (nonzero_sexes == 1) {
    biased_breeds[[breed]] <- breed_1_sex_count
    biased_breed_vector <- c(biased_breed_vector, breed)
  }
}

if (length(biased_breeds) == 0) {
  cat("모든 품종이 성별 분포에 골고루 포함되어 있습니다.\n")
} else {
  cat("성별 분포가 편향된 품종 리스트 및 분포:\n")  
  print(biased_breeds)
  cat("총", length(biased_breed_vector), "개 품종이 편향됨\n")
}

length(biased_breed_vector) # 27

data_01_1e6_30_cov_chrX <- data_01_1e6_30_cov_final[, !..biased_breed_vector, with = FALSE]
dim(data_01_1e6_30_cov_chrX) # 168

colnames(data_01_1e6_30_cov_chrX)
exclude_vars <- c("FID", "IID", paste0("PC", 1:19), "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")
length(exclude_vars) # 25

breed_vars <- setdiff(names(data_01_1e6_30_cov_chrX), exclude_vars)
length(breed_vars)

low_count_breeds <- c()
for (breed in breed_vars) {
  count <- sum(data_01_1e6_30_cov_chrX[[breed]] == 1)
  if (count <= 5) {
    low_count_breeds <- c(low_count_breeds, breed)
  }
}

length(low_count_breeds)
data_01_1e6_30_cov_num <- data_01_1e6_30_cov_chrX[, !..low_count_breeds, with = FALSE]
dim(data_01_1e6_30_cov_num) # 3512  126

fwrite(data_01_1e6_30_cov_num, "../../03_GWAS_Pure/14_Regression_local/04_cov_01_1e6_30.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)
