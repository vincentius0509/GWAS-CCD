getwd()

library(pacman)

p_load(ggplot2, dplyr, data.table)
p_load(RColorBrewer, viridis)

# 03-1_HWE.sh 실행
# 03-2_PCA.sh 실행
# MAF_THRESHOLDS=(0.01 0.05 0.1)
# HWE_THRESHOLDS=(1e-6)
# PCA_DIMENSIONS=(5 10 15 20 25 30)

library(data.table)

# MAF 0.01
data_001_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_5.eigenvec")
data_001_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_10.eigenvec")
data_001_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_15.eigenvec")
data_001_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_20.eigenvec")
data_001_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_25.eigenvec")
data_001_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_30.eigenvec")

colnames(data_001_1e6_5) <- c("FID", "IID", paste0("PC", 1:5))
colnames(data_001_1e6_10) <- c("FID", "IID", paste0("PC", 1:10))
colnames(data_001_1e6_15) <- c("FID", "IID", paste0("PC", 1:15))
colnames(data_001_1e6_20) <- c("FID", "IID", paste0("PC", 1:20))
colnames(data_001_1e6_25) <- c("FID", "IID", paste0("PC", 1:25))
colnames(data_001_1e6_30) <- c("FID", "IID", paste0("PC", 1:30))

# MAF 0.05
data_005_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_5.eigenvec")
data_005_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_10.eigenvec")
data_005_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_15.eigenvec")
data_005_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_20.eigenvec")
data_005_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_25.eigenvec")
data_005_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_30.eigenvec")

colnames(data_005_1e6_5) <- c("FID", "IID", paste0("PC", 1:5))
colnames(data_005_1e6_10) <- c("FID", "IID", paste0("PC", 1:10))    
colnames(data_005_1e6_15) <- c("FID", "IID", paste0("PC", 1:15))
colnames(data_005_1e6_20) <- c("FID", "IID", paste0("PC", 1:20))
colnames(data_005_1e6_25) <- c("FID", "IID", paste0("PC", 1:25))
colnames(data_005_1e6_30) <- c("FID", "IID", paste0("PC", 1:30))

# MAF 0.1
data_01_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_5.eigenvec")
data_01_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_10.eigenvec")
data_01_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_15.eigenvec")
data_01_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_20.eigenvec")
data_01_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_25.eigenvec")
data_01_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_30.eigenvec")

colnames(data_01_1e6_5) <- c("FID", "IID", paste0("PC", 1:5))
colnames(data_01_1e6_10) <- c("FID", "IID", paste0("PC", 1:10))
colnames(data_01_1e6_15) <- c("FID", "IID", paste0("PC", 1:15))
colnames(data_01_1e6_20) <- c("FID", "IID", paste0("PC", 1:20))
colnames(data_01_1e6_25) <- c("FID", "IID", paste0("PC", 1:25))
colnames(data_01_1e6_30) <- c("FID", "IID", paste0("PC", 1:30))

purebred_total <- read.csv("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/00_purebred_total.csv")  
head(purebred_total)
dim(purebred_total)
colnames(purebred_total)

head(data_001_1e6_5)
eigenval_001_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_5.eigenval")
head(eigenval_001_1e6_5)

variance_001_1e6_5 <- eigenval_001_1e6_5 / sum(eigenval_001_1e6_5)
variance_001_1e6_5[1:5]*100

pc1_001_1e6_5 <- round(variance_001_1e6_5[1]*100, 2) # 30.69
pc2_001_1e6_5 <- round(variance_001_1e6_5[2]*100, 2) # 20.42

plot(data_001_1e6_5$PC1, data_001_1e6_5$PC2, xlab = paste0("PC1: ", pc1_001_1e6_5,"%"), 
  ylab = paste0("PC2: ", pc2_001_1e6_5, "%"), main = "PCA plot (MAF 0.01, HWE 1e-6, PCA 5)", col = "blue")

# 디렉토리 경로
output_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/"

# 파일 이름
file_name <- "03_PCA_001_1e6_5.png"

# 전체 경로 생성
output_path <- file.path(output_dir, file_name)
# PNG 파일로 저장
png(filename = output_path, width = 800, height = 600, res = 60)  # res=120 → 해상도 조정 가능
plot(data_001_1e6_5$PC1, data_001_1e6_5$PC2,
     xlab = paste0("PC1: ", pc1_001_1e6_5, "%"),
     ylab = paste0("PC2: ", pc2_001_1e6_5, "%"),
     main = "PCA plot (MAF 0.01, HWE 1e-6, PCA 5)",
     col = "blue")
dev.off()

var_001_1e6_5 <- variance_001_1e6_5[1:5] * 100
cum_001_1e6_5 <- cumsum(var_001_1e6_5)

head(data_001_1e6_10)
eigenval_001_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_10.eigenval")
head(eigenval_001_1e6_10)

variance_001_1e6_10 <- eigenval_001_1e6_10 / sum(eigenval_001_1e6_10)
variance_001_1e6_10[1:10]*100

pc1_001_1e6_10 <- round(variance_001_1e6_10[1]*100, 2) # 19.76
pc2_001_1e6_10 <- round(variance_001_1e6_10[2]*100, 2) # 13.16

output_dir
file_name <- "03_PCA_001_1e6_10.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_001_1e6_10$PC1, data_001_1e6_10$PC2,
     xlab = paste0("PC1: ", pc1_001_1e6_10, "%"),
     ylab = paste0("PC2: ", pc2_001_1e6_10, "%"),
     main = "PCA plot (MAF 0.01, HWE 1e-6, PCA 10)",
     col = "blue")
dev.off()

var_001_1e6_10 <- variance_001_1e6_10[1:10] * 100
cum_001_1e6_10 <- cumsum(var_001_1e6_10)
cum_001_1e6_10

head(data_001_1e6_15)
eigenval_001_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_15.eigenval")
head(eigenval_001_1e6_15)

variance_001_1e6_15 <- eigenval_001_1e6_15 / sum(eigenval_001_1e6_15)
variance_001_1e6_15[1:15]*100

pc1_001_1e6_15 <- round(variance_001_1e6_15[1]*100, 2) # 15.71
pc2_001_1e6_15 <- round(variance_001_1e6_15[2]*100, 2) # 10.45

output_dir
file_name <- "03_PCA_001_1e6_15.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_001_1e6_15$PC1, data_001_1e6_15$PC2,
     xlab = paste0("PC1: ", pc1_001_1e6_15, "%"),
     ylab = paste0("PC2: ", pc2_001_1e6_15, "%"),
     main = "PCA plot (MAF 0.01, HWE 1e-6, PCA 15)",
     col = "blue")
dev.off()

var_001_1e6_15 <- variance_001_1e6_15[1:15] * 100
cum_001_1e6_15 <- cumsum(var_001_1e6_15)
cum_001_1e6_15

head(data_001_1e6_20)
eigen_001_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_20.eigenval")
head(eigen_001_1e6_20)

variance_001_1e6_20 <- eigen_001_1e6_20 / sum(eigen_001_1e6_20)
variance_001_1e6_20[1:20]*100

pc1_001_1e6_20 <- round(variance_001_1e6_20[1]*100, 2) # 13.33
pc2_001_1e6_20 <- round(variance_001_1e6_20[2]*100, 2) # 8.87

output_dir
file_name <- "03_PCA_001_1e6_20.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_001_1e6_20$PC1, data_001_1e6_20$PC2,
     xlab = paste0("PC1: ", pc1_001_1e6_20, "%"),
     ylab = paste0("PC2: ", pc2_001_1e6_20, "%"),
     main = "PCA plot (MAF 0.01, HWE 1e-6, PCA 20)",
     col = "blue")
dev.off()

var_001_1e6_20 <- variance_001_1e6_20[1:20] * 100
cum_001_1e6_20 <- cumsum(var_001_1e6_20)
cum_001_1e6_20

head(data_001_1e6_25)
eigen_001_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_25.eigenval")
head(eigen_001_1e6_25)

variance_001_1e6_25 <- eigen_001_1e6_25 / sum(eigen_001_1e6_25)
variance_001_1e6_25[1:25]*100

pc1_001_1e6_25 <- round(variance_001_1e6_25[1]*100, 2) # 11.69
pc2_001_1e6_25 <- round(variance_001_1e6_25[2]*100, 2) # 7.78

output_dir
file_name <- "03_PCA_001_1e6_25.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_001_1e6_25$PC1, data_001_1e6_25$PC2,
     xlab = paste0("PC1: ", pc1_001_1e6_25, "%"),
     ylab = paste0("PC2: ", pc2_001_1e6_25, "%"),
     main = "PCA plot (MAF 0.01, HWE 1e-6, PCA 25)",
     col = "blue")
dev.off()

var_001_1e6_25 <- variance_001_1e6_25[1:25] * 100
cum_001_1e6_25 <- cumsum(var_001_1e6_25)
cum_001_1e6_25

head(data_001_1e6_30)
eigen_001_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6_PCA_30.eigenval")
head(eigen_001_1e6_30)

variance_001_1e6_30 <- eigen_001_1e6_30 / sum(eigen_001_1e6_30)
variance_001_1e6_30[1:30]*100

pc1_001_1e6_30 <- round(variance_001_1e6_30[1]*100, 2) # 10.52
pc2_001_1e6_30 <- round(variance_001_1e6_30[2]*100, 2) # 7

output_dir
file_name <- "03_PCA_001_1e6_30.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_001_1e6_30$PC1, data_001_1e6_30$PC2,
     xlab = paste0("PC1: ", pc1_001_1e6_30, "%"),
     ylab = paste0("PC2: ", pc2_001_1e6_30, "%"),
     main = "PCA plot (MAF 0.01, HWE 1e-6, PCA 30)",
     col = "blue")
dev.off()

var_001_1e6_30 <- variance_001_1e6_30[1:30] * 100
cum_001_1e6_30 <- cumsum(var_001_1e6_30)
cum_001_1e6_30 

head(data_005_1e6_5)
eigen_005_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_5.eigenval")
head(eigen_005_1e6_5)

variance_005_1e6_5 <- eigen_005_1e6_5 / sum(eigen_005_1e6_5)
variance_005_1e6_5[1:5]*100

pc1_001_1e6_5 <- round(variance_005_1e6_5[1]*100, 2) # 29.82
pc2_001_1e6_5 <- round(variance_005_1e6_5[2]*100, 2) # 28.73

output_dir
file_name <- "03_PCA_005_1e6_5.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_005_1e6_5$PC1, data_005_1e6_5$PC2,
     xlab = paste0("PC1: ", pc1_001_1e6_5, "%"),
     ylab = paste0("PC2: ", pc2_001_1e6_5, "%"),
     main = "PCA plot (MAF 0.05, HWE 1e-6, PCA 5)",
     col = "blue")
dev.off()

var_005_1e6_5 <- variance_005_1e6_5[1:5] * 100
cum_005_1e6_5 <- cumsum(var_005_1e6_5)
cum_005_1e6_5

head(data_005_1e6_10)
eigen_005_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_10.eigenval")
head(eigen_005_1e6_10)

variance_005_1e6_10 <- eigen_005_1e6_10 / sum(eigen_005_1e6_10)
variance_005_1e6_10[1:10]*100

pc1_005_1e6_10 <- round(variance_005_1e6_10[1]*100, 2) # 20.62
pc2_005_1e6_10 <- round(variance_005_1e6_10[2]*100, 2) # 19.85

output_dir
file_name <- "03_PCA_005_1e6_10.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_005_1e6_10$PC1, data_005_1e6_10$PC2,
     xlab = paste0("PC1: ", pc1_005_1e6_10, "%"),
     ylab = paste0("PC2: ", pc2_005_1e6_10, "%"),
     main = "PCA plot (MAF 0.05, HWE 1e-6, PCA 10)",
     col = "blue")
dev.off()

var_005_1e6_10 <- variance_005_1e6_10[1:10] * 100
cum_005_1e6_10 <- cumsum(var_005_1e6_10)
cum_005_1e6_10

head(data_005_1e6_15)
eigen_005_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_15.eigenval")
head(eigen_005_1e6_15)

variance_005_1e6_15 <- eigen_005_1e6_15 / sum(eigen_005_1e6_15)
variance_005_1e6_15[1:15]*100

pc1_005_1e6_15 <- round(variance_005_1e6_15[1]*100, 2) # 17.31
pc2_005_1e6_15 <- round(variance_005_1e6_15[2]*100, 2) # 16.66

output_dir
file_name <- "03_PCA_005_1e6_15.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_005_1e6_15$PC1, data_005_1e6_15$PC2,
     xlab = paste0("PC1: ", pc1_005_1e6_15, "%"),
     ylab = paste0("PC2: ", pc2_005_1e6_15, "%"),
     main = "PCA plot (MAF 0.05, HWE 1e-6, PCA 15)",
     col = "blue")
dev.off()

var_005_1e6_15 <- variance_005_1e6_15[1:15] * 100
cum_005_1e6_15 <- cumsum(var_005_1e6_15)
cum_005_1e6_15

head(data_005_1e6_20)
eigen_005_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_20.eigenval")
head(eigen_005_1e6_20)

variance_005_1e6_20 <- eigen_005_1e6_20 / sum(eigen_005_1e6_20)
variance_005_1e6_20[1:20]*100

pc1_005_1e6_20 <- round(variance_005_1e6_20[1]*100, 2) # 15.72
pc2_005_1e6_20 <- round(variance_005_1e6_20[2]*100, 2) # 15.13

output_dir
file_name <- "03_PCA_005_1e6_20.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_005_1e6_20$PC1, data_005_1e6_20$PC2,
     xlab = paste0("PC1: ", pc1_005_1e6_20, "%"),
     ylab = paste0("PC2: ", pc2_005_1e6_20, "%"),
     main = "PCA plot (MAF 0.05, HWE 1e-6, PCA 20)",
     col = "blue")
dev.off()

var_005_1e6_20 <- variance_005_1e6_20[1:20] * 100
cum_005_1e6_20 <- cumsum(var_005_1e6_20)
cum_005_1e6_20

head(data_005_1e6_25)
eigen_005_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_25.eigenval")
head(eigen_005_1e6_25)

variance_005_1e6_25 <- eigen_005_1e6_25 / sum(eigen_005_1e6_25)
variance_005_1e6_25[1:25]*100

pc1_005_1e6_25 <- round(variance_005_1e6_25[1]*100, 2) # 14.59
pc2_005_1e6_25 <- round(variance_005_1e6_25[2]*100, 2) # 14.04

output_dir
file_name <- "03_PCA_005_1e6_25.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_005_1e6_25$PC1, data_005_1e6_25$PC2,
     xlab = paste0("PC1: ", pc1_005_1e6_25, "%"),
     ylab = paste0("PC2: ", pc2_005_1e6_25, "%"),
     main = "PCA plot (MAF 0.05, HWE 1e-6, PCA 25)",
     col = "blue")
dev.off()

var_005_1e6_25 <- variance_005_1e6_25[1:25] * 100
cum_005_1e6_25 <- cumsum(var_005_1e6_25)
cum_005_1e6_25

head(data_005_1e6_30)
eigen_005_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6_PCA_30.eigenval")
head(eigen_005_1e6_30)

variance_005_1e6_30 <- eigen_005_1e6_30 / sum(eigen_005_1e6_30)
variance_005_1e6_30[1:30]*100

pc1_005_1e6_30 <- round(variance_005_1e6_30[1]*100, 2) # 13.73
pc2_005_1e6_30 <- round(variance_005_1e6_30[2]*100, 2) # 13.22

output_dir
file_name <- "03_PCA_005_1e6_30.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_005_1e6_30$PC1, data_005_1e6_30$PC2,
     xlab = paste0("PC1: ", pc1_005_1e6_30, "%"),
     ylab = paste0("PC2: ", pc2_005_1e6_30, "%"),
     main = "PCA plot (MAF 0.05, HWE 1e-6, PCA 30)",
     col = "blue")
dev.off()

var_005_1e6_30 <- variance_005_1e6_30[1:30] * 100
cum_005_1e6_30 <- cumsum(var_005_1e6_30)
cum_005_1e6_30

head(data_01_1e6_5)
eigen_01_1e6_5 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_5.eigenval")
head(eigen_01_1e6_5)

variance_01_1e6_5 <- eigen_01_1e6_5 / sum(eigen_01_1e6_5)
variance_01_1e6_5[1:5]*100

pc1_01_1e6_5 <- round(variance_01_1e6_5[1]*100, 2) # 43.08
pc2_01_1e6_5 <- round(variance_01_1e6_5[2]*100, 2) # 18.43

output_dir
file_name <- "03_PCA_01_1e6_5.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_01_1e6_5$PC1, data_01_1e6_5$PC2,
     xlab = paste0("PC1: ", pc1_01_1e6_5, "%"),
     ylab = paste0("PC2: ", pc2_01_1e6_5, "%"),
     main = "PCA plot (MAF 0.1, HWE 1e-6, PCA 5)",
     col = "blue")
dev.off()

var_01_1e6_5 <- variance_01_1e6_5[1:5] * 100
cum_01_1e6_5 <- cumsum(var_01_1e6_5)
cum_01_1e6_5

head(data_01_1e6_10)
eigen_01_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_10.eigenval")
head(eigen_01_1e6_10)

variance_01_1e6_10 <- eigen_01_1e6_10 / sum(eigen_01_1e6_10)
variance_01_1e6_10[1:10]*100

pc1_01_1e6_10 <- round(variance_01_1e6_10[1]*100, 2) # 31.61
pc2_01_1e6_10 <- round(variance_01_1e6_10[2]*100, 2) # 13.52

output_dir
file_name <- "03_PCA_01_1e6_10.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_01_1e6_10$PC1, data_01_1e6_10$PC2,
     xlab = paste0("PC1: ", pc1_01_1e6_10, "%"),
     ylab = paste0("PC2: ", pc2_01_1e6_10, "%"),
     main = "PCA plot (MAF 0.1, HWE 1e-6, PCA 10)",
     col = "blue")
dev.off()

var_01_1e6_10 <- variance_01_1e6_10[1:10] * 100
cum_01_1e6_10 <- cumsum(var_01_1e6_10)
cum_01_1e6_10

head(data_01_1e6_15)
eigen_01_1e6_15 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_15.eigenval")
head(eigen_01_1e6_15)

variance_01_1e6_15 <- eigen_01_1e6_15 / sum(eigen_01_1e6_15)
variance_01_1e6_15[1:15]*100

pc1_01_1e6_15 <- round(variance_01_1e6_15[1]*100, 2) # 27.01
pc2_01_1e6_15 <- round(variance_01_1e6_15[2]*100, 2) # 11.55

output_dir
file_name <- "03_PCA_01_1e6_15.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_01_1e6_15$PC1, data_01_1e6_15$PC2,
     xlab = paste0("PC1: ", pc1_01_1e6_15, "%"),
     ylab = paste0("PC2: ", pc2_01_1e6_15, "%"),
     main = "PCA plot (MAF 0.1, HWE 1e-6, PCA 15)",
     col = "blue")
dev.off()

var_01_1e6_15 <- variance_01_1e6_15[1:15] * 100
cum_01_1e6_15 <- cumsum(var_01_1e6_15)
cum_01_1e6_15

head(data_01_1e6_20)
eigen_01_1e6_20 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_20.eigenval")
head(eigen_01_1e6_20)

variance_01_1e6_20 <- eigen_01_1e6_20 / sum(eigen_01_1e6_20)
variance_01_1e6_20[1:20]*100

pc1_01_1e6_20 <- round(variance_01_1e6_20[1]*100, 2) # 24.34
pc2_01_1e6_20 <- round(variance_01_1e6_20[2]*100, 2) # 10.41

output_dir
file_name <- "03_PCA_01_1e6_20.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_01_1e6_20$PC1, data_01_1e6_20$PC2,
     xlab = paste0("PC1: ", pc1_01_1e6_20, "%"),
     ylab = paste0("PC2: ", pc2_01_1e6_20, "%"),
     main = "PCA plot (MAF 0.1, HWE 1e-6, PCA 20)",
     col = "blue")
dev.off()

var_01_1e6_20 <- variance_01_1e6_20[1:20] * 100
cum_01_1e6_20 <- cumsum(var_01_1e6_20)
cum_01_1e6_20

head(data_01_1e6_25)
eigen_01_1e6_25 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_25.eigenval")
head(eigen_01_1e6_25)

variance_01_1e6_25 <- eigen_01_1e6_25 / sum(eigen_01_1e6_25)
variance_01_1e6_25[1:25]*100

pc1_01_1e6_25 <- round(variance_01_1e6_25[1]*100, 2) # 22.53
pc2_01_1e6_25 <- round(variance_01_1e6_25[2]*100, 2) # 9.64

output_dir
file_name <- "03_PCA_01_1e6_25.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_01_1e6_25$PC1, data_01_1e6_25$PC2,
     xlab = paste0("PC1: ", pc1_01_1e6_25, "%"),
     ylab = paste0("PC2: ", pc2_01_1e6_25, "%"),
     main = "PCA plot (MAF 0.1, HWE 1e-6, PCA 25)",
     col = "blue")
dev.off()

var_01_1e6_25 <- variance_01_1e6_25[1:25] * 100
cum_01_1e6_25 <- cumsum(var_01_1e6_25)
cum_01_1e6_25

head(data_01_1e6_30)
eigen_01_1e6_30 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6_PCA_30.eigenval")
head(eigen_01_1e6_30)

variance_01_1e6_30 <- eigen_01_1e6_30 / sum(eigen_01_1e6_30)
variance_01_1e6_30[1:30]*100

pc1_01_1e6_30 <- round(variance_01_1e6_30[1]*100, 2) # 21.22
pc2_01_1e6_30 <- round(variance_01_1e6_30[2]*100, 2) # 9.08

output_dir
file_name <- "03_PCA_01_1e6_30.png"
output_path <- file.path(output_dir, file_name)

png(filename = output_path, width = 800, height = 600, res = 60)
plot(data_01_1e6_30$PC1, data_01_1e6_30$PC2,
     xlab = paste0("PC1: ", pc1_01_1e6_30, "%"),
     ylab = paste0("PC2: ", pc2_01_1e6_30, "%"),
     main = "PCA plot (MAF 0.1, HWE 1e-6, PCA 30)",
     col = "blue")
dev.off()

var_01_1e6_30 <- variance_01_1e6_30[1:30] * 100
cum_01_1e6_30 <- cumsum(var_01_1e6_30)
cum_01_1e6_30

merge_001_1e6_5 <- merge(data_001_1e6_5, purebred_total, by="IID", all.x = TRUE)
merge_001_1e6_10 <- merge(data_001_1e6_10, purebred_total, by="IID", all.x = TRUE)
merge_001_1e6_15 <- merge(data_001_1e6_15, purebred_total, by="IID", all.x = TRUE)
merge_001_1e6_20 <- merge(data_001_1e6_20, purebred_total, by="IID", all.x = TRUE)
merge_001_1e6_25 <- merge(data_001_1e6_25, purebred_total, by="IID", all.x = TRUE)
merge_001_1e6_30 <- merge(data_001_1e6_30, purebred_total, by="IID", all.x = TRUE)
merge_005_1e6_5 <- merge(data_005_1e6_5, purebred_total, by="IID", all.x = TRUE)
merge_005_1e6_10 <- merge(data_005_1e6_10, purebred_total, by="IID", all.x = TRUE)
merge_005_1e6_15 <- merge(data_005_1e6_15, purebred_total, by="IID", all.x = TRUE)
merge_005_1e6_20 <- merge(data_005_1e6_20, purebred_total, by="IID", all.x = TRUE)
merge_005_1e6_25 <- merge(data_005_1e6_25, purebred_total, by="IID", all.x = TRUE)
merge_005_1e6_30 <- merge(data_005_1e6_30, purebred_total, by="IID", all.x = TRUE)
merge_01_1e6_5 <- merge(data_01_1e6_5, purebred_total, by="IID", all.x = TRUE)
merge_01_1e6_10 <- merge(data_01_1e6_10, purebred_total, by="IID", all.x = TRUE)
merge_01_1e6_15 <- merge(data_01_1e6_15, purebred_total, by="IID", all.x = TRUE)
merge_01_1e6_20 <- merge(data_01_1e6_20, purebred_total, by="IID", all.x = TRUE)
merge_01_1e6_25 <- merge(data_01_1e6_25, purebred_total, by="IID", all.x = TRUE)
merge_01_1e6_30 <- merge(data_01_1e6_30, purebred_total, by="IID", all.x = TRUE)

# merge_001_1e6_5
ggplot(merge_001_1e6_5, aes(x = PC1, y = PC2, color = Breed)) +
  geom_point() +
  scale_color_manual(values = scales::hue_pal()(231)) +
  theme_minimal() +
  theme(legend.position = "none")
# top_breeds <- names(sort(table(merge_001_1e6_5$Breed), decreasing = TRUE))[1:10]
# merge_001_1e6_5$Breed_highlight <- ifelse(merge_001_1e6_5$Breed %in% top_breeds,
#                                           merge_001_1e6_5$Breed, "Other")

# ggplot(merge_001_1e6_5, aes(x = PC1, y = PC2, color = Breed_highlight)) +
#   geom_point() +
#   scale_color_manual(values = c(scales::hue_pal()(10), "Other" = "gray")) +
#   theme_minimal()

ggplot(merge_001_1e6_5, aes(x = PC1, y = PC2, color = Weight_Class_5KGBin_at_HLES)) +
  geom_point() +
  theme_minimal()
merge_001_1e6_5$LifeStage_Class_at_HLES


table(merge_001_1e6_5$LifeStage_Class_at_HLES)
ggplot(merge_001_1e6_5, aes(x = PC1, y = PC2, color = LifeStage_Class_at_HLES)) +
  geom_point() +
  theme_minimal()


ggplot(merge_005_1e6_5, aes(x = PC1, y = PC2, color = LifeStage_Class_at_HLES)) +
  geom_point() +
  theme_minimal()

ggplot(merge_01_1e6_5, aes(x = PC1, y = PC2, color = LifeStage_Class_at_HLES)) +
  geom_point() +
  theme_minimal()

head(data_001_1e6_5)

save(
    data_001_1e6_5, data_001_1e6_10, data_001_1e6_15, data_001_1e6_20, data_001_1e6_25, data_001_1e6_30,
    data_005_1e6_5, data_005_1e6_10, data_005_1e6_15, data_005_1e6_20, data_005_1e6_25, data_005_1e6_30,
    data_01_1e6_5, data_01_1e6_10, data_01_1e6_15, data_01_1e6_20, data_01_1e6_25, data_01_1e6_30,
    purebred_total,
    file = "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Code/00_RDataFile/03_GWAS_Pure/03_PCA.rdata"
)
