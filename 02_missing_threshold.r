getwd()
library(ggplot2)
library(dplyr)

# 01_MAF.sh 실행
# 02-1_Missing_Rate.sh 실행
# MAF_LIST=(0.01 0.05 0.1)

# # 1. Threshold:  0.001
# imiss_data <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.001.imiss", header = TRUE)
# summary(imiss_data)

# ggplot(imiss_data, aes(x = F_MISS)) +
#   geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
#   geom_vline(xintercept = 0.2639, color = "red", linetype = "dashed") +
#   labs(title = "Distribution of Missing Data Rates in Individuals (imiss file) [MAF = 0.001]",
#        x = "Missing Rate",
#        y = "Number of Individuals") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 20),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 12)
#   )
# dev.off()
# # --mind: 0.280

# lmiss_data <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.001.lmiss", header = TRUE)
# summary(lmiss_data)
# ggplot(lmiss_data, aes(x = F_MISS)) +
#   geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
#   geom_vline(xintercept = 0.68135, color = "red", linetype = "dashed") +
#   labs(title = "Distribution of Missing Data Rates in SNPs (lmiss file) [MAF = 0.001]",
#        x = "Missing Rate",
#        y = "Number of SNPs") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 20),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 12)
#   )
# dev.off()
# # --geno: 0.68135

# # 2. Threshold: 0.005
# imiss_data <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.005.imiss", header = TRUE)
# summary(imiss_data)
# ggplot(imiss_data, aes(x = F_MISS)) +
#   geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
#   geom_vline(xintercept = 0.2396, color = "red", linetype = "dashed") +
#   labs(title = "Distribution of Missing Data Rates in Individuals (imiss file) [MAF = 0.005]",
#        x = "Missing Rate",
#        y = "Number of Individuals") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 20),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 12)
#   )
# dev.off()
# # --mind: 0.2396

# lmiss_data <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.005.lmiss", header = TRUE)
# summary(lmiss_data)
# ggplot(lmiss_data, aes(x = F_MISS)) +
#   geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
#   geom_vline(xintercept = 0.6379, color = "red", linetype = "dashed") +
#   labs(title = "Distribution of Missing Data Rates in SNPs (lmiss file) [MAF = 0.005]",
#        x = "Missing Rate",
#        y = "Number of SNPs") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 20),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 12)
#   )
# dev.off()
# # --geno: 0.6379

# 3. Threshold: 0.01
imiss_data <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.01.imiss", header = TRUE)
summary(imiss_data)
(quantile(imiss_data$F_MISS, 0.75) + max(imiss_data$F_MISS)) / 2 # 0.3911
ggplot(imiss_data, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  geom_vline(xintercept = 0.4, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in Individuals (imiss file) [MAF = 0.01]",
       x = "Missing Rate",
       y = "Number of Individuals") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()
# --mind: 0.400

lmiss_data <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.01.lmiss", header = TRUE)
summary(lmiss_data)
(quantile(lmiss_data$F_MISS, 0.75) + max(lmiss_data$F_MISS))/2 #0.84925
ggplot(lmiss_data, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  geom_vline(xintercept = 0.84925, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in SNPs (lmiss file) [MAF = 0.01]",
       x = "Missing Rate",
       y = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()
# --geno: 0.84925

# 4. Threshold: 0.05
imiss_data <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.05.imiss", header = TRUE)
summary(imiss_data)
(quantile(imiss_data$F_MISS, 0.75) + max(imiss_data$F_MISS))/2 # 0.1753
ggplot(imiss_data, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  geom_vline(xintercept = 0.1753, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in Individuals (imiss file) [MAF = 0.05]",
       x = "Missing Rate",
       y = "Number of Individuals") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()
# --mind: 0.1753

lmiss_data <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.05.lmiss", header = TRUE)
summary(lmiss_data)
(quantile(lmiss_data$F_MISS, 0.75) + max(lmiss_data$F_MISS))/2 # 0.5737
ggplot(lmiss_data, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  geom_vline(xintercept = 0.5737, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in SNPs (lmiss file) [MAF = 0.05]",
       x = "Missing Rate",
       y = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()
# --geno: 0.5737

# 5. Threshold: 0.1
imiss_data <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.1.imiss", header = TRUE)
summary(imiss_data)
(quantile(imiss_data$F_MISS, 0.75) + max(imiss_data$F_MISS))/2 # 0.19925
ggplot(imiss_data, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  geom_vline(xintercept = 0.19925, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in Individuals (imiss file) [MAF = 0.1]",
       x = "Missing Rate",
       y = "Number of Individuals") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()
# --mind: 0.19925

lmiss_data <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.1.lmiss", header = TRUE)
summary(lmiss_data)
(quantile(lmiss_data$F_MISS, 0.75) + max(lmiss_data$F_MISS))/2 # 0.5868
ggplot(lmiss_data, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  geom_vline(xintercept = 0.5868, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in SNPs (lmiss file) [MAF = 0.1]",
       x = "Missing Rate",
       y = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()
# --geno: 0.5868

# 02-2_Missing_Rate_Filtered.sh 실행
# 02-3_MR_Filtered_Check.sh 실행

# MAF Threshold: 0.01
imiss_data_filtered <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.01_Filtered_Check.imiss", header = TRUE)
ggplot(imiss_data_filtered, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "Purple", color = "black") +
  geom_vline(xintercept = NULL, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in Individuals (imiss file) [Filtered]",
       x = "Missing Rate",
       y = "Number of Individuals") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()

lmiss_data_filtered <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.01_Filtered_Check.lmiss", header = TRUE)
ggplot(lmiss_data_filtered, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "Purple", color = "black") +
  geom_vline(xintercept = NULL, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in SNPs (lmiss file) [Filtered]",
       x = "Missing Rate",
       y = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()

# MAF Threshold: 0.05
imiss_data_filtered <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.05_Filtered_Check.imiss", header = TRUE)
ggplot(imiss_data_filtered, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "Purple", color = "black") +
  geom_vline(xintercept = NULL, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in Individuals (imiss file) [Filtered]",
       x = "Missing Rate",
       y = "Number of Individuals") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()

lmiss_data_filtered <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.05_Filtered_Check.lmiss", header = TRUE)
ggplot(lmiss_data_filtered, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "Purple", color = "black") +
  geom_vline(xintercept = NULL, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in SNPs (lmiss file) [Filtered]",
       x = "Missing Rate",
       y = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()

# MAF Threshold: 0.1
imiss_data_filtered <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.1_Filtered_Check.imiss", header = TRUE)
ggplot(imiss_data_filtered, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "Purple", color = "black") +
  geom_vline(xintercept = NULL, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in Individuals (imiss file) [Filtered]",
       x = "Missing Rate",
       y = "Number of Individuals") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()

lmiss_data_filtered <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_Missing_MAF_0.1_Filtered_Check.lmiss", header = TRUE) 
ggplot(lmiss_data_filtered, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "Purple", color = "black") +
  geom_vline(xintercept = NULL, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Missing Data Rates in SNPs (lmiss file) [Filtered]",
       x = "Missing Rate",
       y = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
dev.off()

# 3527
# 3513
# 3512