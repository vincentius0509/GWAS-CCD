load("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Code/00_RDataFile/03_GWAS_Pure/03_PCA.rdata")

library(pacman)

p_load(ggplot2, dplyr, data.table)

fwrite(data_001_1e6_5, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_001_1e6_5.txt", sep = "\t")
fwrite(data_001_1e6_10, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_001_1e6_10.txt", sep = "\t")
fwrite(data_001_1e6_15, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_001_1e6_15.txt", sep = "\t")
fwrite(data_001_1e6_20, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_001_1e6_20.txt", sep = "\t")
fwrite(data_001_1e6_25, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_001_1e6_25.txt", sep = "\t")
fwrite(data_001_1e6_30, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_001_1e6_30.txt", sep = "\t")
fwrite(data_005_1e6_5, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_005_1e6_5.txt", sep = "\t")
fwrite(data_005_1e6_10, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_005_1e6_10.txt", sep = "\t")
fwrite(data_005_1e6_15, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_005_1e6_15.txt", sep = "\t")
fwrite(data_005_1e6_20, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_005_1e6_20.txt", sep = "\t")
fwrite(data_005_1e6_25, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_005_1e6_25.txt", sep = "\t")
fwrite(data_005_1e6_30, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_005_1e6_30.txt", sep = "\t")
fwrite(data_01_1e6_5, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_01_1e6_5.txt", sep = "\t")
fwrite(data_01_1e6_10, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_01_1e6_10.txt", sep = "\t")
fwrite(data_01_1e6_15, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_01_1e6_15.txt", sep = "\t")
fwrite(data_01_1e6_20, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_01_1e6_20.txt", sep = "\t")
fwrite(data_01_1e6_25, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_01_1e6_25.txt", sep = "\t")
fwrite(data_01_1e6_30, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_PCA_01_1e6_30.txt", sep = "\t")

# FID, IID 추출
data_001_1e6_5_ID <- data_001_1e6_5[, c("FID", "IID")]
data_001_1e6_10_ID <- data_001_1e6_10[, c("FID", "IID")]
data_001_1e6_15_ID <- data_001_1e6_15[, c("FID", "IID")]
data_001_1e6_20_ID <- data_001_1e6_20[, c("FID", "IID")]
data_001_1e6_25_ID <- data_001_1e6_25[, c("FID", "IID")]
data_001_1e6_30_ID <- data_001_1e6_30[, c("FID", "IID")]

data_005_1e6_5_ID <- data_005_1e6_5[, c("FID", "IID")]
data_005_1e6_10_ID <- data_005_1e6_10[, c("FID", "IID")]
data_005_1e6_15_ID <- data_005_1e6_15[, c("FID", "IID")]
data_005_1e6_20_ID <- data_005_1e6_20[, c("FID", "IID")]
data_005_1e6_25_ID <- data_005_1e6_25[, c("FID", "IID")]
data_005_1e6_30_ID <- data_005_1e6_30[, c("FID", "IID")]

data_01_1e6_5_ID <- data_01_1e6_5[, c("FID", "IID")]
data_01_1e6_10_ID <- data_01_1e6_10[, c("FID", "IID")]
data_01_1e6_15_ID <- data_01_1e6_15[, c("FID", "IID")]
data_01_1e6_20_ID <- data_01_1e6_20[, c("FID", "IID")]
data_01_1e6_25_ID <- data_01_1e6_25[, c("FID", "IID")]
data_01_1e6_30_ID <- data_01_1e6_30[, c("FID", "IID")]

is_identical <- identical(data_001_1e6_5_ID, data_001_1e6_10_ID) && identical(data_001_1e6_5_ID, data_001_1e6_15_ID) && identical(data_001_1e6_5_ID, data_001_1e6_20_ID) && identical(data_001_1e6_5_ID, data_001_1e6_25_ID) && identical(data_001_1e6_5_ID, data_001_1e6_30_ID)
is_identical <- is_identical && identical(data_005_1e6_5_ID, data_005_1e6_10_ID) && identical(data_005_1e6_5_ID, data_005_1e6_15_ID) && identical(data_005_1e6_5_ID, data_005_1e6_20_ID) && identical(data_005_1e6_5_ID, data_005_1e6_25_ID) && identical(data_005_1e6_5_ID, data_005_1e6_30_ID)
is_identical <- is_identical && identical(data_01_1e6_5_ID, data_01_1e6_10_ID) && identical(data_01_1e6_5_ID, data_01_1e6_15_ID) && identical(data_01_1e6_5_ID, data_01_1e6_20_ID) && identical(data_01_1e6_5_ID, data_01_1e6_25_ID) && identical(data_01_1e6_5_ID, data_01_1e6_30_ID)
if (is_identical){
  print("모든 파일의 FID와 IID가 일치합니다.")
} else {
  print("파일들 간의 FID와 IID가 일치하지 않습니다.")
}

dim(data_001_1e6_5_ID) # 3527 2
dim(data_005_1e6_5_ID) # 3513 2
dim(data_01_1e6_5_ID) # 3512  2

purebred_total <- read.csv("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/00_purebred_total.csv")  
head(purebred_total)
dim(purebred_total) # 3542  49
colnames(purebred_total)

filtered_cslb_001 <- purebred_total[purebred_total$IID %in% data_001_1e6_5_ID$IID,]
length(unique(filtered_cslb_001$IID)) # 3527
dim(filtered_cslb_001) # 3527 49

filtered_cslb_005 <-purebred_total[purebred_total$IID %in% data_005_1e6_5_ID$IID,]
length(unique(filtered_cslb_005$IID)) # 3513
dim(filtered_cslb_005) # 3513 49

filtered_cslb_01 <- purebred_total[purebred_total$IID %in% data_01_1e6_5_ID$IID,]
length(unique(filtered_cslb_01$IID)) # 3512
dim(filtered_cslb_01) # 3512  49

data_001_1e6_5_unique <- unique(data_001_1e6_5_ID$IID)
length(data_001_1e6_5_unique)
data_005_1e6_5_unique <- unique(data_005_1e6_5_ID$IID)
length(unique(data_005_1e6_5_unique))
data_01_1e6_5_unique <- unique(data_01_1e6_5_ID$IID)
length(unique(data_01_1e6_5_unique))

# 포함된 IID
ids_001 <- unique(filtered_cslb_001$IID)
length(ids_001)
ids_005 <- unique(filtered_cslb_005$IID)
length(ids_005)
ids_01  <- unique(filtered_cslb_01$IID)
length(ids_01)

head(filtered_cslb_001)
pheno_001 <- filtered_cslb_001[,c("FID","IID","cslb_score")]
pheno_001$discrete <- ifelse(pheno_001$cslb_score >= 50, 2, 1)
colnames(pheno_001) <- c("FID", "IID", "phe1", "phe2")
head(pheno_001)
length(unique(pheno_001$IID)) # 3527
table(pheno_001$phe2) # 1: 3479 2: 48
fwrite(pheno_001, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_pheno001.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(filtered_cslb_005)
pheno_005 <- filtered_cslb_005[,c("FID","IID","cslb_score")]
pheno_005$discrete <- ifelse(pheno_005$cslb_score >= 50, 2, 1)
colnames(pheno_005) <- c("FID", "IID", "phe1", "phe2")
head(pheno_005)
length(unique(pheno_005$IID)) # 3513
table(pheno_005$phe2) # 1: 3468 2: 45
fwrite(pheno_005, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_pheno005.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(filtered_cslb_01)
pheno_01 <- filtered_cslb_01[,c("FID","IID","cslb_score")]
pheno_01$discrete <- ifelse(pheno_01$cslb_score >= 50, 2, 1)
colnames(pheno_01) <- c("FID", "IID", "phe1", "phe2")
head(pheno_01)
length(unique(pheno_01$IID)) # 3512
table(pheno_01$phe2) # 1: 3467 2: 45
fwrite(pheno_01, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_pheno01.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(purebred_total)
colnames(purebred_total)
dim(purebred_total) # 3542 49

purebred_cov <- purebred_total

# 성별(0: 암컷, 1: 수컷)
table(purebred_cov$Sex)
table(purebred_cov$Sex_edit)
purebred_cov$Sex_cat <- purebred_cov$Sex_edit
table(purebred_cov$Sex_cat)
purebred_cov$Sex_cov <- ifelse(purebred_cov$Sex_cat == "Female", 0, 1)
table(purebred_cov$Sex_cov)
sum(is.na(purebred_cov$Sex_cov)) # 0

# 나이
purebred_cov$Age_cov <- purebred_cov$Estimated_Age_Years_at_HLES
head(purebred_cov$Age_cov)
sum(is.na(purebred_cov$Age_cov)) # 0

# 중성화 여부 (0: Intact, 1: Neutered | Spayed")
table(purebred_cov$Intact)
purebred_cov$Intact_cov <- ifelse(purebred_cov$Intact == "Intact", 0, 1)
table(purebred_cov$Intact_cov)
sum(is.na(purebred_cov$Intact_cov)) # 0

# 활동 수준 (pa_score)
head(purebred_cov$pa_score)
# View(purebred_cov$pa_score)
purebred_cov$pa_score_cov <- purebred_cov$pa_score
sum(is.na(purebred_cov$pa_score_cov)) # 0

table(purebred_cov$Breed)
length(table(purebred_cov$Breed)) # 221

covariates <- purebred_cov[,c("IID", "Breed", "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")]
head(covariates)
dim(covariates) # 3542  6

breed_counts <- table(covariates$Breed)
breed_counts_sorted <- sort(breed_counts, decreasing = TRUE)
top_breed <- names(breed_counts_sorted)[1]
top_count <- breed_counts_sorted[1]
cat("Most common breed:", top_breed, "with", top_count, "occurrences\n") # Labrador Retriever

library(dplyr)
library(tidyr)

breed_dummy <- covariates %>%
  mutate(Breed = as.factor(Breed)) %>%
  select(IID, Breed) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Breed, values_from = value, values_fill = 0)


colnames(breed_dummy) <- gsub(" ", "_", colnames(breed_dummy))
head(breed_dummy)

covariates_final <- cbind(covariates, breed_dummy[ , -1])

dim(covariates_final)  # 3542  227
head(covariates_final)
table(covariates_final$Labrador_Retriever)

dim(purebred_cov[purebred_cov$Breed == "Labrador Retriever",])
colnames(covariates_final)
covariates_final$Breed <- NULL

head(data_001_1e6_5)
data_001_1e6_5_cov <- merge(data_001_1e6_5, covariates_final, by = "IID", all.x = TRUE)
head(data_001_1e6_5_cov)
setcolorder(data_001_1e6_5_cov, c("FID", "IID", paste0("PC", 1:5), setdiff(names(data_001_1e6_5_cov), c("FID", "IID", paste0("PC", 1:5)))))
head(data_001_1e6_5_cov)
dim(data_001_1e6_5_cov) # 3527  232
names(data_001_1e6_5_cov)[sapply(data_001_1e6_5_cov, function(x) all(x==0))]
table(data_001_1e6_5_cov$Scottish_Deerhound)
table(data_001_1e6_5_cov$English_Foxhound)
zero_cols_001 <- names(data_001_1e6_5_cov)[sapply(data_001_1e6_5_cov, function(x) all(x==0))]
data_001_1e6_5_cov_filtered <- data_001_1e6_5_cov[, !..zero_cols_001, with=FALSE]
print(zero_cols_001) # [1] "Scottish_Deerhound" "English_Foxhound"
print(dim(data_001_1e6_5_cov_filtered)) # 3527 230
fwrite(data_001_1e6_5_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_5.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_001_1e6_10)
data_001_1e6_10_cov <- merge(data_001_1e6_10, covariates_final, by = "IID", all.x = TRUE)
head(data_001_1e6_10_cov)
setcolorder(data_001_1e6_10_cov, c("FID", "IID", paste0("PC", 1:10), setdiff(names(data_001_1e6_10_cov), c("FID", "IID", paste0("PC", 1:10)))))
head(data_001_1e6_10_cov)
dim(data_001_1e6_10_cov) # 3527  237
names(data_001_1e6_10_cov)[sapply(data_001_1e6_10_cov, function(x) all(x==0))] # [1] "Scottish_Deerhound" "English_Foxhound"  
data_001_1e6_10_cov_filtered <- data_001_1e6_10_cov [, !..zero_cols_001, with=FALSE]
print(dim(data_001_1e6_10_cov_filtered)) # 3527 235
fwrite(data_001_1e6_10_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_10.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_001_1e6_15)
data_001_1e6_15_cov <- merge(data_001_1e6_15, covariates_final, by = "IID", all.x = TRUE)
head(data_001_1e6_15_cov)
setcolorder(data_001_1e6_15_cov, c("FID", "IID", paste0("PC", 1:15), setdiff(names(data_001_1e6_15_cov), c("FID", "IID", paste0("PC", 1:15)))))
head(data_001_1e6_15_cov)
dim(data_001_1e6_15_cov) # 3527 242
names(data_001_1e6_15_cov)[sapply(data_001_1e6_15_cov, function(x) all(x==0))] # [1] "Scottish_Deerhound" "English_Foxhound" 
data_001_1e6_15_cov_filtered <- data_001_1e6_15_cov[, !..zero_cols_001, with=FALSE]
print(dim(data_001_1e6_15_cov_filtered)) # 3527 240
fwrite(data_001_1e6_15_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_15.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_001_1e6_20)
data_001_1e6_20_cov <- merge(data_001_1e6_20, covariates_final, by = "IID", all.x = TRUE)
head(data_001_1e6_20_cov)
setcolorder(data_001_1e6_20_cov, c("FID", "IID", paste0("PC", 1:20), setdiff(names(data_001_1e6_20_cov), c("FID", "IID", paste0("PC", 1:20)))))
head(data_001_1e6_20_cov)
dim(data_001_1e6_20_cov) # 3527 247
names(data_001_1e6_20_cov)[sapply(data_001_1e6_20_cov, function(x) all(x==0))] # [1] "Scottish_Deerhound" "English_Foxhound" 
data_001_1e6_20_cov_filtered <- data_001_1e6_20_cov[, !..zero_cols_001, with=FALSE]
print(dim(data_001_1e6_20_cov_filtered)) # 3527 245
fwrite(data_001_1e6_20_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_20.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_001_1e6_25)
data_001_1e6_25_cov <- merge(data_001_1e6_25, covariates_final, by = "IID", all.x = TRUE)
head(data_001_1e6_25_cov)
setcolorder(data_001_1e6_25_cov, c("FID", "IID", paste0("PC", 1:25), setdiff(names(data_001_1e6_25_cov), c("FID", "IID", paste0("PC", 1:25)))))
head(data_001_1e6_25_cov)
dim(data_001_1e6_25_cov) # 3527 252
names(data_001_1e6_25_cov)[sapply(data_001_1e6_25_cov, function(x) all(x==0))] # [1] "Scottish_Deerhound" "English_Foxhound"
data_001_1e6_25_cov_filtered <- data_001_1e6_25_cov[, !..zero_cols_001, with=FALSE]
print(dim(data_001_1e6_25_cov_filtered)) # 3527 250
fwrite(data_001_1e6_25_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_25.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_001_1e6_30)
data_001_1e6_30_cov <- merge(data_001_1e6_30, covariates_final, by = "IID", all.x = TRUE)
head(data_001_1e6_30_cov)
setcolorder(data_001_1e6_30_cov, c("FID", "IID", paste0("PC", 1:30), setdiff(names(data_001_1e6_30_cov), c("FID", "IID", paste0("PC", 1:30)))))
head(data_001_1e6_30_cov)
dim(data_001_1e6_30_cov) # 3527 257
names(data_001_1e6_30_cov)[sapply(data_001_1e6_30_cov, function(x) all(x==0))] # [1] "Scottish_Deerhound" "English_Foxhound"
data_001_1e6_30_cov_filtered <- data_001_1e6_30_cov[, !..zero_cols_001, with=FALSE]
print(dim(data_001_1e6_30_cov_filtered)) # 3527 255
fwrite(data_001_1e6_30_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_30.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_005_1e6_5)
data_005_1e6_5_cov <- merge(data_005_1e6_5, covariates_final, by = "IID", all.x = TRUE)
head(data_005_1e6_5_cov)
setcolorder(data_005_1e6_5_cov, c("FID", "IID", paste0("PC", 1:5), setdiff(names(data_005_1e6_5_cov), c("FID", "IID", paste0("PC", 1:5)))))
head(data_005_1e6_5_cov)
dim(data_005_1e6_5_cov) # 3513  232
names(data_005_1e6_5_cov)[sapply(data_005_1e6_5_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
zero_cols_005 <- names(data_005_1e6_5_cov)[sapply(data_005_1e6_5_cov, function(x) all(x==0))] 
data_005_1e6_5_cov_filtered <- data_005_1e6_5_cov[, !..zero_cols_005, with=FALSE]
print(dim(data_005_1e6_5_cov_filtered)) # 3513  231
fwrite(data_005_1e6_5_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_5.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_005_1e6_10)
data_005_1e6_10_cov <- merge(data_005_1e6_10, covariates_final, by = "IID", all.x = TRUE)
head(data_005_1e6_10_cov)
setcolorder(data_005_1e6_10_cov, c("FID", "IID", paste0("PC", 1:10), setdiff(names(data_005_1e6_10_cov), c("FID", "IID", paste0("PC", 1:10)))))
head(data_005_1e6_10_cov)
dim(data_005_1e6_10_cov) # 3513 237
names(data_005_1e6_10_cov)[sapply(data_005_1e6_10_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
data_005_1e6_10_cov_filtered <- data_005_1e6_10_cov[, !..zero_cols_005, with=FALSE]
print(dim(data_005_1e6_10_cov_filtered)) # 3513 236
fwrite(data_005_1e6_10_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_10.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_005_1e6_15)
data_005_1e6_15_cov <- merge(data_005_1e6_15, covariates_final, by = "IID", all.x = TRUE)
head(data_005_1e6_15_cov)
setcolorder(data_005_1e6_15_cov, c("FID", "IID", paste0("PC", 1:15), setdiff(names(data_005_1e6_15_cov), c("FID", "IID", paste0("PC", 1:15)))))
head(data_005_1e6_15_cov)
dim(data_005_1e6_15_cov) # 3513 242
names(data_005_1e6_15_cov)[sapply(data_005_1e6_15_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
data_005_1e6_15_cov_filtered <- data_005_1e6_15_cov[, !..zero_cols_005, with=FALSE]
print(dim(data_005_1e6_15_cov_filtered)) # 3513 241
fwrite(data_005_1e6_15_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_15.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_005_1e6_20)
data_005_1e6_20_cov <- merge(data_005_1e6_20, covariates_final, by = "IID", all.x = TRUE)
head(data_005_1e6_20_cov)
setcolorder(data_005_1e6_20_cov, c("FID", "IID", paste0("PC", 1:20), setdiff(names(data_005_1e6_20_cov), c("FID", "IID", paste0("PC", 1:20)))))
head(data_005_1e6_20_cov)
dim(data_005_1e6_20_cov) # 3513  247
names(data_005_1e6_20_cov)[sapply(data_005_1e6_20_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
data_005_1e6_20_cov_filtered <- data_005_1e6_20_cov[, !..zero_cols_005, with=FALSE]
print(dim(data_005_1e6_20_cov_filtered)) # 3513 246
fwrite(data_005_1e6_20_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_20.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_005_1e6_25)
data_005_1e6_25_cov <- merge(data_005_1e6_25, covariates_final, by = "IID", all.x = TRUE)
head(data_005_1e6_25_cov)
setcolorder(data_005_1e6_25_cov, c("FID", "IID", paste0("PC", 1:25), setdiff(names(data_005_1e6_25_cov), c("FID", "IID", paste0("PC", 1:25)))))
head(data_005_1e6_25_cov)
dim(data_005_1e6_25_cov) # 3513 252
names(data_005_1e6_25_cov)[sapply(data_005_1e6_25_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
data_005_1e6_25_cov_filtered <- data_005_1e6_25_cov[, !..zero_cols_005, with=FALSE]
print(dim(data_005_1e6_25_cov_filtered)) # 3513  251
fwrite(data_005_1e6_25_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_25.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_005_1e6_30)
data_005_1e6_30_cov <- merge(data_005_1e6_30, covariates_final, by = "IID", all.x = TRUE)
head(data_005_1e6_30_cov)
setcolorder(data_005_1e6_30_cov, c("FID", "IID", paste0("PC", 1:30), setdiff(names(data_005_1e6_30_cov), c("FID", "IID", paste0("PC", 1:30)))))
head(data_005_1e6_30_cov)
dim(data_005_1e6_30_cov) # 3513  257
names(data_005_1e6_30_cov)[sapply(data_005_1e6_30_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
data_005_1e6_30_cov_filtered <- data_005_1e6_30_cov[, !..zero_cols_005, with=FALSE]
print(dim(data_005_1e6_30_cov_filtered)) # 3513 256
fwrite(data_005_1e6_30_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_30.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_01_1e6_5)
data_01_1e6_5_cov <- merge(data_01_1e6_5, covariates_final, by = "IID", all.x = TRUE)
head(data_01_1e6_5_cov)
setcolorder(data_01_1e6_5_cov, c("FID", "IID", paste0("PC", 1:5), setdiff(names(data_01_1e6_5_cov), c("FID", "IID", paste0("PC", 1:5)))))
head(data_01_1e6_5_cov)
dim(data_01_1e6_5_cov) # 3512 232
names(data_01_1e6_5_cov)[sapply(data_01_1e6_5_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
zero_cols_01 <- names(data_01_1e6_5_cov)[sapply(data_01_1e6_5_cov, function(x) all(x==0))]
zero_cols_01
data_01_1e6_5_cov_filtered <- data_01_1e6_5_cov[, !..zero_cols_01, with=FALSE]
dim(data_01_1e6_5_cov_filtered) # 3512  231
fwrite(data_01_1e6_5_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_5.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_01_1e6_10)
data_01_1e6_10_cov <- merge(data_01_1e6_10, covariates_final, by = "IID", all.x = TRUE)
head(data_01_1e6_10_cov)
setcolorder(data_01_1e6_10_cov, c("FID", "IID", paste0("PC", 1:10), setdiff(names(data_01_1e6_10_cov), c("FID", "IID", paste0("PC", 1:10)))))
head(data_01_1e6_10_cov)
dim(data_01_1e6_10_cov) # 3512  237
names(data_01_1e6_10_cov)[sapply(data_01_1e6_10_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
data_01_1e6_10_cov_filtered <- data_01_1e6_10_cov[, !..zero_cols_01, with=FALSE]
dim(data_01_1e6_10_cov_filtered) # 3512 236
fwrite(data_01_1e6_10_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_10.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_01_1e6_15)
data_01_1e6_15_cov <- merge(data_01_1e6_15, covariates_final, by = "IID", all.x = TRUE)
head(data_01_1e6_15_cov)
setcolorder(data_01_1e6_15_cov, c("FID", "IID", paste0("PC", 1:15), setdiff(names(data_01_1e6_15_cov), c("FID", "IID", paste0("PC", 1:15)))))
head(data_01_1e6_15_cov)
dim(data_01_1e6_15_cov) # 3512  242
names(data_01_1e6_15_cov)[sapply(data_01_1e6_15_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
data_01_1e6_15_cov_filtered <- data_01_1e6_15_cov[, !..zero_cols_01, with=FALSE]
dim(data_01_1e6_15_cov_filtered) # 3512 241
fwrite(data_01_1e6_15_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_15.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_01_1e6_20)
data_01_1e6_20_cov <- merge(data_01_1e6_20, covariates_final, by = "IID", all.x = TRUE)
head(data_01_1e6_20_cov)
setcolorder(data_01_1e6_20_cov, c("FID", "IID", paste0("PC", 1:20), setdiff(names(data_01_1e6_20_cov), c("FID", "IID", paste0("PC", 1:20)))))
head(data_01_1e6_20_cov)
dim(data_01_1e6_20_cov) # 3512  247
names(data_01_1e6_20_cov)[sapply(data_01_1e6_20_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
data_01_1e6_20_cov_filtered <- data_01_1e6_20_cov[, !..zero_cols_01, with=FALSE]
dim(data_01_1e6_20_cov_filtered) # 3512 246
fwrite(data_01_1e6_20_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_20.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_01_1e6_25)
data_01_1e6_25_cov <- merge(data_01_1e6_25, covariates_final, by = "IID", all.x = TRUE)
head(data_01_1e6_25_cov)
setcolorder(data_01_1e6_25_cov, c("FID", "IID", paste0("PC", 1:25), setdiff(names(data_01_1e6_25_cov), c("FID", "IID", paste0("PC", 1:25)))))
head(data_01_1e6_25_cov)
dim(data_01_1e6_25_cov) # 3512  252
names(data_01_1e6_25_cov)[sapply(data_01_1e6_25_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
data_01_1e6_25_cov_filtered <- data_01_1e6_25_cov[, !..zero_cols_01, with=FALSE]
dim(data_01_1e6_25_cov_filtered) # 3512 251
fwrite(data_01_1e6_25_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_25.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

head(data_01_1e6_30)
data_01_1e6_30_cov <- merge(data_01_1e6_30, covariates_final, by = "IID", all.x = TRUE)
head(data_01_1e6_30_cov)
setcolorder(data_01_1e6_30_cov, c("FID", "IID", paste0("PC", 1:30), setdiff(names(data_01_1e6_30_cov), c("FID", "IID", paste0("PC", 1:30)))))
head(data_01_1e6_30_cov)
dim(data_01_1e6_30_cov) # 3512  257
names(data_01_1e6_30_cov)[sapply(data_01_1e6_30_cov, function(x) all(x==0))] # [1] "Japanese_Chin"
data_01_1e6_30_cov_filtered <- data_01_1e6_30_cov[, !..zero_cols_01, with=FALSE]
dim(data_01_1e6_30_cov_filtered) # 3512 256
fwrite(data_01_1e6_30_cov_filtered, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_30.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

excluded_breeds_001 <- c(
  "Puli","Finnish_Lapphund","Australian_Labradoodle","Prague_Ratter","Bavarian_Mountain_Scent_Hound",
  "Glen_of_Imaal_Terrier","French_Spaniel","Broholmer","Irish_Red_and_White_Setter","Black_Mouth_Cur",
  "Taiwan_Dog","Portuguese_Podengo_Pequeno","Treeing_Tennessee_Brindle","Perro_de_Presa_Canario",
  "Komondor","Labradoodle","Romanian_Mioritic_Shepherd_Dog","German_Spitz","Akbash","Miki",
  "Dandie_Dinmont_Terrier","Japanese_Chin","Barbet","Basset_Fauve_de_Bretagne","Xoloitzcuintli",
  "Skye_Terrier","Greater_Swiss_Mountain_Dog","Cesky_Terrier","Belgian_Laekenois","Norwegian_Buhund",
  "American_Foxhound","American_Bulldog","Hokkaido","Pharaoh_Hound","Windsprite","Persian_Shepherd",
  "Cobberdog","Kromfohrlander","Working_Kelpie","Mudi","Dogue_de_Bordeaux","Saluki","Kai_Ken",
  "English_Toy_Spaniel","White_Swiss_Shepherd","Lowchen"
)
length(excluded_breeds_001) # 46

excluded_breeds_005 <- c(
  "Puli","Finnish_Lapphund","Australian_Labradoodle","Prague_Ratter","Bavarian_Mountain_Scent_Hound",
  "Glen_of_Imaal_Terrier","Neapolitan_Mastiff","French_Spaniel","Broholmer","Alaskan_Klee_Kai",
  "Irish_Red_and_White_Setter","Black_Mouth_Cur",  "Taiwan_Dog","Portuguese_Podengo_Pequeno","Treeing_Tennessee_Brindle","Perro_de_Presa_Canario",
  "Komondor","American_Water_Spaniel","Labradoodle","Romanian_Mioritic_Shepherd_Dog","German_Spitz","Karelian_Bear_Dog","Akbash","Miki",
  "Dandie_Dinmont_Terrier","Barbet","Basset_Fauve_de_Bretagne","Xoloitzcuintli",
  "Skye_Terrier","Greater_Swiss_Mountain_Dog","Cesky_Terrier","Belgian_Laekenois","Norwegian_Buhund","English_Foxhound",
  "American_Foxhound","American_Bulldog","Hokkaido","Pharaoh_Hound","Windsprite","Persian_Shepherd",
  "Cobberdog","Kromfohrlander","Working_Kelpie","Mudi","Dogue_de_Bordeaux","Saluki","Kai_Ken",
  "English_Toy_Spaniel","White_Swiss_Shepherd","Lowchen"  
)
length(excluded_breeds_005) # 50

excluded_breeds_01 <- c(
  "Puli","Finnish_Lapphund","Australian_Labradoodle","Prague_Ratter","Bavarian_Mountain_Scent_Hound",
  "Glen_of_Imaal_Terrier","Neapolitan_Mastiff","French_Spaniel","Broholmer","Alaskan_Klee_Kai",
  "Irish_Red_and_White_Setter","Black_Mouth_Cur",  "Taiwan_Dog","Portuguese_Podengo_Pequeno","Treeing_Tennessee_Brindle","Perro_de_Presa_Canario",
  "Komondor","American_Water_Spaniel","Labradoodle","Romanian_Mioritic_Shepherd_Dog","German_Spitz","Karelian_Bear_Dog","Akbash","Miki",
  "Dandie_Dinmont_Terrier","Barbet","Basset_Fauve_de_Bretagne","Xoloitzcuintli",
  "Skye_Terrier","Greater_Swiss_Mountain_Dog","Cesky_Terrier","Belgian_Laekenois","Norwegian_Buhund","English_Foxhound",
  "American_Foxhound","American_Bulldog","Hokkaido","Pharaoh_Hound","Windsprite","Persian_Shepherd",
  "Cobberdog","Kromfohrlander","Working_Kelpie","Mudi","Dogue_de_Bordeaux","Saluki","Kai_Ken",
  "English_Toy_Spaniel","White_Swiss_Shepherd","Lowchen"  
)
length(excluded_breeds_01) # 50

##### 001 => 데이터가 1마리 있는 품종들에서 문제가 발생
cols_check <- intersect(excluded_breeds_001,colnames(data_001_1e6_5_cov_filtered))
length(cols_check)

breeds_1_counts <- sapply(cols_check, function(breed) {
  sum(data_001_1e6_5_cov_filtered[[breed]] == 1, na.rm = T)
})

sum(breeds_1_counts)
names(breeds_1_counts[breeds_1_counts == 1]) # 46
breeds_1_counts[breeds_1_counts != 1] # 0

breed_columns <- colnames(data_001_1e6_5_cov_filtered)[
  sapply(data_001_1e6_5_cov_filtered, function(x) all(x %in% c(0,1)))
]

breed_counts <- sapply(breed_columns, function(breed) {
  sum(data_001_1e6_5_cov_filtered[[breed]] == 1, na.rm = TRUE)
})

# 1마리인 품종만 필터링
breeds_with_1_individual <- names(breed_counts[breed_counts == 1])

# 결과 확인
breeds_with_1_individual
length(breeds_with_1_individual)  # 몇 개인지도 같이 확인

##### 005 => 1마리인 품종에서만 문제 발생
cols_check <- intersect(excluded_breeds_005, colnames(data_005_1e6_5_cov_filtered))
length(cols_check)

breeds_1_counts <- sapply(cols_check, function(breed) {
  sum(data_005_1e6_5_cov_filtered[[breed]]==1, na.rm = T)
})

sum(breeds_1_counts)
names(breeds_1_counts[breeds_1_counts == 1]) # 50
breeds_1_counts[breeds_1_counts != 1] # 0

breed_columns <- colnames(data_005_1e6_5_cov_filtered)[
  sapply(data_005_1e6_5_cov_filtered, function(x) all(x %in% c(0,1)))
]

breed_counts <- sapply(breed_columns, function(breed) {
  sum(data_005_1e6_5_cov_filtered[[breed]] == 1, na.rm = TRUE)
})

# 1마리인 품종만 필터링
breeds_with_1_individual <- names(breed_counts[breed_counts == 1])

# 결과 확인
breeds_with_1_individual
length(breeds_with_1_individual)  # 몇 개인지도 같이 확인

##### 01 => 한 마리인 데이터들에서 문제 발생
cols_check <- intersect(excluded_breeds_01, colnames(data_01_1e6_5_cov_filtered))
length(cols_check)

breeds_1_counts <- sapply(cols_check, function(breed) {
  sum(data_01_1e6_5_cov_filtered[[breed]]==1, na.rm = T)
})

sum(breeds_1_counts)
names(breeds_1_counts[breeds_1_counts == 1]) # 50
breeds_1_counts[breeds_1_counts != 1] # 0

breed_columns <- colnames(data_01_1e6_5_cov_filtered)[
  sapply(data_01_1e6_5_cov_filtered, function(x) all(x %in% c(0,1)))
]

breed_counts <- sapply(breed_columns, function(breed) {
  sum(data_01_1e6_5_cov_filtered[[breed]] == 1, na.rm = TRUE)
})

# 1마리인 품종만 필터링
breeds_with_1_individual <- names(breed_counts[breed_counts == 1])

# 결과 확인
breeds_with_1_individual
length(breeds_with_1_individual)  # 몇 개인지도 같이 확인

# 각 covariate 별로 1마리인 품종들은 제거

dim(data_001_1e6_5_cov_filtered) # 3527 230
exclude_columns <- intersect(excluded_breeds_001, names(data_001_1e6_5_cov_filtered))
length(exclude_columns) # 46
data_001_1e6_5_cov_filtered2 <- data_001_1e6_5_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_001_1e6_5_cov_filtered2) # 3527  184
fwrite(data_001_1e6_5_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_5.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_001_1e6_10_cov_filtered) # 3527  235
exclude_columns <- intersect(excluded_breeds_001, names(data_001_1e6_10_cov_filtered))
length(exclude_columns) # 46
data_001_1e6_10_cov_filtered2 <- data_001_1e6_10_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_001_1e6_10_cov_filtered2) # 3527 189
fwrite(data_001_1e6_10_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_10.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_001_1e6_15_cov_filtered) # 3527  240
exclude_columns <- intersect(excluded_breeds_001, names(data_001_1e6_15_cov_filtered))
length(exclude_columns) # 46
data_001_1e6_15_cov_filtered2 <- data_001_1e6_15_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_001_1e6_15_cov_filtered2) # 3527 194
fwrite(data_001_1e6_15_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_15.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_001_1e6_20_cov_filtered) # 3527  245
exclude_columns <- intersect(excluded_breeds_001, names(data_001_1e6_20_cov_filtered))
length(exclude_columns) # 46
data_001_1e6_20_cov_filtered2 <- data_001_1e6_20_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_001_1e6_20_cov_filtered2) # 3527 199
fwrite(data_001_1e6_20_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_20.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_001_1e6_25_cov_filtered) # 3527 250
exclude_columns <- intersect(excluded_breeds_001, names(data_001_1e6_25_cov_filtered))
length(exclude_columns) # 46
data_001_1e6_25_cov_filtered2 <- data_001_1e6_25_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_001_1e6_25_cov_filtered2) # 3527 204
fwrite(data_001_1e6_25_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_25.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_001_1e6_30_cov_filtered) # 3527  255
exclude_columns <- intersect(excluded_breeds_001, names(data_001_1e6_30_cov_filtered))
length(exclude_columns) # 46
data_001_1e6_30_cov_filtered2 <- data_001_1e6_30_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_001_1e6_30_cov_filtered2) # 3527 209
fwrite(data_001_1e6_30_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_30.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_005_1e6_5_cov_filtered) # 3513 231
exclude_columns <- intersect(excluded_breeds_005, names(data_005_1e6_5_cov_filtered))
length(exclude_columns) # 50
data_005_1e6_5_cov_filtered2 <- data_005_1e6_5_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_005_1e6_5_cov_filtered2) # 3513  181
fwrite(data_005_1e6_5_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_5.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_005_1e6_10_cov_filtered) # 3513 236
exclude_columns <- intersect(excluded_breeds_005, names(data_005_1e6_10_cov_filtered))
length(exclude_columns) # 50
data_005_1e6_10_cov_filtered2 <- data_005_1e6_10_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_005_1e6_10_cov_filtered2) # 3513 186
fwrite(data_005_1e6_10_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_10.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_005_1e6_15_cov_filtered) # 3513  241
exclude_columns <- intersect(excluded_breeds_005, names(data_005_1e6_15_cov_filtered))
length(exclude_columns) # 50
data_005_1e6_15_cov_filtered2 <- data_005_1e6_15_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_005_1e6_15_cov_filtered2) # 3513 191
fwrite(data_005_1e6_15_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_15.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_005_1e6_20_cov_filtered) # 3513  246
exclude_columns <- intersect(excluded_breeds_005, names(data_005_1e6_20_cov_filtered))
length(exclude_columns) # 50
data_005_1e6_20_cov_filtered2 <- data_005_1e6_20_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_005_1e6_20_cov_filtered2) # 3513 196
fwrite(data_005_1e6_20_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_20.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_005_1e6_25_cov_filtered) # 3513  251
exclude_columns <- intersect(excluded_breeds_005, names(data_005_1e6_25_cov_filtered))
length(exclude_columns) # 50
data_005_1e6_25_cov_filtered2 <- data_005_1e6_25_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_005_1e6_25_cov_filtered2) # 3513 201
fwrite(data_005_1e6_25_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_25.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_005_1e6_30_cov_filtered) # 3513  256
exclude_columns <- intersect(excluded_breeds_005, names(data_005_1e6_30_cov_filtered))
length(exclude_columns) # 50
data_005_1e6_30_cov_filtered2 <- data_005_1e6_30_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_005_1e6_30_cov_filtered2) # 3513 206
fwrite(data_005_1e6_30_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_005_1e6_30.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_01_1e6_5_cov_filtered) # 3512  231
exclude_columns <- intersect(excluded_breeds_01, names(data_01_1e6_5_cov_filtered))
length(exclude_columns) # 50
data_01_1e6_5_cov_filtered2 <- data_01_1e6_5_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_01_1e6_5_cov_filtered2) # 3512 181
fwrite(data_01_1e6_5_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_5.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_01_1e6_10_cov_filtered) # 3512 236
exclude_columns <- intersect(excluded_breeds_01, names(data_01_1e6_10_cov_filtered))
length(exclude_columns) # 50
data_01_1e6_10_cov_filtered2 <- data_01_1e6_10_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_01_1e6_10_cov_filtered2) # 3512  186
fwrite(data_01_1e6_10_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_10.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_01_1e6_15_cov_filtered) # 3512 241
exclude_columns <- intersect(excluded_breeds_01, names(data_01_1e6_15_cov_filtered))
length(exclude_columns) # 50
data_01_1e6_15_cov_filtered2 <- data_01_1e6_15_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_01_1e6_15_cov_filtered2) # 3512  191
fwrite(data_01_1e6_15_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_15.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_01_1e6_20_cov_filtered) # 3512 246
exclude_columns <- intersect(excluded_breeds_01, names(data_01_1e6_20_cov_filtered))
length(exclude_columns) # 50
data_01_1e6_20_cov_filtered2 <- data_01_1e6_20_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_01_1e6_20_cov_filtered2) # 3512  196
fwrite(data_01_1e6_20_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_20.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_01_1e6_25_cov_filtered) # 3512 251
exclude_columns <- intersect(excluded_breeds_01, names(data_01_1e6_25_cov_filtered))
length(exclude_columns) # 50
data_01_1e6_25_cov_filtered2 <- data_01_1e6_25_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_01_1e6_25_cov_filtered2) # 3512  201
fwrite(data_01_1e6_25_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_25.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

dim(data_01_1e6_30_cov_filtered) # 3512 256
exclude_columns <- intersect(excluded_breeds_01, names(data_01_1e6_30_cov_filtered))
length(exclude_columns) # 50
data_01_1e6_30_cov_filtered2 <- data_01_1e6_30_cov_filtered[, !..exclude_columns, with=FALSE]
dim(data_01_1e6_30_cov_filtered2) # 3512  206
fwrite(data_01_1e6_30_cov_filtered2, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_01_1e6_30.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)


save(
  data_001_1e6_5_cov_filtered2, data_001_1e6_10_cov_filtered2, data_001_1e6_15_cov_filtered2, 
  data_001_1e6_20_cov_filtered2, data_001_1e6_25_cov_filtered2, data_001_1e6_30_cov_filtered2,
  data_005_1e6_5_cov_filtered2, data_005_1e6_10_cov_filtered2, data_005_1e6_15_cov_filtered2,
  data_005_1e6_20_cov_filtered2, data_005_1e6_25_cov_filtered2, data_005_1e6_30_cov_filtered2,
  data_01_1e6_5_cov_filtered2, data_01_1e6_10_cov_filtered2, data_01_1e6_15_cov_filtered2,
  data_01_1e6_20_cov_filtered2, data_01_1e6_25_cov_filtered2, data_01_1e6_30_cov_filtered2,
  file = "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Code/00_RDataFile/03_GWAS_Pure/03-2_Add_Phenotype.rdata"
)
