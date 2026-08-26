getwd()
setwd("/data/02_ByProject/DogAgingProject")

library(data.table)
library(ggplot2)
library(dplyr)

# 1. DAP_2022_DogOverview_v1.0.csv
dap_overview <- read.csv("/data/02_ByProject/DogAgingProject/DogOverview/v2022/DAP_2022_DogOverview_v1.0.csv")
names(dap_overview)
dim(dap_overview) 
head(dap_overview$Breed_Status)
head(dap_overview$Breed)

# 2. Genetic Data
dap_geno <- fread("/data/02_ByProject/DogAgingProject/GeneticData/v2022/DogAgingProject_GeneticData_CuratedRelease_2022.fam", header = FALSE)
head(dap_geno)

colnames(dap_geno) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")
head(dap_geno)
dim(dap_geno)

dap_genoOverview <- merge(dap_geno, dap_overview, by.x = "IID", by.y = "dog_id", all.x = TRUE)

dim(dap_genoOverview) # 7425   40

head(dap_genoOverview)
table(dap_genoOverview$Breed_Status)
table(dap_genoOverview$Breed)
length(table(dap_genoOverview$Breed))

colnames(dap_genoOverview)

table(dap_genoOverview$Breed_Status)

purebred_data <- dap_genoOverview %>% filter(Breed_Status == "Purebred")
head(purebred_data)
dim(purebred_data)

write.csv(purebred_data, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/00_Purebred_Data.csv", row.names = FALSE)

table(purebred_data$Breed)

# 품종 확인
length(table(purebred_data$Breed)) # 총 232개 품종, 3864마리
sum(is.na(purebred_data$Breed)) # 0

colnames(purebred_data)

# 성별 (1 = Female, 2 = Male)
table(purebred_data$Sex)
#    1    2 
# 1805 2059 

table(purebred_data$Sex_Class_at_HLES)
# Female, intact Female, spayed   Male, intact Male, neutered 
#            446           1360            666           1392 

# Female: 1806, Male: 2058

# View(purebred_data[,c("Sex", "Sex_Class_at_HLES")])

purebred_data2 <- purebred_data

purebred_data2$Sex_label <- ifelse(purebred_data$Sex == 1, "Female", "Male")
purebred_data2$Sex_Class_at_HLES <- ifelse(grepl("^Female", purebred_data$Sex_Class_at_HLES), "Female", "Male")

mismatch_rows <- purebred_data2[purebred_data2$Sex_label != purebred_data2$Sex_Class_at_HLES, ]
nrow(mismatch_rows) # 27
mismatch_rows[,c("Sex", "Sex_Class_at_HLES")]

hles_do <- read.csv("/home/yhSong/data/02_ByProject/DogAgingProject/SurveyData/v2022/DAP_2022_HLES_dog_owner_v1.0.csv")
dim(hles_do) # 43517   933
length(unique(hles_do$dog_id)) # 43517

dap_physical_level <- hles_do[, c("dog_id", "pa_activity_level", "pa_avg_daily_active_hours", "pa_avg_daily_active_minutes", "pa_avg_activity_intensity")]
dim(dap_physical_level) # 43517  5

sum(is.na(dap_physical_level$pa_activity_level)) # 0
sum(is.na(dap_physical_level$pa_avg_daily_active_hours)) # 0
sum(is.na(dap_physical_level$pa_avg_daily_active_minutes)) # 0
sum(is.na(dap_physical_level$pa_avg_activity_intensity)) # 0  

dap_physical_level$Total_Active_Minutes <- dap_physical_level$pa_avg_daily_active_hours * 60 + dap_physical_level$pa_avg_daily_active_minutes
head(dap_physical_level)

sum(is.na(dap_physical_level$Total_Active_Minutes)) # 0

colnames(dap_physical_level)

pa_pca_input <- dap_physical_level[, c("pa_activity_level", "pa_avg_activity_intensity", "Total_Active_Minutes")]
pa_pca <- prcomp(pa_pca_input, center = TRUE, scale. = TRUE)
summary(pa_pca)

# pa_pca$x[,1]
dap_physical_level$pa_score <- pa_pca$x[,1]
head(dap_physical_level)

purebred_total <- merge(purebred_data, dap_physical_level, by.x = "IID", by.y = "dog_id", all.x = TRUE)
dim(purebred_total) # 3864   44
head(purebred_total)
colnames(purebred_total)

# cslb 파일
cslb <- fread("/data/02_ByProject/DogAgingProject/SurveyData/v2022/DAP_2022_CSLB_v1.0.csv")
colnames(cslb)
View(cslb)
head(cslb)
dim(cslb)

# 날짜형으로 변환
cslb[, cslb_date := as.Date(cslb_date)]

# dog_id별로 최근 날짜 선택
cslb_latest <- cslb[order(dog_id, -cslb_date), .SD[1], by = dog_id]

# 결과 확인
head(cslb_latest)
dim(cslb_latest) # 30872 18

# 병합 (IID와 dog_id 기준)
merged_data <- merge(purebred_total, cslb_latest[, .(dog_id, cslb_score)],
                     by.x = "IID", by.y = "dog_id",
                     all.x = FALSE, all.y = FALSE)  # 교집합만 남김

head(merged_data)
dim(merged_data) # 3543 47

colSums(is.na(merged_data)) # cslb_score에 NA가 있는지 확인 => 사용할 컬럼에 대해서는 0

# cslb_score가 존재하는 행만 남김 (NA 제거)
purebred_final <- merged_data
head(purebred_final)
table(purebred_final$Breed)

# 품종 정리 1st
View(data.frame(Breed = sort(unique(purebred_final$Breed))))

purebred_final_br <- purebred_final

# Akbash (NON AKC)
# American bulldog -> American Bulldog
nrow(purebred_final[purebred_final$Breed == "American bulldog",]) # 1
purebred_final_br$Breed[purebred_final_br$Breed == "American bulldog"] <- "American Bulldog"
nrow(purebred_final_br[purebred_final_br$Breed == "American Bulldog"]) # 1
nrow(purebred_final_br[purebred_final_br$Breed == "American bulldog"])

# American Pitbull Terrier (NON AKC)
# Australian Labradoodle (NON AKC)
# Black Mouth Cur (NON AKC)

# Chinese Shar-Pei -> Chinese Shar_Pei
nrow(purebred_final[purebred_final$Breed == "Chinese Shar-Pei",]) # 10
purebred_final_br$Breed[purebred_final_br$Breed == "Chinese Shar-Pei"] <- "Chinese Shar_Pei"
nrow(purebred_final_br[purebred_final_br$Breed == "Chinese Shar_Pei"]) # 10
nrow(purebred_final_br[purebred_final_br$Breed == "Chinese Shar-Pei"])

# Cobberdog(NON AKC)

# Finish Spitz -> Finnish Spitz
nrow(purebred_final[purebred_final$Breed == "Finish Spitz",]) # 2
purebred_final_br$Breed[purebred_final_br$Breed == "Finish Spitz"] <- "Finnish Spitz"
nrow(purebred_final_br[purebred_final_br$Breed == "Finnish Spitz"]) # 2
nrow(purebred_final_br[purebred_final_br$Breed == "Finish Spitz"])

# Flat-Coated Retriever -> Flat_Coated Retriever
nrow(purebred_final[purebred_final$Breed == "Flat-Coated Retriever"]) # 19
purebred_final_br$Breed[purebred_final_br$Breed == "Flat-Coated Retriever"] <- "Flat_Coated Retriever"
nrow(purebred_final_br[purebred_final_br$Breed == "Flat_Coated Retriever"]) # 19
nrow(purebred_final_br[purebred_final_br$Breed == "Flat-Coated Retriever"])

# Icelandic Sheepdog/Windsprite outcross (I thought this would be more accurate since Windsprite is not listed in the drop-down) => 제거
nrow(purebred_final[purebred_final$Breed == "Icelandic Sheepdog/Windsprite outcross (I thought this would be more accurate since Windsprite is not listed in the drop-down)"])
nrow(purebred_final_br[purebred_final_br$Breed == "Icelandic Sheepdog/Windsprite outcross (I thought this would be more accurate since Windsprite is not listed in the drop-down)"])
purebred_final_br <- purebred_final_br[
  purebred_final_br$Breed != "Icelandic Sheepdog/Windsprite outcross (I thought this would be more accurate since Windsprite is not listed in the drop-down)",
]
nrow(purebred_final_br[purebred_final_br$Breed == "Icelandic Sheepdog/Windsprite outcross (I thought this would be more accurate since Windsprite is not listed in the drop-down)"])

# Jack Russel Terrier (NON AKC)

# Klee Kai -> Alaskan Klee Kai
nrow(purebred_final[purebred_final$Breed == "Alaskan Klee Kai"]) # 1
nrow(purebred_final[purebred_final$Breed == "Klee Kai"]) # 1
purebred_final_br$Breed[purebred_final_br$Breed == "Klee Kai"] <- "Alaskan Klee Kai"
nrow(purebred_final_br[purebred_final_br$Breed == "Alaskan Klee Kai"]) # 2
nrow(purebred_final_br[purebred_final_br$Breed == "Klee Kai"])

# labradoodle -> Labradoodle (NON AKC)
nrow(purebred_final[purebred_final$Breed == "labradoodle"]) # 1
nrow(purebred_final_br[purebred_final_br$Breed == "labradoodle"]) # 1
purebred_final_br$Breed[purebred_final_br$Breed == "labradoodle"] <- "Labradoodle"
nrow(purebred_final_br[purebred_final_br$Breed == "Labradoodle"]) # 1
nrow(purebred_final_br[purebred_final_br$Breed == "labradoodle"])

# Maremma -> Maremma Sheepdog
nrow(purebred_final[purebred_final$Breed == "Maremma"]) # 1
nrow(purebred_final_br[purebred_final_br$Breed == "Maremma"]) # 1
nrow(purebred_final[purebred_final$Breed == "Maremma Sheepdog"]) # 1
nrow(purebred_final_br[purebred_final_br$Breed == "Maremma Sheepdog"]) # 1
purebred_final_br$Breed[purebred_final_br$Breed == "Maremma"] <- "Maremma Sheepdog"
nrow(purebred_final_br[purebred_final_br$Breed == "Maremma Sheepdog"]) # 2
nrow(purebred_final_br[purebred_final_br$Breed == "Maremma"])

# Mini Australian Shepherd -> Miniature Australian Shepherd (NON AKC)
nrow(purebred_final[purebred_final$Breed == "Mini Australian Shepherd"]) # 1
nrow(purebred_final_br[purebred_final_br$Breed == "Mini Australian Shepherd"]) # 1
nrow(purebred_final[purebred_final$Breed == "Miniature Australian Shepherd"]) # 1
nrow(purebred_final_br[purebred_final_br$Breed == "Miniature Australian Shepherd"]) # 1
purebred_final_br$Breed[purebred_final_br$Breed == "Mini Australian Shepherd"] <- "Miniature Australian Shepherd"
nrow(purebred_final_br[purebred_final_br$Breed == "Miniature Australian Shepherd"]) # 2
nrow(purebred_final_br[purebred_final_br$Breed == "Mini Australian Shepherd"])

# Miniature (Toy) Australian Shepherd -> Miniature Australian Shepherd (NON AKC)
nrow(purebred_final_br[purebred_final_br$Breed == "Miniature (Toy) Australian Shepherd"]) # 1
nrow(purebred_final_br[purebred_final_br$Breed == "Miniature Australian Shepherd"]) # 2
purebred_final_br$Breed[purebred_final_br$Breed == "Miniature (Toy) Australian Shepherd"] <- "Miniature Australian Shepherd"
nrow(purebred_final_br[purebred_final_br$Breed == "Miniature Australian Shepherd"]) # 3
nrow(purebred_final_br[purebred_final_br$Breed == "Miniature (Toy) Australian Shepherd"])

# Perro de Presa Canario (NON AKC) [Presa Canario]
# Persian Shepherd (NON AKC) [Sarabi Dog]

# PODENCO -> Ibizan Hound
nrow(purebred_final_br[purebred_final_br$Breed == "PODENCO"]) # 1
nrow(purebred_final_br[purebred_final_br$Breed == "Ibizan Hound"]) # 1
purebred_final_br$Breed[purebred_final_br$Breed == "PODENCO"] <- "Ibizan Hound"
nrow(purebred_final_br[purebred_final_br$Breed == "Ibizan Hound"]) # 2
nrow(purebred_final_br[purebred_final_br$Breed == "PODENCO"])

# Prague Ratter (NON AKC)
# Shiloh Shepherd (NON AKC)
# Silken Windhound (NON AKC)
# White Swiss Shepherd (NON AKC, Berger Blanc Suisse)
# Windsprite (NON AKC)

View(data.frame(Breed = sort(unique(purebred_final_br$Breed))))
length(unique(purebred_final_br$Breed)) # 221
dim(purebred_final_br) # 3542 49

# FID, IID 컬럼 생성 후 저장
id_file_path <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/00_purebred_id.txt"

# FID, IID 컬럼 생성
IDs <- data.frame(FID = purebred_final_br$IID, IID = purebred_final_br$IID)
write.table(IDs, id_file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 중간 저장
save(purebred_final_br, file = "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Code/00_RDataFile/03_GWAS_Pure/00_datacheck.RData")

load("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Code/00_RDataFile/03_GWAS_Pure/00_datacheck.RData")
head(purebred_final_br)
# 00_Make_Bed.sh 수행 =>  bed, bim, fam 파일 생성

table(purebred_final_br$Sex)
table(purebred_final_br$Sex_Class_at_HLES)

purebred_final_br$Sex_edit <- ifelse(purebred_final_br$Sex == 1, "Female", "Male")
table(purebred_final_br$Sex_edit)
# Female   Male 
#   1648   1894 

# Intact 컬럼 생성
purebred_final_br$Intact <- ifelse(grepl("intact", purebred_final_br$Sex_Class_at_HLES), "Intact", "Neutered | Spayed")
table(purebred_final_br$Intact)
# Intact Neutered | Spayed 
#   963              2579 

dim(purebred_final_br) # 3542 49
head(sort(table(purebred_final_br$Breed), decreasing = TRUE), 5)
length(unique(purebred_final_br$Breed)) # 221
nrow(purebred_final_br) # 3542

write.csv(purebred_final_br, "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/00_purebred_total.csv", row.names = FALSE)


colnames(purebred_final_br)
table(purebred_final_br$Breed_Status)
table(purebred_final_br$LifeStage_Class_at_HLES)




############################################################################################################################################################################

###### 진돗개 확인
getwd()
setwd("/data/02_ByProject/DogAgingProject")

# 1. DAP_2022_DogOverview_v1.0.csv
dap_overview <- read.csv("DogOverview/DAP_2022_DogOverview_v1.0.csv")
names(dap_overview)
dim(dap_overview) 
head(dap_overview$Breed_Status)
head(dap_overview$Breed)


dim(dap_overview[dap_overview$Breed == "Jindo",])

# 2. Genetic Data
dap_geno <- fread("GeneticData/DogAgingProject_GeneticData_CuratedRelease_2022.fam", header = FALSE)
head(dap_geno)


jindo <- dap_overview[dap_overview$Breed == "Jindo",]
length(intersect(jindo$dog_id, dap_geno$V1)) # 0

jindo$dog_id
dap_geno[dap_geno$V1 == 127608,]

table(dap_overview$Breed_Status)

# Mixed breed    Purebred 
#       21607       21910 

mixed <- dap_overview[dap_overview$Breed_Status == "Mixed breed",]

table(mixed$Breed)
mixed_jindo <- mixed[grepl("Jindo", mixed$Breed, ignore.case = TRUE), ]
dim(mixed_jindo)
head(mixed_jindo)

mixed_jindo$Breed # 47
length(intersect(mixed_jindo$dog_id, dap_geno$V1)) # 9

# 9112  10932  22151  24899  30049  32762  91106  91566 106659
dap_geno[dap_geno$V1 == 9112,]
dap_geno[dap_geno$V1 == 10932,]
dap_geno[dap_geno$V1 == 22151,]
dap_geno[dap_geno$V1 == 24899,]
dap_geno[dap_geno$V1 == 30049,]
dap_geno[dap_geno$V1 == 32762,]
dap_geno[dap_geno$V1 == 91106,]
dap_geno[dap_geno$V1 == 91566,]
dap_geno[dap_geno$V1 == 106659,]

dap_overview[dap_overview$dog_id == 9112, "Breed"]
dap_overview[dap_overview$dog_id == 10932, "Breed"]
dap_overview[dap_overview$dog_id == 22151, "Breed"]
dap_overview[dap_overview$dog_id == 24899, "Breed"]
dap_overview[dap_overview$dog_id == 30049, "Breed"]
dap_overview[dap_overview$dog_id == 32762, "Breed"]
dap_overview[dap_overview$dog_id == 91106, "Breed"]
dap_overview[dap_overview$dog_id == 91566, "Breed"]
dap_overview[dap_overview$dog_id == 106659, "Breed"]
