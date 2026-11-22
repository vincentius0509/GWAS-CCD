#!/bin/bash

# --BFILE
# /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-4
# /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6
# /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-4
# /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6
# /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-4
# /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6

# OUTPUT LINK = "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_{MAF THRESHOLD}_{HWE THRESHOLD}_{PC}"

PLINK2_PATH="/data/00_tools/plink2_linux_avx2_20231123/plink2"

HWE_THRESHOLDS=("1e-6")
PCA_DIMENSIONS=(5 10 15 20 25 30)
MAF_VALUES=("001" "005" "01")

# Phenotype 파일 매핑
declare -A PHENO_FILES
PHENO_FILES["001"]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_pheno001.txt"
PHENO_FILES["005"]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_pheno005.txt"
PHENO_FILES["01"]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_pheno01.txt"

# bfile 경로 매핑
declare -A BFILE_PATHS
BFILE_PATHS["001_1e-6"]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6"
BFILE_PATHS["005_1e-6"]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6"
BFILE_PATHS["01_1e-6"]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.1_HWE_1e-6"

for HWE in "${HWE_THRESHOLDS[@]}"; do
  HWE_TAG="${HWE//-/}"   # 1e-6 → 1e6
  for MAF in "${MAF_VALUES[@]}"; do
    PHENO_FILE="${PHENO_FILES[$MAF]}"
    BFILE="${BFILE_PATHS["${MAF}_${HWE}"]}"

    for PC in "${PCA_DIMENSIONS[@]}"; do

      # Covariate 파일 경로
      COVAR_FILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_${MAF}_${HWE_TAG}_${PC}.txt"

      # 헤더에서 covariate 이름 자동 추출 (FID, IID 제외) ✅ Labrador_Retriever 제외
      COVAR_NAMES=$(head -1 "$COVAR_FILE" | tr '\t' '\n' | grep -v -E '^FID$|^IID$|^Labrador_Retriever$' | paste -sd "," -)
      # COVAR_NAMES=$(head -1 "$COVAR_FILE" | tr '\t' '\n' | grep -v -E '^FID$|^IID$' | paste -sd "," -)

      # 출력 prefix
      OUT_PREFIX_LINEAR="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_${MAF}_${HWE_TAG}_${PC}"
      OUT_PREFIX_LOGISTIC="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE2_${MAF}_${HWE_TAG}_${PC}"

      echo "▶ MAF=${MAF}, HWE=${HWE}, PC=${PC}"
      echo "📂 BFILE = $BFILE"
      echo "📂 COVAR = $COVAR_FILE"
      echo "📂 PHENO = $PHENO_FILE"
      echo "📂 Covariates: $COVAR_NAMES"

      # Linear regression (phe1)
      $PLINK2_PATH \
        --bfile "$BFILE" \
        --pheno "$PHENO_FILE" \
        --pheno-name phe1 \
        --covar "$COVAR_FILE" \
        --covar-name $COVAR_NAMES \
        --covar-variance-standardize \
        --dog \
        --glm cols=+a1freq hide-covar \
        --out "${OUT_PREFIX_LINEAR}_linear" \
        --threads 20

      # # Logistic regression (phe2)
      # $PLINK2_PATH \
      #   --bfile "$BFILE" \
      #   --pheno "$PHENO_FILE" \
      #   --pheno-name phe2 \
      #   --covar "$COVAR_FILE" \
      #   --covar-name $COVAR_NAMES \
      #   --covar-variance-standardize \
      #   --dog \
      #   --glm cols=+a1freq hide-covar \
      #   --out "${OUT_PREFIX_LOGISTIC}_logistic" \
      #   --threads 20

    done
  done
done
