#!/bin/bash

# 경로 설정
PLINK2_PATH="/data/00_tools/plink2_linux_avx2_20231123/plink2"

# 고정 조건
MAF="005"
HWE="1e-6"
HWE_TAG="1e6"
BFILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.05_HWE_1e-6"
PHENO_FILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_pheno005.txt"

# 반복할 PC 목록
PCA_DIMENSIONS=(5 10 15 20 25 30)

for PC in "${PCA_DIMENSIONS[@]}"; do
  # covariate 파일 경로
  COVAR_FILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_${MAF}_${HWE_TAG}_${PC}.txt"

  # Covariate 이름 자동 추출 (FID, IID)
  COVAR_NAMES=$(head -1 "$COVAR_FILE" | tr '\t' '\n' | grep -v -E '^FID$|^IID$' | paste -sd "," -)

  # 출력 prefix
  OUT_PREFIX_LINEAR="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_${MAF}_${HWE_TAG}_${PC}_linear"

  # 실행 정보 출력
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
    --glm no-x-sex cols=+a1freq hide-covar \
    --out "$OUT_PREFIX_LINEAR" \
    --threads 20

done
