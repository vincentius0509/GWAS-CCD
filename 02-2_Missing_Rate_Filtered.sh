#!/bin/bash

PLINK_PATH=/data/00_tools/plink_linux_x86_64_20231018/plink

# MAF_THRESHOLDS=(0.001 0.005 0.01 0.05 0.1)

# # 위에서 만든 bed,bim,fam 파일을 이용하여 missing rate 계산
# $PLINK_PATH \
#   --bfile ./analysis/01_GWAS/2.missing_maf_0.05_geno_mind \
#   --dog \
#   --missing \
#   --out ./analysis/01_GWAS/2.missing_maf_0.05_geno_mind_after \
#   --threads 1

# Threshold 설정 후 bed, bim, fam 파일 생성
# $PLINK_PATH \
#   --bfile ./analysis/03_GWAS_Pure/01_MAF_0.001 \
#   --dog \
#   --mind 0.280 \
#   --geno 0.68135 \
#   --make-bed \
#   --out ./analysis/03_GWAS_Pure/02_MAF_0.001_Filtered \
#   --threads 1

# $PLINK_PATH \
#   --bfile ./analysis/03_GWAS_Pure/01_MAF_0.005 \
#   --dog \
#   --mind 0.2396 \
#   --geno 0.6379 \
#   --make-bed \
#   --out ./analysis/03_GWAS_Pure/02_MAF_0.005_Filtered \
#   --threads 1

$PLINK_PATH \
  --bfile ./analysis/03_GWAS_Pure/01_MAF_0.01 \
  --dog \
  --mind 0.400 \
  --geno 0.84925 \
  --make-bed \
  --out ./analysis/03_GWAS_Pure/02_MAF_0.01_Filtered \
  --threads 1

$PLINK_PATH \
  --bfile ./analysis/03_GWAS_Pure/01_MAF_0.05 \
  --dog \
  --mind 0.1753 \
  --geno 0.5737 \
  --make-bed \
  --out ./analysis/03_GWAS_Pure/02_MAF_0.05_Filtered \
  --threads 1

$PLINK_PATH \
  --bfile ./analysis/03_GWAS_Pure/01_MAF_0.1 \
  --dog \
  --mind 0.19925 \
  --geno 0.5868 \
  --make-bed \
  --out ./analysis/03_GWAS_Pure/02_MAF_0.1_Filtered \
  --threads 1