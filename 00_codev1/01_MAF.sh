#!/bin/bash

PLINK_PATH=/data/00_tools/plink_linux_x86_64_20231018/plink
INPUT_PREFIX=./analysis/03_GWAS_Pure/00_purebred_data

# 각각의 THRESHOLD에 대해 PLINK 실행
# MAF_LIST=(0.001 0.005 0.01 0.05 0.1)
MAF_LIST=(0.01 0.05 0.1)

for MAF in "${MAF_LIST[@]}"; do
  OUT_PREFIX="./analysis/03_GWAS_Pure/01_MAF_${MAF}"
  echo "Running PLINK with MAF=${MAF} -> Output: ${OUT_PREFIX}"

  $PLINK_PATH \
    --bfile "$INPUT_PREFIX" \
    --dog \
    --freq \
    --maf "$MAF" \
    --make-bed \
    --out "$OUT_PREFIX" \
    --threads 1

done
