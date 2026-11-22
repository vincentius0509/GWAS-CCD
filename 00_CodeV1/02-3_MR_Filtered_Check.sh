#!/bin/bash

# PLINK 실행 경로
PLINK_PATH=/data/00_tools/plink_linux_x86_64_20231018/plink

# 분석 대상 MAF 파일 경로
INPUT_DIR=./analysis/03_GWAS_Pure

# MAF Threshold 리스트
MAF_THRESHOLDS=(0.01 0.05 0.1)

# 각 MAF threshold별 missing rate 계산
for MAF_THRESHOLD in "${MAF_THRESHOLDS[@]}"; do
  echo "Processing missing rate for MAF threshold: ${MAF_THRESHOLD}"

  INPUT_BFILE=${INPUT_DIR}/02_MAF_${MAF_THRESHOLD}_Filtered
  OUTPUT_PREFIX=${INPUT_DIR}/02_Missing_MAF_${MAF_THRESHOLD}_Filtered_Check

  $PLINK_PATH \
    --bfile "$INPUT_BFILE" \
    --dog \
    --missing \
    --out "$OUTPUT_PREFIX" \
    --threads 1

done
