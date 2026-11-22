#!/bin/bash

PLINK_PATH=/data/00_tools/plink_linux_x86_64_20231018/plink

# HWE_THRESHOLDS=(1e-3 1e-4 1e-5 1e-6 1e-7 1e-8)
HWE_THRESHOLDS=1e-6
PCA_DIMENSIONS=(5 10 15 20 25 30)
MAF_VALUES=(0.01 0.05 0.1)

for MAF in "${MAF_VALUES[@]}"; do
  for HWE in "${HWE_THRESHOLDS[@]}"; do
    INPUT_BFILE=./analysis/03_GWAS_Pure/03_MAF_${MAF}_HWE_${HWE}

    for PC in "${PCA_DIMENSIONS[@]}"; do
      PCA_OUT=${INPUT_BFILE}_PCA_${PC}
      $PLINK_PATH \
        --bfile "$INPUT_BFILE" \
        --dog \
        --pca "$PC" \
        --out "$PCA_OUT" \
        --threads 20
    done
  done
  echo "PCA done for MAF=${MAF}"
done
