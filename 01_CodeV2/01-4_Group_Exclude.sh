#!/bin/bash
# 01-4_Group_Exclude.sh - Removes the group-specific missingness SNPs and recomputes missingness. See README.md.
set -euo pipefail

PLINK2=/data/00_tools/plink2_linux_avx2_20231123/plink2
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023

IN_PREFIX=${OUT_DIR}/01-2_geno0.05
EXCLUDE_FILE=${OUT_DIR}/01-3_group_miss_exclude.txt
OUT_PREFIX=${OUT_DIR}/01-4_groupmiss
MISS_PREFIX=${OUT_PREFIX}_missing
THREADS=32
MEM_MB=200000

for f in ${IN_PREFIX}.pgen ${IN_PREFIX}.pvar ${IN_PREFIX}.psam ${EXCLUDE_FILE}; do
  if [ ! -f "$f" ]; then
    echo "파일 없음: $f"
    exit 1
  fi
done
echo "제외 SNP 수: $(wc -l < ${EXCLUDE_FILE})"   # 44595

$PLINK2 --pfile $IN_PREFIX \
        --dog \
        --exclude $EXCLUDE_FILE \
        --make-pgen \
        --threads $THREADS --memory $MEM_MB \
        --out $OUT_PREFIX

$PLINK2 --pfile $OUT_PREFIX \
        --dog \
        --missing \
        --threads $THREADS --memory $MEM_MB \
        --out $MISS_PREFIX

echo "남은 SNP 수:  $(grep -vc '^#' ${OUT_PREFIX}.pvar)"
echo "남은 개체 수: $(( $(wc -l < ${OUT_PREFIX}.psam) - 1 ))"   # 7166
grep -E "variants? remaining|samples? remaining" ${OUT_PREFIX}.log   || true
