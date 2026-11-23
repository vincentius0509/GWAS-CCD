#!/bin/bash
# 01-4_Mind_Filter.sh - Applies sample missingness < 5% (51 dogs removed). See README.md.
set -euo pipefail

#   01-4_mind${MIND}_missing.{vmiss,smiss}

PLINK2=/data/00_tools/plink2_linux_avx2_20231123/plink2
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023

IN_PREFIX=${OUT_DIR}/01-4_groupmiss
MIND=0.05
OUT_PREFIX=${OUT_DIR}/01-4_mind${MIND}
MISS_PREFIX=${OUT_PREFIX}_missing
THREADS=32
MEM_MB=200000

for f in ${IN_PREFIX}.pgen ${IN_PREFIX}.pvar ${IN_PREFIX}.psam; do
  if [ ! -f "$f" ]; then
    echo "파일 없음: $f"
    exit 1
  fi
done

$PLINK2 --pfile $IN_PREFIX \
        --dog \
        --mind $MIND \
        --make-pgen \
        --threads $THREADS --memory $MEM_MB \
        --out $OUT_PREFIX

$PLINK2 --pfile $OUT_PREFIX \
        --dog \
        --missing \
        --threads $THREADS --memory $MEM_MB \
        --out $MISS_PREFIX

echo "남은 개체 수: $(( $(wc -l < ${OUT_PREFIX}.psam) - 1 ))"
echo "남은 SNP 수:  $(grep -vc '^#' ${OUT_PREFIX}.pvar)"
echo "SNP 결측률 > 0.05: $(awk 'NR > 1 && $NF > 0.05' ${MISS_PREFIX}.vmiss | wc -l)"
grep -E "samples? removed due to missing|samples? remaining|variants? remaining" ${OUT_PREFIX}.log   || true
