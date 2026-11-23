#!/bin/bash
# 01-2_Geno_Filter.sh - Applies variant missingness < 5% and recomputes missingness. See README.md.
set -euo pipefail

PLINK2=/data/00_tools/plink2_linux_avx2_20231123/plink2
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023

IN_PREFIX=${OUT_DIR}/01-1_raw
GENO=0.05
OUT_PREFIX=${OUT_DIR}/01-2_geno${GENO}
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
        --geno $GENO \
        --make-pgen \
        --threads $THREADS --memory $MEM_MB \
        --out $OUT_PREFIX

$PLINK2 --pfile $OUT_PREFIX \
        --dog \
        --missing \
        --threads $THREADS --memory $MEM_MB \
        --out $MISS_PREFIX

echo "남은 SNP 수:  $(grep -vc '^#' ${OUT_PREFIX}.pvar)"
echo "남은 개체 수: $(( $(wc -l < ${OUT_PREFIX}.psam) - 1 ))"
echo "SNP 결측률 최댓값: $(awk 'NR > 1 {if ($NF > m) m = $NF} END {print m}' ${MISS_PREFIX}.vmiss)"
echo "개체 결측률 요약(최소/중앙/최대):"
awk 'NR > 1 {print $NF}' ${MISS_PREFIX}.smiss | sort -g | awk '{a[NR] = $1} END {m = (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2; print a[1], m, a[NR]}'
grep -E "variants? removed due to missing|variants? remaining|samples? remaining" ${OUT_PREFIX}.log   || true
