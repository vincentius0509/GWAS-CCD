#!/bin/bash
# 01-1_Missing_Rate.sh - Extracts the 7,166 dogs with phenotype and covariates, corrects sex, and computes unfiltered missingness and MAF. See README.md.
set -euo pipefail

PLINK2=/data/00_tools/plink2_linux_avx2_20231123/plink2
GENO_DIR=/data/02_ByProject/DogAgingProject/GeneticData
OUT_DIR=${GENO_DIR}/analysis/04_GWAS_2023

INPUT_BFILE=${GENO_DIR}/v2023/DogAgingProject_GeneticData_CuratedRelease_2023_ed
KEEP_FILE=${OUT_DIR}/00_keep_all.txt
COVAR_FILE=${OUT_DIR}/00_covar.txt
SEX_FILE=${OUT_DIR}/01-1_sex_update.txt
DUP_FILE=${OUT_DIR}/01-1_dup_pos_exclude.txt
OUT_PREFIX=${OUT_DIR}/01-1_raw
THREADS=32
MEM_MB=200000

for f in ${INPUT_BFILE}.bed ${INPUT_BFILE}.bim ${INPUT_BFILE}.fam ${KEEP_FILE} ${COVAR_FILE}; do
  if [ ! -f "$f" ]; then
    echo "파일 없음: $f"
    exit 1
  fi
done
echo "KEEP 개체 수: $(wc -l < ${KEEP_FILE})"   # 7166

awk 'BEGIN{OFS="\t"} {s = ($3 == "Male") ? 1 : (($3 == "Female") ? 2 : 0); print $1, $2, s}' ${COVAR_FILE} > ${SEX_FILE}
echo "sex_meta 기준 성별 (1 = Male, 2 = Female):"; awk '{print $3}' ${SEX_FILE} | sort | uniq -c
awk 'NR == FNR {s[$2] = $3; next} ($2 in s) {f = ($5 == 1) ? 2 : (($5 == 2) ? 1 : 0); if (f != s[$2]) n++} END {print ".fam과 다른 개체:", n + 0}' \
    ${SEX_FILE} ${INPUT_BFILE}.fam

awk '{k = $1 ":" $4}
     k == pk { if (!pp) print pid; print $2; pp = 1; next }
     { pp = 0; pk = k; pid = $2 }' ${INPUT_BFILE}.bim > ${DUP_FILE}
echo "중복 위치 변이 수: $(wc -l < ${DUP_FILE})"
echo "중복 위치 수:     $(awk '{print $1":"$4}' ${INPUT_BFILE}.bim | uniq -d | wc -l)"

$PLINK2 --bfile $INPUT_BFILE \
        --dog \
        --keep $KEEP_FILE \
        --update-sex $SEX_FILE \
        --exclude $DUP_FILE \
        --chr 1-38 \
        --make-pgen \
        --threads $THREADS --memory $MEM_MB \
        --out $OUT_PREFIX

$PLINK2 --pfile $OUT_PREFIX \
        --dog \
        --missing \
        --freq \
        --threads $THREADS --memory $MEM_MB \
        --out $OUT_PREFIX

echo "출력 개체 수: $(( $(wc -l < ${OUT_PREFIX}.psam) - 1 ))"   # 7166
echo "출력 SNP 수:  $(grep -vc '^#' ${OUT_PREFIX}.pvar)"
echo "psam 성별 (1 = Male, 2 = Female):"; awk 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "SEX") c = i; next} {print $c}' ${OUT_PREFIX}.psam | sort | uniq -c
grep -E "loaded from|samples? remaining|variants? remaining|--update-sex|missing call" ${OUT_PREFIX}.log   || true
