#!/bin/bash
# 03-3_HWE_Exclude.sh - Removes HWE failures; final QC genotype set (4,842,736 autosomal SNPs, 7,115 dogs). See README.md.
set -euo pipefail

PLINK2=/data/00_tools/plink2_linux_avx2_20231123/plink2
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023

MAF_LIST=(0.01 0.05 0.1)
CHOSEN=W2
THREADS=32
MEM_MB=200000

for MAF in "${MAF_LIST[@]}"; do
  IN_PREFIX=${OUT_DIR}/02-1_maf${MAF}
  EXCLUDE_FILE=${OUT_DIR}/03-2_maf${MAF}_hwe_${CHOSEN}_exclude.txt
  OUT_PREFIX=${OUT_DIR}/03-3_maf${MAF}_hwe
  for f in ${IN_PREFIX}.pgen ${IN_PREFIX}.pvar ${IN_PREFIX}.psam ${EXCLUDE_FILE}; do
    if [ ! -f "$f" ]; then
      echo "파일 없음: $f"
      exit 1
    fi
  done
  echo "MAF ${MAF}: 제외 SNP $(wc -l < ${EXCLUDE_FILE})"

  $PLINK2 --pfile $IN_PREFIX \
          --dog \
          --exclude $EXCLUDE_FILE \
          --make-pgen \
          --threads $THREADS --memory $MEM_MB \
          --out $OUT_PREFIX
  $PLINK2 --pfile $OUT_PREFIX \
          --dog \
          --make-bed \
          --threads $THREADS --memory $MEM_MB \
          --out ${OUT_PREFIX}_bed_tmp
  for e in bed bim fam; do mv ${OUT_PREFIX}_bed_tmp.$e ${OUT_PREFIX}.$e; done
  mv ${OUT_PREFIX}_bed_tmp.log ${OUT_PREFIX}_bed.log

  N_PVAR=$(grep -vc '^#' ${OUT_PREFIX}.pvar); N_BIM=$(wc -l < ${OUT_PREFIX}.bim)
  N_PSAM=$(( $(wc -l < ${OUT_PREFIX}.psam) - 1 )); N_FAM=$(wc -l < ${OUT_PREFIX}.fam)
  echo "  남은 SNP: ${N_PVAR} (bim ${N_BIM}), 개체: ${N_PSAM} (fam ${N_FAM})"
  if [ "${N_PVAR}" -ne "${N_BIM}" ] || [ "${N_PSAM}" -ne "${N_FAM}" ]; then
    echo "pgen과 bed의 SNP·개체 수가 다름: ${OUT_PREFIX}"; exit 1
  fi
done
echo ".fam 성별 (1 = Male, 2 = Female):"; awk '{print $5}' ${OUT_DIR}/03-3_maf0.01_hwe.fam | sort | uniq -c
