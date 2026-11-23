#!/bin/bash
# 02-1_MAF_Filter.sh - Re-applies variant missingness and builds MAF 0.01 / 0.05 / 0.1 sets with allele counts per ancestry group. See README.md.
set -euo pipefail

PLINK2=/data/00_tools/plink2_linux_avx2_20231123/plink2
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023

IN_PREFIX=${OUT_DIR}/01-4_mind0.05
INFO_FILE=${OUT_DIR}/00_all_total.csv
GENO=0.05
MAF_LIST=(0.01 0.05 0.1)
BASE_PREFIX=${OUT_DIR}/02-1_maf${MAF_LIST[0]}
WORK_DIR=${OUT_DIR}/02-1_group_freq
MIN_N=10
THREADS=32
MEM_MB=200000

for f in ${IN_PREFIX}.pgen ${IN_PREFIX}.pvar ${IN_PREFIX}.psam ${INFO_FILE}; do
  if [ ! -f "$f" ]; then
    echo "파일 없음: $f"
    exit 1
  fi
done
mkdir -p ${WORK_DIR}
shopt -s nullglob
rm -f ${WORK_DIR}/G.*

$PLINK2 --pfile $IN_PREFIX \
        --dog \
        --geno $GENO \
        --freq counts \
        --threads $THREADS --memory $MEM_MB \
        --out ${OUT_DIR}/02-1_geno${GENO}

for MAF in "${MAF_LIST[@]}"; do
  echo "MAF = ${MAF}"
  $PLINK2 --pfile $IN_PREFIX \
          --dog \
          --geno $GENO \
          --maf $MAF \
          --make-pgen \
          --threads $THREADS --memory $MEM_MB \
          --out ${OUT_DIR}/02-1_maf${MAF}
done

Rscript - "${BASE_PREFIX}.psam" "${INFO_FILE}" "${WORK_DIR}" "${MIN_N}" <<'EOF'
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
psam <- fread(args[1], select = 1:2, col.names = c("FID", "IID"))
info <- fread(args[2], select = c("IID", "admix_class", "top_pop"))
d <- merge(psam, info, by = "IID", all.x = TRUE, sort = FALSE)
stopifnot(nrow(d) == nrow(psam), !anyNA(d$admix_class))
d[, GRP := fifelse(admix_class == "single_breed", top_pop, admix_class)]
glist <- d[, .N, by = GRP][order(-N, GRP)]
small <- glist[N < as.integer(args[4]), GRP]
d[GRP %in% small, GRP := "NONE"]
fwrite(d[, .(FID, IID, GRP)], file.path(args[3], "group_pheno.txt"), sep = "\t")
fwrite(glist[N >= as.integer(args[4])], file.path(args[3], "group_list.tsv"), sep = "\t")
cat("검정 집단 수:", glist[N >= as.integer(args[4]), .N], " 개체 수:", glist[N >= as.integer(args[4]), sum(N)],
    " / 10마리 미만 집단:", length(small), "개", glist[N < as.integer(args[4]), sum(N)], "마리 / 전체:", nrow(psam), "\n")
EOF

$PLINK2 --pfile $BASE_PREFIX \
        --dog \
        --pheno ${WORK_DIR}/group_pheno.txt \
        --loop-cats GRP \
        --freq counts \
        --threads $THREADS --memory $MEM_MB \
        --out ${WORK_DIR}/G > ${WORK_DIR}/loop.log 2>&1 \
  || { echo "plink2 --loop-cats 실패 (loop.log 마지막 줄):"; tail -5 ${WORK_DIR}/loop.log; exit 1; }
for f in ${WORK_DIR}/G.*.acount; do
  b=$(basename $f .acount)
  awk 'NR == 1 {for (i = 1; i <= NF; i++) {if ($i == "ALT_CTS") a = i; if ($i == "OBS_CT") o = i}; next} {print $a, $o}' $f | gzip > ${WORK_DIR}/${b}.counts.gz
  rm -f $f
done

echo "MAF 필터 전(--geno 재적용 후) SNP 수: $(( $(wc -l < ${OUT_DIR}/02-1_geno${GENO}.acount) - 1 ))"
for MAF in "${MAF_LIST[@]}"; do
  echo "MAF ${MAF} SNP 수: $(grep -vc '^#' ${OUT_DIR}/02-1_maf${MAF}.pvar)"
  grep -E "removed due to|variants? remaining" ${OUT_DIR}/02-1_maf${MAF}.log || true
done
N_SNP=$(grep -vc '^#' ${BASE_PREFIX}.pvar)
for f in ${WORK_DIR}/G.*.counts.gz; do
  n=$(zcat $f | wc -l)
  if [ "$n" -ne "$N_SNP" ]; then echo "SNP 수 불일치: $f ($n)"; fi
done
echo "완료된 집단 수: $(ls ${WORK_DIR}/G.*.counts.gz | wc -l) / $(( $(wc -l < ${WORK_DIR}/group_list.tsv) - 1 ))"
