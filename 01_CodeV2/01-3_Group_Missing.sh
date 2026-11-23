#!/bin/bash
# 01-3_Group_Missing.sh - Counts missing genotypes per SNP within each ancestry group. See README.md.
set -euo pipefail

PLINK2=/data/00_tools/plink2_linux_avx2_20231123/plink2
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023

IN_PREFIX=${OUT_DIR}/01-2_geno0.05
INFO_FILE=${OUT_DIR}/00_all_total.csv
WORK_DIR=${OUT_DIR}/01-3_group_missing
MIN_N=10
THREADS=32
MEM_MB=100000

for f in ${IN_PREFIX}.pgen ${IN_PREFIX}.pvar ${IN_PREFIX}.psam ${INFO_FILE}; do
  if [ ! -f "$f" ]; then
    echo "파일 없음: $f"
    exit 1
  fi
done
mkdir -p ${WORK_DIR}
shopt -s nullglob
rm -f ${WORK_DIR}/G.*

Rscript - "${IN_PREFIX}.psam" "${INFO_FILE}" "${WORK_DIR}" "${MIN_N}" <<'EOF'
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
psam <- fread(args[1], select = 1:2, col.names = c("FID", "IID"))
info <- fread(args[2], select = c("IID", "admix_class", "top_pop"))
d <- merge(psam, info, by = "IID", all.x = TRUE, sort = FALSE)
stopifnot(nrow(d) == nrow(psam), !anyNA(d$admix_class))
d[, GRP := fifelse(admix_class == "single_breed", top_pop, admix_class)]
glist <- d[, .N, by = GRP][order(-N, GRP)]
small <- glist[N < as.integer(args[4]), GRP]
d[GRP %in% small, GRP := "NONE"]   # plink2 범주형 표현형의 결측 이름 (--missing-catname 기본값)
fwrite(d[, .(FID, IID, GRP)], file.path(args[3], "group_pheno.txt"), sep = "\t")
fwrite(glist[N >= as.integer(args[4])], file.path(args[3], "group_list.tsv"), sep = "\t")
cat("검정 집단 수:", glist[N >= as.integer(args[4]), .N], " 개체 수:", glist[N >= as.integer(args[4]), sum(N)],
    " / 10마리 미만 집단:", length(small), "개", glist[N < as.integer(args[4]), sum(N)], "마리\n")   # 67, 6993 / 37, 173
EOF

$PLINK2 --pfile $IN_PREFIX \
        --dog \
        --pheno ${WORK_DIR}/group_pheno.txt \
        --loop-cats GRP \
        --missing variant-only \
        --threads $THREADS --memory $MEM_MB \
        --out ${WORK_DIR}/G > ${WORK_DIR}/loop.log 2>&1 \
  || { echo "plink2 --loop-cats 실패 (loop.log 마지막 줄):"; tail -5 ${WORK_DIR}/loop.log; exit 1; }
for f in ${WORK_DIR}/G.*.vmiss; do
  b=$(basename $f .vmiss)
  awk 'NR > 1 {print $3}' $f | gzip > ${WORK_DIR}/${b}.nmiss.gz
  rm -f $f
done

N_SNP=$(grep -vc '^#' ${IN_PREFIX}.pvar)
echo "SNP 수 (.pvar): ${N_SNP}"   # 4901806
for f in ${WORK_DIR}/G.*.nmiss.gz; do
  n=$(zcat $f | wc -l)
  if [ "$n" -ne "$N_SNP" ]; then echo "SNP 수 불일치: $f ($n)"; fi
done
echo "완료된 집단 수: $(ls ${WORK_DIR}/G.*.nmiss.gz | wc -l) / $(( $(wc -l < ${WORK_DIR}/group_list.tsv) - 1 ))"
grep -E "loop-cats|samples? remaining|Error" ${WORK_DIR}/loop.log | head -5 || true
