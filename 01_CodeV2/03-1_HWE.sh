#!/bin/bash
# 03-1_HWE.sh - Hardy-Weinberg tests in the pooled sample and within ancestry groups. See README.md.
set -euo pipefail

PLINK2=/data/00_tools/plink2_linux_avx2_20231123/plink2
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023

BASE_PREFIX=${OUT_DIR}/02-1_maf0.01
INFO_FILE=${OUT_DIR}/00_all_total.csv
WORK_DIR=${OUT_DIR}/03-1_group_hwe
MIN_N=20
THREADS=32
MEM_MB=200000

for f in ${BASE_PREFIX}.pgen ${BASE_PREFIX}.pvar ${BASE_PREFIX}.psam ${INFO_FILE}; do
  if [ ! -f "$f" ]; then
    echo "파일 없음: $f"
    exit 1
  fi
done
mkdir -p ${WORK_DIR}
shopt -s nullglob
rm -f ${WORK_DIR}/G.*

$PLINK2 --pfile $BASE_PREFIX \
        --dog \
        --hardy \
        --threads $THREADS --memory $MEM_MB \
        --out ${OUT_DIR}/03-1_hwe_all

Rscript - "${BASE_PREFIX}.psam" "${INFO_FILE}" "${WORK_DIR}" "${MIN_N}" <<'EOF'
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
psam <- fread(args[1], select = 1:2, col.names = c("FID", "IID"))
info <- fread(args[2], select = c("IID", "admix_class", "top_pop"))
d <- merge(psam, info, by = "IID", all.x = TRUE, sort = FALSE)
stopifnot(nrow(d) == nrow(psam), !anyNA(d$admix_class))
d[, GRP := fifelse(admix_class == "single_breed", top_pop, "NONE")]
glist <- d[GRP != "NONE", .N, by = GRP][order(-N, GRP)]
small <- glist[N < as.integer(args[4]), GRP]
d[GRP %in% small, GRP := "NONE"]
fwrite(d[, .(FID, IID, GRP)], file.path(args[3], "group_pheno.txt"), sep = "\t")
fwrite(glist[N >= as.integer(args[4])], file.path(args[3], "group_list.tsv"), sep = "\t")
cat("검정 집단 수:", glist[N >= as.integer(args[4]), .N], " 개체 수:", glist[N >= as.integer(args[4]), sum(N)],
    " / 검정 제외: 잡종 집단", d[admix_class != "single_breed", .N], "마리, 20마리 미만 순종형 집단", length(small), "개",
    glist[N < as.integer(args[4]), sum(N)], "마리\n")
EOF

$PLINK2 --pfile $BASE_PREFIX \
        --dog \
        --pheno ${WORK_DIR}/group_pheno.txt \
        --loop-cats GRP \
        --hardy \
        --threads $THREADS --memory $MEM_MB \
        --out ${WORK_DIR}/G > ${WORK_DIR}/loop.log 2>&1 \
  || { echo "plink2 --loop-cats 실패 (loop.log 마지막 줄):"; tail -5 ${WORK_DIR}/loop.log; exit 1; }
for f in ${WORK_DIR}/G.*.hardy; do
  b=$(basename $f .hardy)
  awk 'NR == 1 {for (i = 1; i <= NF; i++) {if ($i == "P") p = i; if ($i == "O(HET_A1)") o = i; if ($i == "E(HET_A1)") e = i}
                if (!p || !o || !e) {print "열 탐색 실패: P / O(HET_A1) / E(HET_A1)" > "/dev/stderr"; exit 1}; next}
       {d = ($o == "NA" || $e == "NA") ? "NA" : $o - $e; print $p, d}' $f | gzip > ${WORK_DIR}/${b}.hwe.gz
  rm -f $f
done

N_SNP=$(grep -vc '^#' ${BASE_PREFIX}.pvar)
echo "SNP 수: ${N_SNP}"   # 4851571
echo "전체 HWE p < 1e-6: $(awk 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "P") p = i; if (!p) {print "열 탐색 실패: P" > "/dev/stderr"; exit 1}; next} $p < 1e-6' ${OUT_DIR}/03-1_hwe_all.hardy | wc -l)"
for f in ${WORK_DIR}/G.*.hwe.gz; do
  n=$(zcat $f | wc -l)
  if [ "$n" -ne "$N_SNP" ]; then echo "SNP 수 불일치: $f ($n)"; fi
done
echo "완료된 집단 수: $(ls ${WORK_DIR}/G.*.hwe.gz | wc -l) / $(( $(wc -l < ${WORK_DIR}/group_list.tsv) - 1 ))"
