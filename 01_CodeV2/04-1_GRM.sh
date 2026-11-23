#!/bin/bash
# 04-1_GRM.sh - Builds the genetic relationship matrix (whole genome and per chromosome for LOCO). See README.md.
set -euo pipefail

GCTA=/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023

BFILE=${OUT_DIR}/03-3_maf0.01_hwe
PRUNE=${OUT_DIR}/03-4_pca/prune.prune.in
QCOVAR=${OUT_DIR}/03-5_qcovar_PC10.txt
GRM_DIR=${OUT_DIR}/04-1_grm
THREADS=32

for f in ${GCTA} ${BFILE}.bed ${BFILE}.bim ${BFILE}.fam ${PRUNE} ${QCOVAR}; do
  if [ ! -f "$f" ]; then
    echo "파일 없음: $f"
    exit 1
  fi
done
mkdir -p ${GRM_DIR}

${GCTA} 2>&1 | grep -m1 "version" || true
echo "bed: $(wc -l < ${BFILE}.fam)마리, $(wc -l < ${BFILE}.bim) SNP / pruned SNP $(wc -l < ${PRUNE})"
N_PRUNE_IN=$(awk 'NR == FNR {p[$1]; next} $2 in p' ${PRUNE} ${BFILE}.bim | wc -l)
echo "pruned SNP 중 bed에 있는 SNP: ${N_PRUNE_IN}"
if [ "${N_PRUNE_IN}" -ne "$(wc -l < ${PRUNE})" ]; then echo "pruned SNP 일부가 bed에 없음"; exit 1; fi

${GCTA} --bfile ${BFILE} --autosome-num 38 --autosome --extract ${PRUNE} \
        --make-grm --thread-num ${THREADS} --out ${GRM_DIR}/grm_all > /dev/null
grep -E "SNPs? (are|is) included|individuals|GRM of" ${GRM_DIR}/grm_all.log | head -5 || true

for CHR in $(seq 1 38); do
  ${GCTA} --bfile ${BFILE} --autosome-num 38 --chr ${CHR} --extract ${PRUNE} \
          --make-grm --thread-num ${THREADS} --out ${GRM_DIR}/grm_chr${CHR} > /dev/null
done
echo "염색체별 GRM: $(ls ${GRM_DIR}/grm_chr*.grm.bin | wc -l) / 38"

Rscript - "${GRM_DIR}" "${BFILE}.fam" "${QCOVAR}" <<'EOF'
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
ids <- fread(file.path(args[1], "grm_all.grm.id"), header = FALSE)
fam <- fread(args[2], header = FALSE)
qc <- fread(args[3], header = FALSE)
# 순서까지 같아야 통과 (GCTA는 FID·IID로 맞추므로 순서가 달라도 분석은 맞지만, 03-5가 fam 순서로 썼으므로 순서까지 확인)
stopifnot(identical(as.character(ids$V2), as.character(fam$V2)), identical(as.character(ids$V2), as.character(qc$V2)))
n <- nrow(ids)
g <- readBin(file.path(args[1], "grm_all.grm.bin"), what = "numeric", n = n * (n + 1) / 2, size = 4)   # 하삼각(대각 포함), 행 순서
idx_diag <- cumsum(seq_len(n))                                                                        # 각 행의 마지막 원소 = 대각
d <- g[idx_diag]; o <- g[-idx_diag]
nsnp <- readBin(file.path(args[1], "grm_all.grm.N.bin"), what = "numeric", n = 1, size = 4)
out <- c(sprintf("개체 %d마리, GRM SNP 수(첫 쌍) %.0f", n, nsnp),
         sprintf("대각(1 + F) 평균 %.3f, 중앙값 %.3f, 범위 %.3f ~ %.3f", mean(d), median(d), min(d), max(d)),
         sprintf("비대각 평균 %.4f, 99%% %.3f, 99.9%% %.3f, 최대 %.3f", mean(o), quantile(o, 0.99), quantile(o, 0.999), max(o)),
         sprintf("비대각 > 0.354 (중복 수준) 쌍: %d, > 0.177 (1촌) 쌍: %d", sum(o > 0.354), sum(o > 0.177)))
# LOCO 전제 확인: 04-5의 --mlma-subtract-grm은 (A_all·N_all - A_c·N_c) / (N_all - N_c)로 LOCO GRM을 만듦
#   → 염색체별 GRM 38개가 전체 GRM과 같은 개체·같은 SNP로 만들어져야 함: 개체 순서 동일, 쌍마다 SNP 수(N) 합 = 전체 GRM의 N
nsum <- numeric(n * (n + 1) / 2)
for (c in 1:38) {
  pc <- file.path(args[1], paste0("grm_chr", c))
  stopifnot(identical(as.character(fread(paste0(pc, ".grm.id"), header = FALSE)$V2), as.character(ids$V2)))
  nc <- readBin(paste0(pc, ".grm.N.bin"), what = "numeric", n = n * (n + 1) / 2, size = 4)
  stopifnot(length(nc) == n * (n + 1) / 2)
  nsum <- nsum + nc
}
nall <- readBin(file.path(args[1], "grm_all.grm.N.bin"), what = "numeric", n = n * (n + 1) / 2, size = 4)
dmax <- max(abs(nsum - nall))
out <- c(out, sprintf("염색체별 GRM 38개 SNP 수 합 vs 전체 GRM SNP 수: 최대 차이 %.1f (0이어야 LOCO 뺄셈이 정확), 쌍별 SNP 수 범위 %.0f ~ %.0f",
                      dmax, min(nall), max(nall)))
writeLines(out); writeLines(out, file.path(args[1], "grm_summary.txt"))
if (dmax > 0.5) stop("염색체별 GRM과 전체 GRM의 SNP 구성이 다름 → 04-1 다시 실행")
EOF
