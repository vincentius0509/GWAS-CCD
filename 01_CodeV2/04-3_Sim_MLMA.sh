#!/bin/bash
# 04-3_Sim_MLMA.sh - Runs GCTA MLMA on 20 simulated breed-structure phenotypes with 0-30 PCs. See README.md.
set -euo pipefail

GCTA=/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023

BFILE=${OUT_DIR}/03-3_maf0.01_hwe
PRUNE=${OUT_DIR}/03-4_pca/prune.prune.in
GRM=${OUT_DIR}/04-1_grm/grm_all
SIM_PHENO=${OUT_DIR}/04-2_sim_pheno.txt
QCOVAR_PREFIX=${OUT_DIR}/04-2_sim_qcovar_PC     # + {k}.txt
WORK_DIR=${OUT_DIR}/04-3_sim
PC_SETS=(0 5 10 15 20 25 30)
N_SIM=20
N_JOBS=16
THREADS_PER_JOB=4
RESUME=${RESUME:-0}

for f in ${GCTA} ${BFILE}.bed ${BFILE}.bim ${BFILE}.fam ${PRUNE} ${GRM}.grm.bin ${GRM}.grm.N.bin ${GRM}.grm.id ${SIM_PHENO}; do
  if [ ! -f "$f" ]; then
    echo "파일 없음: $f"
    exit 1
  fi
done
for k in "${PC_SETS[@]}"; do
  if [ "${k}" -gt 0 ] && [ ! -f "${QCOVAR_PREFIX}${k}.txt" ]; then echo "파일 없음: ${QCOVAR_PREFIX}${k}.txt"; exit 1; fi
done
N_IND=$(wc -l < ${BFILE}.fam)
N_SIM_COL=$(( $(awk 'NR == 1 {print NF}' ${SIM_PHENO}) - 2 ))
if [ "${N_SIM_COL}" -lt "${N_SIM}" ]; then echo "가상 표현형 열 ${N_SIM_COL}개 < N_SIM ${N_SIM}"; exit 1; fi
export N_IND

if [ "${RESUME}" != "1" ]; then rm -rf ${WORK_DIR}/res ${WORK_DIR}/tmp; fi
mkdir -p ${WORK_DIR}/log ${WORK_DIR}/tmp ${WORK_DIR}/res
awk 'NR % 4 == 1' ${PRUNE} > ${WORK_DIR}/test_snps.txt
echo "λ 계산용 검정 SNP: $(wc -l < ${WORK_DIR}/test_snps.txt) / pruned $(wc -l < ${PRUNE}), 개체 ${N_IND}, RESUME=${RESUME}"

run_one() {
  local k=$1 s=$2
  local out=${WORK_DIR}/tmp/PC${k}_SIM${s}
  local res=${WORK_DIR}/res/PC${k}_SIM${s}.tsv
  if [ -s "${res}" ]; then return 0; fi
  local qc=""
  if [ "${k}" -gt 0 ]; then qc="--qcovar ${QCOVAR_PREFIX}${k}.txt"; fi
  ${GCTA} --bfile ${BFILE} --autosome-num 38 --extract ${WORK_DIR}/test_snps.txt --grm ${GRM} \
          --mlma --pheno ${SIM_PHENO} --mpheno ${s} ${qc} \
          --thread-num ${THREADS_PER_JOB} --out ${out} > /dev/null 2>&1 \
    || { echo "GCTA 실패: PC${k} SIM${s} (로그: ${out}.log)"; return 1; }
  if [ ! -s "${out}.mlma" ]; then echo "결과 없음: ${out}.mlma"; return 1; fi
  local n_common
  n_common=$(awk '/individuals are in common/ {print $1}' ${out}.log | tail -1)
  if [ "${n_common:-0}" != "${N_IND}" ]; then echo "분석 개체 수 ${n_common:-?} ≠ ${N_IND}: PC${k} SIM${s}"; return 1; fi
  local lam n
  read lam n < <(awk 'NR == 1 {if ($7 != "b" || $8 != "se") {print "열 확인 실패" > "/dev/stderr"; exit 1}; next}
                      $8 != "nan" && $8 + 0 > 0 && $7 != "nan" {print ($7 / $8)^2}' ${out}.mlma \
                 | sort -g | awk '{a[NR] = $1} END {if (NR == 0) {print "NA 0"; exit}
                                    m = (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2; printf "%.4f %d\n", m / 0.4549364, NR}')
  if ! [[ "${lam}" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then echo "λ 계산 실패 (${lam}): PC${k} SIM${s}"; return 1; fi
  local vg="NA"
  if grep -q "Log-likelihood ratio converged" ${out}.log; then
    vg=$(awk '/^Iter\./ {f = 1; next} f && $1 ~ /^[0-9]+$/ && NF >= 4 {g = $3; e = $4} END {if (g != "" && g + e > 0) printf "%.4f", g / (g + e)}' ${out}.log)
  else
    echo "경고: REML 수렴 문구 없음 → vg_vp = NA: PC${k} SIM${s}"
  fi
  if grep -qi "constrained" ${out}.log; then echo "경고: 분산 성분 경계값 고정(constrained): PC${k} SIM${s}"; fi
  echo -e "${k}\t${s}\t${lam}\t${vg:-NA}\t${n}" > ${res}.part && mv ${res}.part ${res}
  mv ${out}.log ${WORK_DIR}/log/PC${k}_SIM${s}.log
  rm -f ${out}.mlma
}
export -f run_one
export GCTA BFILE PRUNE GRM SIM_PHENO QCOVAR_PREFIX WORK_DIR THREADS_PER_JOB

if ! for k in "${PC_SETS[@]}"; do for s in $(seq 1 ${N_SIM}); do echo "${k} ${s}"; done; done \
     | xargs -P ${N_JOBS} -n 2 bash -c 'run_one "$0" "$1"'; then
  echo "일부 실행 실패 → 위 메시지 확인 후 RESUME=1 bash 04-3_Sim_MLMA.sh 로 실패한 조합만 다시 실행"
fi
rmdir ${WORK_DIR}/tmp 2>/dev/null || true

shopt -s nullglob
res_files=(${WORK_DIR}/res/PC*_SIM*.tsv)
{ echo -e "n_pc\tsim\tlambda\tvg_vp\tn_snp"; if [ ${#res_files[@]} -gt 0 ]; then cat "${res_files[@]}" | sort -k1,1n -k2,2n; fi; } > ${WORK_DIR}/lambda.tsv
N_DONE=${#res_files[@]}
N_EXPECT=$(( ${#PC_SETS[@]} * N_SIM ))
echo "완료: ${N_DONE} / ${N_EXPECT}"
if [ "${N_DONE}" -ne "${N_EXPECT}" ]; then echo "일부 실행 실패 → RESUME=1로 다시 실행"; exit 1; fi
N_BAD=$(awk -F'\t' 'NR > 1 && NF != 5' ${WORK_DIR}/lambda.tsv | wc -l)
if [ "${N_BAD}" -gt 0 ]; then echo "lambda.tsv 형식 오류 ${N_BAD}줄"; exit 1; fi
