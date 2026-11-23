#!/bin/bash
# 04-5_MLMA_LOCO.sh - Primary GWAS: GCTA MLMA-LOCO of CSLB scores with 0-30 PCs. See README.md.
set -euo pipefail

#       QCOVAR=03-5_qcovar_qmat_PC10.txt bash 04-5_MLMA_LOCO.sh / KEEP=00_keep_age8.txt bash 04-5_MLMA_LOCO.sh

GCTA=/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023

MAF=${MAF:-0.01}
MPHENO=${MPHENO:-1}
case "${MPHENO}" in
  1) LABEL=${LABEL:-rint} ;;
  2) LABEL=${LABEL:-raw} ;;
  *) echo "MPHENO는 1(rint) 또는 2(raw): ${MPHENO}"; exit 1 ;;
esac
COVAR=${OUT_DIR}/03-5_covar.txt
BFILE=${OUT_DIR}/03-3_maf${MAF}_hwe
PHENO=${OUT_DIR}/04-2_pheno.txt
GRM_DIR=${OUT_DIR}/04-1_grm
WORK_DIR=${OUT_DIR}/04-5_gwas
N_CONT=5
N_JOBS=${N_JOBS:-38}
THREADS_PER_JOB=${THREADS_PER_JOB:-2}
MODEL_JOBS=${MODEL_JOBS:-3}
RESUME=${RESUME:-0}
mkdir -p ${WORK_DIR}/log

if [ -n "${PC_LIST:-}" ]; then
  if [ -n "${QCOVAR:-}" ] || [ -n "${N_PC:-}" ]; then echo "PC_LIST와 QCOVAR·N_PC는 같이 쓸 수 없음"; exit 1; fi
  for k in ${PC_LIST}; do if ! [[ "${k}" =~ ^[0-9]+$ ]]; then echo "PC_LIST에 숫자가 아닌 값: ${k}"; exit 1; fi; done
  echo "PC_LIST: ${PC_LIST} / 동시 모델 ${MODEL_JOBS}개 × 염색체 ${N_JOBS}개 × 스레드 ${THREADS_PER_JOB} / 시작 $(date '+%F %T')"
  if ! echo ${PC_LIST} | tr ' ' '\n' | xargs -P ${MODEL_JOBS} -I{} bash -c \
       'if [ "${DRYRUN:-0}" = "1" ]; then env PC_LIST= N_PC={} bash "$0"; exit $?; fi
        env PC_LIST= N_PC={} bash "$0" 2>&1 | tee "$1/log/PC{}_run_$(date +%Y%m%d_%H%M%S).log"; exit ${PIPESTATUS[0]}' \
       "$0" "${WORK_DIR}"; then
    echo "일부 모델 실패 → 위 메시지 확인 후 실패한 PC 개수만 RESUME=1 N_PC=k bash $(basename $0) 로 다시 실행"
    exit 1
  fi
  if [ "${DRYRUN:-0}" = "1" ]; then exit 0; fi
  echo "PC_LIST 전체 완료 $(date '+%F %T') — 모델별 λ (자세한 내용은 각 _summary.txt):"
  for k in ${PC_LIST}; do f=$(ls -t ${WORK_DIR}/*_PC${k}_summary.txt 2>/dev/null | head -1 || true); if [ -n "${f}" ]; then echo "PC ${k}: $(grep '^λ' ${f} | tr '\n' ' ')"; fi; done
  exit 0
fi

resolve() { if [ -z "$1" ] || [[ "$1" == /* ]]; then echo "$1"; else echo "${OUT_DIR}/$1"; fi; }

if [ -z "${QCOVAR:-}" ]; then
  if [ -z "${N_PC:-}" ]; then
    if [ ! -f ${OUT_DIR}/04-4_pc_choice.txt ]; then echo "파일 없음: ${OUT_DIR}/04-4_pc_choice.txt (04-4 먼저 실행 또는 N_PC 지정)"; exit 1; fi
    N_PC=$(tr -d '[:space:]' < ${OUT_DIR}/04-4_pc_choice.txt)
  fi
  if ! [[ "${N_PC}" =~ ^[0-9]+$ ]]; then echo "N_PC가 숫자가 아님: '${N_PC}'"; exit 1; fi
  if [ "${N_PC}" -eq 0 ]; then
    QCOVAR=${OUT_DIR}/04-5_qcovar_PC0.txt
    awk -v n=$(( 2 + N_CONT )) '{line = $1; for (i = 2; i <= n; i++) line = line " " $i; print line}' ${OUT_DIR}/03-5_qcovar_PC5.txt > ${QCOVAR}
  else
    QCOVAR=${OUT_DIR}/03-5_qcovar_PC${N_PC}.txt
  fi
else
  QCOVAR=$(resolve "${QCOVAR}")
fi
KEEP=$(resolve "${KEEP:-}")
QTAG=$(basename ${QCOVAR} .txt); QTAG=${QTAG#*qcovar_}          # 03-5_qcovar_PC10 → PC10 
KTAG=""
if [ -n "${KEEP}" ]; then KTAG=$(basename ${KEEP} .txt); KTAG="_keep-${KTAG#*keep_}"; fi     # 00_keep_age8 → _keep-age8
RUN=${LABEL}_maf${MAF}_${QTAG}${KTAG}
TMP=${WORK_DIR}/tmp_${RUN}

for f in ${GCTA} ${BFILE}.bed ${BFILE}.bim ${BFILE}.fam ${PHENO} ${QCOVAR} ${COVAR} ${GRM_DIR}/grm_all.grm.bin ${GRM_DIR}/grm_all.grm.N.bin ${GRM_DIR}/grm_all.grm.id; do
  if [ ! -f "$f" ]; then echo "파일 없음: $f"; exit 1; fi
done
for c in $(seq 1 38); do
  for s in bin N.bin id; do if [ ! -f ${GRM_DIR}/grm_chr${c}.grm.${s} ]; then echo "파일 없음: ${GRM_DIR}/grm_chr${c}.grm.${s}"; exit 1; fi; done
done
if [ -n "${KEEP}" ] && [ ! -f "${KEEP}" ]; then echo "파일 없음: ${KEEP}"; exit 1; fi

N_EXPECT=$(awk -v mp=${MPHENO} -v keep="${KEEP}" -v nq=0 '
  function die(m) {print "입력 검사 실패: " m > "/dev/stderr"; bad = 1; exit 1}
  FILENAME == ARGV[1] {fam[$2]; nfam++; next}
  keep != "" && FILENAME == keep {if (NF < 2) die("keep 파일은 FID IID 두 열"); if (!($2 in fam)) {nk_out++; next}; if ($2 in kp) die("keep ID 중복 " $2); kp[$2]; next}
  {
    f = FILENAME
    if (!(f in nf)) nf[f] = NF
    if (NF != nf[f]) die(f " 열 수가 일정하지 않음 (" FNR "번째 줄)")
    if ((f, $2) in seen) die(f " ID 중복 " $2)
    seen[f, $2]
    if (f == ARGV[2]) {if (NF < 2 + mp) die("표현형 파일에 " mp "번째 표현형 열 없음"); ph[$2] = $(2 + mp)}
    else {for (i = 3; i <= NF; i++) if ($i == "NA" || $i == "nan" || $i == "-9" || $i == "") die(f " 결측 공변량 (" $2 ", " i "번째 열)"); cv[f, $2]}
  }
  END {
    if (bad) exit 1
    n = 0
    for (id in fam) {
      if (keep != "" && !(id in kp)) continue
      if (!(id in ph)) die("표현형 파일에 없는 개체 " id)
      if (ph[id] == "NA" || ph[id] == "-9") continue
      for (j = 3; j < ARGC; j++) if (ARGV[j] != keep && !((ARGV[j], id) in cv)) die(ARGV[j] " 에 없는 개체 " id)
      n++
    }
    if (nk_out) print "참고: keep 개체 중 " nk_out "마리는 최종 표본(fam)에 없어 제외 (00 단계 7,166마리 기준 목록 → 01에서 QC로 빠진 개체)" > "/dev/stderr"
    print n
  }' ${BFILE}.fam ${PHENO} ${QCOVAR} ${COVAR} ${KEEP})
NQ=$(( $(awk 'NR == 1 {print NF}' ${QCOVAR}) - 2 ))
if [[ "${QTAG}" =~ ^PC([0-9]+)$ ]]; then
  k=${BASH_REMATCH[1]}
  if [ "${NQ}" -ne $(( N_CONT + k )) ]; then echo "qcovar 열 수 ${NQ} ≠ 연속형 ${N_CONT} + PC ${k}"; exit 1; fi
fi
echo "실행: ${RUN} / 표현형 열 ${MPHENO} / qcovar $(basename ${QCOVAR}) (공변량 ${NQ}개) / keep ${KEEP:-전체} / 예상 분석 개체 ${N_EXPECT}"
if [ "${N_EXPECT}" -lt 100 ]; then echo "분석 개체가 너무 적음: ${N_EXPECT}"; exit 1; fi

if [ "${DRYRUN:-0}" = "1" ]; then echo "DRYRUN=1: 입력 검사만 하고 종료"; exit 0; fi
if [ "${RESUME}" != "1" ]; then rm -rf ${TMP}; rm -f ${WORK_DIR}/log/${RUN}_chr*.log; fi
mkdir -p ${TMP}

run_chr() {
  local c=$1
  local out=${TMP}/chr${c}
  if [ -s ${out}.mlma ] && [ -f ${WORK_DIR}/log/${RUN}_chr${c}.log ]; then return 0; fi
  rm -f ${out}.mlma
  local keep_opt=""
  if [ -n "${KEEP}" ]; then keep_opt="--keep ${KEEP}"; fi
  ${GCTA} --bfile ${BFILE} --autosome-num 38 --chr ${c} ${keep_opt} \
          --grm ${GRM_DIR}/grm_all --mlma-subtract-grm ${GRM_DIR}/grm_chr${c} --mlma \
          --pheno ${PHENO} --mpheno ${MPHENO} --qcovar ${QCOVAR} --covar ${COVAR} \
          --thread-num ${THREADS_PER_JOB} --out ${out} > /dev/null 2>&1 \
    || { echo "GCTA 실패: chr${c} (로그: ${out}.log)"; rm -f ${out}.mlma; return 1; }
  local n_common
  n_common=$(awk '/individuals are in common/ {print $1}' ${out}.log | tail -1)
  if [ "${n_common:-0}" != "${N_EXPECT}" ]; then echo "chr${c}: 분석 개체 ${n_common:-?} ≠ 예상 ${N_EXPECT} (로그: ${out}.log)"; rm -f ${out}.mlma; return 1; fi
  if ! grep -q "Log-likelihood ratio converged" ${out}.log; then echo "경고: chr${c} REML 수렴 문구 없음 (로그 확인)"; fi
  if grep -qi "constrained" ${out}.log; then echo "경고: chr${c} 분산 성분 경계값 고정(constrained)"; fi
  mv ${out}.log ${WORK_DIR}/log/${RUN}_chr${c}.log
}
export -f run_chr
export GCTA BFILE KEEP GRM_DIR PHENO MPHENO QCOVAR COVAR WORK_DIR RUN TMP THREADS_PER_JOB N_EXPECT

if ! seq 1 38 | xargs -P ${N_JOBS} -I{} bash -c 'run_chr {}'; then
  echo "일부 염색체 실패 → 위 메시지 확인 후 같은 설정 앞에 RESUME=1을 붙여 다시 실행하면 끝난 염색체는 건너뜀"
  exit 1
fi

shopt -s nullglob
res_files=(${TMP}/chr*.mlma)
if [ "${#res_files[@]}" -ne 38 ]; then echo "염색체 결과 ${#res_files[@]} / 38 → RESUME=1로 다시 실행"; exit 1; fi
for c in $(seq 1 38); do
  if [ "$(head -1 ${TMP}/chr${c}.mlma | cut -f1-9)" != "$(printf 'Chr\tSNP\tbp\tA1\tA2\tFreq\tb\tse\tp')" ]; then echo "chr${c} 결과 헤더가 예상과 다름"; exit 1; fi
done
{ head -1 ${TMP}/chr1.mlma; for c in $(seq 1 38); do tail -n +2 ${TMP}/chr${c}.mlma; done; } | gzip > ${WORK_DIR}/${RUN}.mlma.gz
N_RES=$(( $(zcat ${WORK_DIR}/${RUN}.mlma.gz | wc -l) - 1 ))
N_BIM=$(wc -l < ${BFILE}.bim)
awk '{n[$1]++} END {for (c in n) print c, n[c]}' ${BFILE}.bim > ${TMP}/bim_count.txt
CHR_DIFF=$(for c in $(seq 1 38); do
  nb=$(awk -v c=${c} '$1 == c {print $2}' ${TMP}/bim_count.txt); nb=${nb:-0}; nr=$(( $(wc -l < ${TMP}/chr${c}.mlma) - 1 ))
  if [ "${nb}" -ne "${nr}" ]; then echo -n "chr${c} ${nr}/${nb}  "; fi
done)
N_NAN=$(zcat ${WORK_DIR}/${RUN}.mlma.gz | awk 'NR > 1 && ($9 == "nan" || $8 == "nan" || $8 + 0 <= 0)' | wc -l)
DIFF=$(( N_BIM - N_RES ))
echo "결과 SNP 수: ${N_RES} / bim ${N_BIM} (차이 ${DIFF}), p·se = nan 또는 se ≤ 0: ${N_NAN}"
if [ -n "${CHR_DIFF}" ]; then echo "SNP 수가 다른 염색체: ${CHR_DIFF}"; fi
if [ -z "${KEEP}" ] && [ "${DIFF}" -ne 0 ]; then echo "경고: 전체 표본인데 결과 SNP가 bim과 다름 (${DIFF}개) → 로그 확인"; fi
if [ "${DIFF}" -gt $(( N_BIM / 100 )) ]; then echo "경고: 결과 SNP가 bim보다 1% 넘게 적음 (${DIFF}개)"; fi
if [ "${DIFF}" -lt 0 ]; then echo "오류: 결과 SNP가 bim보다 많음 (중복 합치기 의심)"; exit 1; fi
rm -rf ${TMP}

SUMMARY=${WORK_DIR}/${RUN}_summary.txt
{
  echo "실행: ${RUN}, 표현형 열 ${MPHENO}, qcovar $(basename ${QCOVAR}) (공변량 ${NQ}개), covar $(basename ${COVAR}), keep ${KEEP:-전체}"
  echo "분석 개체 ${N_EXPECT}, 결과 SNP ${N_RES} / bim ${N_BIM}, nan·se ≤ 0 SNP ${N_NAN}${CHR_DIFF:+, SNP 수 차이: ${CHR_DIFF}}"
  zcat ${WORK_DIR}/${RUN}.mlma.gz | awk '$8 != "nan" && $8 + 0 > 0 && $7 != "nan" && NR > 1 {print ($7 / $8)^2}' | sort -g \
    | awk '{a[NR] = $1} END {if (NR == 0) {print "λ 계산 불가 (유효 SNP 없음)"; exit}
                            m = (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2; l = m / 0.4549364
                            printf "λ (median χ² / 0.4549): %.4f (SNP %d)\n", l, NR
                            printf "λ 범위 판정 (기준 0.95 ~ 1.05): %s\n", (l >= 0.95 && l <= 1.05) ? "통과" : "불통과 → 다른 기준과 함께 검토"}'
  zcat ${WORK_DIR}/${RUN}.mlma.gz | awk 'NR > 1 && $9 != "nan" {if ($9 < 5e-8) a++; if ($9 < 1e-6) b++; if ($9 < 1e-5) c++}
       END {printf "p < 5e-8: %d / p < 1e-6: %d / p < 1e-5: %d\n", a, b, c}'
  echo "상위 10개 SNP (Chr SNP bp A1 A2 Freq b se p):"
  zcat ${WORK_DIR}/${RUN}.mlma.gz | awk 'NR > 1 && $9 != "nan"' | sort -k9,9g | head -10 || true
  echo "염색체별 V(G)/Vp (LOCO GRM이 설명하는 표현형 분산 비율, 공변량 보정 후):"
  for c in $(seq 1 38); do
    lg=${WORK_DIR}/log/${RUN}_chr${c}.log
    if grep -q "Log-likelihood ratio converged" ${lg}; then
      v=$(awk '/^Iter\./ {f = 1; next} f && $1 ~ /^[0-9]+$/ && NF >= 4 {g = $3; e = $4} END {if (g != "" && g + e > 0) printf "%.4f", g / (g + e); else printf "NA"}' ${lg})
    else v="NA(미수렴)"; fi
    printf "chr%s %s  " ${c} "${v}"
  done; echo
} | tee ${SUMMARY}
