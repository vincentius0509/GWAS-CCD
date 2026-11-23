#!/bin/bash
# 05-2_Sens_MLMA.sh - Runs the GCTA MLMA-LOCO models for strata and alternative covariates. See README.md.
set -euo pipefail

#   CHECK=1 bash 00_Code/05_GWAS_v2023/05-2_Sens_MLMA.sh
#   DRYRUN=1 bash 00_Code/05_GWAS_v2023/05-2_Sens_MLMA.sh
#   nohup bash 00_Code/05_GWAS_v2023/05-2_Sens_MLMA.sh > 04_GWAS_2023/05-2_run.log 2>&1 &
#   DRYRUN=1 RUNS="young8 sex_m sex_f wt_t1 wt_t2 wt_t3" bash 00_Code/05_GWAS_v2023/05-2_Sens_MLMA.sh
#   RUNS="young8 sex_m sex_f wt_t1 wt_t2 wt_t3" nohup bash 00_Code/05_GWAS_v2023/05-2_Sens_MLMA.sh > 04_GWAS_2023/05-2_run_strata.log 2>&1 &

GCTA=/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023/
GWAS_045=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023/04-5_gwas
RUNS_TSV=${OUT_DIR}/05-1_runs.tsv
BFILE=${OUT_DIR}/03-3_maf0.01_hwe
GRM_DIR=${OUT_DIR}/04-1_grm
WORK_DIR=${OUT_DIR}/05-2_gwas
N_JOBS=${N_JOBS:-110}
THREADS_PER_JOB=${THREADS_PER_JOB:-2}
RESUME=${RESUME:-0}
DRYRUN=${DRYRUN:-0}
CHECK=${CHECK:-0}
mkdir -p ${WORK_DIR}/log ${WORK_DIR}/cfg ${WORK_DIR}/done

for f in ${GCTA} ${RUNS_TSV} ${BFILE}.bed ${BFILE}.bim ${BFILE}.fam ${GRM_DIR}/grm_all.grm.bin ${GRM_DIR}/grm_all.grm.N.bin ${GRM_DIR}/grm_all.grm.id; do
  if [ ! -f "$f" ]; then echo "파일 없음: $f"; exit 1; fi
done
for c in $(seq 1 38); do for s in bin N.bin id; do
  if [ ! -f ${GRM_DIR}/grm_chr${c}.grm.${s} ]; then echo "파일 없음: ${GRM_DIR}/grm_chr${c}.grm.${s}"; exit 1; fi
done; done

row_of() {
  awk -F'\t' -v id="$1" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
    $h["id"] == id {cv = $h["covar_file"]; if (cv == "") cv = "-"; print $h["id"], $h["n"], $h["pheno_file"], $h["qcovar_file"], cv, $h["keep_file"], $h["engine"]}' ${RUNS_TSV}
}
if [ "${CHECK}" = "1" ]; then
  SEL="main_check"
elif [ -n "${RUNS:-}" ]; then
  SEL="${RUNS}"
else
  SEL=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["engine"] == "gcta" {print $h["id"]}' ${RUNS_TSV} | tr '\n' ' ')
fi
N_RUN=0
for id in ${SEL}; do
  read -r rid n ph qc cv kp eng <<< "$(row_of ${id})" || true
  if [ -z "${rid:-}" ]; then echo "실행 표에 없는 id: ${id}"; exit 1; fi
  if [ "${CHECK}" != "1" ] && [ "${eng}" != "gcta" ]; then echo "${id}: engine = ${eng} → 05-2 대상 아님 (gcta만)"; exit 1; fi
  for f in ${ph} ${qc} ${kp}; do if [ ! -f "$f" ]; then echo "${id}: 파일 없음 $f"; exit 1; fi; done
  if [ "${cv}" != "-" ] && [ ! -f "${cv}" ]; then echo "${id}: 파일 없음 ${cv}"; exit 1; fi
  awk -v n=${n} -v id=${id} '
    function die(m) {print id ": 입력 검사 실패 - " m > "/dev/stderr"; bad = 1; exit 1}
    function isnum(v) {return v ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/}
    FILENAME == ARGV[1] {fam[$2]; next}
    {f = FILENAME; if (!(f in nf)) nf[f] = NF; if (NF != nf[f]) die(f " 열 수 일정하지 않음 (" FNR "줄)")
     for (i = 1; i <= NF; i++) if ($i == "NA" || $i == "nan" || $i == "") die(f " 결측 (" FNR "줄 " i "열)")
     if (f == ARGV[3] || f == ARGV[4]) for (i = 3; i <= NF; i++) if (!isnum($i)) die(f " 숫자가 아닌 값 (" FNR "줄 " i "열: " $i ") — 헤더가 있거나 형식 오류")
     if (f == ARGV[2]) {if (!($2 in fam)) die("fam에 없는 개체 " $2); ord[FNR] = $2; nk = FNR}
     else if (ord[FNR] != $2) die(f " 개체 순서·구성이 keep과 다름 (" FNR "줄)")
     cnt[f] = FNR}
    END {if (bad) exit 1; if (nk != n) die("keep 개체 " nk " ≠ 실행 표 n " n)
         for (f in cnt) if (cnt[f] != nk) die(f " 개체 수 " cnt[f] " ≠ " nk)}' \
    ${BFILE}.fam ${kp} ${ph} ${qc} $([ "${cv}" != "-" ] && echo ${cv})
  {
    echo "RID=${rid}"; echo "N_EXPECT=${n}"; echo "PHENO=${ph}"; echo "QCOVAR=${qc}"; echo "KEEP=${kp}"
    if [ "${cv}" != "-" ]; then echo "COVAR=${cv}"; else echo "COVAR="; fi
  } > ${WORK_DIR}/cfg/${id}.cfg
  echo "모형 ${id}: 개체 ${n}, qcovar $(( $(awk 'NR == 1 {print NF}' ${qc}) - 2 ))개, covar $([ "${cv}" != "-" ] && echo $(( $(awk 'NR == 1 {print NF}' ${cv}) - 2 )) || echo 0)개 — 입력 검사 통과"
  N_RUN=$(( N_RUN + 1 ))
done

run_job() {
  local id=$1 c=$2 extra=${3:-}
  source ${WORK_DIR}/cfg/${id}.cfg
  local tmp=${WORK_DIR}/tmp_${id}
  local out=${tmp}/chr${c}
  mkdir -p ${tmp}
  if [ -s ${out}.mlma ] && [ -f ${WORK_DIR}/log/${id}_chr${c}.log ]; then return 0; fi
  rm -f ${out}.mlma
  local covar_opt=""; if [ -n "${COVAR}" ]; then covar_opt="--covar ${COVAR}"; fi
  ${GCTA} --bfile ${BFILE} --autosome-num 38 --chr ${c} --keep ${KEEP} ${extra} \
          --grm ${GRM_DIR}/grm_all --mlma-subtract-grm ${GRM_DIR}/grm_chr${c} --mlma \
          --pheno ${PHENO} --qcovar ${QCOVAR} ${covar_opt} \
          --thread-num ${THREADS_PER_JOB} --out ${out} > /dev/null 2>&1 \
    || { echo "GCTA 실패: ${id} chr${c} (로그: ${out}.log)"; rm -f ${out}.mlma; return 1; }
  local n_common
  n_common=$(awk '/individuals are in common/ {print $1}' ${out}.log | tail -1)
  if [ "${n_common:-0}" != "${N_EXPECT}" ]; then echo "${id} chr${c}: 분석 개체 ${n_common:-?} ≠ 예상 ${N_EXPECT} (로그: ${out}.log)"; rm -f ${out}.mlma; return 1; fi
  if ! grep -q "Log-likelihood ratio converged" ${out}.log; then echo "경고: ${id} chr${c} REML 수렴 문구 없음"; fi
  if grep -qi "constrained" ${out}.log; then echo "경고: ${id} chr${c} 분산 성분 경계값 고정(constrained)"; fi
  mv ${out}.log ${WORK_DIR}/log/${id}_chr${c}.log
  touch ${WORK_DIR}/done/${id}_chr${c}
}

merge_run() {
  local id=$1
  source ${WORK_DIR}/cfg/${id}.cfg
  local tmp=${WORK_DIR}/tmp_${id}
  local nres=0
  for c in $(seq 1 38); do if [ -s ${tmp}/chr${c}.mlma ] && [ -f ${WORK_DIR}/log/${id}_chr${c}.log ]; then nres=$(( nres + 1 )); fi; done
  if [ ${nres} -ne 38 ]; then echo "${id}: 염색체 결과 ${nres} / 38 → 합치지 않음 (RESUME=1 RUNS=\"${id}\"로 다시 실행)"; return 1; fi
  for c in $(seq 1 38); do
    if [ "$(head -1 ${tmp}/chr${c}.mlma | cut -f1-9)" != "$(printf 'Chr\tSNP\tbp\tA1\tA2\tFreq\tb\tse\tp')" ]; then echo "${id} chr${c} 결과 헤더가 예상과 다름"; return 1; fi
  done
  { head -1 ${tmp}/chr1.mlma; for c in $(seq 1 38); do tail -n +2 ${tmp}/chr${c}.mlma; done; } | gzip > ${WORK_DIR}/${id}.mlma.gz
  local N_RES N_BIM N_NAN DIFF CHR_DIFF
  N_RES=$(( $(zcat ${WORK_DIR}/${id}.mlma.gz | wc -l) - 1 ))
  N_BIM=$(wc -l < ${BFILE}.bim)
  CHR_DIFF=$(awk 'FILENAME == ARGV[1] {nb[$1]++; next} FNR > 1 {nr[FILENAME]++}
                  END {for (c = 1; c <= 38; c++) {f = ARGV[1 + c]; if (nb[c] != nr[f]) printf "chr%d %d/%d  ", c, nr[f], nb[c]}}' \
             ${BFILE}.bim $(for c in $(seq 1 38); do echo ${tmp}/chr${c}.mlma; done))
  N_NAN=$(zcat ${WORK_DIR}/${id}.mlma.gz | awk 'NR > 1 && ($9 == "nan" || $8 == "nan" || $8 + 0 <= 0)' | wc -l)
  DIFF=$(( N_BIM - N_RES ))
  if [ "${DIFF}" -lt 0 ]; then echo "${id}: 오류 - 결과 SNP가 bim보다 많음"; return 1; fi
  {
    echo "실행: ${id}, 개체 ${N_EXPECT}, qcovar $(basename ${QCOVAR}) ($(( $(awk 'NR == 1 {print NF}' ${QCOVAR}) - 2 ))개), covar $(if [ -n "${COVAR}" ]; then basename ${COVAR}; else echo 없음; fi), keep $(basename ${KEEP})"
    echo "결과 SNP ${N_RES} / bim ${N_BIM} (차이 ${DIFF}: 하위집합에서 단형성이 된 SNP 등), nan·se ≤ 0 SNP ${N_NAN}${CHR_DIFF:+, SNP 수 다른 염색체: ${CHR_DIFF}}"
    if [ "${DIFF}" -gt $(( N_BIM / 100 )) ]; then echo "경고: 결과 SNP가 bim보다 1% 넘게 적음"; fi
    zcat ${WORK_DIR}/${id}.mlma.gz | awk '$8 != "nan" && $8 + 0 > 0 && $7 != "nan" && NR > 1 {print ($7 / $8)^2}' | sort -g -S 2G \
      | awk '{a[NR] = $1} END {if (NR == 0) {print "λ 계산 불가"; exit}
                               m = (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2; l = m / 0.4549364
                               printf "λ (median χ² / 0.4549): %.4f (SNP %d)\n", l, NR
                               printf "λ 범위 판정 (기준 0.95 ~ 1.05): %s\n", (l >= 0.95 && l <= 1.05) ? "통과" : "불통과"}'
    zcat ${WORK_DIR}/${id}.mlma.gz | awk 'NR > 1 && $9 != "nan" {if ($9 < 5e-8) a++; if ($9 < 1e-6) b++; if ($9 < 1e-5) c++; if ($9 < 1e-4) d++}
         END {printf "p < 5e-8: %d / p < 1e-6: %d / p < 1e-5: %d / p < 1e-4: %d\n", a, b, c, d}'
    echo "상위 10개 SNP (Chr SNP bp A1 A2 Freq b se p):"
    zcat ${WORK_DIR}/${id}.mlma.gz | awk 'NR > 1 && $9 != "nan"' | sort -k9,9g -S 2G | head -10 || true
    echo "염색체별 V(G)/Vp:"
    for c in $(seq 1 38); do
      lg=${WORK_DIR}/log/${id}_chr${c}.log
      if grep -q "Log-likelihood ratio converged" ${lg}; then
        v=$(awk '/^Iter\./ {f = 1; next} f && $1 ~ /^[0-9]+$/ && NF >= 4 {g = $3; e = $4} END {if (g != "" && g + e > 0) printf "%.4f", g / (g + e); else printf "NA"}' ${lg})
      else v="NA(미수렴)"; fi
      printf "chr%s %s  " ${c} "${v}"
    done; echo
  } > ${WORK_DIR}/${id}_summary.txt
  rm -rf ${tmp}
  echo "합침: ${id} ($(grep '^λ (median' ${WORK_DIR}/${id}_summary.txt | cut -d' ' -f5), $(grep '^p < 5e-8' ${WORK_DIR}/${id}_summary.txt)) $(date '+%F %T')"
}
export -f run_job merge_run
export GCTA BFILE GRM_DIR WORK_DIR THREADS_PER_JOB

if [ "${CHECK}" = "1" ]; then
  MAIN_PC=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["id"] == "main" {print $h["ref_result"]}' ${RUNS_TSV})
  if [ ! -f "${MAIN_PC}" ]; then echo "주 분석 결과 없음: ${MAIN_PC}"; exit 1; fi
  mkdir -p ${WORK_DIR}/check
  awk '$1 == 38 && ++n <= 2000 {print $2}' ${BFILE}.bim > ${WORK_DIR}/check/snps.txt
  rm -rf ${WORK_DIR}/tmp_main_check; rm -f ${WORK_DIR}/log/main_check_chr38.log
  echo "CHECK: main_check chr38 SNP 2,000개 실행 $(date '+%F %T') (REML 포함 약 10~20분)"
  THREADS_PER_JOB=16 run_job main_check 38 "--extract ${WORK_DIR}/check/snps.txt"
  zcat ${MAIN_PC} | awk 'NR == FNR {s[$1]; next} ($2 in s)' ${WORK_DIR}/check/snps.txt - > ${WORK_DIR}/check/main_045.txt
  awk 'NR == FNR {b[$2] = $7; se[$2] = $8; p[$2] = $9; next}
       FNR > 1 && ($2 in b) {n++; d = $7 - b[$2]; if (d < 0) d = -d; if (d > mb) mb = d
                             d = $8 - se[$2]; if (d < 0) d = -d; if (d > ms) ms = d
                             if (p[$2] > 0 && $9 > 0) {d = log($9 / p[$2]) / log(10); if (d < 0) d = -d; if (d > mp) mp = d}}
       END {printf "CHECK 비교 SNP %d개: |b 차이| 최대 %.2e, |se 차이| 최대 %.2e, |log10 p 차이| 최대 %.2e\n", n, mb, ms, mp
            if (n < 1900 || mb > 1e-4 || ms > 1e-4 || mp > 1e-3) {print "CHECK 실패: 04-5 주 분석과 다름 → 05-1 입력 확인"; exit 1}
            print "CHECK 통과: 05-1 입력 + 05-2 실행기 = 04-5 주 분석과 같은 결과"}' \
      ${WORK_DIR}/check/main_045.txt ${WORK_DIR}/tmp_main_check/chr38.mlma | tee ${WORK_DIR}/check/check_result.txt
  rm -rf ${WORK_DIR}/tmp_main_check; rm -f ${WORK_DIR}/log/main_check_chr38.log ${WORK_DIR}/done/main_check_chr*
  exit 0
fi

if [ "${RESUME}" = "1" ]; then
  KEEP_SEL=""
  for id in ${SEL}; do
    if [ -s ${WORK_DIR}/${id}.mlma.gz ] && [ -s ${WORK_DIR}/${id}_summary.txt ]; then echo "이미 완료: ${id} (건너뜀)"; else KEEP_SEL="${KEEP_SEL} ${id}"; fi
  done
  SEL="${KEEP_SEL}"; N_RUN=$(echo ${SEL} | wc -w)
  if [ ${N_RUN} -eq 0 ]; then echo "선택한 모형이 모두 완료됨"; exit 0; fi
fi
awk '{n[$1]++} END {for (c in n) print c, n[c]}' ${BFILE}.bim > ${WORK_DIR}/chr_nsnp.txt
for id in ${SEL}; do
  n=$(grep '^N_EXPECT=' ${WORK_DIR}/cfg/${id}.cfg | cut -d= -f2)
  awk -v id=${id} -v n=${n} '{printf "%s %s %.0f\n", id, $1, n * $2}' ${WORK_DIR}/chr_nsnp.txt
done | sort -k3,3gr | cut -d' ' -f1,2 > ${WORK_DIR}/jobs.txt
N_TOTAL=$(wc -l < ${WORK_DIR}/jobs.txt)
echo "작업: 모형 ${N_RUN}개 × 염색체 38 = ${N_TOTAL}개, 동시 ${N_JOBS}개 × 스레드 ${THREADS_PER_JOB}"
if [ "${DRYRUN}" = "1" ]; then echo "DRYRUN=1: 입력 검사·작업 목록(${WORK_DIR}/jobs.txt)만 만들고 종료"; exit 0; fi

if [ "${RESUME}" != "1" ]; then
  for id in ${SEL}; do rm -rf ${WORK_DIR}/tmp_${id}; rm -f ${WORK_DIR}/log/${id}_chr*.log ${WORK_DIR}/${id}.mlma.gz ${WORK_DIR}/${id}_summary.txt ${WORK_DIR}/done/${id}_chr*; done
fi
echo "시작 $(date '+%F %T')"

merge_watch() {
  local merged=" "
  while true; do
    for id in ${SEL}; do
      case "${merged}" in *" ${id} "*) continue ;; esac
      shopt -s nullglob; local marks=(${WORK_DIR}/done/${id}_chr*); shopt -u nullglob
      if [ ${#marks[@]} -eq 38 ]; then
        if merge_run ${id}; then merged="${merged}${id} "; fi
      fi
    done
    if [ -f ${WORK_DIR}/.jobs_finished ]; then break; fi
    sleep 20
  done
}
rm -f ${WORK_DIR}/.jobs_finished
merge_watch &
WATCH_PID=$!
JOB_STATUS=0
xargs -P ${N_JOBS} -L 1 bash -c 'run_job "$0" "$1"' < ${WORK_DIR}/jobs.txt || JOB_STATUS=$?
touch ${WORK_DIR}/.jobs_finished
wait ${WATCH_PID} || true
rm -f ${WORK_DIR}/.jobs_finished

FAILED=""
for id in ${SEL}; do
  if [ ! -f ${WORK_DIR}/${id}.mlma.gz ]; then
    if [ -d ${WORK_DIR}/tmp_${id} ] && merge_run ${id}; then continue; fi
    FAILED="${FAILED} ${id}"
  fi
done
echo "끝 $(date '+%F %T') / 작업 종료 상태 ${JOB_STATUS}"
if [ -n "${FAILED}" ]; then
  echo "완료되지 않은 모형:${FAILED} → 위 메시지(GCTA 실패·개체 수) 확인 후"
  echo "  RESUME=1 RUNS=\"${FAILED# }\" bash $(basename $0)  로 끝나지 않은 염색체만 다시 실행"
  exit 1
fi
echo "모든 모형 완료 → 05-4로 비교"
