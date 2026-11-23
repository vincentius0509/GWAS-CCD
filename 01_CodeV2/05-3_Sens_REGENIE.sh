#!/bin/bash
# 05-3_Sens_REGENIE.sh - REGENIE models: SNP x age interaction and CSLB >= 50 case-control (all dogs and dogs >= 8 years). See README.md.
set -euo pipefail

#   DRYRUN=1 bash 00_Code/05_GWAS_v2023/05-3_Sens_REGENIE.sh
#   nohup bash 00_Code/05_GWAS_v2023/05-3_Sens_REGENIE.sh > 04_GWAS_2023/05-3_run.log 2>&1 &
#   step1/, log/

REGENIE=${REGENIE:-/data/00_tools/regenie/regenie}
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023
RUNS_TSV=${OUT_DIR}/05-1_runs.tsv
BFILE=${OUT_DIR}/03-3_maf0.01_hwe
PRUNE=${OUT_DIR}/03-4_pca/prune.prune.in
WORK_DIR=${OUT_DIR}/05-3_regenie
RUNS=${RUNS:-"bin_all bin_age8 int_age"}
N_JOBS=${N_JOBS:-64}
STEP1_THREADS=${STEP1_THREADS:-32}
DRYRUN=${DRYRUN:-0}
RESUME=${RESUME:-0}
mkdir -p ${WORK_DIR}/step1 ${WORK_DIR}/log ${WORK_DIR}/tmp ${WORK_DIR}/cfg

for f in ${RUNS_TSV} ${BFILE}.bed ${BFILE}.bim ${BFILE}.fam ${PRUNE}; do if [ ! -f "$f" ]; then echo "파일 없음: $f"; exit 1; fi; done
if [ ! -x "${REGENIE}" ]; then
  echo "REGENIE 실행 파일 없음: ${REGENIE} (머리말의 설치 방법 참고, 경로가 다르면 REGENIE=경로 를 앞에 붙여 실행)"
  if [ "${DRYRUN}" != "1" ]; then exit 1; fi
else
  echo "REGENIE: $(${REGENIE} --version 2>&1 | awk 'NR == 1' || true)"
fi
N_PRUNE=$(awk 'NR == FNR {s[$1]; next} ($2 in s)' ${PRUNE} ${BFILE}.bim | wc -l)
echo "step 1 SNP (pruned, bed에 있는 것): ${N_PRUNE}"

for id in ${RUNS}; do
  line=$(awk -F'\t' -v id=${id} 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
           $h["id"] == id {print $h["engine"] "\t" $h["n"] "\t" $h["pheno_file"] "\t" $h["qcovar_file"] "\t" $h["keep_file"] "\t" $h["covar"] "\t" $h["n_case"]}' ${RUNS_TSV})
  if [ -z "${line}" ]; then echo "실행 표에 없는 id: ${id}"; exit 1; fi
  IFS=$'\t' read -r eng n ph cvf kp catv ncase <<< "${line}"
  case "${eng}" in regenie_bt|regenie_int) ;; *) echo "${id}: engine ${eng} → 05-3 대상 아님"; exit 1 ;; esac
  for f in ${ph} ${cvf} ${kp}; do if [ ! -f "$f" ]; then echo "${id}: 파일 없음 $f"; exit 1; fi; done
  awk -v n=${n} -v id=${id} -v bt=$([ "${eng}" = "regenie_bt" ] && echo 1 || echo 0) '
    function die(m) {print id ": 입력 검사 실패 - " m > "/dev/stderr"; bad = 1; exit 1}
    FILENAME == ARGV[1] {fam[$2]; next}
    FILENAME == ARGV[2] {if (!($2 in fam)) die("fam에 없는 개체 " $2); ord[FNR] = $2; nk = FNR; next}
    FNR == 1 {nf[FILENAME] = NF; if ($1 != "FID" || $2 != "IID") die(FILENAME " 헤더가 FID IID로 시작하지 않음"); next}
    {if (NF != nf[FILENAME]) die(FILENAME " 열 수 일정하지 않음"); if (ord[FNR - 1] != $2) die(FILENAME " 개체 순서·구성이 keep과 다름 (" FNR "줄)")
     for (i = 3; i <= NF; i++) if ($i == "NA" || $i == "") die(FILENAME " 결측 (" FNR "줄)")
     if (FILENAME == ARGV[3] && bt && $3 != "0" && $3 != "1") die("이분형 표현형이 0/1이 아님: " $3)
     if (FILENAME == ARGV[3] && bt) cs += $3
     cnt[FILENAME] = FNR - 1}
    END {if (bad) exit 1; if (nk != n) die("keep 개체 " nk " ≠ 실행 표 n " n); for (f in cnt) if (cnt[f] != n) die(f " 개체 수 " cnt[f] " ≠ " n)
         if (bt) printf "%d", cs}' ${BFILE}.fam ${kp} ${ph} ${cvf} > ${WORK_DIR}/cfg/${id}.ncase
  pname=$(head -1 ${ph} | awk '{print $3}')
  {
    echo "ENG=${eng}"; echo "N_EXPECT=${n}"; echo "PHENO=${ph}"; echo "COVARF=${cvf}"; echo "KEEP=${kp}"; echo "CATV=${catv}"; echo "PNAME=${pname}"
  } > ${WORK_DIR}/cfg/${id}.cfg
  if [ "${eng}" = "regenie_bt" ]; then
    echo "모형 ${id}: 이분형 ${pname}, 개체 ${n}, case $(cat ${WORK_DIR}/cfg/${id}.ncase) (실행 표 ${ncase}), 범주형 공변량 ${catv} — 입력 검사 통과"
  else
    if ! awk 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "age_cslb") f = 1; exit} END {exit !f}' ${cvf}; then echo "${id}: 공변량 파일에 age_cslb 없음 (상호작용 변수)"; exit 1; fi
    echo "모형 ${id}: 연속형 ${pname} × age_cslb 상호작용, 개체 ${n}, 범주형 공변량 ${catv} — 입력 검사 통과"
  fi
done
if [ "${DRYRUN}" = "1" ]; then echo "DRYRUN=1: 입력 검사만 하고 종료"; exit 0; fi

common_opts() {
  local o="--bed ${BFILE} --keep ${KEEP} --phenoFile ${PHENO} --phenoCol ${PNAME} --covarFile ${COVARF} --nauto 38"
  if [ -n "${CATV}" ]; then o="${o} --catCovarList ${CATV}"; fi
  if [ "${ENG}" = "regenie_bt" ]; then o="${o} --bt"; fi
  echo "${o}"
}

step1() {
  local id=$1
  source ${WORK_DIR}/cfg/${id}.cfg
  if [ "${RESUME}" = "1" ] && [ -s ${WORK_DIR}/step1/${id}_pred.list ]; then echo "step 1 건너뜀 (완료): ${id}"; return 0; fi
  rm -f ${WORK_DIR}/step1/${id}_*
  ${REGENIE} --step 1 $(common_opts) --extract ${PRUNE} --bsize 1000 --lowmem --lowmem-prefix ${WORK_DIR}/tmp/${id}_l0 \
             --threads ${STEP1_THREADS} --out ${WORK_DIR}/step1/${id} > ${WORK_DIR}/log/${id}_step1.out 2>&1 \
    || { echo "step 1 실패: ${id} (로그: ${WORK_DIR}/step1/${id}.log, ${WORK_DIR}/log/${id}_step1.out)"; return 1; }
  if [ ! -s ${WORK_DIR}/step1/${id}_pred.list ]; then echo "step 1 결과 없음: ${id}"; return 1; fi
  if grep -qE "ERROR" ${WORK_DIR}/log/${id}_step1.out ${WORK_DIR}/step1/${id}.log 2>/dev/null; then
    echo "step 1 로그에 오류 문구: ${id} (로그: ${WORK_DIR}/step1/${id}.log)"; return 1
  fi
  echo "step 1 완료: ${id} $(date '+%F %T')"
}

step2() {
  local id=$1 c=$2
  source ${WORK_DIR}/cfg/${id}.cfg
  local out=${WORK_DIR}/tmp/${id}_chr${c}
  if [ -s ${out}_${PNAME}.regenie ] && [ -f ${WORK_DIR}/log/${id}_chr${c}.log ]; then return 0; fi
  rm -f ${out}_*
  local extra="--bsize 400"
  if [ "${ENG}" = "regenie_bt" ]; then extra="${extra} --firth --approx --pThresh 0.01"; fi
  if [ "${ENG}" = "regenie_int" ]; then extra="${extra} --interaction age_cslb"; fi
  ${REGENIE} --step 2 $(common_opts) --chr ${c} --pred ${WORK_DIR}/step1/${id}_pred.list ${extra} \
             --threads 2 --out ${out} > /dev/null 2>&1 \
    || { echo "step 2 실패: ${id} chr${c} (로그: ${out}.log)"; return 1; }
  if [ ! -s ${out}_${PNAME}.regenie ]; then echo "step 2 결과 없음: ${id} chr${c}"; return 1; fi
  local nused
  nused=$(awk '/number of individuals used in analysis/ {n = $NF} END {print n}' ${out}.log)
  if [ "${nused:-0}" != "${N_EXPECT}" ]; then echo "${id} chr${c}: 분석 개체 ${nused:-?} ≠ 예상 ${N_EXPECT} (로그: ${out}.log)"; return 1; fi
  mv ${out}.log ${WORK_DIR}/log/${id}_chr${c}.log
}
export -f step1 step2 common_opts
export REGENIE BFILE PRUNE WORK_DIR STEP1_THREADS RESUME

echo "step 1 시작 $(date '+%F %T')"
if ! echo ${RUNS} | tr ' ' '\n' | xargs -P 3 -I{} bash -c 'step1 {}'; then echo "step 1 일부 실패 → 위 메시지 확인 후 RESUME=1로 다시 실행"; exit 1; fi

awk '{n[$1]++} END {for (c in n) print c, n[c]}' ${BFILE}.bim > ${WORK_DIR}/chr_nsnp.txt
for id in ${RUNS}; do
  source ${WORK_DIR}/cfg/${id}.cfg
  w=1; if [ "${ENG}" = "regenie_int" ]; then w=3; fi
  awk -v id=${id} -v n=${N_EXPECT} -v w=${w} '{printf "%s %s %.0f\n", id, $1, $2 * n * w}' ${WORK_DIR}/chr_nsnp.txt
done | sort -k3,3gr | cut -d' ' -f1,2 > ${WORK_DIR}/jobs.txt
if [ "${RESUME}" != "1" ]; then for id in ${RUNS}; do rm -f ${WORK_DIR}/tmp/${id}_chr* ${WORK_DIR}/log/${id}_chr*.log; done; fi
echo "step 2 시작: 작업 $(wc -l < ${WORK_DIR}/jobs.txt)개, 동시 ${N_JOBS}개 $(date '+%F %T')"
STATUS=0
xargs -P ${N_JOBS} -L 1 bash -c 'step2 "$0" "$1"' < ${WORK_DIR}/jobs.txt || STATUS=$?

FAILED=""
for id in ${RUNS}; do
  source ${WORK_DIR}/cfg/${id}.cfg
  ndone=0; for c in $(seq 1 38); do if [ -s ${WORK_DIR}/tmp/${id}_chr${c}_${PNAME}.regenie ] && [ -f ${WORK_DIR}/log/${id}_chr${c}.log ]; then ndone=$(( ndone + 1 )); fi; done
  if [ ${ndone} -ne 38 ]; then FAILED="${FAILED} ${id}"; echo "${id}: 염색체 결과 ${ndone} / 38"; continue; fi
  { head -1 ${WORK_DIR}/tmp/${id}_chr1_${PNAME}.regenie; for c in $(seq 1 38); do tail -n +2 ${WORK_DIR}/tmp/${id}_chr${c}_${PNAME}.regenie; done; } \
    | gzip > ${WORK_DIR}/${id}.regenie.gz
  KP=$(awk 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "LOG10P") print i; exit}' ${WORK_DIR}/tmp/${id}_chr1_${PNAME}.regenie)
  {
    echo "실행: ${id} (${ENG}), 표현형 ${PNAME}, 개체 ${N_EXPECT}, 공변량 $(basename ${COVARF})"
    for t in $(zcat ${WORK_DIR}/${id}.regenie.gz | awk 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "TEST") k = i; next} {print $k}' | sort -u); do
      zcat ${WORK_DIR}/${id}.regenie.gz | awk -v t=${t} 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
        $h["TEST"] == t && $h["CHISQ"] != "NA" {print $h["CHISQ"]}' | sort -g -S 2G \
        | awk -v t=${t} '{a[NR] = $1} END {if (NR == 0) {print t ": 유효 SNP 없음"; exit}
             m = (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2; printf "TEST %s: λ %.4f (SNP %d)\n", t, m / 0.4549364, NR}'
      zcat ${WORK_DIR}/${id}.regenie.gz | awk -v t=${t} 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
        $h["TEST"] == t && $h["LOG10P"] != "NA" {l = $h["LOG10P"]; if (l > 7.30103) a++; if (l > 5) b++; if (l > 4) c++}
        END {printf "  p < 5e-8: %d / p < 1e-5: %d / p < 1e-4: %d\n", a, b, c}'
      echo "  상위 10개 (CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P):"
      zcat ${WORK_DIR}/${id}.regenie.gz | awk -v t=${t} 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["TEST"] == t && $h["LOG10P"] != "NA"' \
        | sort -k${KP},${KP}gr -S 2G | awk 'NR <= 10'
    done
  } > ${WORK_DIR}/${id}_summary.txt
  rm -f ${WORK_DIR}/tmp/${id}_chr*
  echo "합침: ${id} $(date '+%F %T')"; grep "^TEST" ${WORK_DIR}/${id}_summary.txt
done
if [ -n "${FAILED}" ] || [ ${STATUS} -ne 0 ]; then
  echo "완료되지 않은 모형:${FAILED} → 위 메시지 확인 후 RESUME=1 RUNS=\"${FAILED# }\" bash $(basename $0)"
  exit 1
fi
echo "REGENIE 모형 완료 $(date '+%F %T') → 05-4로 비교"
