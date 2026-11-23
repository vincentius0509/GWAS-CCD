#!/bin/bash
# 13-6_int_wt_REGENIE.sh - Genome-wide SNP x body weight interaction GWAS (REGENIE). See README.md.
set -euo pipefail

REGENIE=${REGENIE:-/data/00_tools/regenie/regenie}
D0=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023
OUT_DIR=/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023
BFILE=${D0}/03-3_maf0.01_hwe
PRUNE=${D0}/03-4_pca/prune.prune.in
PHENO=${D0}/05-1_inputs/int_age_regenie_pheno.txt
COV_IN=${D0}/05-1_inputs/int_age_regenie_covar.txt
KEEP=${D0}/05-1_inputs/int_age_keep.txt
CATV=Sex,Intact,Vision,Cohort_BHS
IVAR=log2_wt
N_EXPECT=7115
TEST=${TEST:-0}
ID=int_wt; if [ "${TEST}" = "1" ]; then ID=int_wt_test; fi
WORK_DIR=${OUT_DIR}/13-6_regenie
N_JOBS=${N_JOBS:-38}
STEP1_THREADS=${STEP1_THREADS:-32}
RESUME=${RESUME:-0}
DRYRUN=${DRYRUN:-0}
mkdir -p ${WORK_DIR}/step1 ${WORK_DIR}/log ${WORK_DIR}/tmp
for f in ${BFILE}.bed ${BFILE}.bim ${BFILE}.fam ${PRUNE} ${PHENO} ${COV_IN} ${KEEP}; do if [ ! -f "$f" ]; then echo "파일 없음: $f"; exit 1; fi; done
if [ ! -x "${REGENIE}" ]; then echo "REGENIE 실행 파일 없음: ${REGENIE}"; exit 1; fi
echo "REGENIE: $(${REGENIE} --version 2>&1 | awk 'NR == 1' || true) · 실행 ${ID} $(date '+%F %T')"

COVARF=${WORK_DIR}/${ID}_covar.txt
awk 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "wt_cov") w = i; if (!w) {print "wt_cov 열 없음" > "/dev/stderr"; exit 1}; print $0, "log2_wt"; next}
     {if ($w <= 0) {print "체중 ≤ 0: " $2 > "/dev/stderr"; exit 1}; printf "%s %.6f\n", $0, log($w) / log(2)}' ${COV_IN} > ${COVARF}
awk -v n=${N_EXPECT} '
  function die(m) {print "입력 검사 실패 - " m > "/dev/stderr"; bad = 1; exit 1}
  FILENAME == ARGV[1] {fam[$2]; next}
  FILENAME == ARGV[2] {if (!($2 in fam)) die("fam에 없는 개체 " $2); ord[FNR] = $2; nk = FNR; next}
  FNR == 1 {nf[FILENAME] = NF; if ($1 != "FID" || $2 != "IID") die(FILENAME " 헤더"); next}
  {if (NF != nf[FILENAME]) die(FILENAME " 열 수"); if (ord[FNR - 1] != $2) die(FILENAME " 개체 순서가 keep과 다름 (" FNR "줄)")
   for (i = 3; i <= NF; i++) if ($i == "NA" || $i == "") die(FILENAME " 결측 (" FNR "줄)"); cnt[FILENAME] = FNR - 1}
  END {if (bad) exit 1; if (nk != n) die("keep 개체 " nk " ≠ " n); for (f in cnt) if (cnt[f] != n) die(f " 개체 수 " cnt[f])}' ${BFILE}.fam ${KEEP} ${PHENO} ${COVARF}
PNAME=$(head -1 ${PHENO} | awk '{print $3}')
echo "입력 검사 통과: ${N_EXPECT}마리, 표현형 ${PNAME}, 상호작용 ${IVAR} (= log2 wt_cov), 범주형 ${CATV}"
awk 'NR > 1 {v = $NF; s += v; ss += v * v; if (NR == 2 || v < mn) mn = v; if (NR == 2 || v > mx) mx = v} END {m = s / (NR - 1); printf "log2_wt: 평균 %.2f (%.1f kg) · SD %.2f · 범위 %.2f ~ %.2f\n", m, 2 ^ m, sqrt(ss / (NR - 1) - m * m), mn, mx}' ${COVARF}
if [ "${DRYRUN}" = "1" ]; then echo "DRYRUN=1: 입력 검사만 하고 종료"; exit 0; fi
COMMON="--bed ${BFILE} --keep ${KEEP} --phenoFile ${PHENO} --phenoCol ${PNAME} --covarFile ${COVARF} --catCovarList ${CATV} --nauto 38"

EXTRACT1=${PRUNE}
if [ "${TEST}" = "1" ]; then awk 'NR == FNR {s[$1]; next} ($1 >= 35) && ($2 in s) {print $2}' ${PRUNE} ${BFILE}.bim > ${WORK_DIR}/tmp/prune_test.txt; EXTRACT1=${WORK_DIR}/tmp/prune_test.txt; fi
if [ "${RESUME}" = "1" ] && [ -s ${WORK_DIR}/step1/${ID}_pred.list ]; then echo "step 1 건너뜀 (완료)"
else
  rm -f ${WORK_DIR}/step1/${ID}_*
  echo "step 1 시작 $(date '+%F %T') (SNP $(wc -l < ${EXTRACT1}))"
  ${REGENIE} --step 1 ${COMMON} --extract ${EXTRACT1} --bsize 1000 --lowmem --lowmem-prefix ${WORK_DIR}/tmp/${ID}_l0 \
             --threads ${STEP1_THREADS} --out ${WORK_DIR}/step1/${ID} > ${WORK_DIR}/log/${ID}_step1.out 2>&1 \
    || { echo "step 1 실패 (로그: ${WORK_DIR}/log/${ID}_step1.out)"; exit 1; }
  if [ ! -s ${WORK_DIR}/step1/${ID}_pred.list ] || grep -qE "ERROR" ${WORK_DIR}/log/${ID}_step1.out; then echo "step 1 결과 없음 또는 오류 문구"; exit 1; fi
  echo "step 1 완료 $(date '+%F %T')"
fi

step2() {
  local c=$1 out=${WORK_DIR}/tmp/${ID}_chr$1
  if [ -s ${out}_${PNAME}.regenie ] && [ -f ${WORK_DIR}/log/${ID}_chr${c}.log ]; then return 0; fi
  rm -f ${out}_*
  ${REGENIE} --step 2 ${COMMON} --chr ${c} --pred ${WORK_DIR}/step1/${ID}_pred.list --bsize 400 --interaction ${IVAR} \
             --threads 2 --out ${out} > /dev/null 2>&1 || { echo "step 2 실패: chr${c} (로그: ${out}.log)"; return 1; }
  if [ ! -s ${out}_${PNAME}.regenie ]; then echo "step 2 결과 없음: chr${c}"; return 1; fi
  local nused; nused=$(awk '/number of individuals used in analysis/ {n = $NF} END {print n}' ${out}.log)
  if [ "${nused:-0}" -lt 7100 ]; then echo "chr${c}: 분석 개체 ${nused:-?} (기대 ${N_EXPECT})"; return 1; fi
  if [ "${nused}" != "${N_EXPECT}" ]; then echo "참고 chr${c}: 분석 개체 ${nused} (기대 ${N_EXPECT})"; fi
  mv ${out}.log ${WORK_DIR}/log/${ID}_chr${c}.log
}
export -f step2; export REGENIE WORK_DIR ID PNAME COMMON IVAR N_EXPECT
CHRS=$(awk '{n[$1]++} END {for (c in n) print n[c], c}' ${BFILE}.bim | sort -k1,1gr | awk '{print $2}')
if [ "${TEST}" = "1" ]; then CHRS=38; fi
if [ "${RESUME}" != "1" ]; then rm -f ${WORK_DIR}/tmp/${ID}_chr* ${WORK_DIR}/log/${ID}_chr*.log; fi
NCHR=$(echo ${CHRS} | wc -w)
echo "step 2 시작: 염색체 ${NCHR}개, 동시 ${N_JOBS}개 $(date '+%F %T')"
STATUS=0
echo ${CHRS} | tr ' ' '\n' | xargs -P ${N_JOBS} -I{} bash -c 'step2 {}' || STATUS=$?

ndone=0; for c in ${CHRS}; do if [ -s ${WORK_DIR}/tmp/${ID}_chr${c}_${PNAME}.regenie ] && [ -f ${WORK_DIR}/log/${ID}_chr${c}.log ]; then ndone=$(( ndone + 1 )); fi; done
if [ ${ndone} -ne ${NCHR} ] || [ ${STATUS} -ne 0 ]; then echo "염색체 결과 ${ndone} / ${NCHR} → 메시지 확인 후 RESUME=1로 다시 실행"; exit 1; fi
FIRST=$(echo ${CHRS} | awk '{print $1}')
{ head -1 ${WORK_DIR}/tmp/${ID}_chr${FIRST}_${PNAME}.regenie; for c in $(echo ${CHRS} | tr ' ' '\n' | sort -n); do tail -n +2 ${WORK_DIR}/tmp/${ID}_chr${c}_${PNAME}.regenie; done; } | gzip > ${WORK_DIR}/${ID}.regenie.gz
KP=$(awk 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "LOG10P") print i; exit}' ${WORK_DIR}/tmp/${ID}_chr${FIRST}_${PNAME}.regenie)
{
  echo "실행: ${ID} (REGENIE, CSLB_RINT × ${IVAR} 상호작용), 개체 ${N_EXPECT}, 공변량 $(basename ${COVARF})"
  for t in $(zcat ${WORK_DIR}/${ID}.regenie.gz | awk 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "TEST") k = i; next} {print $k}' | sort -u); do
    zcat ${WORK_DIR}/${ID}.regenie.gz | awk -v t=${t} 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["TEST"] == t && $h["CHISQ"] != "NA" {print $h["CHISQ"]}' | sort -g -S 2G \
      | awk -v t=${t} '{a[NR] = $1} END {if (NR == 0) {print t ": 유효 SNP 없음"; exit}
           m = (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2; d0 = (t ~ /2DF/) ? 1.386294 : 0.4549364   # 2-자유도 검정은 χ²(2) 중앙값
           printf "TEST %s: λ %.4f (SNP %d)\n", t, m / d0, NR}'
    zcat ${WORK_DIR}/${ID}.regenie.gz | awk -v t=${t} 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
      $h["TEST"] == t && $h["LOG10P"] != "NA" {l = $h["LOG10P"]; if (l > 7.30103) a++; if (l > 5) b++; if (l > 4) c++}
      END {printf "  p < 5e-8: %d / p < 1e-5: %d / p < 1e-4: %d\n", a, b, c}'
    echo "  상위 10개:"
    zcat ${WORK_DIR}/${ID}.regenie.gz | awk -v t=${t} 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["TEST"] == t && $h["LOG10P"] != "NA"' \
      | sort -k${KP},${KP}gr -S 2G | awk 'NR <= 10'
  done
} > ${WORK_DIR}/${ID}_summary.txt
rm -f ${WORK_DIR}/tmp/${ID}_chr*
cat ${WORK_DIR}/${ID}_summary.txt | grep -E "^TEST|p < 5e-8"
echo "완료 ${ID} $(date '+%F %T') → 13-7 (MAGMA)"
