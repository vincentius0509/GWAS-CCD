#!/bin/bash
set -euo pipefail

# ============================================================
# [수정 필요] 아래 변수들을 본인 환경에 맞게 확인/수정하세요
# ============================================================

PLINK2_PATH="$HOME/tools/plink2/plink2"

# 외장하드(J드라이브) 기준 데이터 경로
DATA_ROOT="/mnt/j/01_lab_lee/03_GWAS_Pure"

# covariate 파일 및 regression 출력 파일이 저장되는 경로
OUT_DIR="${DATA_ROOT}/14_Regression_local"
COVAR_DIR="${OUT_DIR}"

# 로컬 CPU(2코어/4스레드) 기준
THREADS=3
MEMORY_MB=12000

# ============================================================
# 고정 조건 (MAF 0.1, HWE 1e-6)
# ============================================================

MAF="01"
HWE="1e-6"
HWE_TAG="1e6"
BFILE="${DATA_ROOT}/03_MAF_0.1_HWE_1e-6"
PHENO_FILE="${DATA_ROOT}/04_pheno01.txt"

PCA_DIMENSIONS=(5 10 15 20 25 30)

mkdir -p "$OUT_DIR"

RUN_LOG="${OUT_DIR}/local_run_timing_01_1e6.log"
echo "=== 로컬 실행 시작: $(date) ===" > "$RUN_LOG"

for PC in "${PCA_DIMENSIONS[@]}"; do
  COVAR_FILE="${COVAR_DIR}/04_cov_${MAF}_${HWE_TAG}_${PC}.txt"

  if [[ ! -f "$COVAR_FILE" ]]; then
    echo "❌ COVAR 파일을 찾을 수 없습니다: $COVAR_FILE"
    echo "   COVAR_DIR 경로를 확인해주세요."
    exit 1
  fi

  COVAR_NAMES=$(head -1 "$COVAR_FILE" | tr '\t' '\n' | grep -v -E '^FID$|^IID$' | paste -sd "," -)

  OUT_PREFIX_LINEAR="${OUT_DIR}/05_GWAS_PHE1_${MAF}_${HWE_TAG}_${PC}_linear"

  echo "▶ MAF=${MAF}, HWE=${HWE}, PC=${PC}"
  echo "📂 BFILE = $BFILE"
  echo "📂 COVAR = $COVAR_FILE"
  echo "📂 PHENO = $PHENO_FILE"
  echo "📂 Covariates: $COVAR_NAMES"

  RUN_START=$(date +%s)

  "$PLINK2_PATH" \
    --bfile "$BFILE" \
    --pheno "$PHENO_FILE" \
    --pheno-name phe1 \
    --covar "$COVAR_FILE" \
    --covar-name $COVAR_NAMES \
    --covar-variance-standardize \
    --dog \
    --glm no-x-sex cols=+a1freq hide-covar \
    --out "$OUT_PREFIX_LINEAR" \
    --threads "$THREADS" \
    --memory "$MEMORY_MB"

  RUN_END=$(date +%s)
  RUN_ELAPSED=$((RUN_END - RUN_START))
  echo "PC=${PC} : $((RUN_ELAPSED / 60))분 $((RUN_ELAPSED % 60))초" | tee -a "$RUN_LOG"

done

echo "=== 로컬 실행 종료: $(date) ===" >> "$RUN_LOG"
echo "🎉 MAF=0.1, HWE=1e-6, PC 6종 전체 완료. 소요시간은 ${RUN_LOG} 참고"