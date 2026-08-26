#!/bin/bash
set -euo pipefail

# ============================================================
# [수정 필요] 아래 변수들을 본인 환경에 맞게 확인/수정하세요
# ============================================================

PLINK_PATH="$HOME/tools/plink/plink"

# BFILE(원본 bed/bim/fam)이 있는 경로
BFILE="/mnt/j/01_lab_lee/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6"
BIM="${BFILE}.bim"

# 10_Breed_Keep_Files.R에서 생성한 keep_*.txt 파일들이 있는 폴더
KEEP_DIR="/mnt/j/01_lab_lee/03_GWAS_Pure/14_Regression_local/breed_keep"

# 결과가 저장될 경로 (breed별로 하위 폴더 생성)
OUTDIR="/mnt/j/01_lab_lee/03_GWAS_Pure/14_Regression_local/09_BreedStratified_LD/09_LD_breed_stratified"
mkdir -p "$OUTDIR"

# 로컬 CPU(2코어/4스레드) 기준
THREADS=3
MEMORY_MB=12000

#############################################
# INDEX SNP 목록 — CHR:BP:A1:A2 형식
# circadian entrainment 관련 5개 SNP
#############################################
INDEX_SNPS=(
  "5:25529184:C:T"
  "24:12787767:A:G"
  "24:12787779:T:C"
  "27:44543118:A:T"
  "26:13787793:A:G"
)

# LD window 설정
LD_WINDOW_KB=1000
LD_WINDOW_R2_MIN=0.0

echo "▶ BFILE: $BFILE"
echo "▶ KEEP_DIR: $KEEP_DIR"
echo "▶ WINDOW: ${LD_WINDOW_KB}kb, r² ≥ ${LD_WINDOW_R2_MIN}"
echo "▶ OUTDIR: $OUTDIR"

# ------------------------------------------------------------------------------
# breed_keep 폴더의 keep_*.txt 파일들을 자동으로 스캔
# ------------------------------------------------------------------------------
shopt -s nullglob
KEEP_FILES=("${KEEP_DIR}"/keep_*.txt)
shopt -u nullglob

if [[ ${#KEEP_FILES[@]} -eq 0 ]]; then
  echo "❌ ${KEEP_DIR}에서 keep_*.txt 파일을 찾지 못했습니다."
  echo "   10_Breed_Keep_Files.R을 먼저 실행했는지 확인해주세요."
  exit 1
fi

echo "▶ 발견된 breed keep 파일: ${#KEEP_FILES[@]}개"
for kf in "${KEEP_FILES[@]}"; do
  echo "   - $(basename "$kf")"
done

TMP_SNP_LIST="$(mktemp)"
trap 'rm -f "$TMP_SNP_LIST"' EXIT

RUN_LOG="${OUTDIR}/breed_stratified_run.log"
echo "=== breed-stratified LD 실행 시작: $(date) ===" > "$RUN_LOG"

# ------------------------------------------------------------------------------
# breed x index SNP 이중 루프
# ------------------------------------------------------------------------------
for KEEP_FILE in "${KEEP_FILES[@]}"; do

  # 파일명에서 breed명 추출: keep_Labrador_Retriever.txt -> Labrador_Retriever
  BREED=$(basename "$KEEP_FILE" .txt)
  BREED=${BREED#keep_}

  N_INDIV=$(wc -l < "$KEEP_FILE")

  echo "==============================================="
  echo "▶ BREED: ${BREED} (N=${N_INDIV})"
  echo "==============================================="

  BREED_OUTDIR="${OUTDIR}/${BREED}"
  mkdir -p "$BREED_OUTDIR"

  for SNP in "${INDEX_SNPS[@]}"; do
    IFS=':' read -r CHR BP A1 A2 <<< "$SNP"

    OUT_PREFIX="${BREED_OUTDIR}/09_LD_${CHR}_${BP}_${A1}_${A2}_${BREED}"

    echo "-----------------------------------------------"
    echo "▶ INDEX SNP: chr=${CHR}, bp=${BP}, alleles=${A1}/${A2}, breed=${BREED}"

    IDX_ID=$(awk -v c="$CHR" -v p="$BP" -v a1="$A1" -v a2="$A2" '
        ($1==c && $4==p) &&
        (($5==a1 && $6==a2) || ($5==a2 && $6==a1)) {
          print $2;
          exit;
        }' "$BIM")

    if [[ -z "$IDX_ID" ]]; then
      echo "❌ 정확히 일치하는 allele 조합의 SNP을 BIM 파일에서 찾지 못했습니다!"
      echo "   BIM에서 직접 확인하세요: CHR=$CHR, BP=$BP, Alleles=$A1/$A2"
      continue
    fi

    echo "✅ 선택된 SNP ID: ${IDX_ID}"
    printf "%s\n" "$IDX_ID" > "$TMP_SNP_LIST"

    "$PLINK_PATH" \
      --bfile "$BFILE" \
      --keep "$KEEP_FILE" \
      --dog \
      --r2 dprime with-freqs \
      --ld-snp-list "$TMP_SNP_LIST" \
      --ld-window-kb "$LD_WINDOW_KB" \
      --ld-window 999999 \
      --ld-window-r2 "$LD_WINDOW_R2_MIN" \
      --threads "$THREADS" \
      --memory "$MEMORY_MB" \
      --out "$OUT_PREFIX" \
      || { echo "⚠️ ${BREED} / ${SNP} 실행 중 오류 (표본 부족 등) - 건너뜀"; continue; }

    gzip -f "${OUT_PREFIX}.ld"
    echo "✅ 완료: ${OUT_PREFIX}.ld.gz"
    echo "${BREED},${SNP},${N_INDIV}" >> "$RUN_LOG"

  done
done

echo "=== breed-stratified LD 실행 종료: $(date) ===" >> "$RUN_LOG"
echo "🎉 모든 breed x SNP 조합 처리 완료. 로그: ${RUN_LOG}"