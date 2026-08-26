#!/bin/bash
set -euo pipefail

# ============================================================
# [수정 필요] 아래 변수들을 본인 환경에 맞게 확인/수정하세요
# ============================================================

PLINK_PATH="$HOME/tools/plink/plink"

# BFILE(원본 bed/bim/fam)이 있는 경로
BFILE="/mnt/j/01_lab_lee/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6"
BIM="${BFILE}.bim"

# 결과가 저장될 경로
OUTDIR="/mnt/j/01_lab_lee/03_GWAS_Pure/14_Regression_local"
mkdir -p "$OUTDIR"

# 로컬 CPU(2코어/4스레드) 기준
THREADS=3
MEMORY_MB=12000

#############################################
# INDEX SNP 목록 — CHR:BP:A1:A2 형식
# circadian entrainment 관련 5개 SNP
#   5:25529184  (GUCY1A2)
#   24:12787767 (PLCB4)
#   24:12787779 (PLCB4, 2nd variant)
#   27:44543118 (CACNA1C)
#   26:13787793 (NOS1)
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
echo "▶ WINDOW: ${LD_WINDOW_KB}kb, r² ≥ ${LD_WINDOW_R2_MIN}"
echo "▶ OUTDIR: $OUTDIR"

TMP_SNP_LIST="$(mktemp)"
trap 'rm -f "$TMP_SNP_LIST"' EXIT


###############################################################
# 예: 24:12787779:T:C → chr=24, bp=12787779, a1=T, a2=C로 파싱
###############################################################
for SNP in "${INDEX_SNPS[@]}"; do
  IFS=':' read -r CHR BP A1 A2 <<< "$SNP"

  OUT_PREFIX="${OUTDIR}/09_LD_${CHR}_${BP}_${A1}_${A2}"

  echo "---------------------------------------------"
  echo "▶ INDEX SNP: chr=${CHR}, bp=${BP}, alleles=${A1}/${A2}"

  #########################################################################
  # BIM 파일에서 정확히 (CHR, BP, allele1, allele2) 모두 일치하는 라인을 찾기
  # BIM 포맷: CHR SNPID CM BP ALLELE1 ALLELE2
  #
  # allele 순서는 PLINK 내부 기준이라 T C 또는 C T 모두 허용하게 처리
  #########################################################################
  IDX_ID=$(awk -v c="$CHR" -v p="$BP" -v a1="$A1" -v a2="$A2" '
      ($1==c && $4==p) &&
      (($5==a1 && $6==a2) || ($5==a2 && $6==a1)) {
        print $2;
        exit;
      }' "$BIM")

  if [[ -z "$IDX_ID" ]]; then
    echo "❌ 정확히 일치하는 allele 조합의 SNP을 BIM 파일에서 찾지 못했습니다!"
    echo "   BIM에서 직접 확인하세요: CHR=$CHR, BP=$BP, Alleles=$A1/$A2"
    exit 1
  fi

  echo "✅ 선택된 SNP ID: ${IDX_ID}"

  printf "%s\n" "$IDX_ID" > "$TMP_SNP_LIST"
  echo "▶ OUT: ${OUT_PREFIX}.ld.gz"


  #########################################################################
  # PLINK LD 계산
  #########################################################################
  "$PLINK_PATH" \
    --bfile "$BFILE" \
    --dog \
    --r2 dprime with-freqs \
    --ld-snp-list "$TMP_SNP_LIST" \
    --ld-window-kb "$LD_WINDOW_KB" \
    --ld-window 999999 \
    --ld-window-r2 "$LD_WINDOW_R2_MIN" \
    --threads "$THREADS" \
    --memory "$MEMORY_MB" \
    --out "$OUT_PREFIX"

  gzip -f "${OUT_PREFIX}.ld"
  echo "✅ 완료: ${OUT_PREFIX}.ld.gz"

done

echo "🎉 모든 SNP 처리 완료"