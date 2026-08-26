#!/bin/bash
set -euo pipefail

PLINK_PATH="/data/00_tools/plink_linux_x86_64_20231018/plink"

# 단일 bfile
BFILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6"
BIM="${BFILE}.bim"

#############################################
# INDEX SNP 목록 — 반드시 CHR:BP:A1:A2 형식
# 예: 24:12787779:T:C
#############################################
INDEX_SNPS=("24:12787779:T:C")

# LD window 설정
LD_WINDOW_KB=1000
LD_WINDOW_R2_MIN=0.0
THREADS=20

OUTDIR="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure"

echo "▶ BFILE: $BFILE"
echo "▶ WINDOW: ${LD_WINDOW_KB}kb, r² ≥ ${LD_WINDOW_R2_MIN}"
echo "▶ OUTDIR: $OUTDIR"

TMP_SNP_LIST="$(mktemp)"
trap 'rm -f "$TMP_SNP_LIST"' EXIT


###############################################################
# 24:12787779:T:C → chr=24, bp=12787779, a1=T, a2=C로 파싱
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
    --out "$OUT_PREFIX"

  gzip -f "${OUT_PREFIX}.ld"
  echo "✅ 완료: ${OUT_PREFIX}.ld.gz"

done

echo "🎉 모든 SNP 처리 완료"
