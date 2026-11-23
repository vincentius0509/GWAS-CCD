#!/bin/bash
set -euo pipefail

# ============================================================
# [수정 필요] 아래 변수들을 본인 환경에 맞게 확인/수정하세요
# ============================================================

PLINK2_PATH="$HOME/tools/plink2/plink2"

# raw curated release (필터링 이전 원본 genotype)
INPUT_BFILE="/mnt/j/01_lab_lee/05_Data/genetic/v2022/DogAgingProject_GeneticData_CuratedRelease_2022"

# phenotype / covariate 파일
PHENO_FILE="/mnt/j/01_lab_lee/03_GWAS_Pure/04_pheno001.txt"
COVAR_FILE="/mnt/j/01_lab_lee/03_GWAS_Pure/14_Regression_local/04_cov_001_1e6_10.txt"

# 결과 저장 경로
OUTDIR="/mnt/j/01_lab_lee/03_GWAS_Pure/14_Regression_local/12_Geno_Beta"
mkdir -p "$OUTDIR"

# 로컬 CPU(2코어/4스레드) 기준
THREADS=3
MEMORY_MB=12000

# 각 index SNP 기준 beta 계산 window (±bp)
WINDOW_BP=50000

echo "▶ INPUT_BFILE: $INPUT_BFILE"
echo "▶ OUTDIR:      $OUTDIR"

# ------------------------------------------------------------------------------
# 5개 index SNP (CHR:BP:A1:A2, Gene 라벨) — beta 계산 시 region 중심으로 사용
# ------------------------------------------------------------------------------
INDEX_LOCI=(
  "5:25529184:C:T:GUCY1A2"
  "24:12787767:A:G:PLCB4"
  "24:12787779:T:C:PLCB4_2"
  "27:44543118:A:T:CACNA1C"
  "26:13787793:A:G:NOS1"
)

# LD_Regional_Plot.R이 만든 locus별 hits CSV(r2>=0.2)가 있는 폴더
LD_CSV_DIR="/mnt/j/01_lab_lee/00_Picture/04_GWAS_pure/02_LD_local"
echo "▶ LD_CSV_DIR:  $LD_CSV_DIR"

# ------------------------------------------------------------------------------
# LD_CSV_DIR의 "*_LD_hits_r2_02.csv" 파일들을 전부 읽어서 SNP 컬럼(1번째 컬럼)만
# 추출 -> 하나의 전체 SNP 목록으로 통합 (PRS에 쓸 전체 candidate SNP)
# ------------------------------------------------------------------------------
SNP_LIST_FILE="${OUTDIR}/candidate_snp_list.txt"

shopt -s nullglob
HIT_CSVS=("${LD_CSV_DIR}"/*_LD_hits_r2_02.csv)
shopt -u nullglob

if [[ ${#HIT_CSVS[@]} -eq 0 ]]; then
  echo "❌ ${LD_CSV_DIR}에서 *_LD_hits_r2_02.csv 파일을 찾지 못했습니다."
  exit 1
fi

echo "▶ 발견된 LD hits CSV: ${#HIT_CSVS[@]}개"
for f in "${HIT_CSVS[@]}"; do
  echo "   - $(basename "$f")"
done

: > "$SNP_LIST_FILE"
for f in "${HIT_CSVS[@]}"; do
  awk -F',' 'NR>1 {print $1}' "$f" >> "$SNP_LIST_FILE"
done

# 5개 index SNP 자체도 hits 목록에서 빠졌을 경우를 대비해 명시적으로 추가
for LOCUS in "${INDEX_LOCI[@]}"; do
  IFS=':' read -r CHR BP A1 A2 GENE <<< "$LOCUS"
  echo "${CHR}:${BP}:${A1}:${A2}" >> "$SNP_LIST_FILE"
done

sort -u "$SNP_LIST_FILE" -o "$SNP_LIST_FILE"
N_SNPS=$(wc -l < "$SNP_LIST_FILE")
echo "▶ 통합된 candidate SNP 개수: ${N_SNPS}개 -> ${SNP_LIST_FILE}"

# ------------------------------------------------------------------------------
# 2) SNP별로 개별 유전형 추출 (recode A) - PRS 파이프라인에서 SNP 단위 파일이 필요
# ------------------------------------------------------------------------------
echo "-----------------------------------------------"
echo "▶ [1/2] SNP별 유전형 추출 (recode A)"

while IFS= read -r SNP; do
  [[ -z "$SNP" ]] && continue

  # 파일명에 쓰기 위해 콜론을 언더스코어로 변환 (예: 24:12787779:T:C -> 24_12787779_T_C)
  SAFE_NAME=$(echo "$SNP" | tr ':' '_')

  echo "▶ 추출 중: ${SNP}"

  "$PLINK2_PATH" \
    --bfile "$INPUT_BFILE" \
    --snp "$SNP" \
    --recode A \
    --dog \
    --threads "$THREADS" \
    --memory "$MEMORY_MB" \
    --out "${OUTDIR}/snp_${SAFE_NAME}_genotype" \
    || { echo "⚠️ ${SNP} 추출 실패 - 건너뜀"; continue; }

done < "$SNP_LIST_FILE"

echo "✅ SNP별 유전형 추출 완료 (${OUTDIR}/snp_*_genotype.raw)"

# ------------------------------------------------------------------------------
# 3) index SNP별로 ±WINDOW_BP region에 대해 beta 계산 (--glm)
# ------------------------------------------------------------------------------
echo "-----------------------------------------------"
echo "▶ [2/2] region별 beta 계산 (±${WINDOW_BP}bp)"

for LOCUS in "${INDEX_LOCI[@]}"; do
  IFS=':' read -r CHR BP A1 A2 GENE <<< "$LOCUS"

  FROM_BP=$((BP - WINDOW_BP))
  TO_BP=$((BP + WINDOW_BP))

  OUT_PREFIX="${OUTDIR}/Beta_chr${CHR}_${BP}_${GENE}"

  echo "▶ ${GENE}: chr${CHR}:${FROM_BP}-${TO_BP} (center=${BP})"

  "$PLINK2_PATH" \
    --bfile "$INPUT_BFILE" \
    --chr "$CHR" \
    --from-bp "$FROM_BP" \
    --to-bp "$TO_BP" \
    --glm cols=+a1freq hide-covar \
    --pheno "$PHENO_FILE" \
    --pheno-name phe1 \
    --covar "$COVAR_FILE" \
    --covar-variance-standardize \
    --dog \
    --threads "$THREADS" \
    --memory "$MEMORY_MB" \
    --out "$OUT_PREFIX"

  echo "✅ 완료: ${OUT_PREFIX}.phe1.glm.linear"
done

echo "🎉 모든 candidate SNP 유전형 추출 + region beta 계산 완료 -> ${OUTDIR}"
