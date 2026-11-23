# ==============================================================================
# LD_SNP.R
#
# 목적: circadian entrainment 관련 5개 SNP(GUCY1A2, PLCB4 x2, CACNA1C, NOS1)의
#       GWAS 결과(REF, ALT, Effect Allele(A1), β, p-value)를 뽑아
#       Locus (Chr:Pos) | REF | ALT | Effect Allele (A1) | Gene Symbol | β | p-value
#       형식의 CSV로 저장. Gene Symbol은 직접 채워 넣을 것이므로 빈 컬럼으로 둠.
# ==============================================================================

library(data.table)

# ------------------------------------------------------------------------------
# GWAS 결과(.glm.linear) 파일 경로
# ------------------------------------------------------------------------------
gwas_path <- "J:/01_lab_lee/03_GWAS_Pure/14_Regression_local/05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear"

out_dir <- "J:/01_lab_lee/03_GWAS_Pure/14_Regression_local"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# circadian entrainment 관련 5개 SNP (ID 형식: CHR:POS:REF:ALT)
# ------------------------------------------------------------------------------
target_snps <- c(
  "5:25529184:C:T",    # GUCY1A2
  "24:12787767:A:G",   # PLCB4
  "24:12787779:T:C",   # PLCB4 (2nd variant)
  "27:44543118:A:T",   # CACNA1C
  "26:13787793:A:G"    # NOS1
)

# ------------------------------------------------------------------------------
# GWAS 결과 로드 및 필터링
# ------------------------------------------------------------------------------
gwas_dt <- fread(gwas_path)

matched <- gwas_dt[ID %in% target_snps]

if (nrow(matched) < length(target_snps)) {
  missing_snps <- setdiff(target_snps, matched$ID)
  cat("⚠️ GWAS 결과에서 찾지 못한 SNP:\n")
  print(missing_snps)
  cat("   -> 파일 경로/조합(MAF·PC)이 맞는지 확인해주세요.\n")
}

head(matched)
# ------------------------------------------------------------------------------
# 요청하신 표 형식으로 정리
# ------------------------------------------------------------------------------
ld_snp_table <- matched[, .(
  `Locus (Chr:Pos)` = paste0(`#CHROM`, ":", POS),
  REF,
  ALT,
  `Effect Allele (A1)` = A1,
  `Gene Symbol` = "",          # 직접 채워 넣을 컬럼
  `β (Effect Size)` = round(BETA, 4),
  `p-value` = P
)]

# target_snps 순서(GUCY1A2 -> PLCB4 -> PLCB4 -> CACNA1C -> NOS1)대로 정렬
match_order <- match(paste0(matched$`#CHROM`, ":", matched$POS),
                     sapply(strsplit(target_snps, ":"), function(x) paste0(x[1], ":", x[2])))
ld_snp_table <- ld_snp_table[order(match_order)]

print(ld_snp_table)

# ------------------------------------------------------------------------------
# 저장
# ------------------------------------------------------------------------------
out_path <- file.path(out_dir, "08_LD_SNP_CircadianEntrainment.csv")
fwrite(ld_snp_table, out_path)
cat("✅ 저장 완료:", out_path, "\n")
