# ==============================================================================
# 10_Breed_Keep_Files.R
#
# 목적: 03-2_Add_Phenotype.rdata의 원본(필터링 이전) breed one-hot 컬럼에서
#       개체별 단일 breed 라벨을 재구성하고, 개체수가 충분한 상위 breed에 대해
#       PLINK --keep 형식(FID, IID, 헤더 없음, 탭 구분)의 파일을 생성합니다.
#
#       이 --keep 파일들은 09_Linkage_Disequilibrium_local.sh에서
#       breed-stratified LD를 돌릴 때 사용됩니다.
# ==============================================================================

library(data.table)
getwd()
# ------------------------------------------------------------------------------
# [수정 필요] 경로 확인
# ------------------------------------------------------------------------------
load("J:/01_lab_lee/00_Code/00_RDataFile/03_GWAS_Pure/03-2_Add_Phenotype.rdata")

out_dir <- "J:/01_lab_lee/03_GWAS_Pure/14_Regression_local/breed_keep"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# [확인 필요] 사용할 원본 covariate 객체 이름
# MAF=0.01, PC=10(1e6_10) 조합의 "필터링 이전" 버전을 씁니다.
# (Golden_Retriever, German_Shepherd_Dog 등이 아직 빠지지 않은 상태의 객체)
# 정확한 객체명이 다르면 아래 한 줄만 바꿔주세요.
# ------------------------------------------------------------------------------
cat("사용 가능한 객체 목록:\n")
print(ls())

cov_raw <- data_001_1e6_10_cov_filtered2  # [필요시 이 객체명 수정]

# ------------------------------------------------------------------------------
# breed one-hot 컬럼 식별
#   FID, IID, PC1~PCn, Sex_cov, Age_cov, Intact_cov, pa_score_cov을 제외한
#   나머지 0/1 컬럼을 전부 breed dummy로 간주
# ------------------------------------------------------------------------------
non_breed_cols <- c("FID", "IID",
                    grep("^PC[0-9]+$", names(cov_raw), value = TRUE),
                    "Sex_cov", "Age_cov", "Intact_cov", "pa_score_cov")

breed_cols <- setdiff(names(cov_raw), non_breed_cols)
cat("breed 컬럼 개수:", length(breed_cols), "\n")

# ------------------------------------------------------------------------------
# 개체별 단일 breed 라벨 재구성 (one-hot -> 라벨)
#   purebred 전제이므로 각 행에 breed 컬럼 중 정확히 1개만 1이어야 함
# ------------------------------------------------------------------------------
breed_matrix <- as.matrix(cov_raw[, ..breed_cols])
n_breeds_per_row <- rowSums(breed_matrix)

if (any(n_breeds_per_row != 1)) {
  cat("⚠️ breed 컬럼 합이 1이 아닌 개체가", sum(n_breeds_per_row != 1), "마리 있습니다.\n")
  cat("   (0개: breed 컬럼에 없는 개체 / 2개 이상: 중복 표시 개체 - 확인 필요)\n")
}

cov_raw[, breed := breed_cols[apply(breed_matrix, 1, function(r) {
  idx <- which(r == 1)
  if (length(idx) == 1) idx else NA
})]]

# ------------------------------------------------------------------------------
# breed별 개체수 집계 및 상위 breed 선정
# ------------------------------------------------------------------------------
breed_counts <- cov_raw[!is.na(breed), .N, by = breed][order(-N)]
print(breed_counts)

fwrite(breed_counts, file.path(out_dir, "breed_counts_full.csv"))

# [수정 필요] 최소 개체수 기준 (LD 추정이 안정적이려면 어느 정도 표본이 필요)
MIN_N <- 10
top_breeds <- breed_counts[N >= MIN_N]
cat("\n선정된 breed (N >=", MIN_N, "):", nrow(top_breeds), "개\n")
print(top_breeds)

# ------------------------------------------------------------------------------
# breed별 --keep 파일 생성 (FID, IID, 헤더 없음)
# ------------------------------------------------------------------------------
for (b in top_breeds$breed) {
  keep_dt <- cov_raw[breed == b, .(FID, IID)]
  keep_path <- file.path(out_dir, paste0("keep_", b, ".txt"))
  fwrite(keep_dt, keep_path, sep = "\t", col.names = FALSE)
  cat(sprintf("✅ %s (N=%d) -> %s\n", b, nrow(keep_dt), keep_path))
}

cat("\n모든 breed --keep 파일 생성 완료:", out_dir, "\n")
