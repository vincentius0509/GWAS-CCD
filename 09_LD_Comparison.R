# ==============================================================================
# 09_LD_Breed_Comparison.R
#
# 목적: 09_Linkage_Disequilibrium_local.sh(pooled)와
#       09_Linkage_Disequilibrium_local_BreedStratified.sh(breed별)로 각각
#       생성된 .ld.gz 파일들을 전부 읽어서 하나의 테이블로 합치고,
#       index SNP별로 pooled vs breed 간 LD 구조(r², D')를 비교/요약합니다.
#
#       -> "pooled LD가 breed 구조 때문에 왜곡된 건 아닌지" 판단하는 게 목적
# ==============================================================================

library(data.table)

# ------------------------------------------------------------------------------
# [수정 필요] 경로 확인
# ------------------------------------------------------------------------------
pooled_dir <- "J:/01_lab_lee/03_GWAS_Pure/14_Regression_local"
breed_dir  <- "J:/01_lab_lee/03_GWAS_Pure/14_Regression_local/09_BreedStratified_LD/09_LD_breed_stratified"
run_log    <- file.path(breed_dir, "breed_stratified_run.log")

out_dir <- "J:/01_lab_lee/03_GWAS_Pure/14_Regression_local"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# circadian entrainment 관련 5개 index SNP
# ------------------------------------------------------------------------------
index_snps <- data.table(
  CHR = c(5, 24, 24, 27, 26),
  BP  = c(25529184, 12787767, 12787779, 44543118, 13787793),
  A1  = c("C", "A", "T", "A", "A"),
  A2  = c("T", "G", "C", "T", "G"),
  Gene = c("GUCY1A2", "PLCB4", "PLCB4_2", "CACNA1C", "NOS1")
)

# ------------------------------------------------------------------------------
# breed별 표본수(N) 파싱 (bash 스크립트가 남긴 run log에서)
# 로그 형식: "breed,CHR:BP:A1:A2,N" 형태의 줄만 추출
# ------------------------------------------------------------------------------
breed_n <- data.table()
if (file.exists(run_log)) {
  log_lines <- readLines(run_log)
  data_lines <- log_lines[grepl("^[^=].*,.*,[0-9]+$", log_lines)]
  if (length(data_lines) > 0) {
    parsed <- tstrsplit(data_lines, ",")
    breed_n <- unique(data.table(Breed = parsed[[1]], N = as.integer(parsed[[3]])))
  }
}
print(breed_n)

# ------------------------------------------------------------------------------
# .ld.gz 파일 하나를 읽어 표준 컬럼으로 정리하는 함수
# PLINK1.9 --r2 dprime with-freqs 출력 컬럼:
#   CHR_A BP_A SNP_A MAF_A CHR_B BP_B SNP_B MAF_B R2 DP
# ------------------------------------------------------------------------------
read_ld_file <- function(path, group_label, n_indiv = NA_integer_) {
  if (!file.exists(path)) {
    cat("⚠️ 파일 없음:", path, "\n")
    return(NULL)
  }
  dt <- fread(path)
  if (nrow(dt) == 0) return(NULL)
  dt[, Group := group_label]
  dt[, N := n_indiv]
  dt
}

# ------------------------------------------------------------------------------
# 메인 루프: index SNP마다 pooled + 전체 breed 결과를 읽어서 합침
# ------------------------------------------------------------------------------
all_ld <- data.table()

breed_folders <- list.dirs(breed_dir, full.names = FALSE, recursive = FALSE)
cat("발견된 breed 폴더:", length(breed_folders), "개\n")
print(breed_folders)

for (i in seq_len(nrow(index_snps))) {
  chr <- index_snps$CHR[i]; bp <- index_snps$BP[i]
  a1  <- index_snps$A1[i];  a2 <- index_snps$A2[i]
  gene <- index_snps$Gene[i]
  
  # --- pooled ---
  pooled_path <- file.path(pooled_dir, sprintf("09_LD_%s_%s_%s_%s.ld.gz", chr, bp, a1, a2))
  pooled_dt <- read_ld_file(pooled_path, "Pooled", NA_integer_)
  if (!is.null(pooled_dt)) {
    pooled_dt[, Gene := gene]
    all_ld <- rbind(all_ld, pooled_dt, fill = TRUE)
  }
  
  # --- breed별 ---
  for (breed in breed_folders) {
    breed_path <- file.path(breed_dir, breed, sprintf("09_LD_%s_%s_%s_%s_%s.ld.gz", chr, bp, a1, a2, breed))
    n_indiv <- if (nrow(breed_n[Breed == breed]) > 0) breed_n[Breed == breed]$N[1] else NA_integer_
    breed_dt <- read_ld_file(breed_path, breed, n_indiv)
    if (!is.null(breed_dt)) {
      breed_dt[, Gene := gene]
      all_ld <- rbind(all_ld, breed_dt, fill = TRUE)
    }
  }
}

cat("\n총", nrow(all_ld), "개 LD pair 로드 완료\n")

# ------------------------------------------------------------------------------
# 전체 결합 테이블 저장 (long format: 모든 pair, 모든 group)
# ------------------------------------------------------------------------------
fwrite(all_ld, file.path(out_dir, "09_LD_Combined_Pooled_vs_Breed.csv"))
cat("✅ 저장 완료:", file.path(out_dir, "09_LD_Combined_Pooled_vs_Breed.csv"), "\n")

# ------------------------------------------------------------------------------
# 요약 테이블: Group x Index SNP(Gene) 별로
#   - window 내 SNP 개수
#   - r2 >= 0.8(강한 LD)인 SNP 개수
#   - mean r2, max r2 (자기 자신과의 쌍(R2=1)은 제외)
# ------------------------------------------------------------------------------
summary_dt <- all_ld[SNP_A != SNP_B, .(
  N_SNPs_in_window = .N,
  N_SNPs_R2_ge_0.8 = sum(R2 >= 0.8, na.rm = TRUE),
  Mean_R2 = round(mean(R2, na.rm = TRUE), 4),
  Max_R2  = round(max(R2, na.rm = TRUE), 4),
  N_indiv = N[1]
), by = .(Group, Gene)]

setorder(summary_dt, Gene, -N_indiv)
print(summary_dt)

fwrite(summary_dt, file.path(out_dir, "09_LD_Summary_Pooled_vs_Breed.csv"))
cat("✅ 요약 테이블 저장:", file.path(out_dir, "09_LD_Summary_Pooled_vs_Breed.csv"), "\n")

# ------------------------------------------------------------------------------
# index SNP별 비교 그래프: Group(Pooled + breed들)별 R2 분포 boxplot
# ------------------------------------------------------------------------------
for (i in seq_len(nrow(index_snps))) {
  gene <- index_snps$Gene[i]
  sub <- all_ld[Gene == gene & SNP_A != SNP_B]
  if (nrow(sub) == 0) next
  
  # Pooled을 맨 왼쪽에 오도록 순서 지정
  group_order <- c("Pooled", setdiff(unique(sub$Group), "Pooled"))
  sub$Group <- factor(sub$Group, levels = group_order)
  
  png(file.path(out_dir, sprintf("09_LD_Boxplot_%s.png", gene)),
      width = max(2400, 300 * length(group_order)), height = 1200, res = 200)
  par(mar = c(8, 4.5, 3, 1))
  boxplot(R2 ~ Group, data = sub, las = 2, col = "#2A4F9E20", border = "#2A4F9E",
          main = sprintf("Pooled vs Breed-stratified LD (r^2) around %s", gene),
          ylab = expression(r^2), xlab = "")
  dev.off()
  cat("  -> 저장:", sprintf("09_LD_Boxplot_%s.png", gene), "\n")
}

cat("\n모든 index SNP에 대한 pooled vs breed LD 비교 완료.\n")
