# ==============================================================================
# Manhattan / QQ plot 재생성 - 제목의 PC 개수를 실제 covariate 파일 기준으로 표기
#
# 배경: 파일명(예: "05_GWAS_PHE1_001_1e6_10_linear")의 마지막 숫자(10)는
#       "PCA를 10개까지 계산한 소스"를 가리키는 명명 규칙일 뿐이며,
#       실제 GWAS 회귀분석에 covariate로 들어간 PC 개수는 collinearity 등의
#       이유로 그보다 적을 수 있습니다(예: 실제로는 PC1-8만 사용).
#       이 코드는 04_cov_*.txt 파일에서 실제 PC 컬럼 개수를 세어, 그 값을
#       그림 제목에 정확히 반영합니다. GWAS 결과(p-value, λGC)는 재계산하지
#       않고 기존 값을 그대로 사용합니다.
# ==============================================================================

library(qqman)
library(data.table)
getwd()

base_dir <- getwd()
out_dir  <- "J:/01_lab_lee/00_Picture/04_GWAS_pure/01_LR"

# 조합 정의 및 MAF 텍스트 라벨 <-> 실제 숫자값 매핑 (매칭 버그 방지의 핵심)
# ------------------------------------------------------------------------------
maf_vals   <- c("001", "005", "01")
maf_numeric_map <- c("001" = 0.01, "005" = 0.05, "01" = 0.1)   # 실제 MAF 값
maf_title_map   <- c("001" = "0.01", "005" = "0.05", "01" = "0.1")  # 제목 표시용
hwe_val    <- "1e6"
hwe_numeric <- 1e6
pc_vals    <- c(5, 10, 15, 20, 25, 30)

genomewide_threshold <- 5e-8
suggestive_threshold <- 1e-5
low_threshold <- 1e-4

# ------------------------------------------------------------------------------
# 기존에 계산해두신 λGC 결과표(results)를 그대로 재사용
# (results <- fread("05_Genomic_Inflation_Factor.csv") 로 다시 불러오셔도 됩니다)
# ------------------------------------------------------------------------------
results <- fread(file.path("05_Genomic_Inflation_Factor.csv"))
head(results)

calculate_lambda_gc <- function(p_values) {
  p_values <- p_values[is.finite(p_values) & p_values > 0 & p_values < 1]
  if (length(p_values) == 0) return(NA_real_)
  chisq <- qchisq(1 - p_values, df = 1)
  median(chisq) / qchisq(0.5, df = 1)
}

# ------------------------------------------------------------------------------
# 실제 covariate 파일에서 PC 컬럼 개수를 세는 함수
# ------------------------------------------------------------------------------
get_actual_pc_count <- function(maf, pc_label) {
  cov_path <- file.path(base_dir, sprintf("04_cov_%s_1e6_%d.txt", maf, pc_label))
  if (!file.exists(cov_path)) {
    cat("⚠️ covariate 파일을 찾을 수 없습니다:", cov_path, "\n")
    return(NA_integer_)
  }
  cov_data <- fread(cov_path, nrows = 1)
  pc_cols <- grep("^PC[0-9]+$", colnames(cov_data), value = TRUE)
  length(pc_cols)
}

# ------------------------------------------------------------------------------
# GWAS 결과 컬럼 정리 (qqman 요구 형식: CHR, BP, SNP, P)
# ------------------------------------------------------------------------------
prep_gwas_df <- function(df) {
  df <- copy(df)
  setnames(df, old = "#CHROM", new = "CHR", skip_absent = TRUE)
  setnames(df, old = "POS",    new = "BP",  skip_absent = TRUE)
  setnames(df, old = "ID",     new = "SNP", skip_absent = TRUE)
  df$CHR <- as.numeric(gsub("X", "39", df$CHR))
  df
}

# ------------------------------------------------------------------------------
# [핵심] Manhattan + QQ를 한 이미지(가로 1x2)로 결합 출력하는 함수
# ------------------------------------------------------------------------------
plot_combined_manhattan_qq_horizontal <- function(df,
                                                  title_suffix,
                                                  out_prefix,
                                                  lambda_gc_value,
                                                  thr_gw = 5e-8, thr_sg = 1e-5, thr_ex = 1e-4,
                                                  ymax = NULL) {
  df <- prep_gwas_df(df)
  if (is.null(ymax)) {
    ymax <- ceiling(max(-log10(df$P), na.rm = TRUE) * 1.05)
    ymax <- max(ymax, 10)
  }
  
  ## ------------------ PDF (가로 1x2 패널) ------------------
  pdf(paste0(out_prefix, "_horizontal.pdf"), width = 12, height = 5, useDingbats = FALSE)
  par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3.0, 1.0), mgp = c(2.2, 0.8, 0))
  
  manhattan(df, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
            col = c("#2A4F9E", "#FF9D00"),
            suggestiveline = -log10(thr_sg), genomewideline = -log10(thr_gw),
            ylim = c(0, ymax), cex = 0.42, cex.axis = 0.85,
            main = paste("Manhattan plot of the GWAS", "\n", title_suffix),
            ylab = expression(-log[10](italic(p))), xlab = "Chromosome")
  abline(h = -log10(thr_ex), lty = 2, col = "#666666")
  
  qq(df$P,
     main = paste0("QQ plot of the GWAS ", title_suffix, "\n",
                   "lambdaGC = ", sprintf("%.3f", lambda_gc_value)),
     xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
     col = "#2A4F9E", cex = 0.6)
  abline(0, 1, col = "#666666", lwd = 1.2)
  dev.off()
  
  ## ------------------ PNG (고해상도, 한 파일에 두 패널) ------------------
  png(paste0(out_prefix, "_horizontal.png"), width = 3600, height = 1500, res = 300)
  par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3.0, 1.0), mgp = c(2.2, 0.8, 0))
  
  manhattan(df, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
            col = c("#2A4F9E", "#FF9D00"),
            suggestiveline = -log10(thr_sg), genomewideline = -log10(thr_gw),
            ylim = c(0, ymax), cex = 0.42, cex.axis = 0.85,
            main = paste("Manhattan plot of the GWAS", "\n", title_suffix),
            ylab = expression(-log[10](italic(p))), xlab = "Chromosome")
  abline(h = -log10(thr_ex), lty = 2, col = "#666666")
  
  qq(df$P,
     main = paste0("QQ plot of the GWAS ", title_suffix, "\n",
                   "lambdaGC = ", sprintf("%.3f", lambda_gc_value)),
     xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
     col = "#2A4F9E", cex = 0.6)
  abline(0, 1, col = "#666666", lwd = 1.2)
  dev.off()
  
  cat("  -> 저장 완료:", basename(paste0(out_prefix, "_horizontal.png/.pdf")), "\n")
}

# ------------------------------------------------------------------------------
# 메인 루프: 모든 MAF x PC 조합에 대해 결합 이미지 재생성
# ------------------------------------------------------------------------------
for (maf in maf_vals) {
  for (pc_label in pc_vals) {
    
    # (1) 실제 covariate PC 개수 (제목 표시용)
    actual_pc <- get_actual_pc_count(maf, pc_label)
    if (is.na(actual_pc)) next
    pc_title <- if (actual_pc == 1) "PC1" else sprintf("PC1-%d", actual_pc)
    
    cat(sprintf("MAF=%s, 파일명 PC=%d -> 실제 covariate PC 개수=%d (%s)\n",
                maf, pc_label, actual_pc, pc_title))
    
    # (2) GWAS 결과 로딩
    linear_path <- file.path(base_dir, sprintf("05_GWAS_PHE1_%s_%s_%d_linear.phe1.glm.linear", maf, hwe_val, pc_label))
    if (!file.exists(linear_path)) {
      cat("⚠️ GWAS 결과 파일을 찾을 수 없습니다:", linear_path, "\n")
      next
    }
    plink_dt <- fread(linear_path)
    
    # (3) lambdaGC 즉석 재계산 (results 테이블 매칭 버그를 우회, 항상 정확함)
    lambda_val <- calculate_lambda_gc(plink_dt$P)
    
    # (4) 결합 이미지 생성 (제목: 실제 MAF 숫자값 + 실제 PC 개수)
    title_suffix <- sprintf("(MAF=%s, HWE=1e-6, %s)", maf_title_map[maf], pc_title)
    out_prefix <- file.path(out_dir, sprintf("01_Combined_%s_1e6_%d", maf, pc_label))
    
    plot_combined_manhattan_qq_horizontal(
      df = plink_dt,
      title_suffix = title_suffix,
      out_prefix = out_prefix,
      lambda_gc_value = lambda_val
    )
    
    rm(plink_dt); gc()
  }
}

cat("\n모든 조합에 대한 결합 이미지 생성 완료.\n")
