# ==============================================================================
# 05_Regression_Result3.r
#
# 목적: MAF(0.01/0.05/0.1) x PC(5~30) x Model(linear / linear_noActivity)
#       전체 36개 조합에 대해:
#         1) genomic inflation factor(λGC) 계산 (05_Regression_Result.r 방식)
#         2) horizontal Manhattan+QQ 결합 플롯 생성 (05_Regression_Result2.r 방식)
#         3) λGC 요약 테이블 저장 + linear vs noActivity 비교 그래프 생성
#
#       결과 이미지/표는 전부 01_LR_local 폴더에 저장됩니다.
#       (05_Regression_Result2.r의 get_actual_pc_count()처럼, 그림 제목의
#        PC 개수는 파일명 숫자가 아니라 실제 covariate 파일의 PC 컬럼 수를
#        기준으로 표기합니다.)
# ==============================================================================

library(qqman)
library(data.table)

# ------------------------------------------------------------------------------
# [수정 필요] 경로 확인
# ------------------------------------------------------------------------------
# regression 결과(.glm.linear)와 covariate 파일(04_cov_*.txt)이 있는 폴더
base_dir <- "J:/01_lab_lee/03_GWAS_Pure/14_Regression_local"

# 결합 플롯 및 요약 결과가 저장될 폴더
out_dir  <- "J:/01_lab_lee/00_Picture/04_GWAS_pure/01_LR_local"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 조합 정의
# ------------------------------------------------------------------------------
maf_vals        <- c("001", "005", "01")
maf_title_map   <- c("001" = "0.01", "005" = "0.05", "01" = "0.1")
hwe_val         <- "1e6"
pc_vals         <- c(5, 10, 15, 20, 25, 30)

# 모델 변형: linear(activity 포함, 기존) / linear_noActivity(activity 제외, sensitivity)
model_vals <- c("linear", "linear_noActivity")
model_label_map <- c("linear" = "With Activity", "linear_noActivity" = "Without Activity")

genomewide_threshold <- 5e-8
suggestive_threshold <- 1e-5
low_threshold <- 1e-4

# ------------------------------------------------------------------------------
# λGC 계산 함수 (05_Regression_Result.r 방식)
# ------------------------------------------------------------------------------
calculate_lambda_gc <- function(p_values) {
  p_values <- p_values[is.finite(p_values) & p_values > 0 & p_values < 1]
  if (length(p_values) == 0) return(NA_real_)
  chisq <- qchisq(1 - p_values, df = 1)
  median(chisq) / qchisq(0.5, df = 1)
}

# ------------------------------------------------------------------------------
# 실제 covariate 파일에서 PC 컬럼 개수를 세는 함수 (05_Regression_Result2.r 방식)
# noActivity 모델도 covariate 파일 자체는 activity 포함 원본 파일을 그대로 쓰고,
# --covar-name에서 pa_score_cov만 빼고 실행한 것이므로 PC 개수 파악 로직은 동일함
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
# horizontal Manhattan+QQ 결합 플롯 함수 (05_Regression_Result2.r 방식)
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
  
  manhattan_title <- paste0("Manhattan plot of the GWAS\n", title_suffix)
  qq_title <- paste0("QQ plot of the GWAS ", title_suffix,
                     "\nlambdaGC = ", sprintf("%.3f", lambda_gc_value))
  
  ## ------------------ PDF (가로 1x2 패널) ------------------
  pdf(paste0(out_prefix, "_horizontal.pdf"), width = 14, height = 5.5, useDingbats = FALSE)
  par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3.6, 1.0), mgp = c(2.2, 0.8, 0), cex.main = 0.85)
  
  manhattan(df, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
            col = c("#2A4F9E", "#FF9D00"),
            suggestiveline = -log10(thr_sg), genomewideline = -log10(thr_gw),
            ylim = c(0, ymax), cex = 0.42, cex.axis = 0.85,
            main = manhattan_title,
            ylab = expression(-log[10](italic(p))), xlab = "Chromosome")
  abline(h = -log10(thr_ex), lty = 2, col = "#666666")
  
  qq(df$P,
     main = qq_title,
     xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
     col = "#2A4F9E", cex = 0.6)
  abline(0, 1, col = "#666666", lwd = 1.2)
  dev.off()
  
  ## ------------------ PNG (고해상도, 한 파일에 두 패널) ------------------
  png(paste0(out_prefix, "_horizontal.png"), width = 4200, height = 1650, res = 300)
  par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3.6, 1.0), mgp = c(2.2, 0.8, 0), cex.main = 0.85)
  
  manhattan(df, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
            col = c("#2A4F9E", "#FF9D00"),
            suggestiveline = -log10(thr_sg), genomewideline = -log10(thr_gw),
            ylim = c(0, ymax), cex = 0.42, cex.axis = 0.85,
            main = manhattan_title,
            ylab = expression(-log[10](italic(p))), xlab = "Chromosome")
  abline(h = -log10(thr_ex), lty = 2, col = "#666666")
  
  qq(df$P,
     main = qq_title,
     xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
     col = "#2A4F9E", cex = 0.6)
  abline(0, 1, col = "#666666", lwd = 1.2)
  dev.off()
  
  cat("  -> 저장 완료:", basename(paste0(out_prefix, "_horizontal.png/.pdf")), "\n")
}

# ------------------------------------------------------------------------------
# 메인 루프: MAF x PC x Model(linear/noActivity) 전체 조합
#   - λGC 계산 -> 요약 테이블에 누적
#   - 실제 PC 개수 반영한 결합 Manhattan+QQ 플롯 생성 -> 01_LR_local 저장
# ------------------------------------------------------------------------------
lambda_results <- data.table()

for (maf in maf_vals) {
  for (pc_label in pc_vals) {
    
    # 실제 covariate PC 개수 (제목 표시용, 두 모델 공통)
    actual_pc <- get_actual_pc_count(maf, pc_label)
    if (is.na(actual_pc)) next
    pc_title <- if (actual_pc == 1) "PC1" else sprintf("PC1-%d", actual_pc)
    
    for (model in model_vals) {
      
      linear_path <- file.path(
        base_dir,
        sprintf("05_GWAS_PHE1_%s_%s_%d_%s.phe1.glm.linear", maf, hwe_val, pc_label, model)
      )
      
      if (!file.exists(linear_path)) {
        cat("⚠️ GWAS 결과 파일을 찾을 수 없습니다:", linear_path, "\n")
        next
      }
      
      cat(sprintf("▶ MAF=%s, PC=%d(%s), Model=%s\n", maf, pc_label, pc_title, model))
      
      plink_dt <- fread(linear_path)
      
      # λGC 계산
      lambda_val <- calculate_lambda_gc(plink_dt$P)
      
      lambda_results <- rbind(
        lambda_results,
        data.table(
          MAF = maf,
          HWE = hwe_val,
          PC_label = pc_label,
          Actual_PC = actual_pc,
          Model = model,
          Lambda_GC = lambda_val
        )
      )
      
      # 결합 Manhattan+QQ 플롯
      title_suffix <- sprintf("(MAF=%s, HWE=1e-6, %s, %s)",
                              maf_title_map[maf], pc_title, model_label_map[model])
      out_prefix <- file.path(out_dir, sprintf("01_Combined_%s_1e6_%d_%s", maf, pc_label, model))
      
      plot_combined_manhattan_qq_horizontal(
        df = plink_dt,
        title_suffix = title_suffix,
        out_prefix = out_prefix,
        lambda_gc_value = lambda_val
      )
      
      rm(plink_dt); gc()
    }
  }
}

# ------------------------------------------------------------------------------
# λGC 요약 테이블 저장
# ------------------------------------------------------------------------------
fwrite(lambda_results, file.path(out_dir, "05_Genomic_Inflation_Factor_linear_vs_noActivity.csv"))
cat("\n✅ λGC 요약 테이블 저장:", file.path(out_dir, "05_Genomic_Inflation_Factor_linear_vs_noActivity.csv"), "\n")
print(lambda_results)

# ------------------------------------------------------------------------------
# 요약 비교 그래프: MAF별로 PC 개수에 따른 λGC 추이 (Activity 포함 vs 제외)
# ------------------------------------------------------------------------------
png(file.path(out_dir, "05_Lambda_Comparison_linear_vs_noActivity.png"),
    width = 3600, height = 1200, res = 300)
par(mfrow = c(1, 3), mar = c(4.2, 4.6, 3.0, 1.0), mgp = c(2.2, 0.8, 0))

for (maf in maf_vals) {
  sub <- lambda_results[MAF == maf]
  sub_linear <- sub[Model == "linear"][order(PC_label)]
  sub_noAct  <- sub[Model == "linear_noActivity"][order(PC_label)]
  
  y_range <- range(c(sub_linear$Lambda_GC, sub_noAct$Lambda_GC), na.rm = TRUE)
  y_range <- c(min(y_range, 0.95, na.rm = TRUE), max(y_range, 1.05, na.rm = TRUE))
  
  plot(sub_linear$PC_label, sub_linear$Lambda_GC, type = "b", pch = 19, lwd = 2,
       col = "#2A4F9E", ylim = y_range,
       xlab = "PC (file label)", ylab = expression(lambda[GC]),
       main = sprintf("MAF = %s", maf_title_map[maf]))
  lines(sub_noAct$PC_label, sub_noAct$Lambda_GC, type = "b", pch = 17, lwd = 2, col = "#FF9D00")
  abline(h = 1, lty = 2, col = "#999999")
  
  legend("topright", legend = c("With Activity", "Without Activity"),
         col = c("#2A4F9E", "#FF9D00"), pch = c(19, 17), lwd = 2, bty = "n", cex = 0.8)
}

dev.off()
cat("✅ 비교 그래프 저장:", file.path(out_dir, "05_Lambda_Comparison_linear_vs_noActivity.png"), "\n")

cat("\n모든 조합에 대한 λGC 계산 및 결합 이미지 생성 완료.\n")

# ------------------------------------------------------------------------------
#  Activity 포함("linear") 모델의 SNP threshold 개수 집계
# ------------------------------------------------------------------------------
snp_counts <- data.table()

for (maf in maf_vals) {
  for (pc_label in pc_vals) {
    
    linear_path <- file.path(
      base_dir,
      sprintf("05_GWAS_PHE1_%s_%s_%d_linear.phe1.glm.linear", maf, hwe_val, pc_label)
    )
    if (!file.exists(linear_path)) {
      cat("⚠️ GWAS 결과 파일을 찾을 수 없습니다:", linear_path, "\n")
      next
    }
    
    plink_dt <- fread(linear_path, select = "P")
    
    snp_counts <- rbind(
      snp_counts,
      data.table(
        MAF = maf,
        PC_label = pc_label,
        SNPs_5e8 = sum(plink_dt$P < genomewide_threshold, na.rm = TRUE),
        SNPs_1e5 = sum(plink_dt$P < suggestive_threshold, na.rm = TRUE),
        SNPs_1e4 = sum(plink_dt$P < low_threshold, na.rm = TRUE)
      )
    )
    
    rm(plink_dt); gc()
  }
}

# 메인 루프에서 이미 계산해둔 Lambda_GC(Model=="linear")와 SNP count를 merge
snp_summary_activity <- merge(
  lambda_results[Model == "linear"],
  snp_counts,
  by = c("MAF", "PC_label")
)[order(MAF, PC_label)]

snp_summary_table <- data.table(
  `MAF Threshold` = as.numeric(maf_title_map[snp_summary_activity$MAF]),
  `HWE Threshold` = "1e-6",
  `PCs included`  = snp_summary_activity$Actual_PC,
  `Lambda_GC`      = snp_summary_activity$Lambda_GC,
  `SNPs (p < 5e-8)`  = snp_summary_activity$SNPs_5e8,
  `SNPs (p < 1e-5)`  = snp_summary_activity$SNPs_1e5,
  `SNPs (p < 1e-4)`  = snp_summary_activity$SNPs_1e4
)
head(snp_summary_table)
fwrite(snp_summary_table, file.path(out_dir, "05_SNP_Threshold_Summary_WithActivity.csv"))
cat("✅ SNP threshold 요약표(Activity 포함) 저장:",
    file.path(out_dir, "05_SNP_Threshold_Summary_WithActivity.csv"), "\n")
print(snp_summary_table)

cat("\n모든 조합에 대한 λGC 계산 및 결합 이미지 생성 완료.\n")


target_path <- file.path(base_dir, "05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear")

gwas_dt <- fread(target_path)

# 5e-8 미만 SNP만 필터링, p값 오름차순 정렬
sig_snps <- gwas_dt[P < 5e-8][order(P)]

# 확인
print(sig_snps)

sig_snps_table <- sig_snps[, .(
  ID,
  `Effect allele` = A1,
  `Effect allele frequency` = round(A1_FREQ, 4),
  `β` = round(BETA, 4),
  SE = round(SE, 4),
  `95% CI` = sprintf("%.4f, %.4f", BETA - 1.96 * SE, BETA + 1.96 * SE)
)]

print(sig_snps_table)

fwrite(sig_snps_table, file.path(out_dir, "05_SigSNPs_5e8_001_1e6_PC10_WithActivity_table.csv"))
cat("✅ 저장 완료:", file.path(out_dir, "05_SigSNPs_5e8_001_1e6_PC10_WithActivity_table.csv"), "\n")


# ------------------------------------------------------------------------------
# Exploratory threshold(1e-4) 만족 SNP 목록 저장 (CHR, SNP, POS, A1, P, BETA 형식)
# ------------------------------------------------------------------------------
low_thresh_snps <- gwas_dt[P < low_threshold][order(P)]

low_thresh_table <- low_thresh_snps[, .(
  CHR = `#CHROM`,
  SNP = ID,
  POS,
  A1,
  P,
  BETA
)]

low_thresh_path <- file.path(out_dir, "05_Regression_Result_001_1e6_10_LT.csv")
fwrite(low_thresh_table, low_thresh_path)
cat("✅ Exploratory threshold(1e-4) SNP 목록 저장:", low_thresh_path, "\n")
print(low_thresh_table)
dim(low_thresh_table)

