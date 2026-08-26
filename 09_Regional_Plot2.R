# ==============================================================================
# Regional association plot 재출력
#   - 경로: J:/01_lab_lee/03_GWAS_Pure (base) -> J:\01_lab_lee\00_Picture\04_GWAS_pure\02_LD (output)
#   - [수정 1] 라벨 겹침 방지: index label과 neighbor label이 근접할 때
#     (예: chr24 12787767 <-> 12787779, r²=0.89) 서로 밀어내도록 y_lab 오프셋 로직 보강
#   - [수정 2] index variant의 p-value를 콘솔에 명시적으로 출력하여
#     Table 2와의 수치 일치 여부를 매 실행마다 바로 확인 가능하게 함
#   - [수정 3] 24:12787779는 반드시 "_T_C" LD 파일을 사용 (T:A 파일과 혼동 방지)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(scales)
})
getwd()
base_dir <- "J:/01_lab_lee/03_GWAS_Pure/14_Regression_local"
out_dir  <- "J:/01_lab_lee/00_Picture/04_GWAS_pure/02_LD_local"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

make_regional_plot <- function(
    ld_prefix, assoc_path, out_prefix,
    r2_cut = 0.20, top_k = 5, width_in = 10, height_in = 4
){
  stopifnot(file.exists(paste0(ld_prefix, ".ld.gz")))
  stopifnot(file.exists(assoc_path))
  
  ld <- fread(paste0(ld_prefix, ".ld.gz"))
  need_cols <- c("CHR_A","BP_A","SNP_A","CHR_B","BP_B","SNP_B","R2")
  if (!all(need_cols %in% names(ld))) {
    stop("LD file missing columns: ", paste(setdiff(need_cols, names(ld)), collapse=", "))
  }
  
  idx_chr <- as.integer(ld$CHR_A[1])
  idx_bp  <- as.integer(ld$BP_A[1])
  idx_snp <- as.character(ld$SNP_A[1])
  
  ld2 <- ld %>%
    transmute(CHR = as.integer(CHR_B), BP = as.integer(BP_B),
              SNP = as.character(SNP_B), R2 = as.numeric(R2)) %>%
    arrange(BP)
  
  a <- fread(assoc_path)
  pick_col <- function(nn, cands){ x <- intersect(cands, nn); if (length(x)) x[1] else NA_character_ }
  nn <- names(a)
  map <- list(
    CHR = pick_col(nn, c("#CHROM","CHROM","CHR","Chrom","chrom","chromosome")),
    BP  = pick_col(nn, c("POS","BP","position","bp")),
    SNP = pick_col(nn, c("ID","SNP","RSID","variant_id","MarkerName","VARIANT","MARKER")),
    P   = pick_col(nn, c("P","PVAL","P.value","p","pval","PVALUE","P_BOLT_LMM","P_BOLT_LMM_INF","P1df","P2df"))
  )
  if (any(is.na(unlist(map)))) stop("Could not find CHR/BP/SNP/P in assoc header.")
  
  assoc <- a %>%
    rename(CHR_raw = !!map$CHR, POS = !!map$BP, ID = !!map$SNP, P = !!map$P) %>%
    mutate(
      CHR = gsub("^chr","", CHR_raw, ignore.case = TRUE),
      CHR = gsub("^X$", "39", CHR, ignore.case = TRUE),
      CHR = suppressWarnings(as.integer(CHR)),
      POS = suppressWarnings(as.integer(POS)),
      P   = suppressWarnings(as.numeric(P)),
      ID  = as.character(ID)
    ) %>%
    filter(!is.na(CHR), !is.na(POS), !is.na(P))
  
  # ----------------------------------------------------------------------
  # [수정 2] index variant의 assoc 값을 SNP ID로 직접 조회하여 콘솔에 명시적으로 출력
  #          -> 그림에 실제로 반영될 p-value를 매 실행마다 즉시 검증 가능
  # ----------------------------------------------------------------------
  idx_check <- assoc[ID == idx_snp]
  if (nrow(idx_check) == 0) {
    cat("⚠️ 경고: index SNP(", idx_snp, ")이 assoc 파일에서 ID로 직접 매칭되지 않습니다.\n")
    cat("   CHR:BP 기준으로 재검색합니다.\n")
    idx_check <- assoc[CHR == idx_chr & POS == idx_bp]
  }
  cat(sprintf("\n[Index 확인] %s -> P = %.6g, -log10(P) = %.4f (assoc 파일 기준: %s)\n",
              idx_snp,
              ifelse(nrow(idx_check) > 0, idx_check$P[1], NA),
              ifelse(nrow(idx_check) > 0, -log10(idx_check$P[1]), NA),
              basename(assoc_path)))
  if (nrow(idx_check) > 1) {
    cat("⚠️ 경고: index 위치에 매칭되는 행이 2개 이상입니다 (multi-allelic 가능성). 아래 전체 확인:\n")
    print(idx_check)
  }
  
  j1 <- inner_join(ld2, assoc %>% select(CHR, POS, ID, P), by = c("SNP"="ID"))
  left_ld <- anti_join(ld2, j1 %>% distinct(SNP), by = "SNP")
  j2 <- inner_join(left_ld, assoc %>% select(CHR, POS, P), by = c("CHR"="CHR", "BP"="POS"))
  
  df <- bind_rows(mutate(j1, join_key="SNP"), mutate(j2, join_key="CHR+BP")) %>%
    mutate(logP = -log10(pmax(P, .Machine$double.xmin))) %>%
    arrange(BP)
  if (nrow(df)==0) stop("No overlap between LD and assoc; check IDs or build.")
  
  message(sprintf("Join summary: by SNP = %d, by CHR+BP = %d, total = %d",
                  nrow(j1), nrow(j2), nrow(df)))
  
  idx_core <- sub(":[ACGT]+:[ACGT]+$", "", idx_snp, perl = TRUE)
  lab_df <- df %>%
    filter(!is.na(R2)) %>%
    filter(!(SNP == idx_snp |
               sub(":[ACGT]+:[ACGT]+$", "", SNP, perl = TRUE) == idx_core |
               (BP == idx_bp & R2 >= 0.99))) %>%
    arrange(desc(R2), BP) %>%
    slice_head(n = top_k) %>%
    mutate(lbl = sprintf("%s (r²=%.2f)", SNP, R2))
  
  # ----------------------------------------------------------------------
  # [수정 1] 라벨 겹침 방지
  #   - index label은 항상 그림 상단(y=Inf 근처)에 고정되므로,
  #     neighbor label이 index와 같은 BP 부근(±거리 임계값 이내)에 있으면
  #     추가로 더 큰 y 오프셋을 주어 물리적으로 떨어뜨림
  # ----------------------------------------------------------------------
  yr <- range(df$logP, na.rm = TRUE); y_span <- diff(yr)
  if (nrow(lab_df) > 0) {
    bp_range <- diff(range(df$BP, na.rm = TRUE))
    near_index <- abs(lab_df$BP - idx_bp) < (bp_range * 0.02)  # index와 아주 가까운 이웃
    
    base_offset <- rep(c(0.18, -0.18, 0.32, -0.32, 0.46), length.out = nrow(lab_df)) * y_span
    extra_offset <- ifelse(near_index, 0.25 * y_span, 0)  # index 근접 라벨은 추가로 더 밀어냄
    
    lab_df$y_lab <- lab_df$logP + base_offset + extra_offset * sign(base_offset + 1e-9)
  }
  
  cols <- c("#D0D0D0","#46ACC8","#E17C05","#B2182B")
  
  # [수정 4] index 라벨 y위치: Inf 대신 실제 데이터 범위 기반 고정 좌표 사용
  #          (Inf + clip="off" 조합에서 라벨이 렌더링 시 누락되는 문제 해결)
  label_y_top <- max(df$logP, na.rm = TRUE) * 1.12
  plot_ylim_top <- max(df$logP, na.rm = TRUE) * 1.20
  
  p <- ggplot(df, aes(x = BP, y = logP)) +
    geom_point(aes(color = R2), size = 1.4, alpha = 0.95, na.rm = TRUE) +
    scale_color_gradientn(colors = cols, limits = c(0,1), name = expression(r^2)) +
    geom_vline(xintercept = idx_bp, linetype = "dashed", linewidth = 0.5) +
    # index point를 명시적으로 별도 강조 (검게 테두리 친 채움 없는 큰 다이아몬드)
    geom_point(data = df %>% filter(SNP == idx_snp | (BP == idx_bp & R2 >= 0.99)),
               aes(x = BP, y = logP), shape = 23, size = 3.5, stroke = 1.0,
               colour = "black", fill = "red") +
    annotate("label", x = idx_bp, y = label_y_top, vjust = 0, label.size = 0.3,
             label = paste0("Index: ", idx_snp), size = 3.3, fill = "white") +
    coord_cartesian(ylim = c(min(df$logP, na.rm = TRUE) - 0.1, plot_ylim_top), clip = "off") +
    labs(
      title   = sprintf("Regional association plot — chr%s:%s", idx_chr, format(idx_bp, big.mark=",")),
      x       = sprintf("Genomic position on chr%s (bp)", idx_chr),
      y       = expression(-log[10](p)),
      caption = sprintf("LD source: %s | Assoc: %s | Points: %d",
                        basename(paste0(ld_prefix, ".ld.gz")), basename(assoc_path), nrow(df))
    ) +
    scale_x_continuous(labels = label_number(big.mark=",")) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face="bold", size=13, margin=margin(b=6)),
      axis.title = element_text(size=12),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 12, 10, 10),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  
  if (nrow(lab_df) > 0) {
    p <- p +
      geom_point(data = lab_df, aes(x = BP, y = logP),
                 shape = 21, size = 3.0, stroke = 0.7, fill = NA, colour = "black") +
      geom_segment(data = lab_df,
                   aes(x = BP, xend = BP, y = logP, yend = y_lab),
                   linewidth = 0.4, linetype = "dotted", alpha = 0.7) +
      geom_label(data = lab_df, aes(x = BP, y = y_lab, label = lbl),
                 size = 3.1, label.size = 0.25,
                 label.padding = grid::unit(0.10, "lines"),
                 fill = "white")
  }
  
  png_file <- paste0(out_prefix, ".png")
  pdf_file <- paste0(out_prefix, ".pdf")
  ggsave(filename = png_file, plot = p, width = width_in, height = height_in, dpi = 600, type = "cairo")
  ggsave(filename = pdf_file, plot = p, width = width_in, height = height_in, useDingbats = FALSE)
  message("✅ Saved: ", png_file); message("✅ Saved: ", pdf_file)
  
  hits <- df %>%
    mutate(side = ifelse(BP < idx_bp, "left", "right")) %>%
    filter(!is.na(R2), R2 >= r2_cut) %>%
    arrange(desc(R2), BP) %>%
    transmute(SNP, CHR, BP, `-log10P` = round(logP,3), R2 = round(R2,4), side, join_key)
  csv_file <- paste0(out_prefix, "_LD_hits_r2_", gsub("\\.","", format(r2_cut)), ".csv")
  fwrite(hits, csv_file)
  message("📝 Exported: ", csv_file)
  
  invisible(list(plot=p, data=df, hits=hits, index=list(chr=idx_chr, bp=idx_bp, snp=idx_snp)))
}

# ------------------------------------------------------------------------------
# 5개 locus 전체 재출력 (경로는 base_dir/out_dir 기준으로 통일)
# chr24:12787779는 반드시 "_T_C" LD 파일을 사용
# ------------------------------------------------------------------------------
loci <- list(
  list(ld = "09_LD_5_25529184_C_T",     out = "07_Regional_chr5_25529184"),
  list(ld = "09_LD_24_12787767_A_G",    out = "08_Regional_chr24_12787767"),
  list(ld = "09_LD_24_12787779_T_C",    out = "09_Regional_chr24_12787779_T_C"),
  list(ld = "09_LD_26_13787793_A_G",    out = "10_Regional_chr26_13787793"),
  list(ld = "09_LD_27_44543118_A_T",    out = "11_Regional_chr27_44543118")
)

assoc_path <- file.path(base_dir, "05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear")

for (locus in loci) {
  cat("\n==============================\n")
  cat("처리 중:", locus$ld, "\n")
  cat("==============================\n")
  
  ld_prefix  <- file.path(base_dir, locus$ld)
  out_prefix <- file.path(out_dir, locus$out)
  
  result <- tryCatch(
    make_regional_plot(ld_prefix, assoc_path, out_prefix,
                       r2_cut = 0.20, top_k = 5, width_in = 10, height_in = 4),
    error = function(e) {
      cat("❌ 에러 발생:", conditionMessage(e), "\n")
      NULL
    }
  )
}

cat("\n모든 locus 재출력 완료.\n")