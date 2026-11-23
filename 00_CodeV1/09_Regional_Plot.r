# ld_prefix  <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/09_LD_5_25529184"
# image_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/09_LD_5_25529184"

# ld_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/09_LD_27_44543118"
# image_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/09_LD_27_44543118"

# ld_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/09_LD_24_12787767"
# image_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/09_LD_24_12787767"

ld_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/09_LD_24_12787779_T_C"
image_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/09_LD_24_12787779"

# ld_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/09_LD_23_32288824"
# image_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/09_LD_23_32288824"

bim_path   <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6.bim"

assoc_path <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear"

out_png    <- paste0(image_prefix, "_region.png")
out_pdf    <- paste0(image_prefix, "_region.pdf")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

stopifnot(file.exists(paste0(ld_prefix, ".ld.gz")))
ld <- fread(paste0(ld_prefix, ".ld.gz"))

need_cols <- c("CHR_A","BP_A","SNP_A","CHR_B","BP_B","SNP_B","R2")
if (!all(need_cols %in% names(ld))) {
  stop("LD 파일에 필요한 컬럼이 없습니다: ", paste(setdiff(need_cols, names(ld)), collapse=", "))
}

idx_chr <- as.integer(ld$CHR_A[1])
idx_bp  <- as.integer(ld$BP_A[1])
idx_snp <- as.character(ld$SNP_A[1])

bim <- fread(bim_path, col.names = c("CHR","SNP","CM","BP","A1","A2"))

ld2 <- ld %>%
  transmute(
    CHR = as.integer(CHR_B),
    BP  = as.integer(BP_B),
    SNP = as.character(SNP_B),
    R2  = as.numeric(R2)
  ) %>%
  arrange(BP)

stopifnot(file.exists(assoc_path))
a <- fread(assoc_path)

std_cols <- function(dt) {
  nn <- names(dt)
  chr_col <- c("#CHROM","CHROM","CHR","Chrom","chrom","chromosome")
  pos_col <- c("POS","BP","position","bp")
  snp_col <- c("ID","SNP","RSID","variant_id","MarkerName","VARIANT","MARKER")
  p_col   <- c("P","PVAL","P.value","p","pval","PVALUE","P_BOLT_LMM","P_BOLT_LMM_INF","P1df","P2df")
  pick <- function(cands) { x <- intersect(cands, nn); if (length(x)) x[1] else NA_character_ }
  data.frame(
    CHR = pick(chr_col),
    BP  = pick(pos_col),
    SNP = pick(snp_col),
    P   = pick(p_col),
    stringsAsFactors = FALSE
  )
}
m <- std_cols(a)
if (any(is.na(m))) {
  stop("assoc 파일의 표준 컬럼(CHR/BP/SNP/P) 일부를 찾지 못했습니다. 파일 헤더를 확인하세요.")
}

assoc <- a %>%
  rename(CHR_raw = !!m$CHR, POS = !!m$BP, ID = !!m$SNP, P = !!m$P) %>%
  mutate(
    CHR = as.character(CHR_raw),
    CHR = gsub("^chr", "", CHR, ignore.case = TRUE),
    CHR = gsub("^X$", "39", CHR, ignore.case = TRUE),  # 개 유전체에서 X→39
    CHR = suppressWarnings(as.integer(CHR)),
    POS = suppressWarnings(as.integer(POS)),
    P   = suppressWarnings(as.numeric(P)),
    ID  = as.character(ID)
  ) %>%
  filter(!is.na(CHR), !is.na(POS), !is.na(P))

# 1) SNP 매칭
j1 <- ld2 %>% inner_join(assoc %>% select(CHR, POS, ID, P), by = c("SNP" = "ID"))
# 2) SNP로 못 붙은 LD 행은 CHR+BP로 보충
left_ld <- anti_join(ld2, j1 %>% select(SNP) %>% distinct(), by = "SNP")

j2 <- left_ld %>%
  inner_join(assoc %>% select(CHR, POS, P),
             by = c("CHR" = "CHR", "BP" = "POS"))

df <- bind_rows(
        j1 %>% mutate(join_key = "SNP"),
        j2 %>% mutate(join_key = "CHR+BP")
      ) %>%
      mutate(logP = -log10(pmax(P, .Machine$double.xmin))) %>%
      arrange(BP)

message(sprintf("조인 요약: SNP기반 %d개, CHR+BP기반 %d개, 총 %d개",
                nrow(j1), nrow(j2), nrow(df)))

if (nrow(df) == 0) {
  stop("LD와 assoc의 교집합이 없습니다. SNP ID 규칙 또는 좌표계를 확인하세요.")
}

p <- ggplot(df, aes(x = BP, y = logP)) +
  geom_point(aes(color = R2), alpha = 0.9, size = 1.8) +
  scale_color_gradient(low = "grey70", high = "red", name = expression(R^2)) +
  geom_vline(xintercept = idx_bp, linetype = "dashed") +
  annotate("text", x = idx_bp, y = Inf, vjust = 1.5,
           label = paste0("Index: ", idx_snp), size = 3.2) +
  scale_x_continuous(labels = label_number(big.mark=",")) +
  labs(
    title = paste0("Regional plot (Chr ", idx_chr, ":", format(idx_bp, big.mark=","), ")"),
    x = paste0("Position on chr", idx_chr, " (bp)"),
    y = expression(-log[10](P)),
    caption = paste0("LD: ", basename(paste0(ld_prefix, ".ld.gz")),
                     " | Assoc: ", basename(assoc_path),
                     " | n=", nrow(assoc))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    legend.background= element_rect(fill = "white", colour = NA),
    legend.key       = element_rect(fill = "white", colour = NA)
  )

cap <- sprintf("LD: %s | Assoc: %s | plotted=%d",
               basename(paste0(ld_prefix, ".ld.gz")),
               basename(assoc_path),
               nrow(df))
p <- p + labs(caption = cap)

top_k <- 5
y_col <- if ("logP" %in% names(df)) "logP" else "R2"

idx_core <- sub(":[ACGT]+:[ACGT]+$", "", idx_snp, perl = TRUE)

lab_df <- df %>%
  dplyr::filter(!is.na(R2)) %>%
      dplyr::filter(!(SNP == idx_snp |
                  sub(":[ACGT]+:[ACGT]+$", "", SNP, perl = TRUE) == idx_core |
                  (BP == idx_bp & R2 >= 0.99))) %>%   # 추가 안전장치
  dplyr::arrange(dplyr::desc(R2), BP) %>%
  dplyr::slice_head(n = top_k) %>%
  dplyr::mutate(lbl = sprintf("%s(%.2f)", SNP, R2))

# 겹침 줄이기 위한 y 오프셋
yr <- range(df[[y_col]], na.rm = TRUE); y_span <- diff(yr)
offsets <- seq(0.12, 0.28, length.out = nrow(lab_df)) * y_span
signs <- rep(c(1, -1), length.out = nrow(lab_df))
lab_df$y_lab <- lab_df[[y_col]] + signs * offsets

p <- p +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 40, 10, 10)) +
  # 점 강조(링)
  geom_point(
    data = lab_df,
    aes(x = BP, y = .data[[y_col]]),
    shape = 21, size = 3.2, stroke = 0.7, fill = NA, color = "black"
  ) +
  # 점→라벨 연결선
  geom_segment(
    data = lab_df,
    aes(x = BP, xend = BP, y = .data[[y_col]], yend = y_lab),
    linetype = "dotted", linewidth = 0.4, alpha = 0.7
  ) +
  # 라벨
  geom_label(
    data = lab_df,
    aes(x = BP, y = y_lab, label = lbl),
    size = 3.2, label.size = 0.2,
    label.padding = grid::unit(0.12, "lines"),
    fill = "white"
  )

ggsave(out_png, p, width = 9, height = 5.5, dpi = 300)
ggsave(out_pdf, p, width = 9, height = 5.5)
message("✅ Saved: ", out_png)
message("✅ Saved: ", out_pdf)

r2_cut <- 0.20

hits <- df %>%
  dplyr::mutate(side = ifelse(BP < idx_bp, "left", "right")) %>%
  dplyr::filter(!is.na(R2), R2 >= r2_cut) %>%
  dplyr::arrange(dplyr::desc(R2), BP) %>%
  dplyr::select(SNP, CHR, BP, logP, R2, side, join_key)

cat(sprintf("R^2 ≥ %.2f 변이 개수: %d\n", r2_cut, nrow(hits)))
print(utils::head(hits, 20)) 


####################################################################################################################################
#####논문용 figure
####################################################################################################################################


# ---------- Regional association (LD) plot for publication ----------
# Dependencies: data.table, dplyr, ggplot2, scales
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(scales)
})

make_regional_plot <- function(
  ld_prefix,                   # path prefix of PLINK --ld output (".ld.gz" expected)
  assoc_path,                  # PLINK assoc (e.g., *.glm.linear)
  out_prefix,                  # output prefix (no extension)
  r2_cut   = 0.20,             # label/filter threshold for R^2
  top_k    = 5,                # how many high-LD neighbors to annotate
  width_in = 10, height_in = 4 # figure size in inches (landscape, 2-column friendly)
){
  stopifnot(file.exists(paste0(ld_prefix, ".ld.gz")))
  stopifnot(file.exists(assoc_path))

  # --- read LD ---
  ld <- fread(paste0(ld_prefix, ".ld.gz"))
  need_cols <- c("CHR_A","BP_A","SNP_A","CHR_B","BP_B","SNP_B","R2")
  if (!all(need_cols %in% names(ld))) {
    stop("LD file missing columns: ", paste(setdiff(need_cols, names(ld)), collapse=", "))
  }

  # index variant from LD header
  idx_chr <- as.integer(ld$CHR_A[1])
  idx_bp  <- as.integer(ld$BP_A[1])
  idx_snp <- as.character(ld$SNP_A[1])

  ld2 <- ld %>%
    transmute(
      CHR = as.integer(CHR_B),
      BP  = as.integer(BP_B),
      SNP = as.character(SNP_B),
      R2  = as.numeric(R2)
    ) %>% arrange(BP)

  # --- read assoc robustly (auto-detect standard columns) ---
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
      CHR = gsub("^X$", "39", CHR, ignore.case = TRUE),  # dog genome: X->39
      CHR = suppressWarnings(as.integer(CHR)),
      POS = suppressWarnings(as.integer(POS)),
      P   = suppressWarnings(as.numeric(P)),
      ID  = as.character(ID)
    ) %>%
    filter(!is.na(CHR), !is.na(POS), !is.na(P))

  # --- join LD with assoc (by SNP, then by CHR:BP fallback) ---
  j1 <- inner_join(ld2, assoc %>% select(CHR, POS, ID, P), by = c("SNP"="ID"))
  left_ld <- anti_join(ld2, j1 %>% distinct(SNP), by = "SNP")
  j2 <- inner_join(left_ld, assoc %>% select(CHR, POS, P), by = c("CHR"="CHR", "BP"="POS"))

  df <- bind_rows(
          mutate(j1, join_key="SNP"),
          mutate(j2, join_key="CHR+BP")
        ) %>%
        mutate(logP = -log10(pmax(P, .Machine$double.xmin))) %>%
        arrange(BP)
  if (nrow(df)==0) stop("No overlap between LD and assoc; check IDs or build.")

  message(sprintf("Join summary: by SNP = %d, by CHR+BP = %d, total = %d",
                  nrow(j1), nrow(j2), nrow(df)))

  # --- choose top LD neighbors to annotate (avoid index variant & perfect duplicates) ---
  idx_core <- sub(":[ACGT]+:[ACGT]+$", "", idx_snp, perl = TRUE)
  lab_df <- df %>%
    filter(!is.na(R2)) %>%
    filter(!(SNP == idx_snp |
             sub(":[ACGT]+:[ACGT]+$", "", SNP, perl = TRUE) == idx_core |
             (BP == idx_bp & R2 >= 0.99))) %>%
    arrange(desc(R2), BP) %>%
    slice_head(n = top_k) %>%
    mutate(lbl = sprintf("%s (r²=%.2f)", SNP, R2))

  # y-offsets to reduce overlap
  yr <- range(df$logP, na.rm = TRUE); y_span <- diff(yr)
  if (nrow(lab_df) > 0) {
    lab_df$y_lab <- lab_df$logP + rep(c(0.18, -0.18), length.out = nrow(lab_df)) * y_span
  }

  # --- palette (color-blind safe gradient) ---
  cols <- c("#D0D0D0","#46ACC8","#E17C05","#B2182B")  # grey → teal → orange → burgundy

  # --- plot ---
  p <- ggplot(df, aes(x = BP, y = logP)) +
    geom_point(aes(color = R2), size = 1.4, alpha = 0.95, na.rm = TRUE) +
    scale_color_gradientn(colors = cols, limits = c(0,1), name = expression(r^2)) +
    geom_vline(xintercept = idx_bp, linetype = "dashed", linewidth = 0.5) +
    annotate("label", x = idx_bp, y = Inf, vjust = 1.3, label.size = 0,
             label = paste0("Index: ", idx_snp), size = 3.3) +
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

  # ring highlight + leader lines + labels
  if (nrow(lab_df) > 0) {
    p <- p +
      coord_cartesian(clip = "off") +
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

  # --- save high-res files (no svglite required) ---
  png_file <- paste0(out_prefix, ".png")
  pdf_file <- paste0(out_prefix, ".pdf")
  ggsave(filename = png_file, plot = p, width = width_in, height = height_in, dpi = 600, type = "cairo")
  ggsave(filename = pdf_file, plot = p, width = width_in, height = height_in, useDingbats = FALSE)
  message("✅ Saved: ", png_file); message("✅ Saved: ", pdf_file)

  # --- export table of high-LD hits (for Methods/Supplementary) ---
  hits <- df %>%
    mutate(side = ifelse(BP < idx_bp, "left", "right")) %>%
    filter(!is.na(R2), R2 >= r2_cut) %>%
    arrange(desc(R2), BP) %>%
    transmute(SNP, CHR, BP, `-log10P` = round(logP,3), R2 = round(R2,3), side, join_key)
  csv_file <- paste0(out_prefix, "_LD_hits_r2_", gsub("\\.","", format(r2_cut)), ".csv")
  fwrite(hits, csv_file)
  message("📝 Exported: ", csv_file)

  invisible(list(plot=p, data=df, hits=hits,
                 index=list(chr=idx_chr, bp=idx_bp, snp=idx_snp)))
}

# ---------------- Example (your chr23 locus) ----------------
# ld_prefix  <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/09_LD_23_32288824"
# ld_prefix  <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/09_LD_5_25529184"
# ld_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/09_LD_27_44543118"
# ld_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/09_LD_24_12787767"
ld_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/09_LD_24_12787779_T_C"
assoc_path <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear"
out_prefix <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/00_ArticleFigure/09_LD_24_12787779"

make_regional_plot(ld_prefix, assoc_path, out_prefix,
                   r2_cut = 0.20, top_k = 5, width_in = 10, height_in = 4)
