## 1) 경로 지정
tar_file <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/13_eQTL/GTEx_Analysis_v8_eQTL.tar"
out_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/13_eQTL"

setwd(out_dir)

## 2) tar 안의 Brain 관련 파일 목록 얻기
brain_files <- system(
  paste(
    "tar -tf", shQuote(tar_file),
    "| grep -E 'GTEx_Analysis_v8_eQTL/Brain_.*signif_variant_gene_pairs.txt.gz'"
  ),
  intern = TRUE
)

## 확인
print(brain_files)

## 3) Brain 파일들만 선택적으로 추출
if (length(brain_files) > 0) {
  cmd <- paste(
    "tar -xvf",
    shQuote(tar_file),
    paste(shQuote(brain_files), collapse = " ")
  )
  system(cmd)
} else {
  stop("No Brain_*signif_variant_gene_pairs files found in tar")
}

library(data.table)
library(ggplot2)
library(stringr)

# dir.create("plots", showWarnings = FALSE)

files <- list.files(
  "GTEx_Analysis_v8_eQTL",
  pattern = "^Brain_.*signif_variant_gene_pairs.txt.gz$",
  full.names = TRUE
)

# 파일명에서 tissue 이름 추출
get_tissue_name <- function(x) {
  x <- basename(x)
  x <- sub("^Brain_", "", x)
  x <- sub("\\.v8\\.signif_variant_gene_pairs\\.txt\\.gz$", "", x)
  x
}
tissues <- sapply(files, get_tissue_name)
names(files) <- tissues
print(names(files))

files[["Amygdala"]]
fread(files[["Amygdala"]])

# 2) target gene 설정 (GTEx query 결과 기반)
target_gene_label <- "CACNA1C"
ensg_core <- "ENSG00000151067"  # CACNA1C core (버전 제외)


# 3) 각 tissue 파일에서 CACNA1C (gene_id core 매칭)만 로딩
dt_list <- lapply(names(files), function(tissue) {
  dt <- fread(files[[tissue]])

  # GTEx v8 signif 파일에는 gene_name이 없고 gene_id가 있음
  if (!all(c("gene_id", "pval_nominal", "slope", "variant_id") %in% names(dt))) {
    stop("Expected columns missing in: ", files[[tissue]],
         "\nFound columns: ", paste(names(dt), collapse = ", "))
  }

  # 버전(.23 등) 무시: ENSG core로 매칭
  dt <- dt[startsWith(gene_id, paste0(ensg_core, ".")) | gene_id == ensg_core]

  if (nrow(dt) == 0) return(NULL)

  dt[, Tissue := tissue]
  dt
})
dt_list <- Filter(Negate(is.null), dt_list)

gtex <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
cat("Total significant eQTL rows for", target_gene_label, ":", nrow(gtex), "\n")

wrap_subtitle <- function(x, width = 70) {
  str_wrap(x, width = width)
}

save_pub <- function(plot, filename, w=7.5, h=6, dpi=600) {
  ggsave(filename, plot,
         width = w, height = h, units = "in",
         dpi = dpi, bg = "white",
         limitsize = FALSE)
}

dist_df <- gtex[, .(
  Tissue,
  NES = slope,
  pval = pval_nominal
)]

dist_n <- dist_df[, .N, by = Tissue][order(-N)]
min_n <- 1

# tissue 정렬: median NES 기준
ord <- dist_df[, .(med = median(NES)), by = Tissue][order(med)]$Tissue
dist_df[, Tissue := factor(Tissue, levels = ord)]
subtxt <- paste0("Gene: ", target_gene_label, " (", ensg_core, ".*); significant variant–gene pairs (per tissue n ≥ ", min_n, ")")

p_dist <- ggplot(dist_df, aes(x = NES, y = Tissue)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(
    title = "GTEx v8 brain cis-eQTL NES distribution",
    subtitle = wrap_subtitle(subtxt, 75),
    x = "Normalized Effect Size (NES)",
    y = ""
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = margin(t=30, r=20, b=15, l=15)
  )

print(p_dist)
save_pub(p_dist, file.path("plots", paste0("GTEx_v8_", target_gene_label, "_NES_distribution_brain.png")), w=8.5, h=6)

count_tbl <- gtex[, .(
  n_total = .N,
  n_pos = sum(slope > 0),
  n_neg = sum(slope < 0),
  n_zero = sum(slope == 0)
), by = Tissue]

# 비율도 함께 (논문에 유용)
count_tbl[, `:=`(
  frac_pos = n_pos / n_total,
  frac_neg = n_neg / n_total
)]

# 보기 좋게 정렬: (1) n_total 큰 순 or (2) frac_pos 순
setorder(count_tbl, -n_total)      # 표본 많은 tissue 먼저
# setorder(count_tbl, frac_pos)    # NES>0 비율 낮은 순(원하면 이걸로)

print(count_tbl)

target_gene_label <- "PLCB4"
ensg_core <- "ENSG00000101333"

dt_list <- lapply(names(files), function(tissue) {
  dt <- fread(files[[tissue]])

  # GTEx v8 signif 파일에는 gene_name이 없고 gene_id가 있음
  if (!all(c("gene_id", "pval_nominal", "slope", "variant_id") %in% names(dt))) {
    stop("Expected columns missing in: ", files[[tissue]],
         "\nFound columns: ", paste(names(dt), collapse = ", "))
  }

  # 버전(.23 등) 무시: ENSG core로 매칭
  dt <- dt[startsWith(gene_id, paste0(ensg_core, ".")) | gene_id == ensg_core]

  if (nrow(dt) == 0) return(NULL)

  dt[, Tissue := tissue]
  dt
})
dt_list <- Filter(Negate(is.null), dt_list)

gtex <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
cat("Total significant eQTL rows for", target_gene_label, ":", nrow(gtex), "\n")

wrap_subtitle <- function(x, width = 70) {
  str_wrap(x, width = width)
}

save_pub <- function(plot, filename, w=7.5, h=6, dpi=600) {
  ggsave(filename, plot,
         width = w, height = h, units = "in",
         dpi = dpi, bg = "white",
         limitsize = FALSE)
}

dist_df <- gtex[, .(
  Tissue,
  NES = slope,
  pval = pval_nominal
)]

dist_n <- dist_df[, .N, by = Tissue][order(-N)]
min_n <- 1

# tissue 정렬: median NES 기준
ord <- dist_df[, .(med = median(NES)), by = Tissue][order(med)]$Tissue
dist_df[, Tissue := factor(Tissue, levels = ord)]
subtxt <- paste0("Gene: ", target_gene_label, " (", ensg_core, ".*); significant variant–gene pairs (per tissue n ≥ ", min_n, ")")

p_dist <- ggplot(dist_df, aes(x = NES, y = Tissue)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(
    title = "GTEx v8 brain cis-eQTL NES distribution",
    subtitle = wrap_subtitle(subtxt, 75),
    x = "Normalized Effect Size (NES)",
    y = ""
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = margin(t=30, r=20, b=15, l=15)
  )

print(p_dist)
save_pub(p_dist, file.path("plots", paste0("GTEx_v8_", target_gene_label, "_NES_distribution_brain.png")), w=8.5, h=6)

count_tbl <- gtex[, .(
  n_total = .N,
  n_pos = sum(slope > 0),
  n_neg = sum(slope < 0),
  n_zero = sum(slope == 0)
), by = Tissue]

# 비율도 함께 (논문에 유용)
count_tbl[, `:=`(
  frac_pos = n_pos / n_total,
  frac_neg = n_neg / n_total
)]

# 보기 좋게 정렬: (1) n_total 큰 순 or (2) frac_pos 순
setorder(count_tbl, -n_total)      # 표본 많은 tissue 먼저
# setorder(count_tbl, frac_pos)    # NES>0 비율 낮은 순(원하면 이걸로)

print(count_tbl)

target_gene_label <- "GUCY1A2"
ensg_core <- "ENSG00000152402"

dt_list <- lapply(names(files), function(tissue) {
  dt <- fread(files[[tissue]])

  # GTEx v8 signif 파일에는 gene_name이 없고 gene_id가 있음
  if (!all(c("gene_id", "pval_nominal", "slope", "variant_id") %in% names(dt))) {
    stop("Expected columns missing in: ", files[[tissue]],
         "\nFound columns: ", paste(names(dt), collapse = ", "))
  }

  # 버전(.23 등) 무시: ENSG core로 매칭
  dt <- dt[startsWith(gene_id, paste0(ensg_core, ".")) | gene_id == ensg_core]

  if (nrow(dt) == 0) return(NULL)

  dt[, Tissue := tissue]
  dt
})
dt_list <- Filter(Negate(is.null), dt_list)

gtex <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
cat("Total significant eQTL rows for", target_gene_label, ":", nrow(gtex), "\n")

wrap_subtitle <- function(x, width = 70) {
  str_wrap(x, width = width)
}

save_pub <- function(plot, filename, w=7.5, h=6, dpi=600) {
  ggsave(filename, plot,
         width = w, height = h, units = "in",
         dpi = dpi, bg = "white",
         limitsize = FALSE)
}

dist_df <- gtex[, .(
  Tissue,
  NES = slope,
  pval = pval_nominal
)]

dist_n <- dist_df[, .N, by = Tissue][order(-N)]
min_n <- 1

# tissue 정렬: median NES 기준
ord <- dist_df[, .(med = median(NES)), by = Tissue][order(med)]$Tissue
dist_df[, Tissue := factor(Tissue, levels = ord)]
subtxt <- paste0("Gene: ", target_gene_label, " (", ensg_core, ".*); significant variant–gene pairs (per tissue n ≥ ", min_n, ")")

p_dist <- ggplot(dist_df, aes(x = NES, y = Tissue)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(
    title = "GTEx v8 brain cis-eQTL NES distribution",
    subtitle = wrap_subtitle(subtxt, 75),
    x = "Normalized Effect Size (NES)",
    y = ""
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = margin(t=30, r=20, b=15, l=15)
  )

print(p_dist)
save_pub(p_dist, file.path("plots", paste0("GTEx_v8_", target_gene_label, "_NES_distribution_brain.png")), w=8.5, h=6)

count_tbl <- gtex[, .(
  n_total = .N,
  n_pos = sum(slope > 0),
  n_neg = sum(slope < 0),
  n_zero = sum(slope == 0)
), by = Tissue]

# 비율도 함께 (논문에 유용)
count_tbl[, `:=`(
  frac_pos = n_pos / n_total,
  frac_neg = n_neg / n_total
)]

# 보기 좋게 정렬: (1) n_total 큰 순 or (2) frac_pos 순
setorder(count_tbl, -n_total)      # 표본 많은 tissue 먼저
# setorder(count_tbl, frac_pos)    # NES>0 비율 낮은 순(원하면 이걸로)

print(count_tbl)
