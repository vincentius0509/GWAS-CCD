getwd()

library(data.table)

# ------------------------------------------------------------------
# DAVID output (.txt) -> formatted supplementary-table CSV
# 오직 아래 4개 txt 파일의 컬럼(Count, List Total, User Ids 등)만 사용.
# %는 각 파일 내부의 Count / List Total * 100 으로 계산 (외부 참조 없음).
# ------------------------------------------------------------------

base_dir <- "J:/01_lab_lee/03_GWAS_Pure/14_Regression_local"

files <- c(
  GOTERM_BP = "08_001_1e6_10_LT_GOTERM_BP_DIRECT.txt",
  GOTERM_CC = "08_001_1E6_10_LT_GOTERM_CC_DIRECT.txt",
  GOTERM_MF = "08_001_1E6_10_LT_GOTERM_MF_DIRECT.txt",
  KEGG      = "08_001_1E6_10_LT_KEGG_PATHWAY.txt"
)

for (nm in names(files)) {
  
  in_path  <- file.path(base_dir, files[[nm]])
  out_path <- file.path(base_dir, paste0("08_",nm, "_formatted.csv"))
  
  dt <- fread(in_path, sep = "\t", header = TRUE)
  
  # User Ids -> Genes (쉼표+공백으로 재구성)
  dt[, Genes := gsub(",", ", ", `User Ids`)]
  
  # % = Count / List Total * 100 (해당 txt 파일 내부 값만 사용)
  dt[, PercentCol := round(Count / `List Total` * 100, 1)]
  
  out <- dt[, .(
    Term              = `Full Term`,
    Count             = Count,
    `%`               = PercentCol,
    `P-value`         = `P-Value`,
    Genes             = Genes,
    `List Total`      = `List Total`,
    `Pop Hits`        = `Pop Hits`,
    `Pop Total`       = `Pop Total`,
    `Fold Enrichment` = `Fold Enrichment`,
    Bonferroni        = Bonferroni,
    Benjamini         = Benjamini,
    FDR               = FDR
  )]
  
  out <- out[order(`P-value`)]
  
  fwrite(out, out_path)
  cat(sprintf("[%s] %d rows -> %s\n", nm, nrow(out), out_path))
}
