library(data.table)
library(dplyr)
#############################################
# 1. Load GWAS (PLINK GLM) result file
#############################################

glm_file <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear"
gwas_df <- fread(glm_file, header = TRUE, sep = "\t", data.table = FALSE)
head(gwas_df)

gwas_df[gwas_df$ID == "24:12787767:A:G", ]
gwas_df[gwas_df$ID == "24:12787779:T:C", ]

#############################################
# 2. User-provided SNP list
#############################################
snp_list <- c(
  "X5.25529184.C.T_C",
  "X27.44543118.A.T_A",
  "X24.12787767.A.G_A",
  "X24.12787779.T.C_T",
  "X23.32288824.G.A_G"
)

#############################################
# 3. Parse SNP list (convert to chr-pos-ref-alt)
#############################################
parse_snp <- function(x) {
  x_clean <- sub("^X", "", x)
  key <- strsplit(x_clean, "_")[[1]][1]
  parts <- strsplit(key, "\\.")[[1]]
  
  data.frame(
    SNP_raw = x,
    CHROM = parts[1],
    POS = as.numeric(parts[2]),
    REF_user = parts[3],
    ALT_user = parts[4],
    stringsAsFactors = FALSE
  )
}

parsed <- bind_rows(lapply(snp_list, parse_snp))

#############################################
# 4. Merge with GWAS results
#############################################
matched <- parsed %>%
  left_join(gwas_df,
            by = c("CHROM", "POS"),
            multiple = "all") %>%
  # 필터조건: REF/ALT 동일 또는 REF/ALT 뒤집힌 경우 허용
  filter(
    (REF_user == REF & ALT_user == ALT) |
    (REF_user == ALT & ALT_user == REF)
  )

print(matched)

#############################################
# 5. Biomart format output
#############################################
biomart_df <- matched %>%
  mutate(
    Start = POS,
    End = POS,
    Ref = REF_user,
    Alt = ALT_user
  ) %>%
  select(Chr = CHROM, Start, End, Ref, Alt, SNP_ID = SNP_raw)

print(biomart_df)

# Save
write.table(biomart_df,
            file = "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/11_biomart.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
