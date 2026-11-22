glm_file <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear"

# 패키지 로드
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
library(data.table)

gwas_df <- fread(glm_file, header = TRUE, sep = "\t", data.table = FALSE)
head(gwas_df)

cat("Columns:\n"); print(names(gwas_df))
cat("\nFirst rows:\n"); print(head(gwas_df, 6))

table(gwas_df$TEST)
dim(gwas_df)

is_add <- gwas_df$TEST == "ADD"
no_NA <- !is.na(gwas_df$P)

gwas_clean <- gwas_df[is_add & no_NA,]
dim(gwas_clean)

cat(sprintf("\n원본 행수: %d, ADD+P값 존재 행수: %d\n",
            nrow(gwas_df), nrow(gwas_clean)))

thr <- 1e-4
sig <- subset(gwas_clean, P < thr)
dim(sig)

cat(sprintf("P < %.0e 인 SNP 수: %d\n", thr, nrow(sig)))

head(sig)

write.table(sig$ID, file = "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/10_snplist.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

head(sig[order(sig$P), ], 10)

raw <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/10_snp_selected.txt.raw", header = T, stringsAsFactors = F)
head(raw)
dim(raw)
length(raw$FID) # 3527
length(unique(raw$FID))# 3527

raw$X5.25529184.C.T_C
raw$X27.44543118.A.T_A
raw$X24.12787767.A.G_G
raw$X24.12787779.T.C_C
raw$X23.32288824.G.A_G

raw_select <- raw[, c("FID","IID","PAT","MAT","SEX","PHENOTYPE","X5.25529184.C.T_C",
    "X27.44543118.A.T_A","X24.12787767.A.G_G","X24.12787779.T.C_C","X23.32288824.G.A_G")]
dim(raw_select)
head(raw_select)

cslb_2022 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/SurveyData/v2022/DAP_2022_CSLB_v1.0.csv")
cslb_2023 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/SurveyData/v2023/DAP_2023_CSLB_v1.0.csv")

cslb2022_ed <- cslb_2022[,c("dog_id","cslb_score")]
head(cslb2022_ed)
dim(cslb2022_ed)
length(cslb2022_ed$dog_id) # 59977
length(unique(cslb2022_ed$dog_id)) # 30872

cslb2023_ed <- cslb_2023[,c("dog_id","cslb_score")]
head(cslb2023_ed)
dim(cslb2023_ed)
length(cslb2023_ed$dog_id) # 73898
length(unique(cslb2023_ed$dog_id)) # 33190

cslb2022_agg <- cslb2022_ed[, .(
  max_score = max(cslb_score, na.rm = TRUE)  # 해당 개체의 최고점
), by = dog_id]

length(cslb2022_agg$dog_id)
head(cslb2022_agg)
cslb2022_agg[, CCD_Status := ifelse(max_score >= 50, "CCD", "Control")]
head(cslb2022_agg)


cslb2023_agg <- cslb2023_ed[, .(
  max_score = max(cslb_score, na.rm = TRUE)
), by = dog_id]


length(cslb2023_agg$dog_id) # 33190
head(cslb2023_agg)
cslb2023_agg[, CCD_Status := ifelse(max_score >= 50, "CCD", "Control")]
head(cslb2023_agg)

cslb2022_agg[, IID := as.character(dog_id)]
cslb2023_agg[, IID := as.character(dog_id)]

merged2022_df <- merge(raw_select, cslb2022_agg[, .(IID, CCD_Status)],
                   by = "IID", all.x = TRUE)

cat("병합된 개체 수:", nrow(merged2022_df), "\n")
cat("CCD/Control 분포:\n")
print(table(merged2022_df$CCD_Status, useNA = "ifany"))
    # CCD Control 
    #  60    3467 


merged2023_df <- merge(raw_select, cslb2023_agg[, .(IID, CCD_Status)],
                    by = "IID", all.x = TRUE)

cat("병합된 개체 수:", nrow(merged2023_df), "\n")
cat("CCD/Control 분포:\n")
print(table(merged2023_df$CCD_Status, useNA = "ifany"))
head(merged2023_df)

merged2023_df <- merged2023_df[!is.na(merged2023_df$CCD_Status),]
table(merged2023_df$CCD_Status, useNA = "ifany")
    # CCD Control 
    #  79    3430 

snp_cols <- c("X5.25529184.C.T_C",
              "X27.44543118.A.T_A",
              "X24.12787767.A.G_G",
              "X24.12787779.T.C_C",
              "X23.32288824.G.A_G")

# normalize_geno <- function(x) {
#   y <- as.integer(x)
#   factor(ifelse(y == 0, "AA",
#                 ifelse(y == 1, "AB",
#                        ifelse(y == 2, "BB", NA))),
#          levels = c("AA","AB","BB"))
# }

outdir <- "/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/04_GWAS_pure/10_barplots"
# dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

normalize_geno_CT <- function(x) {
  y <- as.integer(x)
  factor(ifelse(y == 0, "CC",
                ifelse(y == 1, "CT",
                       ifelse(y == 2, "TT", NA))),
         levels = c("CC","CT","TT"))
}

normalize_geno_AT <- function(x) {
  y <- as.integer(x)
  factor(ifelse(y == 0, "AA",
                ifelse(y == 1, "AT",
                       ifelse(y == 2, "TT", NA))),
         levels = c("AA","AT","TT"))
}

normalize_geno_AG <- function(x) {
  y <- as.integer(x)
  factor(ifelse(y == 0, "AA",
                ifelse(y == 1, "AG",
                       ifelse(y == 2, "GG", NA))),
         levels = c("AA","AG","GG"))
}

normalize_geno_TC <- function(x) {
  y <- as.integer(x)
  factor(ifelse(y == 0, "TT",
                ifelse(y == 1, "TC",
                       ifelse(y == 2, "CC", NA))),
         levels = c("TT","TC","CC"))
}

normalize_geno_GA <- function(x) {
  y <- as.integer(x)
  factor(ifelse(y == 0, "GG",
                ifelse(y == 1, "GA",
                       ifelse(y == 2, "AA", NA))),
         levels = c("GG","GA","AA"))
}

head(merged2022_df)
head(merged2023_df)

snp <- snp_cols[5]

df_one <- data.frame(
  Genotype   = normalize_geno_GA(merged2023_df[[snp]]),
  CCD_Status = merged2023_df$CCD_Status
)

tab <- table(df_one$Genotype, df_one$CCD_Status)
prop_tab <- prop.table(tab, 2)
heights <- t(prop_tab)

cols <- c("Control" = "#2AB7CA", "CCD" = "#9D4EDD")
ylim_max <- max(heights, na.rm = TRUE) + 0.1

snp

outfile <- file.path(outdir, paste0(snp, "_proportion_GA_grouped_2023.png"))
png(outfile, width = 1400, height = 900, res = 150)
op <- par(mar = c(6, 6, 4, 10), xpd = TRUE)
bp <- barplot(
  heights,
  beside = TRUE,
  ylim = c(0, ylim_max),
  col  = cols[rownames(heights)],
  border = NA,
  names.arg = colnames(heights),  # CC, CT, TT
  xlab = "Genotype",
  ylab = "Proportion within group",
  main = paste0("23:32288824:G:A — genotype representation by group (v2023) [CCD: 79, Control: 3430]")
)

# 퍼센트 라벨 표시
labs <- paste0(round(as.vector(heights) * 100, 1), "%")
text(x = as.vector(bp), y = as.vector(heights), labels = labs, pos = 3, cex = 0.9)

legend("topright",
       inset = c(-0.2, 0),          # inset의 첫 번째 값(음수)으로 오른쪽 바깥으로 이동
       bty = "n",
       legend = rownames(heights),
       fill = cols[rownames(heights)],
       cex = 1)
par(op)
dev.off()

cat("✅ 저장됨: ", outfile, "\n")


raw_total <- read.table("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/10_snp_selected_total.raw", header = T, stringsAsFactors = F)
dim(raw_total)

raw_total_select <- raw_total[, c("FID","IID","PAT","MAT","SEX","PHENOTYPE","X5.25529184.C.T_C",
    "X27.44543118.A.T_A","X24.12787767.A.G_G","X24.12787779.T.C_C","X23.32288824.G.A_G")]
dim(raw_total_select)

dim(cslb2022_agg)
dim(cslb2023_agg)

merged_total2022_df <- merge(raw_total_select, cslb2022_agg[, .(IID, CCD_Status)],
                   by = "IID", all.x = TRUE)
nrow(merged_total2022_df)
table(merged_total2022_df$CCD_Status)

merged_total2023_df <- merge(raw_total_select, cslb2023_agg[, .(IID, CCD_Status)],
                   by = "IID", all.x = TRUE)

snp_cols
outdir

dim(merged_total2022_df)
table(merged_total2022_df$CCD_Status) # 60  3482

dim(merged_total2023_df)
table(merged_total2023_df$CCD_Status) # 79  3445

snp <- snp_cols[5]

df_one <- data.frame(
  Genotype   = normalize_geno_GA(merged_total2023_df[[snp]]),
  CCD_Status = merged_total2023_df$CCD_Status
)

tab <- table(df_one$Genotype, df_one$CCD_Status)

ft <- fisher.test(tab)           # 두 집단(열) x 3 유전자형(행)
p_overall <- ft$p.value
p_overall_lab <- if (p_overall < 1e-4) {
  sprintf("Fisher p = %.2e", p_overall)
} else {
  sprintf("Fisher p = %.4f", p_overall)
}

prop_tab <- prop.table(tab, 2)
heights <- t(prop_tab)

cols <- c("Control" = "#2AB7CA", "CCD" = "#9D4EDD")
ylim_max <- max(heights, na.rm = TRUE) + 0.1

snp

outfile <- file.path(outdir, paste0(snp, "_totalproportion_GA_grouped_2023.png"))
png(outfile, width = 1400, height = 900, res = 150)
op <- par(mar = c(6, 6, 5.5, 12), xpd = TRUE)
bp <- barplot(
  heights,
  beside = TRUE,
  ylim = c(0, ylim_max),
  col  = cols[rownames(heights)],
  border = NA,
  names.arg = colnames(heights),  # CC, CT, TT
  xlab = "Genotype",
  ylab = "Proportion within group",
  main = paste0("23:32288824:T:C — genotype representation by group (v2023) [CCD: 79, Normal: 3445]")
)

# 퍼센트 라벨 표시
labs <- paste0(round(as.vector(heights) * 100, 1), "%")
text(x = as.vector(bp), y = as.vector(heights), labels = labs, pos = 3, cex = 0.9)

legend("topright",
       inset = c(-0.2, 0),          # inset의 첫 번째 값(음수)으로 오른쪽 바깥으로 이동
       bty = "n",
       legend = rownames(heights),
       fill = cols[rownames(heights)],
       cex = 1)

mid_x <- mean(range(bp))  # x축 중앙 계산
text(x = mid_x, y = ylim_max - 0.02,
     labels = sprintf("Fisher p = %.4f", p_overall),
     cex = 1.2, font = 2, adj = 0.5)

par(op)
dev.off()

cat("✅ 저장됨: ", outfile, "\n")


head(merged_total2022_df)
head(merged_total2023_df)

length(merged_total2022_df$IID)
length(merged_total2023_df$IID)

length(unique(merged_total2022_df$IID))
length(unique(merged_total2023_df$IID))

table(merged_total2022_df$SEX)
table(merged_total2023_df$SEX)

table(merged_total2022_df$CCD_Status)
table(merged_total2023_df$CCD_Status)

library(data.table)

# ---------- 1) SNP 선택 ----------
snp <- snp_cols[5]  

# ---------- 2) 데이터 준비 ----------
setDT(merged_total2023_df)
dt <- merged_total2023_df[, .(IID, SEX, CCD_Status, snp_val = get(snp))]
dt <- dt[!is.na(SEX) & !is.na(CCD_Status) & !is.na(snp_val)]

snp

dt[, SexLabel := factor(as.integer(SEX), levels = c(1,2), labels = c("Male","Female"))]
dt[, CCD_Status := ifelse(CCD_Status %in% c("Control","Normal"), "Normal", "CCD")]
dt[, CCD_Status := factor(CCD_Status, levels = c("CCD","Normal"))]
dt[, Genotype   := normalize_geno_GA(snp_val)]  # -----------------------------------------------> 수정

# ---------- 3) 성별×유전자형 Fisher 검정 (그룹별) ----------
fisher_fmt <- function(p) if (is.na(p)) "NA" else if (p < 1e-4) sprintf("%.2e", p) else sprintf("%.4f", p)

tab_ccd     <- table(dt[CCD_Status=="CCD"]$SexLabel,   dt[CCD_Status=="CCD"]$Genotype)
tab_normal  <- table(dt[CCD_Status=="Normal"]$SexLabel, dt[CCD_Status=="Normal"]$Genotype)

p_ccd    <- tryCatch(fisher.test(tab_ccd)$p.value,    error=function(e) NA_real_)
p_normal <- tryCatch(fisher.test(tab_normal)$p.value, error=function(e) NA_real_)

p_line <- paste0("Fisher (Sex × Genotype) — CCD: p=", fisher_fmt(p_ccd),
                 "   |   Normal: p=", fisher_fmt(p_normal))

# ---------- 4) 스택형 barplot용 집계 ----------
cnt <- dt[, .N, by = .(CCD_Status, Genotype, SexLabel)]
tot_group <- dt[, .N, by = .(CCD_Status)]
cnt <- merge(cnt, tot_group, by = "CCD_Status", suffixes = c("", "_group"))
cnt[, PropSeg := N / N_group]                          # 세그먼트 높이(그룹 내 비율)

# 막대 전체 높이(성별 합친 유전자형 비율)
cnt_tot <- cnt[, .(N = sum(N), N_group = first(N_group)), by = .(CCD_Status, Genotype)]
cnt_tot[, PropBar := N / N_group]

# wide 변환
prop_cast  <- dcast(cnt,     CCD_Status + Genotype ~ SexLabel, value.var = "PropSeg", fill = 0)
count_cast <- dcast(cnt,     CCD_Status + Genotype ~ SexLabel, value.var = "N",       fill = 0)
bar_tot    <- dcast(cnt_tot, CCD_Status + Genotype ~ 1,        value.var = "PropBar", fill = 0); setnames(bar_tot, "1", "PropBar")

# 열 순서(원하면 수정 가능)
ord <- with(prop_cast, order(CCD_Status, Genotype))
prop_cast  <- prop_cast[ord]; count_cast <- count_cast[ord]; bar_tot <- bar_tot[ord]

mat_prop <- t(as.matrix(prop_cast[, .(Male, Female)]))
mat_cnt  <- t(as.matrix(count_cast[, .(Male, Female)]))
colnames(mat_prop) <- paste(prop_cast$CCD_Status, prop_cast$Genotype, sep = "_")
colnames(mat_cnt)  <- colnames(mat_prop)

# ---------- 5) 그리기 ----------
cols <- c("Male" = "#4E79A7", "Female" = "#F28E2B")
outfile <- file.path(outdir, paste0(snp, "_Genotype_SexStacked_withCounts_Fisher_2023.png"))

png(outfile, width = 1400, height = 900, res = 150)
op <- par(mar = c(6,6,6.8,12), xpd = TRUE)

snp

ylim_top <- max(colSums(mat_prop), na.rm = TRUE) + 0.12
bp <- barplot(
  mat_prop, ylim = c(0, ylim_top),
  col = cols, border = NA,
  names.arg = gsub("_", "\n", colnames(mat_prop)),
  xlab = "Group and Genotype",
  ylab = "Proportion within group",
  main = "23:32288824:G:A — genotype × sex distribution (v2023) [CCD vs Normal]"
)

# (A) 막대 전체 비율(%)
bar_tops <- colSums(mat_prop)
text(bp, bar_tops, labels = paste0(round(bar_tops*100,1), "%"), pos = 3, cex = 0.95, font = 2)

# (B) 세그먼트 내부 n 라벨(파랑/주황)
cum_prop    <- apply(mat_prop, 2, cumsum)
seg_centers <- rbind(mat_prop[1,]/2, cum_prop[1,] + mat_prop[2,]/2)
lab_male    <- ifelse(mat_cnt[1,] > 0, paste0("n=", mat_cnt[1,]), "")
lab_female  <- ifelse(mat_cnt[2,] > 0, paste0("n=", mat_cnt[2,]), "")
text(bp, seg_centers[1,], labels = lab_male,   col = "white", font = 2, cex = 0.85)
text(bp, seg_centers[2,], labels = lab_female, col = "white", font = 2, cex = 0.85)

# ✅ (C) 제목 아래 중앙에 p-value(들) 표시
mtext(p_line, side = 3, line = 1.3, adj = 0.5, cex = 1.05, font = 2)

legend("topright", inset = c(-0.2, 0), bty = "n",
       legend = names(cols), fill = cols, cex = 1)

par(op); dev.off()

cat("✅ 저장됨: ", outfile, "\n")
cat(sprintf("📊 Fisher p-values — CCD: %.5g | Normal: %.5g\n", p_ccd, p_normal))

head(merged_total2022_df)
dim(merged_total2022_df)

DogOverview <- read.csv("/home/yhSong/data/02_ByProject/DogAgingProject/DogOverview/v2022/DAP_2022_DogOverview_v1.0.csv")
head(DogOverview)

head(DogOverview$Breed)

# 두 데이터셋을 data.table 형태로 변환
setDT(merged_total2022_df)
setDT(DogOverview)

# 병합 실행: IID ↔ dog_id 기준 (유전형 데이터 기준으로 유지)
merged_full2022_df <- merge(
  merged_total2022_df,
  DogOverview,
  by.x = "IID",      # merged_total2022_df의 키
  by.y = "dog_id",   # DogOverview의 키
  all.x = TRUE       # 유전형 데이터 기준으로 병합
)

# 병합 결과 확인
cat("✅ 병합 완료!\n")
cat("행 수:", nrow(merged_full2022_df), "\n") # 3542
cat("열 수:", ncol(merged_full2022_df), "\n") # 46

# 주요 컬럼 미리보기
head(merged_full2022_df[, .(IID, Breed, CCD_Status)])

table(merged_full2022_df$Breed)
table(merged_full2022_df$CCD_Status)

merged_full2022_select <- merged_full2022_df[,c("IID", "FID", "PAT", "MAT", "SEX", "PHENOTYPE", "X5.25529184.C.T_C", "X27.44543118.A.T_A", 
    "X24.12787767.A.G_G", "X24.12787779.T.C_C", "X23.32288824.G.A_G", "CCD_Status", "Estimated_Age_Years_at_HLES", "Breed_Status", "Breed", "Sex_Class_at_HLES",
    "Weight_Class_5KGBin_at_HLES", "Weight_Class_10KGBin_at_HLES", "Breed_Class", "Breed_Size_Class_at_HLES", "Age_Years_Class_at_HLES", "LifeStage_Class_at_HLES")]
head(merged_full2022_select)
dim(merged_full2022_select) # 3542 22

table(merged_full2022_select$CCD_Status)

breed_sorted <- sort(table(merged_full2022_select$Breed), decreasing = TRUE)
breed2022_ordered <- names(breed_sorted)

# 데이터 준비
dt <- as.data.table(merged_full2022_select)
dt[, CCD_Status := ifelse(CCD_Status == "CCD", "CCD", "Normal")]

snp_list <- c("X5.25529184.C.T_C",
              "X27.44543118.A.T_A",
              "X24.12787767.A.G_G",
              "X24.12787779.T.C_C",
              "X23.32288824.G.A_G")

# 품종 순서
breed_sorted <- sort(table(dt$Breed), decreasing = TRUE)
breed_ordered <- names(breed_sorted)

results_breed_snp <- list()

for (breed in breed_ordered) {
  
  sub_breed <- dt[Breed == breed]
  n_ccd    <- sum(sub_breed$CCD_Status == "CCD")
  n_normal <- sum(sub_breed$CCD_Status == "Normal")
  
  # CCD/Normal 숫자 체크
  if (n_ccd < 2 | n_normal < 2) {
    reason <- c()
    if (n_ccd < 2)    reason <- c(reason, paste0("CCD < 2 (", n_ccd, ")"))
    if (n_normal < 2) reason <- c(reason, paste0("Normal < 2 (", n_normal, ")"))
    
    cat(paste0(breed, ": Skip because ", paste(reason, collapse = " & "), "\n"))
    next
  }
  
  cat("\n===============================\n")
  cat("Breed:", breed, "\n")
  cat(" CCD:", n_ccd, " / Normal:", n_normal, "\n")
  cat("===============================\n")
  
  results_breed_snp[[breed]] <- list()
  
  for (snp in snp_list) {
    
    geno <- sub_breed[[snp]]
    
    temp <- data.table(
      CCD = sub_breed$CCD_Status,
      G   = geno
    )[!is.na(G)]
    
    # 유전자형 variation 체크
    if (length(unique(temp$G)) < 2) {
      cat(paste0("  - SNP ", snp, ": Skip because genotype variation < 2\n"))
      next
    }
    
    # 2×3 테이블
    tab <- table(temp$CCD, temp$G)
    
    # ✅ 여기서 행/열 개수 체크 추가
    if (nrow(tab) < 2 | ncol(tab) < 2) {
      reason <- c()
      if (nrow(tab) < 2) {
        reason <- c(reason, paste0("only one CCD status level after NA filtering (",
                                   paste(rownames(tab), collapse = ", "), ")"))
      }
      if (ncol(tab) < 2) {
        reason <- c(reason, paste0("only one genotype level after NA filtering (",
                                   paste(colnames(tab), collapse = ", "), ")"))
      }
      cat(paste0("  - SNP ", snp, ": Skip because ", paste(reason, collapse = " & "), "\n"))
      next
    }
    
    cat("\n--- SNP:", snp, " ---\n")
    print(tab)
    
    # Fisher exact test (이제는 최소 2×2 보장)
    ft <- fisher.test(tab)
    
    cat("  Fisher p:", ft$p.value, "\n")
    
    results_breed_snp[[breed]][[snp]] <- list(
      table = tab,
      p = ft$p.value,
      test = ft
    )
  }
}


sig_list <- list()

for (breed in names(results_breed_snp)) {
  for (snp in names(results_breed_snp[[breed]])) {
    
    pval <- results_breed_snp[[breed]][[snp]]$p
    
    if (!is.null(pval) && !is.na(pval) && pval < 0.05) {
      sig_list[[paste0(breed, "_", snp)]] <- pval
    }
  }
}

sig_list


table(merged_full2022_df$Breed_Size_Class_at_HLES)
x <- merged_full2022_df$Breed_Size_Class_at_HLES

# 새로운 size class 생성
merged_full2022_df[, Size_Class :=
  fifelse(grepl("^Giant", x), "Giant",
  fifelse(grepl("^Large", x), "Large",
  fifelse(grepl("^Medium", x), "Medium",
  fifelse(grepl("^Standard", x), "Standard",
  fifelse(grepl("^Toy and Small", x), "Toy and Small", NA_character_)))))
]

# 결과 확인
table(merged_full2022_df$Size_Class)

table(merged_full2022_df$CCD_Status)

snp_list <- c("X5.25529184.C.T_C",
              "X27.44543118.A.T_A",
              "X24.12787767.A.G_G",
              "X24.12787779.T.C_C",
              "X23.32288824.G.A_G")

library(data.table)

dt <- as.data.table(merged_full2022_df)

# 결과 저장
results_size_snp <- list()

# 분석 대상 Size_Class
size_classes <- unique(dt$Size_Class)

for (sz in size_classes) {
  
  sub_sz <- dt[Size_Class == sz]

  # CCD/Normal 개수 부족 시 skip
  if (sum(sub_sz$CCD_Status == "CCD") < 2 |
      sum(sub_sz$CCD_Status == "Control") < 2) {

    cat("\n[", sz, "] → Skip: CCD 또는 Normal 표본 수 부족\n")
    results_size_snp[[sz]] <- "Skip: sample size insufficient"
    next
  }
  
  cat("\n===============================\n")
  cat("Size Class:", sz, "\n")
  cat("===============================\n")
  
  results_size_snp[[sz]] <- list()
  
  for (snp in snp_list) {
    
    geno <- sub_sz[[snp]]

    temp <- data.table(
      CCD = sub_sz$CCD_Status,
      G   = geno
    )[!is.na(G)]

    # 유효 유전자형이 2종류 미만이면 skip
    if (length(unique(temp$G)) < 2) {
      cat("[", sz, "] SNP:", snp, " → Skip (genotype variation 부족)\n")
      results_size_snp[[sz]][[snp]] <- "Skip: insufficient genotype variation"
      next
    }

    # 2×3 또는 2×2 교차표
    tab <- table(temp$CCD, temp$G)

    cat("\n--- SNP:", snp, " ---\n")
    print(tab)

    # Fisher exact test
    ft <- fisher.test(tab)

    cat("Fisher p:", ft$p.value, "\n")

    results_size_snp[[sz]][[snp]] <- list(
      table = tab,
      p = ft$p.value,
      test = ft
    )
  }
}

sig_results <- lapply(results_size_snp, function(x) {
  if (is.list(x)) {
    Filter(function(y) is.list(y) && !is.null(y$p) && y$p < 0.05, x)
  }
})

sig_results
