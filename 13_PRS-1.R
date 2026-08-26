library(data.table)
library(dplyr)
library(groupdata2)

getwd()
setwd("J:/01_lab_lee/03_GWAS_Pure/14_Regression_local")
set.seed(2025)

# plink_001_1e6_10 <- fread("/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear")
plink_001_1e6_10 <- fread("05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear")
head(plink_001_1e6_10)

############################################################ 1 Data processing - beta #####################################################

# 1. Beta
# 5:25529184:C:T
snp5_25529184_C_T <- plink_001_1e6_10[plink_001_1e6_10$ID == "5:25529184:C:T"]
snp5_25529184_C_T

# 24:12787767:A:G
snp24_12787767_A_G <- plink_001_1e6_10[plink_001_1e6_10$ID == "24:12787767:A:G"]
snp24_12787767_A_G

# 24:12787779:T:C
snp24_12787779_T_C <- plink_001_1e6_10[plink_001_1e6_10$ID == "24:12787779:T:C"]
snp24_12787779_T_C

# 26:13787793:A:G
snp26_13787793_A_G <- plink_001_1e6_10[plink_001_1e6_10$ID == "26:13787793:A:G"]
snp26_13787793_A_G

# 27:44543118:A:T
snp27_44543118_A_T <- plink_001_1e6_10[plink_001_1e6_10$ID == "27:44543118:A:T"]
snp27_44543118_A_T

############################################################ 1 Data processing - genotype #####################################################

# 2. Genotype
# 5:25529184:C:T
geno5_25529184_C_T <- fread("12_Geno_Beta/snp_5_25529184_C_T_genotype.raw")
head(geno5_25529184_C_T)

# 24:12787767:A:G
geno24_12787767_A_G <- fread("12_Geno_Beta/snp_24_12787767_A_G_genotype.raw")
head(geno24_12787767_A_G)

# 24:12787779:T:C
geno24_12787779_T_C <- fread("12_Geno_Beta/snp_24_12787779_T_C_genotype.raw")
head(geno24_12787779_T_C)

# 26:13787793:A:G
geno26_13787793_A_G <- fread("12_Geno_Beta/snp_26_13787793_A_G_genotype.raw")
head(geno26_13787793_A_G)

# 27:44543118:A:T
geno27_44543118_A_T <- fread("12_Geno_Beta/snp_27_44543118_A_T_genotype.raw")
head(geno27_44543118_A_T)

############################################################ 1 Data processing - covariate #####################################################

# 3. PCs,Age,Sex,Breed,Neuter Status, Activity Level, CSLB score
covariates <- fread("04_cov_001_1e6_10.txt")
head(covariates)

cslb <- fread("../04_pheno001.txt")
head(cslb)

############################################################ 1 Data processing - merge data (full data) #####################################################

setnames(cslb, old = "phe1", new = "cslb_score")

merged_df <- merge(
  covariates,
  cslb[, .(FID, IID, cslb_score)],   # FID, IID, cslb_score only
  by = c("FID", "IID"),
  all.x = TRUE
)

head(merged_df)
dim(merged_df)

beta_df <- rbind(
  snp5_25529184_C_T,
  snp24_12787767_A_G,
  snp24_12787779_T_C,
  snp26_13787793_A_G,
  snp27_44543118_A_T
)[, .(ID, A1, BETA)]

beta_df

beta_df[ID %in% c("24:12787767:A:G", "24:12787779:T:C"), .(ID, A1, BETA)]

extract_geno <- function(dt) {
  snp_col <- grep("^\\d", colnames(dt), value = TRUE)  
  dt[, .(FID, IID, snp = get(snp_col))]
}

g5  <- extract_geno(geno5_25529184_C_T)
g24_1 <- extract_geno(geno24_12787767_A_G)
g24_2 <- extract_geno(geno24_12787779_T_C)
g26 <- extract_geno(geno26_13787793_A_G)
g27 <- extract_geno(geno27_44543118_A_T)

# 컬럼명 통일
setnames(g5,  "snp", "snp5_25529184")
setnames(g24_1, "snp", "snp24_12787767")
setnames(g24_2, "snp", "snp24_12787779")
setnames(g26, "snp", "snp26_13787793")
setnames(g27, "snp", "snp27_44543118")

# 3) 하나의 genotype DF로 merge
geno_df <- Reduce(
  function(x,y) merge(x, y, by = c("FID","IID"), all = TRUE),
  list(g5, g24_1, g24_2, g26, g27)
)

head(geno_df)
dim(geno_df)

full_df <- merge(
  merged_df,    # cov + cslb_score
  geno_df,      # SNP genotype
  by = c("FID", "IID"),
  all.x = TRUE
)

head(full_df)
dim(full_df)


# ------------------------------------------------------------------------------
# full_df 생성 직후 - Golden Retriever 컬럼 존재 여부 확인
# ------------------------------------------------------------------------------

# 1. Golden Retriever 관련 컬럼이 full_df에 있는지 확인
grep("Golden", names(full_df), value = TRUE, ignore.case = TRUE)

# 2. breed_cols를 만드는 원래 탐지 로직으로 확인 (아직 exclude_vars 적용 전 상태로)
all_breed_candidate_cols <- names(full_df)[sapply(full_df, function(x) all(x %in% c(0, 1), na.rm = TRUE))]
grep("Golden", all_breed_candidate_cols, value = TRUE, ignore.case = TRUE)

# 3. full_df 전체에서 모든 breed dummy가 0인 개체(=Unknown이 될 개체) 미리 확인
breed_dummy_cols_temp <- setdiff(all_breed_candidate_cols, c("Sex_cov", "Intact_cov"))  # 필요시 조정
row_sums <- rowSums(full_df[, ..breed_dummy_cols_temp], na.rm = TRUE)
sum(row_sums == 0)
unknown_ids_full <- full_df$IID[row_sums == 0]
length(unknown_ids_full)

# 4. 원본 품종 라벨 파일과 대조
purebred_total <- read.csv("../00_Purebred_total.csv", stringsAsFactors = FALSE)
names(purebred_total)  # 실제 컬럼명 먼저 확인

matched_breeds <- purebred_total[purebred_total$IID %in% unknown_ids_full, c("IID", "Breed")]
table(matched_breeds$Breed)

# ------------------------------------------------------------------------------
# 원본 Breed 라벨을 full_df에 직접 병합 (dummy 컬럼 재구성 대신, 진짜 라벨 사용)
# ------------------------------------------------------------------------------

purebred_total <- read.csv("../00_Purebred_total.csv", stringsAsFactors = FALSE)

# 타입 일치 확인 후 병합
full_df <- as.data.frame(full_df)
full_df <- merge(
  full_df,
  purebred_total[, c("IID", "Breed")],
  by = "IID",
  all.x = TRUE
)

# 병합 후 결측 확인 (원본에 없는 IID가 있는지)
sum(is.na(full_df$Breed))

# 컬럼명 공백/특수문자를 breed_cols 표기 방식과 통일 (예: "Golden Retriever" → "Golden_Retriever")
full_df$Breed <- gsub(" ", "_", full_df$Breed)

table(full_df$Breed)[order(-table(full_df$Breed))][1:10]  # 상위 10개 품종 확인
# ------------------------------------------------------------------------------
# STEP 1. 기본 설정값 정의
# ------------------------------------------------------------------------------


snp_cols <- c(
  "snp5_25529184",
  "snp24_12787767",
  "snp24_12787779",
  "snp26_13787793",
  "snp27_44543118"
)

table(full_df$snp5_25529184, useNA = "ifany")
#    0    1    2 <NA> 
#    1  137 2730  659

table(full_df$snp24_12787767, useNA = "ifany")
#    0    1    2 <NA> 
#   34  333 2156 1004 

table(full_df$snp24_12787779, useNA = "ifany")
#    0    1    2 <NA> 
#   32  315 2304  876 

table(full_df$snp26_13787793, useNA = "ifany")
#  0    1    2 <NA> 
#  1   35 1685 1806 

table(full_df$snp27_44543118, useNA = "ifany")
#    0    1    2 <NA> 
#    6  104 2762  655 

# GWAS 원 모델과 동일한 형태로 SNP effect size를 재추정할 때 같이 넣을 조정 공변량
adjust_covariates <- c("Age_cov", "Sex_cov", "Intact_cov", paste0("PC", 1:8), "pa_score_cov")

# 컬럼명에 특수문자가 있는 품종명 정리 (기존 스크립트와 동일)
fix_col <- function(df) {
  names(df)[names(df) == "Poodle_(Toy)"] <- "Poodle_Toy"
  return(df)
}

full_df <- fix_col(as.data.frame(full_df))

head(full_df)
full_df2 <- full_df

dim(full_df)
dim(full_df2)

# breed dummy 컬럼 자동 탐지
breed_cols <- names(full_df)[sapply(full_df, function(x) all(x %in% c(0, 1), na.rm = TRUE))]
exclude_vars <- c(
  "Sex_cov", "Intact_cov", "Age_cov", "pa_score_cov",
  paste0("PC", 1:8), snp_cols, "PRS", "cslb_score"
)
breed_cols <- setdiff(breed_cols, exclude_vars)

cat("Number of breed covariates detected:", length(breed_cols), "\n") # 102

# ------------------------------------------------------------------------------
# STEP 2. 분석 대상 데이터 준비 (CSLB 점수 결측만 제외, genotype 결측은 fold 내에서 처리)
# ------------------------------------------------------------------------------

analysis_df <- full_df[!is.na(full_df$cslb_score), ]
analysis_df$IID <- factor(analysis_df$IID)

cat("Analysis dataset dimensions:", dim(analysis_df)[1], "rows,", dim(analysis_df)[2], "columns\n")

dim(analysis_df)

# CSLB 점수 분포 기반 층화를 위한 구간 나누기
qs <- unique(quantile(analysis_df$cslb_score, probs = seq(0, 1, 0.2), na.rm = TRUE))
analysis_df$cslb_strata <- cut(analysis_df$cslb_score, breaks = qs, include.lowest = TRUE)

K       <- 10
REPEATS <- 5

# ------------------------------------------------------------------------------
# STEP 3. Fold 배정을 한 번만 생성해서 모든 모델이 동일하게 공유
# ------------------------------------------------------------------------------
# fold_assignments[[r]] = 해당 repeat의 fold_id 벡터 (analysis_df와 같은 순서/길이)

fold_assignments <- list()
for (r in 1:REPEATS) {
  set.seed(2025 + r)
  folds <- fold(analysis_df, k = K, cat_col = "cslb_strata", id_col = "IID")
  fold_assignments[[r]] <- folds$.folds
}

# ------------------------------------------------------------------------------
# STEP 4. 공통 헬퍼 함수: 하나의 fold에 대해 결측치 대체 + beta 재추정 + PRS 생성
# ------------------------------------------------------------------------------

prepare_fold_data <- function(train_fold, test_fold) {
  
  # (a) 결측 genotype 대체: train fold의 AF만 사용
  af_train <- sapply(snp_cols, function(s) mean(train_fold[[s]], na.rm = TRUE) / 2)
  for (s in snp_cols) {
    train_fold[[s]][is.na(train_fold[[s]])] <- 2 * af_train[s]
    test_fold[[s]][is.na(test_fold[[s]])]   <- 2 * af_train[s]
  }
  
  # (b) SNP effect size(beta) 재추정: train fold에서만
  beta_refit <- sapply(snp_cols, function(s) {
    f <- as.formula(paste("cslb_score ~", s, "+", paste(adjust_covariates, collapse = " + ")))
    fit <- lm(f, data = train_fold)
    coef(fit)[s]
  })
  names(beta_refit) <- snp_cols
  
  # (c) Formula 1: genome-wide PRS
  geno_mat_train <- as.matrix(train_fold[, snp_cols])
  geno_mat_test  <- as.matrix(test_fold[, snp_cols])
  train_fold$PRS <- as.numeric(geno_mat_train %*% beta_refit)
  test_fold$PRS  <- as.numeric(geno_mat_test  %*% beta_refit)
  
  # (d) Formula 2: chromosome별 개별 PRS
  train_fold$PRS_5    <- train_fold$snp5_25529184  * beta_refit["snp5_25529184"]
  train_fold$PRS_24_1 <- train_fold$snp24_12787767 * beta_refit["snp24_12787767"]
  train_fold$PRS_24_2 <- train_fold$snp24_12787779 * beta_refit["snp24_12787779"]
  train_fold$PRS_26   <- train_fold$snp26_13787793 * beta_refit["snp26_13787793"]
  train_fold$PRS_27   <- train_fold$snp27_44543118 * beta_refit["snp27_44543118"]
  
  test_fold$PRS_5    <- test_fold$snp5_25529184  * beta_refit["snp5_25529184"]
  test_fold$PRS_24_1 <- test_fold$snp24_12787767 * beta_refit["snp24_12787767"]
  test_fold$PRS_24_2 <- test_fold$snp24_12787779 * beta_refit["snp24_12787779"]
  test_fold$PRS_26   <- test_fold$snp26_13787793 * beta_refit["snp26_13787793"]
  test_fold$PRS_27   <- test_fold$snp27_44543118 * beta_refit["snp27_44543118"]
  
  list(train = train_fold, test = test_fold)
}

prepare_fold_data <- function(train_fold, test_fold) {
 
  # (a) 결측 genotype 대체: train fold의 AF만 사용
  af_train <- sapply(snp_cols, function(s) mean(train_fold[[s]], na.rm = TRUE) / 2)
  for (s in snp_cols) {
    train_fold[[s]][is.na(train_fold[[s]])] <- 2 * af_train[s]
    test_fold[[s]][is.na(test_fold[[s]])]   <- 2 * af_train[s]
  }
 
  # (b) SNP effect size(beta) 재추정: train fold에서만
  beta_refit <- sapply(snp_cols, function(s) {
    f <- as.formula(paste("cslb_score ~", s, "+", paste(adjust_covariates, collapse = " + ")))
    fit <- lm(f, data = train_fold)
    coef(fit)[s]
  })
  names(beta_refit) <- snp_cols
 
  # (c) Formula 1: genome-wide PRS
  geno_mat_train <- as.matrix(train_fold[, snp_cols])
  geno_mat_test  <- as.matrix(test_fold[, snp_cols])
  train_fold$PRS <- as.numeric(geno_mat_train %*% beta_refit)
  test_fold$PRS  <- as.numeric(geno_mat_test  %*% beta_refit)
 
  # (d) Formula 2: chromosome별 개별 PRS
  train_fold$PRS_5    <- train_fold$snp5_25529184  * beta_refit["snp5_25529184"]
  train_fold$PRS_24_1 <- train_fold$snp24_12787767 * beta_refit["snp24_12787767"]
  train_fold$PRS_24_2 <- train_fold$snp24_12787779 * beta_refit["snp24_12787779"]
  train_fold$PRS_26   <- train_fold$snp26_13787793 * beta_refit["snp26_13787793"]
  train_fold$PRS_27   <- train_fold$snp27_44543118 * beta_refit["snp27_44543118"]
 
  test_fold$PRS_5    <- test_fold$snp5_25529184  * beta_refit["snp5_25529184"]
  test_fold$PRS_24_1 <- test_fold$snp24_12787767 * beta_refit["snp24_12787767"]
  test_fold$PRS_24_2 <- test_fold$snp24_12787779 * beta_refit["snp24_12787779"]
  test_fold$PRS_26   <- test_fold$snp26_13787793 * beta_refit["snp26_13787793"]
  test_fold$PRS_27   <- test_fold$snp27_44543118 * beta_refit["snp27_44543118"]
 
  list(train = train_fold, test = test_fold)
}
 
compute_r2_rmse <- function(actual, predicted) {
  r2   <- cor(actual, predicted, use = "complete.obs")^2
  rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
  list(R2 = r2, RMSE = rmse)
}

# ------------------------------------------------------------------------------
# STEP 5-A. Baseline 모델 (PRS 없음) - 독립된 for문
# ------------------------------------------------------------------------------

results_baseline <- list()
idx <- 1

for (r in 1:REPEATS) {
  analysis_df$fold_id <- fold_assignments[[r]]
  
  for (k in 1:K) {
    train_fold <- analysis_df[analysis_df$fold_id != k, ]
    test_fold  <- analysis_df[analysis_df$fold_id == k, ]
    
    prepped <- prepare_fold_data(train_fold, test_fold)
    train_fold <- prepped$train
    test_fold  <- prepped$test
    
    f_baseline <- as.formula(paste("cslb_score ~", paste(c(adjust_covariates, breed_cols), collapse = " + ")))
    model_baseline <- lm(f_baseline, data = train_fold)
    pred <- predict(model_baseline, newdata = test_fold)
    metrics <- compute_r2_rmse(test_fold$cslb_score, pred)
    
    results_baseline[[idx]] <- data.frame(
      repeat_num = r, fold_num = k,
      R2_baseline = metrics$R2, RMSE_baseline = metrics$RMSE
    )
    idx <- idx + 1
  }
}
results_baseline_df <- do.call(rbind, results_baseline)
cat("Baseline model: done.\n")

# ------------------------------------------------------------------------------
# STEP 5-B. Model 1 - Formula 1 (genome-wide PRS only) - 독립된 for문
# ------------------------------------------------------------------------------

results_model1_f1 <- list()
idx <- 1

for (r in 1:REPEATS) {
  analysis_df$fold_id <- fold_assignments[[r]]
  
  for (k in 1:K) {
    train_fold <- analysis_df[analysis_df$fold_id != k, ]
    test_fold  <- analysis_df[analysis_df$fold_id == k, ]
    
    prepped <- prepare_fold_data(train_fold, test_fold)
    train_fold <- prepped$train
    test_fold  <- prepped$test
    
    model1_f1 <- lm(cslb_score ~ PRS, data = train_fold)
    pred <- predict(model1_f1, newdata = test_fold)
    metrics <- compute_r2_rmse(test_fold$cslb_score, pred)
    
    results_model1_f1[[idx]] <- data.frame(
      repeat_num = r, fold_num = k,
      R2_model1_formula1 = metrics$R2, RMSE_model1_formula1 = metrics$RMSE
    )
    idx <- idx + 1
  }
}
results_model1_f1_df <- do.call(rbind, results_model1_f1)
cat("Model 1 - Formula 1: done.\n")

# ------------------------------------------------------------------------------
# STEP 5-C. Model 1 - Formula 2 (chromosome-specific PRS only) - 독립된 for문
# ------------------------------------------------------------------------------

results_model1_f2 <- list()
idx <- 1

for (r in 1:REPEATS) {
  analysis_df$fold_id <- fold_assignments[[r]]
  
  for (k in 1:K) {
    train_fold <- analysis_df[analysis_df$fold_id != k, ]
    test_fold  <- analysis_df[analysis_df$fold_id == k, ]
    
    prepped <- prepare_fold_data(train_fold, test_fold)
    train_fold <- prepped$train
    test_fold  <- prepped$test
    
    model1_f2 <- lm(cslb_score ~ PRS_5 + PRS_24_1 + PRS_24_2 + PRS_26 + PRS_27, data = train_fold)
    pred <- predict(model1_f2, newdata = test_fold)
    metrics <- compute_r2_rmse(test_fold$cslb_score, pred)
    
    results_model1_f2[[idx]] <- data.frame(
      repeat_num = r, fold_num = k,
      R2_model1_formula2 = metrics$R2, RMSE_model1_formula2 = metrics$RMSE
    )
    idx <- idx + 1
  }
}
results_model1_f2_df <- do.call(rbind, results_model1_f2)
cat("Model 1 - Formula 2: done.\n")

# ------------------------------------------------------------------------------
# STEP 5-D. Model 2 - Formula 1 (genome-wide PRS + covariates) - 독립된 for문
# ------------------------------------------------------------------------------

results_model2_f1 <- list()
idx <- 1

for (r in 1:REPEATS) {
  analysis_df$fold_id <- fold_assignments[[r]]
  
  for (k in 1:K) {
    train_fold <- analysis_df[analysis_df$fold_id != k, ]
    test_fold  <- analysis_df[analysis_df$fold_id == k, ]
    
    prepped <- prepare_fold_data(train_fold, test_fold)
    train_fold <- prepped$train
    test_fold  <- prepped$test
    
    f_model2_f1 <- as.formula(paste("cslb_score ~", paste(c("PRS", adjust_covariates, breed_cols), collapse = " + ")))
    model2_f1 <- lm(f_model2_f1, data = train_fold)
    pred <- predict(model2_f1, newdata = test_fold)
    metrics <- compute_r2_rmse(test_fold$cslb_score, pred)
    
    results_model2_f1[[idx]] <- data.frame(
      repeat_num = r, fold_num = k,
      R2_model2_formula1 = metrics$R2, RMSE_model2_formula1 = metrics$RMSE
    )
    idx <- idx + 1
  }
}
results_model2_f1_df <- do.call(rbind, results_model2_f1)
cat("Model 2 - Formula 1: done.\n")

# ------------------------------------------------------------------------------
# STEP 5-E. Model 2 - Formula 2 (chromosome-specific PRS + covariates) - 독립된 for문
# ------------------------------------------------------------------------------

results_model2_f2 <- list()
idx <- 1

for (r in 1:REPEATS) {
  analysis_df$fold_id <- fold_assignments[[r]]
  
  for (k in 1:K) {
    train_fold <- analysis_df[analysis_df$fold_id != k, ]
    test_fold  <- analysis_df[analysis_df$fold_id == k, ]
    
    prepped <- prepare_fold_data(train_fold, test_fold)
    train_fold <- prepped$train
    test_fold  <- prepped$test
    
    f_model2_f2 <- as.formula(paste(
      "cslb_score ~",
      paste(c("PRS_5", "PRS_24_1", "PRS_24_2", "PRS_26", "PRS_27", adjust_covariates, breed_cols), collapse = " + ")
    ))
    model2_f2 <- lm(f_model2_f2, data = train_fold)
    pred <- predict(model2_f2, newdata = test_fold)
    metrics <- compute_r2_rmse(test_fold$cslb_score, pred)
    
    results_model2_f2[[idx]] <- data.frame(
      repeat_num = r, fold_num = k,
      R2_model2_formula2 = metrics$R2, RMSE_model2_formula2 = metrics$RMSE
    )
    idx <- idx + 1
  }
}
results_model2_f2_df <- do.call(rbind, results_model2_f2)
cat("Model 2 - Formula 2: done.\n")

# ==============================================================================
# Model 3 ~ 5 (mixed-effects models)
# ==============================================================================

library(lme4)
library(lmerTest)
library(mclust)
library(performance)

# ------------------------------------------------------------------------------
# 공통 헬퍼: breed dummy 컬럼들로부터 Breed factor 생성 (train 기준으로 level 고정)
# ------------------------------------------------------------------------------

# make_breed_factor <- function(df, breed_dummy_cols) {
#   apply(df[, breed_dummy_cols], 1, function(x) {
#     col <- breed_dummy_cols[which(x == 1)]
#     if (length(col) == 0) return("Unknown")
#     return(col)
#   })
# }

# 공통 헬퍼: mixed model 결과에서 out-of-sample R²/RMSE + marginal/conditional R² 추출
evaluate_mixed_model <- function(model, test_fold, actual_col = "cslb_score") {
  pred <- predict(model, newdata = test_fold, allow.new.levels = TRUE)
  oos  <- compute_r2_rmse(test_fold[[actual_col]], pred)
  
  r2_decomp <- tryCatch(performance::r2(model), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
  
  list(
    R2_out_of_sample   = oos$R2,
    RMSE_out_of_sample = oos$RMSE,
    R2_marginal        = as.numeric(r2_decomp$R2_marginal),
    R2_conditional      = as.numeric(r2_decomp$R2_conditional)
  )
}

# ------------------------------------------------------------------------------
# STEP 5-F. Model 3: PRS + covariates + (1|Breed)  -- Formula 1 & 2
# ------------------------------------------------------------------------------

results_model3_f1 <- list()
results_model3_f2 <- list()
idx1 <- 1; idx2 <- 1

for (r in 1:REPEATS) {
  analysis_df$fold_id <- fold_assignments[[r]]
  
  for (k in 1:K) {
    train_fold <- analysis_df[analysis_df$fold_id != k, ]
    test_fold  <- analysis_df[analysis_df$fold_id == k, ]
    
    prepped <- prepare_fold_data(train_fold, test_fold)
    train_fold <- prepped$train
    test_fold  <- prepped$test
    
    # Breed factor 생성 (train 기준 level 고정, test는 train level에 맞춤)
    train_fold$Breed <- factor(train_fold$Breed)
    test_fold$Breed  <- factor(test_fold$Breed, levels = levels(train_fold$Breed))
    
    # --- Formula 1: genome-wide PRS ---
    f3_f1 <- as.formula(paste("cslb_score ~ PRS +", paste(adjust_covariates, collapse = " + "), "+ (1 | Breed)"))
    model3_f1 <- lmer(f3_f1, data = train_fold, REML = TRUE)
    m3f1 <- evaluate_mixed_model(model3_f1, test_fold)
    
    results_model3_f1[[idx1]] <- data.frame(
      repeat_num = r, fold_num = k,
      is_singular = lme4::isSingular(model3_f1),
      R2_model3_formula1_oos = m3f1$R2_out_of_sample,
      RMSE_model3_formula1   = m3f1$RMSE_out_of_sample,
      R2_model3_formula1_marginal    = m3f1$R2_marginal,
      R2_model3_formula1_conditional = m3f1$R2_conditional
    )
    idx1 <- idx1 + 1
    
    # --- Formula 2: chromosome-specific PRS ---
    f3_f2 <- as.formula(paste(
      "cslb_score ~ PRS_5 + PRS_24_1 + PRS_24_2 + PRS_26 + PRS_27 +",
      paste(adjust_covariates, collapse = " + "), "+ (1 | Breed)"
    ))
    model3_f2 <- lmer(f3_f2, data = train_fold, REML = TRUE)
    m3f2 <- evaluate_mixed_model(model3_f2, test_fold)
    
    results_model3_f2[[idx2]] <- data.frame(
      repeat_num = r, fold_num = k,
      is_singular = lme4::isSingular(model3_f2),
      R2_model3_formula2_oos = m3f2$R2_out_of_sample,
      RMSE_model3_formula2   = m3f2$RMSE_out_of_sample,
      R2_model3_formula2_marginal    = m3f2$R2_marginal,
      R2_model3_formula2_conditional = m3f2$R2_conditional
    )
    idx2 <- idx2 + 1
  }
}
results_model3_f1_df <- do.call(rbind, results_model3_f1)
results_model3_f2_df <- do.call(rbind, results_model3_f2)
cat("Model 3 (Breed random effect): done.\n")

dim(results_model3_f1_df)  # 50행(REPEATS×K)이 다 채워졌는지 확인
sum(is.na(results_model3_f1_df$R2_model3_formula1_conditional))  # NA가 몇 개인지
sum(is.na(results_model3_f1_df$R2_model3_formula1_marginal))

# ------------------------------------------------------------------------------
# STEP 5-G. Model 4: PRS + covariates + (1|Breed) + (1|PC cluster) -- Formula 1 & 2
# ------------------------------------------------------------------------------

pc_cols <- paste0("PC", 1:8)
results_model4_f1 <- list()   # 빈 리스트로 시작
results_model4_f2 <- list()
idx1 <- 1    
idx2 <- 1

for (r in 1:REPEATS) {
  analysis_df$fold_id <- fold_assignments[[r]]
  
  for (k in 1:K) {
    train_fold <- analysis_df[analysis_df$fold_id != k, ]
    test_fold  <- analysis_df[analysis_df$fold_id == k, ]
    
    prepped <- prepare_fold_data(train_fold, test_fold)
    train_fold <- prepped$train
    test_fold  <- prepped$test
    
    train_fold$Breed <- factor(train_fold$Breed)
    test_fold$Breed  <- factor(test_fold$Breed, levels = levels(train_fold$Breed))
    
    # PC cluster: mclust을 train fold PC1-8로만 적합 (test는 predict로 배정)
    # G(클러스터 개수)는 train fold의 BIC 최적값을 매 fold마다 자동 선택 (test 성능을 보고 고르지 않음)
    pc_train_mat <- as.matrix(train_fold[, pc_cols])
    pc_test_mat  <- as.matrix(test_fold[, pc_cols])
    
    set.seed(2025)
    bic_result <- mclustBIC(pc_train_mat, verbose = FALSE)
    best_idx <- which(bic_result == max(bic_result, na.rm = TRUE), arr.ind = TRUE)
    best_G <- as.numeric(rownames(bic_result)[best_idx[1, "row"]])
    if (is.na(best_G) || best_G < 2) best_G <- 2
    
    pc_mclust <- Mclust(pc_train_mat, G = best_G, verbose = FALSE)
    
    train_fold$PC_cluster <- factor(pc_mclust$classification)
    test_cluster <- predict(pc_mclust, newdata = pc_test_mat)$classification
    test_fold$PC_cluster <- factor(test_cluster, levels = levels(train_fold$PC_cluster))
    
    # --- Formula 1 ---
    f4_f1 <- as.formula(paste(
      "cslb_score ~ PRS +", paste(adjust_covariates, collapse = " + "),
      "+ (1 | Breed) + (1 | PC_cluster)"
    ))
    model4_f1 <- lmer(f4_f1, data = train_fold, REML = TRUE)
    m4f1 <- evaluate_mixed_model(model4_f1, test_fold)
    
    results_model4_f1[[idx1]] <- data.frame(
      repeat_num = r, fold_num = k,
      is_singular = lme4::isSingular(model4_f1),
      selected_G = best_G,
      R2_model4_formula1_oos = m4f1$R2_out_of_sample,
      RMSE_model4_formula1   = m4f1$RMSE_out_of_sample,
      R2_model4_formula1_marginal    = m4f1$R2_marginal,
      R2_model4_formula1_conditional = m4f1$R2_conditional
    )
    idx1 <- idx1 + 1
    
    # --- Formula 2 ---
    f4_f2 <- as.formula(paste(
      "cslb_score ~ PRS_5 + PRS_24_1 + PRS_24_2 + PRS_26 + PRS_27 +",
      paste(adjust_covariates, collapse = " + "),
      "+ (1 | Breed) + (1 | PC_cluster)"
    ))
    model4_f2 <- lmer(f4_f2, data = train_fold, REML = TRUE)
    m4f2 <- evaluate_mixed_model(model4_f2, test_fold)
    
    results_model4_f2[[idx2]] <- data.frame(
      repeat_num = r, fold_num = k,
      is_singular = lme4::isSingular(model4_f2),
      selected_G = best_G,
      R2_model4_formula2_oos = m4f2$R2_out_of_sample,
      RMSE_model4_formula2   = m4f2$RMSE_out_of_sample,
      R2_model4_formula2_marginal    = m4f2$R2_marginal,
      R2_model4_formula2_conditional = m4f2$R2_conditional
    )
    idx2 <- idx2 + 1
  }
}
results_model4_f1_df <- do.call(rbind, results_model4_f1)
results_model4_f2_df <- do.call(rbind, results_model4_f2)
cat("Model 4 (Breed + PC cluster random effects): done.\n")

dim(results_model4_f1_df)  # 50행(REPEATS×K)이 다 채워졌는지 확인
sum(is.na(results_model4_f1_df$R2_model4_formula1_conditional))  # NA가 몇 개인지
sum(is.na(results_model4_f1_df$R2_model4_formula1_marginal))

table(results_model4_f1_df$selected_G)
table(train_fold$Breed, train_fold$PC_cluster)

getwd()
write.csv(table(train_fold$Breed, train_fold$PC_cluster),
          file = "10_breed_cluster_table_method1model4.csv",
          row.names = TRUE)

# ------------------------------------------------------------------------------
# STEP 5-H. Model 5: PRS + covariates + (1|Breed) + (1|GMM cluster) -- Formula 1 & 2
# ------------------------------------------------------------------------------
# GMM cluster: 원 스크립트의 pseudo-admixture 방식(PC1,PC2 기준 breed centroid와의 거리)을
# train fold 안에서만 계산하고, test fold는 train에서 만든 centroid/모델로 투영합니다.

# # 단일 fold 테스트 (MAX_DIMS를 100으로 늘려서 시간 측정)
# MAX_DIMS_TEST <- 100
# 
# r_test <- 1; k_test <- 1
# df_r <- analysis_df
# df_r$fold_id <- fold_assignments[[r_test]]
# train_fold <- df_r[df_r$fold_id != k_test, ]
# test_fold  <- df_r[df_r$fold_id == k_test, ]
# 
# prepped <- prepare_fold_data(train_fold, test_fold)
# train_fold <- prepped$train
# test_fold  <- prepped$test
# 
# train_fold$Breed <- factor(train_fold$Breed)
# breed_centroids <- aggregate(train_fold[, gbc_pc_cols], by = list(Breed = train_fold$Breed), FUN = mean)
# cent_mat <- as.matrix(breed_centroids[, gbc_pc_cols])
# rownames(cent_mat) <- breed_centroids$Breed
# pc_train_mat <- as.matrix(train_fold[, gbc_pc_cols])
# train_GBC <- calc_pseudo_admix(pc_train_mat, cent_mat)
# gbc_train_trim <- train_GBC[, -ncol(train_GBC), drop = FALSE]
# gbc_pca <- prcomp(gbc_train_trim, scale. = TRUE)
# 
# cum_var <- summary(gbc_pca)$importance[3, ]
# n_dims_needed <- which(cum_var >= 0.85)[1]
# cat("85%를 달성하려면 필요한 차원 수:", n_dims_needed, "\n")  # ← 이게 핵심 확인 포인트
# 
# # 실제 mclust 소요시간 측정
# train_GBC_reduced <- gbc_pca$x[, 1:min(n_dims_needed, MAX_DIMS_TEST), drop = FALSE]
# system.time({
#   bic_result <- mclustBIC(train_GBC_reduced, verbose = FALSE,
#                           modelNames = c("EII", "VII", "EEI", "VVI", "VVV"))
# })


library(foreach)
library(doSNOW)
library(proxy)
library(mclust)
library(lme4)


gbc_pc_cols <- c("PC1", "PC2")
TARGET_VAR  <- 0.7   # 목표 누적 분산 설명 비율
MAX_DIMS    <- 80     # 계산 시간 제어용 상한

calc_pseudo_admix <- function(pc_mat, cent_mat) {
  # 유클리드 거리를 행렬 연산으로 한 번에 계산
  # dist_mat[i, k] = ||pc_mat[i,] - cent_mat[k,]||
  dist_mat <- as.matrix(proxy::dist(pc_mat, cent_mat, method = "Euclidean"))
  
  out <- 1 / (dist_mat + 1e-8)
  out <- out / rowSums(out)
  colnames(out) <- rownames(cent_mat)
  out
}

results_model5_f1 <- list()
results_model5_f2 <- list()
idx1 <- 1; idx2 <- 1

# (r, k) 조합을 미리 그리드로 만들기
grid <- expand.grid(r = 1:REPEATS, k = 1:K)

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeSOCKcluster(n_cores)      
registerDoSNOW(cl)                 
cat("Model 5: using", n_cores, "cores for", nrow(grid), "folds.\n")

pb <- txtProgressBar(max = nrow(grid), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

model5_results <- foreach(
  i = 1:nrow(grid),
  .packages = c("lme4", "mclust", "proxy"),
  .export = c("analysis_df", "fold_assignments", "gbc_pc_cols",
              "TARGET_VAR", "MAX_DIMS",
              "calc_pseudo_admix", "prepare_fold_data", "evaluate_mixed_model",
              "adjust_covariates"),
  .options.snow = opts
) %dopar% {
  
  r <- grid$r[i]
  k <- grid$k[i]
  
  df_r <- analysis_df
  df_r$fold_id <- fold_assignments[[r]]
  
  train_fold <- df_r[df_r$fold_id != k, ]
  test_fold  <- df_r[df_r$fold_id == k, ]
  
  prepped <- prepare_fold_data(train_fold, test_fold)
  train_fold <- prepped$train
  test_fold  <- prepped$test
  
  train_fold$Breed <- factor(train_fold$Breed)
  test_fold$Breed  <- factor(test_fold$Breed, levels = levels(train_fold$Breed))
  
  # Breed별 PC1/PC2 중심점 계산 (train fold 기준)
  breed_centroids <- aggregate(train_fold[, gbc_pc_cols], by = list(Breed = train_fold$Breed), FUN = mean)
  cent_mat <- as.matrix(breed_centroids[, gbc_pc_cols])
  rownames(cent_mat) <- breed_centroids$Breed
  
  pc_train_mat <- as.matrix(train_fold[, gbc_pc_cols])
  pc_test_mat  <- as.matrix(test_fold[, gbc_pc_cols])
  
  train_GBC <- calc_pseudo_admix(pc_train_mat, cent_mat)
  test_GBC  <- calc_pseudo_admix(pc_test_mat, cent_mat)
  
  # --- 차원 축소: compositional 데이터(행 합=1)의 선형종속 컬럼 하나 제거 후 PCA ---
  gbc_train_trim <- train_GBC[, -ncol(train_GBC), drop = FALSE]
  gbc_test_trim  <- test_GBC[,  -ncol(test_GBC),  drop = FALSE]
  
  gbc_pca <- prcomp(gbc_train_trim, scale. = TRUE)
  
  # ------------------------------------------------------------------------
  # [진단 + 차원 선택] 누적 분산 TARGET_VAR를 처음 넘기는 시점까지 차원 수 자동 결정
  # ------------------------------------------------------------------------
  cum_var <- summary(gbc_pca)$importance[3, ]
  n_dims  <- which(cum_var >= TARGET_VAR)[1]
  if (is.na(n_dims)) n_dims <- length(cum_var)      # 목표치를 아예 못 넘기면 전체 사용
  n_dims  <- min(n_dims, MAX_DIMS, ncol(gbc_pca$x))
  
  var_explained <- cum_var[n_dims]                   # 실제 사용된 차원에서의 누적 분산 설명 비율
  
  train_GBC_reduced <- gbc_pca$x[, 1:n_dims, drop = FALSE]
  test_GBC_reduced  <- predict(gbc_pca, newdata = gbc_test_trim)[, 1:n_dims, drop = FALSE]
  
  # GMM (mclust)을 축소된 pseudo-admixture 공간에 적합, G는 train BIC로 자동 선택
  set.seed(2025)
  bic_result <- mclustBIC(train_GBC_reduced, verbose = FALSE,
                          modelNames = c("EII", "VII", "EEI", "VVI", "VVV"))
  best_idx <- which(bic_result == max(bic_result, na.rm = TRUE), arr.ind = TRUE)
  best_G <- as.numeric(rownames(bic_result)[best_idx[1, "row"]])
  if (is.na(best_G) || best_G < 2) best_G <- 2
  
  gbc_mclust <- Mclust(train_GBC_reduced, G = best_G, verbose = FALSE)
  
  train_fold$GMM_cluster <- factor(gbc_mclust$classification)
  test_cluster <- predict(gbc_mclust, newdata = test_GBC_reduced)$classification
  test_fold$GMM_cluster <- factor(test_cluster, levels = levels(train_fold$GMM_cluster))
  
  # --- Formula 1 ---
  f5_f1 <- as.formula(paste(
    "cslb_score ~ PRS +", paste(adjust_covariates, collapse = " + "),
    "+ (1 | Breed) + (1 | GMM_cluster)"
  ))
  model5_f1 <- lmer(f5_f1, data = train_fold, REML = TRUE,
                    control = lmerControl(optimizer = "bobyqa"))
  m5f1 <- evaluate_mixed_model(model5_f1, test_fold)
  
  row_f1 <- data.frame(
    repeat_num = r, fold_num = k,
    is_singular = lme4::isSingular(model5_f1),
    selected_G = best_G,
    n_dims_used = n_dims,                # ← 추가: 실제 사용된 PCA 차원 수
    pca_var_explained = var_explained,   # ← 그 차원에서의 누적 분산 설명 비율
    R2_model5_formula1_oos = m5f1$R2_out_of_sample,
    RMSE_model5_formula1   = m5f1$RMSE_out_of_sample,
    R2_model5_formula1_marginal    = m5f1$R2_marginal,
    R2_model5_formula1_conditional = m5f1$R2_conditional
  )
  
  # --- Formula 2 ---
  f5_f2 <- as.formula(paste(
    "cslb_score ~ PRS_5 + PRS_24_1 + PRS_24_2 + PRS_26 + PRS_27 +",
    paste(adjust_covariates, collapse = " + "),
    "+ (1 | Breed) + (1 | GMM_cluster)"
  ))
  model5_f2 <- lmer(f5_f2, data = train_fold, REML = TRUE,
                    control = lmerControl(optimizer = "bobyqa"))
  m5f2 <- evaluate_mixed_model(model5_f2, test_fold)
  
  row_f2 <- data.frame(
    repeat_num = r, fold_num = k,
    is_singular = lme4::isSingular(model5_f2),
    selected_G = best_G,
    n_dims_used = n_dims,                # ← 추가
    pca_var_explained = var_explained,   # ← 추가
    R2_model5_formula2_oos = m5f2$R2_out_of_sample,
    RMSE_model5_formula2   = m5f2$RMSE_out_of_sample,
    R2_model5_formula2_marginal    = m5f2$R2_marginal,
    R2_model5_formula2_conditional = m5f2$R2_conditional
  )
  
  list(f1 = row_f1, f2 = row_f2,
       breed_cluster_table = table(train_fold$Breed, train_fold$GMM_cluster))
}

stopCluster(cl)
close(pb)

results_model5_f1_df <- do.call(rbind, lapply(model5_results, `[[`, "f1"))
results_model5_f2_df <- do.call(rbind, lapply(model5_results, `[[`, "f2"))

cat("Model 5 (Breed + GMM cluster random effects, PCA-reduced + parallel): done.\n")

dim(results_model5_f1_df)  # 50행(REPEATS×K)이 다 채워졌는지 확인
sum(is.na(results_model5_f1_df$R2_model5_formula1_conditional))  # NA가 몇 개인지
sum(is.na(results_model5_f1_df$R2_model5_formula1_marginal))

# singular fit 비율 확인 (예전엔 100%였던 부분)
cat("Formula1 singular fit 비율:", mean(results_model5_f1_df$is_singular), "\n")
cat("Formula2 singular fit 비율:", mean(results_model5_f2_df$is_singular), "\n")

table(results_model5_f1_df$n_dims_used)

table(results_model5_f1_df$selected_G)

# ------------------------------------------------------------------------------
# [진단 결과 요약] 실제 사용된 차원 수(n_dims_used)와 그때의 분산 설명 비율
# ------------------------------------------------------------------------------
cat("\n=== 사용된 PCA 차원 수 (목표 누적분산 =", TARGET_VAR, ") ===\n")
cat("평균 n_dims_used:", round(mean(results_model5_f1_df$n_dims_used), 2), "\n") # 평균 n_dims_used: 68.14 
cat("범위:", min(results_model5_f1_df$n_dims_used), "~", max(results_model5_f1_df$n_dims_used), "\n") # 범위: 64 ~ 72 
print(table(results_model5_f1_df$n_dims_used)) 

cat("\n=== 실제 달성된 누적 분산 설명 비율 ===\n")
cat("평균:", round(mean(results_model5_f1_df$pca_var_explained), 4), "\n") # 평균: 0.7022 
cat("범위:", round(min(results_model5_f1_df$pca_var_explained), 4), "~",
    round(max(results_model5_f1_df$pca_var_explained), 4), "\n") # 범위: 0.7 ~ 0.7044 

# 마지막 fold 기준 breed x GMM_cluster 교차표
model5_results[[length(model5_results)]]$breed_cluster_table

write.csv(model5_results[[length(model5_results)]]$breed_cluster_table,
          file = "10_breed_cluster_table_method1model5.csv",
          row.names = TRUE)

# ==============================================================================
# STEP 6. 모든 결과(Baseline, Model 1~5)를 repeat_num, fold_num 기준으로 병합
# ==============================================================================


library(dplyr)

# ------------------------------------------------------------------------------
# 6-1. 충돌 가능한 컬럼명 정리 (원본 results_model*_df 는 그대로 두고 새 이름으로 저장)
# ------------------------------------------------------------------------------

results_model3_f1_df <- results_model3_f1_df %>%
  rename(is_singular_model3_f1 = is_singular)
head(results_model3_f1_df)

results_model3_f2_df <- results_model3_f2_df %>%
  rename(is_singular_model3_f2 = is_singular)
head(results_model3_f2_df)

results_model4_f1_df <- results_model4_f1_df %>%
  rename(is_singular_model4_f1 = is_singular,
         selected_G_model4_f1  = selected_G)
head(results_model4_f1_df)

results_model4_f2_df <- results_model4_f2_df %>%
  rename(is_singular_model4_f2 = is_singular,
         selected_G_model4_f2  = selected_G)
head(results_model4_f2_df)

results_model5_f1_df <- results_model5_f1_df %>%
  rename(is_singular_model5_f1       = is_singular,
         selected_G_model5_f1        = selected_G,
         n_dims_used_model5_f1       = n_dims_used,
         pca_var_explained_model5_f1 = pca_var_explained)
head(results_model5_f1_df)

results_model5_f2_df <- results_model5_f2_df %>%
  rename(is_singular_model5_f2       = is_singular,
         selected_G_model5_f2        = selected_G,
         n_dims_used_model5_f2       = n_dims_used,
         pca_var_explained_model5_f2 = pca_var_explained)
head(results_model5_f2_df)

# ------------------------------------------------------------------------------
# 6-2. Baseline + Model 1 + Model 2 (기존) + Model 3~5 (신규) 전체 병합
# ------------------------------------------------------------------------------

results_df <- results_baseline_df %>%
  left_join(results_model1_f1_df, by = c("repeat_num", "fold_num")) %>%
  left_join(results_model1_f2_df, by = c("repeat_num", "fold_num")) %>%
  left_join(results_model2_f1_df, by = c("repeat_num", "fold_num")) %>%
  left_join(results_model2_f2_df, by = c("repeat_num", "fold_num")) %>%
  left_join(results_model3_f1_df, by = c("repeat_num", "fold_num")) %>%
  left_join(results_model3_f2_df, by = c("repeat_num", "fold_num")) %>%
  left_join(results_model4_f1_df, by = c("repeat_num", "fold_num")) %>%
  left_join(results_model4_f2_df, by = c("repeat_num", "fold_num")) %>%
  left_join(results_model5_f1_df, by = c("repeat_num", "fold_num")) %>%
  left_join(results_model5_f2_df, by = c("repeat_num", "fold_num"))

head(results_df)

# ------------------------------------------------------------------------------
# 6-3. Incremental R² 계산 (병합 후 한 번에, Baseline 대비 증분)
# ------------------------------------------------------------------------------

# Model 1/2 (기존)
results_df$incremental_R2_formula1 <- results_df$R2_model2_formula1 - results_df$R2_baseline
results_df$incremental_R2_formula2 <- results_df$R2_model2_formula2 - results_df$R2_baseline

# Model 3 (out-of-sample R² 기준)
results_df$incremental_R2_model3_formula1 <- results_df$R2_model3_formula1_oos - results_df$R2_baseline
results_df$incremental_R2_model3_formula2 <- results_df$R2_model3_formula2_oos - results_df$R2_baseline

# Model 4
results_df$incremental_R2_model4_formula1 <- results_df$R2_model4_formula1_oos - results_df$R2_baseline
results_df$incremental_R2_model4_formula2 <- results_df$R2_model4_formula2_oos - results_df$R2_baseline

# Model 5
results_df$incremental_R2_model5_formula1 <- results_df$R2_model5_formula1_oos - results_df$R2_baseline
results_df$incremental_R2_model5_formula2 <- results_df$R2_model5_formula2_oos - results_df$R2_baseline

cat("\nAll models merged. Final results_df dimensions:", dim(results_df)[1], "rows,", dim(results_df)[2], "columns\n")

# ==============================================================================
# STEP 7. 결과 집계: 평균 + 95% CI
# ==============================================================================

summarize_metric <- function(x) {
  m  <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  data.frame(mean = m, se = se, ci_lower = m - 1.96 * se, ci_upper = m + 1.96 * se,
             n_valid = sum(!is.na(x)), n_na = sum(is.na(x)))
}

# ------------------------------------------------------------------------------
# 7-1. 연속형 성능 지표 (R², RMSE, incremental R², marginal/conditional R² 등)
# ------------------------------------------------------------------------------

metric_cols <- c(
  # Baseline
  "R2_baseline", "RMSE_baseline",
  
  # Model 1
  "R2_model1_formula1", "RMSE_model1_formula1",
  "R2_model1_formula2", "RMSE_model1_formula2",
  
  # Model 2
  "R2_model2_formula1", "RMSE_model2_formula1", "incremental_R2_formula1",
  "R2_model2_formula2", "RMSE_model2_formula2", "incremental_R2_formula2",
  
  # Model 3
  "R2_model3_formula1_oos", "RMSE_model3_formula1",
  "R2_model3_formula1_marginal", "R2_model3_formula1_conditional", "incremental_R2_model3_formula1",
  "R2_model3_formula2_oos", "RMSE_model3_formula2",
  "R2_model3_formula2_marginal", "R2_model3_formula2_conditional", "incremental_R2_model3_formula2",
  
  # Model 4
  "R2_model4_formula1_oos", "RMSE_model4_formula1",
  "R2_model4_formula1_marginal", "R2_model4_formula1_conditional", "incremental_R2_model4_formula1",
  "R2_model4_formula2_oos", "RMSE_model4_formula2",
  "R2_model4_formula2_marginal", "R2_model4_formula2_conditional", "incremental_R2_model4_formula2",
  
  # Model 5
  "R2_model5_formula1_oos", "RMSE_model5_formula1",
  "R2_model5_formula1_marginal", "R2_model5_formula1_conditional", "incremental_R2_model5_formula1",
  "n_dims_used_model5_f1", "pca_var_explained_model5_f1",
  "R2_model5_formula2_oos", "RMSE_model5_formula2",
  "R2_model5_formula2_marginal", "R2_model5_formula2_conditional", "incremental_R2_model5_formula2",
  "n_dims_used_model5_f2", "pca_var_explained_model5_f2"
)

summary_table <- do.call(rbind, lapply(metric_cols, function(col) {
  s <- summarize_metric(results_df[[col]])
  s$metric <- col
  s
}))
summary_table <- summary_table[, c("metric", "mean", "se", "ci_lower", "ci_upper", "n_valid", "n_na")]

cat("\n=== Nested CV Summary (Baseline / Model 1-5, Formula 1 & 2) ===\n")
print(summary_table, row.names = FALSE)

# ------------------------------------------------------------------------------
# 7-2. Singular fit 비율 요약 (Model 3~5)
# ------------------------------------------------------------------------------

singular_cols <- c(
  "is_singular_model3_f1", "is_singular_model3_f2",
  "is_singular_model4_f1", "is_singular_model4_f2",
  "is_singular_model5_f1", "is_singular_model5_f2"
)

singular_summary <- data.frame(
  metric = singular_cols,
  singular_rate = sapply(singular_cols, function(col) mean(results_df[[col]], na.rm = TRUE)),
  n_singular    = sapply(singular_cols, function(col) sum(results_df[[col]], na.rm = TRUE)),
  n_total       = sapply(singular_cols, function(col) sum(!is.na(results_df[[col]])))
)

cat("\n=== Singular Fit Rate Summary (Model 3-5) ===\n")
print(singular_summary, row.names = FALSE)

# ------------------------------------------------------------------------------
# 7-3. Model 4, 5의 selected_G 분포 요약
# ------------------------------------------------------------------------------

cat("\n=== selected_G distribution ===\n")
cat("Model 4 - Formula 1:\n"); print(table(results_df$selected_G_model4_f1))
cat("Model 4 - Formula 2:\n"); print(table(results_df$selected_G_model4_f2))
cat("Model 5 - Formula 1:\n"); print(table(results_df$selected_G_model5_f1))
cat("Model 5 - Formula 2:\n"); print(table(results_df$selected_G_model5_f2))

# ==============================================================================
# STEP 8. 결과 저장
# ==============================================================================

write.csv(results_df,        "10_PRS_all_models_nestedCV_perfold_Method1.csv",   row.names = FALSE)
write.csv(summary_table,     "10_PRS_all_models_nestedCV_summary_Method1.csv",   row.names = FALSE)
write.csv(singular_summary,  "10_PRS_all_models_singular_fit_summary_Method1.csv", row.names = FALSE)

cat("\n결과 저장 완료:\n")
cat(" - PRS_all_models_nestedCV_perfold.csv   (fold별 전체 원자료, ", nrow(results_df), "행)\n")
cat(" - PRS_all_models_nestedCV_summary.csv   (지표별 평균/CI 요약)\n")
cat(" - PRS_all_models_singular_fit_summary.csv (Model 3-5 singular fit 비율)\n")

