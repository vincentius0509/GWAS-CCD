library(data.table)
library(dplyr)
library(groupdata2)
library(lme4)
library(lmerTest)
library(mclust)
library(performance)
library(proxy)

set.seed(2025)
# ------------------------------------------------------------------------------
# STEP  0. 스크립트 골격과 SNP/블록 정의
# ------------------------------------------------------------------------------

# 4개 블록, 44개 SNP 전체 목록
block_map <- list(
  chr5  = c("5:25529184:C:T", "5:25529170:A:T"),
  chr24 = c("24:12787767:A:G","24:12787779:T:C","24:12787773:G:A",
            "24:12787779:T:A","24:12785743:C:T","24:12785679:G:A",
            "24:12782819:C:T","24:12782528:G:A","24:12780920:G:A",
            "24:12780924:T:C","24:12752625:A:T","24:12724858:A:T",
            "24:12724860:A:T","24:12818735:A:G","24:12818736:T:A",
            "24:12753225:G:A","24:12754140:A:G"),
  chr26 = c("26:13787793:A:G","26:13787261:C:T","26:13770897:T:G",
            "26:13783025:T:C","26:13393883:C:T","26:13259438:T:C",
            "26:13435785:G:T","26:14169952:T:A","26:13235107:G:A",
            "26:13185590:C:T","26:13196151:G:A","26:14504206:G:A",
            "26:13097561:T:C","26:13099805:G:A","26:13104092:C:T"),
  chr27 = c("27:44543118:A:T","27:44801518:G:A","27:44847773:C:T",
            "27:44860521:T:C","27:44790338:C:T","27:45065570:A:G",
            "27:44235256:G:C","27:44941531:C:T","27:44520608:T:G",
            "27:44929991:C:T")
)
block_map
snp_cols_raw <- unlist(block_map, use.names = FALSE)  # 44개 SNP ID (원래 표기: "5:25529184:C:T")

# 대립유전자(REF:ALT)까지 포함해서 컬럼명을 만들어 좌표 충돌 방지
snp_to_colname <- function(snp_id) {
  paste0("snp", gsub(":", "_", snp_id))   # ":" 제거하는 sub() 부분을 뺌
}
snp_cols <- sapply(snp_cols_raw, snp_to_colname)

# 확인: 32개가 전부 고유한 이름인지
cat("중복 컬럼명 여부:", any(duplicated(snp_cols)), "\n")
print(snp_cols)


# ------------------------------------------------------------------------------
# STEP 1. 원본 파일 로딩 및 병합
# ------------------------------------------------------------------------------

plink_001_1e6_10 <- fread("05_GWAS_PHE1_001_1e6_10_linear.phe1.glm.linear")
head(plink_001_1e6_10)

############################################################ 1 Data processing - genotype #####################################################

# ------------------------------------------------------------------------------
# 1. genotype .raw 파일 일괄 로딩 함수
#      기존 코드(geno5_25529184_C_T <- fread("snp_5_25529184_C_T_genotype.raw") 등)를
#      보면 파일명 규칙이 두 가지로 섞여 있습니다:
#        - "snp_{chr}_{pos}_{ref}_{alt}_genotype.raw" (일부 index SNP)
#        - "snp_{chr}_{pos}_genotype.raw" (대부분의 LD-tagged SNP)
#      두 패턴을 순서대로 시도하는 방식으로 통일합니다.
# ------------------------------------------------------------------------------

geno_dir <- "12_Geno_Beta/"   # 기존 코드의 실제 경로로 교체

load_one_snp_geno <- function(snp_id, snp_colname) {
  parts <- strsplit(snp_id, ":")[[1]]
  chr <- parts[1]; pos <- parts[2]; ref <- parts[3]; alt <- parts[4]
  
  fname_with_allele <- paste0(geno_dir, "snp_", chr, "_", pos, "_", ref, "_", alt, "_genotype.raw")
  fname_no_allele    <- paste0(geno_dir, "snp_", chr, "_", pos, "_genotype.raw")
  
  fname <- if (file.exists(fname_with_allele)) fname_with_allele else fname_no_allele
  if (!file.exists(fname)) stop("파일을 찾을 수 없습니다: ", snp_id)
  
  dt <- fread(fname)
  snp_col <- grep("^\\d", colnames(dt), value = TRUE)[1]
  out <- dt[, .(FID, IID, dosage = get(snp_col))]
  setnames(out, "dosage", snp_colname)
  out
}

geno_list <- Map(load_one_snp_geno, snp_cols_raw, snp_cols)
geno_df <- Reduce(function(x, y) merge(x, y, by = c("FID","IID"), all = TRUE), geno_list)

dim(geno_df)   # 확인: 행 수 = 표본 수, 열 수 = 2(FID,IID) + 32(SNP)
head(geno_df)

# ------------------------------------------------------------------------------
# 1. covariates, cslb 파일 로딩
# ------------------------------------------------------------------------------

# 3. PCs,Age,Sex,Breed,Neuter Status, Activity Level, CSLB score
covariates <- fread("04_cov_001_1e6_10.txt")
head(covariates)

cslb <- fread("../04_pheno001.txt")
head(cslb)

setnames(cslb, old = "phe1", new = "cslb_score")

base_df <- merge(
  covariates,
  cslb[, .(FID, IID, cslb_score)],
  by = c("FID", "IID"),
  all.x = TRUE
)
head(base_df)

# ------------------------------------------------------------------------------
# 1. Breed 원본 라벨 병합 
# ------------------------------------------------------------------------------

purebred_total <- read.csv("../00_Purebred_total.csv", stringsAsFactors = FALSE)

base_df <- merge(
  as.data.frame(base_df),
  purebred_total[, c("IID", "Breed")],
  by = "IID",
  all.x = TRUE
)
base_df$Breed <- gsub(" ", "_", base_df$Breed)

sum(is.na(base_df$Breed))   # 확인: 매칭 안 된 개체가 있는지

# ------------------------------------------------------------------------------
# 1-6. 최종 full_df 병합 (covariates + cslb + Breed + genotype)
# ------------------------------------------------------------------------------
full_df <- merge(
  base_df,
  geno_df,
  by = c("FID", "IID"),
  all.x = TRUE
)

dim(full_df)
head(full_df)

# ------------------------------------------------------------------------------
# 2-1. [PRS-1.R에서 그대로 복사] 조정 공변량 정의
# ------------------------------------------------------------------------------
adjust_covariates <- c("Age_cov", "Sex_cov", "Intact_cov", paste0("PC", 1:8), "pa_score_cov")

# ------------------------------------------------------------------------------
# 2-2. R²/RMSE 계산 함수
# ------------------------------------------------------------------------------
compute_r2_rmse <- function(actual, predicted) {
  r2   <- cor(actual, predicted, use = "complete.obs")^2
  rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
  list(R2 = r2, RMSE = rmse)
}

# ------------------------------------------------------------------------------
# 2-3. mixed model 평가 함수
#      (out-of-sample R²/RMSE + marginal/conditional R² 동시 산출)
# ------------------------------------------------------------------------------
evaluate_mixed_model <- function(model, test_fold, actual_col = "cslb_score") {
  pred <- predict(model, newdata = test_fold, allow.new.levels = TRUE)
  oos  <- compute_r2_rmse(test_fold[[actual_col]], pred)
  
  r2_decomp <- tryCatch(performance::r2(model), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
  
  list(
    R2_out_of_sample    = oos$R2,
    RMSE_out_of_sample  = oos$RMSE,
    R2_marginal         = as.numeric(r2_decomp$R2_marginal),
    R2_conditional       = as.numeric(r2_decomp$R2_conditional)
  )
}

# ------------------------------------------------------------------------------
# 2-4. pseudo-admixture 벡터화 계산 함수 (Model 5용)
# ------------------------------------------------------------------------------
calc_pseudo_admix <- function(pc_mat, cent_mat) {
  dist_mat <- as.matrix(proxy::dist(pc_mat, cent_mat, method = "Euclidean"))
  out <- 1 / (dist_mat + 1e-8)
  out <- out / rowSums(out)
  colnames(out) <- rownames(cent_mat)
  out
}

# ------------------------------------------------------------------------------
# 2-5. prepare_fold_data() — 32개 SNP / 4개 block 버전
# ------------------------------------------------------------------------------

prepare_fold_data <- function(train_fold, test_fold) {
  
  # (a) 결측치 대체: train fold의 allele frequency만 사용
  af_train <- sapply(snp_cols, function(s) mean(train_fold[[s]], na.rm = TRUE) / 2)
  names(af_train) <- snp_cols
  for (s in snp_cols) {
    train_fold[[s]][is.na(train_fold[[s]])] <- 2 * af_train[s]
    test_fold[[s]][is.na(test_fold[[s]])]   <- 2 * af_train[s]
  }
  
  # (b) 32개 SNP 각각 effect size(beta) 재추정: train fold에서만
  beta_refit <- sapply(snp_cols, function(s) {
    f <- as.formula(paste("cslb_score ~", s, "+", paste(adjust_covariates, collapse = " + ")))
    coef(lm(f, data = train_fold))[s]
  })
  names(beta_refit) <- snp_cols
  
  # (c) Formula 1: genome-wide 단일 PRS (32개 SNP 전체 가중합)
  train_fold$PRS <- as.numeric(as.matrix(train_fold[, snp_cols]) %*% beta_refit)
  test_fold$PRS  <- as.numeric(as.matrix(test_fold[, snp_cols])  %*% beta_refit)
  
  # (d) Formula 2: 4개 block(chr5, chr24,chr26, chr27)별 locus-level PRS
  for (blk in names(block_map)) {
    blk_snp_cols <- sapply(block_map[[blk]], snp_to_colname)
    train_fold[[paste0("PRS_", blk)]] <- as.numeric(as.matrix(train_fold[, blk_snp_cols]) %*% beta_refit[blk_snp_cols])
    test_fold[[paste0("PRS_", blk)]]  <- as.numeric(as.matrix(test_fold[, blk_snp_cols])  %*% beta_refit[blk_snp_cols])
  }
  
  list(train = train_fold, test = test_fold, beta_refit = beta_refit, af_train = af_train)
}

# ==============================================================================
# STEP 3. analysis_df, fold_assignments 준비
# ==============================================================================

# ------------------------------------------------------------------------------
# 3-1. 분석 대상 데이터 준비 (CSLB 점수 결측만 제외, genotype 결측은 fold 내에서 처리)
# ------------------------------------------------------------------------------

analysis_df <- full_df[!is.na(full_df$cslb_score), ]
analysis_df <- as.data.frame(analysis_df)

analysis_df$IID   <- factor(analysis_df$IID)
analysis_df$Breed <- factor(analysis_df$Breed)   # 문자형 → factor (Model 2의 고정효과, Model 3~5의 랜덤효과 모두에 사용)

row.names(analysis_df) <- NULL   # 인덱스를 1:N으로 재정렬 (fold 인덱스 관리를 위해 필수)

cat("Analysis dataset dimensions:", dim(analysis_df)[1], "rows,", dim(analysis_df)[2], "columns\n")
cat("Number of Breed levels:", nlevels(analysis_df$Breed), "\n")

# ------------------------------------------------------------------------------
# 3-2. CSLB 점수 분포 기반 층화를 위한 구간 나누기
# ------------------------------------------------------------------------------

qs <- unique(quantile(analysis_df$cslb_score, probs = seq(0, 1, 0.2), na.rm = TRUE))
analysis_df$cslb_strata <- cut(analysis_df$cslb_score, breaks = qs, include.lowest = TRUE)

# ------------------------------------------------------------------------------
# 3-3. Repeated K-fold 설정
# ------------------------------------------------------------------------------

K       <- 10   # fold 수
REPEATS <- 5    # 반복 횟수 (계산 시간 보고 조정 가능)

# ------------------------------------------------------------------------------
# 3-4. Fold 배정을 한 번만 생성해서 모든 모델이 동일하게 공유
#      fold_assignments[[r]] = 해당 repeat의 fold_id 벡터 (analysis_df와 같은 순서/길이)
# ------------------------------------------------------------------------------

fold_assignments <- list()
for (r in 1:REPEATS) {
  set.seed(2025 + r)
  folds <- groupdata2::fold(analysis_df, k = K, cat_col = "cslb_strata", id_col = "IID")
  fold_assignments[[r]] <- folds$.folds
}

cat("Fold assignments created for", REPEATS, "repeats x", K, "folds.\n")

analysis_df <- merge(
  base_df,
  geno_df,
  by = c("FID", "IID"),
  all.x = TRUE
)

dim(analysis_df)
head(analysis_df)

# ------------------------------------------------------------------------------
# 희귀 품종을 "Other_rare_breed"로 그룹핑 (고정효과 모델 전용)
#   - 표본 크기만 기준으로 하는 규칙이라 outcome(cslb_score)과 무관 → data leakage 아님
#   - Model 3~5의 랜덤효과에는 이 그룹핑을 적용하지 않고, 원본 개별 품종 라벨을 그대로 사용
# ------------------------------------------------------------------------------

BREED_MIN_N <- 5   # 이 값 미만이면 "Other_rare_breed"로 그룹핑 (기준값은 조정 가능)

breed_counts <- table(analysis_df$Breed)
rare_breeds  <- names(breed_counts[breed_counts < BREED_MIN_N])

cat("전체 품종 수:", length(breed_counts), "\n")
cat(BREED_MIN_N, "마리 미만인 희귀 품종 수:", length(rare_breeds), "\n")
cat("그 희귀 품종에 속한 개체 수 합계:", sum(breed_counts[rare_breeds]), "\n")

# 고정효과 모델(Baseline/Model1/Model2) 전용 변수 생성
analysis_df$Breed_grouped <- ifelse(
  analysis_df$Breed %in% rare_breeds,
  "Other_rare_breed",
  as.character(analysis_df$Breed)
)
analysis_df$Breed_grouped <- factor(analysis_df$Breed_grouped)

cat("그룹핑 후 최종 레벨 수:", nlevels(analysis_df$Breed_grouped), "\n")

# ------------------------------------------------------------------------------
# STEP 4-A. Baseline 모델 (PRS 없음) - 독립된 for문
# ------------------------------------------------------------------------------

str(geno_df)  # 전체 컬럼 타입을 한 번에 확인

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
    
    # ------------------------------------------------------------------------
    # 모델 적합 (train_fold의 Breed_grouped factor 레벨을 그대로 사용)
    # ------------------------------------------------------------------------
    f_baseline <- as.formula(paste("cslb_score ~", paste(c(adjust_covariates, "Breed_grouped"), collapse = " + ")))
    model_baseline <- lm(f_baseline, data = train_fold)
    
    # ------------------------------------------------------------------------
    # [핵심] 모델이 실제로 기억하는 레벨(xlevels)을 test_fold에 그대로 적용
    #        -> 수동 재계산으로 인한 레벨 불일치 문제를 원천적으로 방지
    # ------------------------------------------------------------------------
    model_levels <- model_baseline$xlevels$Breed_grouped
    test_fold$Breed_grouped <- factor(as.character(test_fold$Breed_grouped), levels = model_levels)
    
    n_before <- nrow(test_fold)
    test_fold <- test_fold[!is.na(test_fold$Breed_grouped), ]
    n_excluded <- n_before - nrow(test_fold)
    if (n_excluded > 0) {
      cat("경고: fold", k, "(repeat", r, ") -", n_excluded,
          "개체 제외 (train 모델이 모르는 Breed_grouped 레벨)\n")
    }
    
    # ------------------------------------------------------------------------
    # 예측 및 성능 평가
    # ------------------------------------------------------------------------
    pred <- predict(model_baseline, newdata = test_fold)
    metrics <- compute_r2_rmse(test_fold$cslb_score, pred)
    
    results_baseline[[idx]] <- data.frame(
      repeat_num = r, fold_num = k,
      n_excluded = n_excluded,
      R2_baseline = metrics$R2, RMSE_baseline = metrics$RMSE
    )
    idx <- idx + 1
  }
}

results_baseline_df <- do.call(rbind, results_baseline)
cat("Baseline model: done.\n")

dim(results_baseline_df)
summary(results_baseline_df$R2_baseline)
cat("\n총 제외된 개체 수 (전체 50 fold 합계):", sum(results_baseline_df$n_excluded), "\n")

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
    
    # Breed_grouped를 쓰지 않는 회귀식이므로 레벨 매칭 불필요
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

summary(results_model1_f1_df$R2_model1_formula1)

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
    
    model1_f2 <- lm(cslb_score ~ PRS_chr5 + PRS_chr24 + PRS_chr26 + PRS_chr27, data = train_fold)
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
summary(results_model1_f2_df$R2_model1_formula2)

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
    
    f_model2_f1 <- as.formula(paste("cslb_score ~", paste(c("PRS", adjust_covariates, "Breed_grouped"), collapse = " + ")))
    model2_f1 <- lm(f_model2_f1, data = train_fold)
    
    # [핵심] 모델이 실제로 학습한 레벨(xlevels)을 test_fold에 그대로 적용
    model_levels <- model2_f1$xlevels$Breed_grouped
    test_fold$Breed_grouped <- factor(as.character(test_fold$Breed_grouped), levels = model_levels)
    
    n_before <- nrow(test_fold)
    test_fold <- test_fold[!is.na(test_fold$Breed_grouped), ]
    n_excluded <- n_before - nrow(test_fold)
    if (n_excluded > 0) {
      cat("경고: fold", k, "(repeat", r, ") -", n_excluded,
          "개체 제외 (train 모델이 모르는 Breed_grouped 레벨)\n")
    }
    
    pred <- predict(model2_f1, newdata = test_fold)
    metrics <- compute_r2_rmse(test_fold$cslb_score, pred)
    
    results_model2_f1[[idx]] <- data.frame(
      repeat_num = r, fold_num = k,
      n_excluded = n_excluded,
      R2_model2_formula1 = metrics$R2, RMSE_model2_formula1 = metrics$RMSE
    )
    idx <- idx + 1
  }
}
results_model2_f1_df <- do.call(rbind, results_model2_f1)
cat("Model 2 - Formula 1: done.\n")
summary(results_model2_f1_df$R2_model2_formula1)

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
      paste(c("PRS_chr5", "PRS_chr24", "PRS_chr26", "PRS_chr27", adjust_covariates, "Breed_grouped"), collapse = " + ")
    ))
    model2_f2 <- lm(f_model2_f2, data = train_fold)
    
    # [핵심] 이전에는 이 부분이 통째로 빠져 있었음 -> 추가
    model_levels <- model2_f2$xlevels$Breed_grouped
    test_fold$Breed_grouped <- factor(as.character(test_fold$Breed_grouped), levels = model_levels)
    
    n_before <- nrow(test_fold)
    test_fold <- test_fold[!is.na(test_fold$Breed_grouped), ]
    n_excluded <- n_before - nrow(test_fold)
    if (n_excluded > 0) {
      cat("경고: fold", k, "(repeat", r, ") -", n_excluded,
          "개체 제외 (train 모델이 모르는 Breed_grouped 레벨)\n")
    }
    
    pred <- predict(model2_f2, newdata = test_fold)
    metrics <- compute_r2_rmse(test_fold$cslb_score, pred)
    
    results_model2_f2[[idx]] <- data.frame(
      repeat_num = r, fold_num = k,
      n_excluded = n_excluded,
      R2_model2_formula2 = metrics$R2, RMSE_model2_formula2 = metrics$RMSE
    )
    idx <- idx + 1
  }
}
results_model2_f2_df <- do.call(rbind, results_model2_f2)
cat("Model 2 - Formula 2: done.\n")
summary(results_model2_f2_df$R2_model2_formula2)

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


# ==============================================================================
# 단일 fold 진단 - PRS 분산, block별 분산/상관, beta 이상치 확인
#   Model 3 for문 (STEP 5-F) 돌리기 전에 실행
# ==============================================================================

cat("========================================\n")
cat("진단: 단일 fold(repeat=1, fold=1)로 PRS/beta 상태 점검\n")
cat("========================================\n\n")

diag_r <- 1
diag_k <- 1

analysis_df$fold_id <- fold_assignments[[diag_r]]
diag_train <- analysis_df[analysis_df$fold_id != diag_k, ]
diag_test  <- analysis_df[analysis_df$fold_id == diag_k, ]

diag_prepped <- prepare_fold_data(diag_train, diag_test)
diag_train   <- diag_prepped$train

# ------------------------------------------------------------------------------
# 1. PRS(genome-wide, 32개 SNP 합산) 분산 확인
# ------------------------------------------------------------------------------
cat("--- 1. PRS (genome-wide) 분산 확인 ---\n")
cat("sd(PRS):", sd(diag_train$PRS), "\n")
cat("range(PRS):", paste(round(range(diag_train$PRS), 4), collapse = " ~ "), "\n")
cat("mean(PRS):", mean(diag_train$PRS), "\n\n")
 
# ------------------------------------------------------------------------------
# 2. Block별 PRS 분산 및 서로 간 상관관계 확인
# ------------------------------------------------------------------------------
cat("--- 2. Block별 PRS 분산 ---\n")
block_prs_cols <- paste0("PRS_", names(block_map))
for (col in block_prs_cols) {
  cat(sprintf("  sd(%s): %.5f, range: %.4f ~ %.4f\n",
              col, sd(diag_train[[col]]), min(diag_train[[col]]), max(diag_train[[col]])))
}
 
cat("\n--- Block 간 상관관계 행렬 ---\n")
print(cor(diag_train[, block_prs_cols]))
cat("\n")
 
# ------------------------------------------------------------------------------
# 3. 재추정된 beta 값 중 이상치(비정상적으로 크거나 작은 값) 확인
# ------------------------------------------------------------------------------
cat("--- 3. 재추정된 beta (32개 SNP) ---\n")
beta_check <- diag_prepped$beta_refit
print(sort(beta_check))
 
cat("\nbeta 절댓값 상위 5개 (이상치 의심):\n")
print(sort(abs(beta_check), decreasing = TRUE)[1:5])
 
cat("\n해당 SNP들의 train fold 내 allele frequency (극단적으로 낮은/높은지 확인):\n")
extreme_snps <- names(sort(abs(beta_check), decreasing = TRUE)[1:5])
print(diag_prepped$af_train[extreme_snps])
 
cat("\n========================================\n")
cat("진단 완료. 위 결과를 확인한 뒤 Model 3(STEP 5-F)으로 진행하세요.\n")
cat("========================================\n")

# ------------------------------------------------------------------------------
# 4. [신규] PRS 전체 NA의 원인 SNP 찾기
# ------------------------------------------------------------------------------
cat("--- 4. train_fold 내 SNP별 결측 개수 확인 ---\n")

na_counts <- colSums(is.na(diag_train[, snp_cols]))
print(sort(na_counts, decreasing = TRUE))


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
      "cslb_score ~ PRS_chr5 + PRS_chr24 + PRS_chr26 + PRS_chr27 +",
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
      "cslb_score ~ PRS_chr5 + PRS_chr24 + PRS_chr26 + PRS_chr27 +",
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

write.csv(table(train_fold$Breed, train_fold$PC_cluster),
          file = "10_breed_cluster_table_method2model4.csv",
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
    "cslb_score ~ PRS_chr5 + PRS_chr24 + PRS_chr26 + PRS_chr27 +",
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
cat("평균 n_dims_used:", round(mean(results_model5_f1_df$n_dims_used), 2), "\n") # 평균 n_dims_used: 67.94 
cat("범위:", min(results_model5_f1_df$n_dims_used), "~", max(results_model5_f1_df$n_dims_used), "\n") # 범위: 61 ~ 74 
print(table(results_model5_f1_df$n_dims_used))

cat("\n=== 실제 달성된 누적 분산 설명 비율 ===\n")
cat("평균:", round(mean(results_model5_f1_df$pca_var_explained), 4), "\n") # 평균: 0.7024  
cat("범위:", round(min(results_model5_f1_df$pca_var_explained), 4), "~",
    round(max(results_model5_f1_df$pca_var_explained), 4), "\n") # 범위: 0.7 ~ 0.7043 

# 마지막 fold 기준 breed x GMM_cluster 교차표
model5_results[[length(model5_results)]]$breed_cluster_table

write.csv(model5_results[[length(model5_results)]]$breed_cluster_table,
          file = "10_breed_cluster_table_method2model5.csv",
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

write.csv(results_df,        "10_PRS2_all_models_nestedCV_perfold.csv",   row.names = FALSE)
write.csv(summary_table,     "10_PRS2_all_models_nestedCV_summary.csv",   row.names = FALSE)
write.csv(singular_summary,  "10_PRS2_all_models_singular_fit_summary.csv", row.names = FALSE)

cat("\n결과 저장 완료:\n")
cat(" - PRS2_all_models_nestedCV_perfold.csv   (fold별 전체 원자료, ", nrow(results_df), "행)\n")
cat(" - PRS2_all_models_nestedCV_summary.csv   (지표별 평균/CI 요약)\n")
cat(" - PRS2_all_models_singular_fit_summary.csv (Model 3-5 singular fit 비율)\n")

