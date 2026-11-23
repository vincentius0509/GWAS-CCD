# 00_datacheck.R - Checks the DAP 2023 release files and defines the analysis sample, covariates and phenotype (most recent CSLB score). See README.md.
# ==============================================================================
#   00_Picture/05_GWAS_2023/00_*.png

suppressMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "/home/yhSong/data/02_ByProject/DogAgingProject"
geno_dir <- file.path(base_dir, "GeneticData/v2023")
geno_pre <- file.path(geno_dir, "DogAgingProject_GeneticData_CuratedRelease_2023")
out_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")   
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_prefix <- file.path(out_dir, "00")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
dir.create(pic_dir, recursive = TRUE, showWarnings = FALSE)

f_overview <- file.path(base_dir, "DogOverview/v2023/DAP_2023_DogOverview_v1.1.csv")
f_cslb <- file.path(base_dir, "SurveyData/v2023/DAP_2023_CSLB_v1.0.csv")
f_hles <- file.path(base_dir, "SurveyData/v2023/DAP_2023_HLES_dog_owner_v1.1.csv")
f_afus <- file.path(base_dir, "SurveyData/v2023/DAP_2023_AFUS_dog_owner_v1.1.csv")
f_fam22 <- file.path(base_dir, "GeneticData/v2022/DogAgingProject_GeneticData_CuratedRelease_2022.fam")
f_hles_hs <- file.path(base_dir, "SurveyData/v2023/DAP_2023_HLES_health_condition_v1.0.csv")
f_afus_hs <- file.path(base_dir, "SurveyData/v2023/DAP_2023_AFUS_health_condition_v1.0.csv")
f_trt <- file.path(base_dir, "AdditionalStudies/v2023/DAP_2023_CognitiveGames_123Treat_v1.0.csv")

cb  <- fread("/data/02_ByProject/DogAgingProject/DAP_2023_CODEBOOK_v1.1.csv")
row <- cb[DataFile == "HLES_health_condition" & Variable == "hs_condition"]
map <- data.table(code  = trimws(strsplit(row$Values, "\\|")[[1]]),
                  label = trimws(strsplit(row$ValueLabels, "\\|")[[1]]))
map[code %in% c(101, 102, 103, 115, 118, 119, 201, 205, 1113, 1201, 1202, 1205, 1215, 1216)]

vision_codes <- c(101, 102, 115, 103, 118, 119)
vision_codes_narrow <- c(101, 102, 115)
hearing_codes <- c(201, 205)
seizure_codes <- 1215
vestibular_codes <- 1216
dementia_codes <- 1205
oa_codes <- 1113
brain_codes <- c(1201, 1202)
other_neuro_recode <- list(
  dementia   = "cognitive|dementia|senil|\\bccd\\b|\\bcds\\b",
  vestibular = "vestibular|vertigo|nystagmus",
  seizure    = "seizure|epilep",
  brain      = "infarct|\\btia\\b|stroke|batten|\\bncl\\b|meningitis|encephal|cerebel|brain|hydroceph")

own_cut <- c(discordant = 0.5, confirmed = 0.85, admixed = 0.45)
age_floor_main <- 8
age_floor_cand <- c(0, 6, 7, 8, 9, 10, 11, 12)
lbs2kg <- 0.4536

belgian <- c("belgian_malinois", "belgian_tervuren", "belgian_sheepdog")
aussie <- "australian_shepherd"
alias <- list(
  "Poodle"           = c("poodle", "toy_poodle"),
  "Poodle (Toy)"     = c("toy_poodle", "poodle"),
  "Bulldog"          = "english_bulldog",
  "St. Bernard"      = "saint_bernard",
  "Chinese Shar-Pei" = "shar_pei",
  "Cocker Spaniel"   = "american_cocker_spaniel",
  "Entlebucher Mountain Dog" = "entlebucher",
  "Finnish Spitz"    = "finish_spitz",
  "Belgian Malinois" = belgian, "Belgian Tervuren" = belgian, "Belgian Sheepdog" = belgian, "Belgian Laekenois" = belgian,
  "Miniature American Shepherd" = aussie, "Miniature Australian Shepherd" = aussie,
  "Mini Australian Shepherd" = aussie, "Miniature (Toy) Australian Shepherd" = aussie,
  "Russell Terrier"  = "jack_russell_terrier",
  "Parson Russell Terrier" = "jack_russell_terrier",
  "Miniature Bull Terrier" = "bull_terrier"
)
norm <- function(s) gsub("[^a-z]", "", tolower(s))

cat("1) release 점검\n")
fam <- fread(paste0(geno_pre, "_ed.fam"), header = FALSE, col.names = c("FID", "IID", "PID", "MID", "sex_fam", "pheno"))
fam_raw <- fread(paste0(geno_pre, ".fam"), header = FALSE)
stopifnot(identical(fam$IID, fam_raw$V2), identical(fam$sex_fam, fam_raw$V5))
cat(sprintf("   genotype %d마리, 원본·수정 .fam 개체 순서·성별 동일, 수정본 FID = IID: %s\n",
            nrow(fam), all(fam$FID == fam$IID)))

meta <- fread(paste0(geno_pre, ".tsv"))[, .(IID = dog, sex_meta = sex, coi, size_score = size)]
head(meta)

cat(sprintf("   메타데이터(.tsv) %d마리, genotype과 일치: %s\n", nrow(meta), setequal(meta$IID, fam$IID)))

sex_tab <- merge(fam[, .(IID, sex_fam)], meta[, .(IID, sex_meta)], by = "IID")[, .N, by = .(sex_fam, sex_meta)]
head(sex_tab)
#    sex_fam sex_meta     N
#      <int>   <char> <int>
# 1:       2        M  3892
# 2:       1        F  3723
# 3:       1        M    12

print(dcast(sex_tab, sex_fam ~ sex_meta, value.var = "N", fill = 0))
code_F <- sex_tab[sex_meta == "F"][which.max(N), sex_fam]
fam[, Sex_fam := fifelse(sex_fam == code_F, "Female", "Male")]
cat(sprintf("   .fam 성별 코드 %d = Female (DAP 성별 확인 결과 기준)\n", code_F))

sex_x <- meta$sex_meta[match(fam$IID, meta$IID)]
stopifnot(!anyNA(sex_x), all(sex_x %in% c("F", "M")))
fam[, Sex := fifelse(sex_x == "F", "Female", "Male")]
cat(sprintf("   성별은 sex_meta(X 염색체 추정치) 사용, .fam 성별과 다른 개체: %d마리 (제외하지 않고 정정)\n", fam[Sex != Sex_fam, .N]))

fam22 <- fread(f_fam22, header = FALSE)$V2
head(fam22)

fam[, in_v2022 := IID %in% fam22]
release_tab <- data.table(item = c("v2022 genotype", "v2023 genotype", "both", "v2022 only (v2023에서 제외)", "v2023 only (신규)"),
                          N = c(length(fam22), nrow(fam), sum(fam$in_v2022), length(setdiff(fam22, fam$IID)), sum(!fam$in_v2022)))
print(release_tab)
#                           item     N
#                         <char> <int>
# 1:              v2022 genotype  7425
# 2:              v2023 genotype  7627
# 3:                        both  7378

cat("\n2) DogOverview\n")
ov <- fread(f_overview, select = c("dog_id", "Breed_Status", "Breed", "Estimated_DOB", "Sex_Class_at_HLES",
                                   "Breed_Size_Class_at_HLES", "Status", "Cohort", "AgeDOB_Estimation_Method"))
dt <- merge(fam[, .(FID, IID, Sex, Sex_fam, in_v2022)], ov, by.x = "IID", by.y = "dog_id", all.x = TRUE)
dt <- merge(dt, meta, by = "IID", all.x = TRUE)
dt[, `:=`(Estimated_DOB = as.Date(Estimated_DOB),
          Intact_cov = fifelse(grepl("intact", Sex_Class_at_HLES, ignore.case = TRUE), 0L, 1L),
          sex_owner = fifelse(grepl("^Female", Sex_Class_at_HLES), "Female", "Male"))]
print(dt[, .N, by = Breed_Status])
cat(sprintf("   보호자 응답 성별과 불일치: Sex(sex_meta) %d마리, .fam %d마리 (성별은 sex_meta 기준 사용)\n",
            dt[!is.na(Sex_Class_at_HLES) & Sex != sex_owner, .N], dt[!is.na(Sex_Class_at_HLES) & Sex_fam != sex_owner, .N]))

head(dt)

# 3. ancestry
cat("\n3) ancestry\n")
anc <- fread(paste0(geno_pre, ".ancestry.tsv"))
head(anc)

anc[, `:=`(dog = as.character(dog), pct = as.numeric(pct))]
pops <- unique(anc$pop)
setorder(anc, dog, -pct)
anc_sum <- anc[, .(top_pop = pop[1], top_pct = pct[1], second_pop = pop[2], second_pct = pct[2],
                   n_pop_05 = sum(pct >= 0.05)), by = dog]
head(anc_sum)

dt <- merge(dt[, IID_chr := as.character(IID)], anc_sum, by.x = "IID_chr", by.y = "dog", all.x = TRUE)
head(dt)

breeds <- dt[Breed_Status == "Purebred", .N, by = Breed]
breeds[, panel_pops := lapply(Breed, function(b) {
  if (b %in% names(alias)) return(alias[[b]])
  hit <- pops[norm(pops) == norm(b)]
  if (length(hit)) hit else character(0) })]
breeds[, in_panel := lengths(panel_pops) > 0]
fwrite(breeds[order(-N), .(Breed, N, in_panel, panel_pops = vapply(panel_pops, paste, "", collapse = ";"))],
       paste0(out_prefix, "_breed_panel_map.csv"))
cat(sprintf("   순종 보고 품종 %d개 중 참조 패널에 있음 %d개 (개체 %d / %d)\n",
            nrow(breeds), sum(breeds$in_panel), breeds[in_panel == TRUE, sum(N)], sum(breeds$N)))

pop_of <- setNames(breeds$panel_pops, breeds$Breed)
dt[, own_pct := NA_real_]
dt[Breed_Status == "Purebred" & Breed %in% breeds[in_panel == TRUE, Breed],
   own_pct := mapply(function(id, b) anc[dog == id & pop %in% pop_of[[b]], sum(pct)], IID_chr, Breed)]

dt[, top_is_own := mapply(function(b, tp) !is.na(tp) && b %in% names(pop_of) && tp %in% pop_of[[b]], Breed, top_pop)]

yj <- rbindlist(lapply(seq(0.5, 0.99, 0.01), function(t) 
dt[!is.na(top_pct), .(t = t,
        sens = mean(top_pct[Breed_Status == "Purebred"] >= t), spec = mean(top_pct[Breed_Status != "Purebred"] < t))]))[, J := sens + spec - 1]

cat(sprintf("   순종 경계 확인: 문헌 0.85의 J = %.3f (민감도 %.3f, 특이도 %.3f) / 데이터 최적(Youden) %.2f의 J = %.3f\n",
            yj[abs(t - 0.85) < 1e-9, J], yj[abs(t - 0.85) < 1e-9, sens], yj[abs(t - 0.85) < 1e-9, spec], yj[which.max(J), t], max(yj$J)))

dt[, admix_class := fcase(top_pct >= own_cut[["confirmed"]], "single_breed",
                          top_pct >= own_cut[["admixed"]], "intermediate",
                          !is.na(top_pct), "highly_admixed")]
print(dcast(dt[, .N, by = .(Breed_Status, admix_class)], Breed_Status ~ admix_class, value.var = "N", fill = 0))
# Key: <Breed_Status>
#    Breed_Status highly_admixed intermediate single_breed
#          <char>          <int>        <int>        <int>
# 1:  Mixed breed           2021         1464          157
# 2:     Purebred            402          464         3119

dt[, genetic_class := fcase(
  Breed_Status != "Purebred", "mixed",
  is.na(own_pct), "pure_unverifiable",
  own_pct >= own_cut[["confirmed"]], "pure_confirmed",
  own_pct >= own_cut[["discordant"]] | top_is_own, "pure_likely",
  default = "pure_mismatch")]

em1 <- function(x, iter = 5000, tol = 1e-10) {
  km <- kmeans(x, 2, nstart = 20)
  mu <- as.numeric(tapply(x, km$cluster, mean)); s <- as.numeric(tapply(x, km$cluster, sd))
  w <- as.numeric(table(km$cluster)) / length(x); ll_old <- -Inf
  for (it in seq_len(iter)) {
    dn <- sapply(1:2, function(k) w[k] * dnorm(x, mu[k], s[k]))
    ll <- sum(log(rowSums(dn))); r <- dn / rowSums(dn)
    w <- colMeans(r); mu <- colSums(r * x) / colSums(r)
    s <- sqrt(colSums(r * outer(x, mu, "-")^2) / colSums(r))
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(r = r, mu = mu, s = s, w = w, iter = it)
}

ok <- !is.na(dt$coi) & dt$coi > 0
fit_em <- em1(log(dt$coi[ok]))
pure_k <- which.max(fit_em$mu)
dt[ok, purelike_prob := fit_em$r[, pure_k]]
dt[, genetic_type := fifelse(purelike_prob >= 0.5, "purelike", "mixlike")]
cross <- dt[abs(purelike_prob - 0.5) < 0.05, round(range(coi), 3)]
cat(sprintf("   log(COI) 혼합분포(EM %d회): 순종형 비율 %.2f, 성분 평균 COI %.3f / %.3f, 경계(사후확률 0.5) COI ≈ %.3f~%.3f\n",
            fit_em$iter, fit_em$w[pure_k], exp(fit_em$mu[pure_k]), exp(fit_em$mu[-pure_k]), cross[1], cross[2]))

anc_tab <- dt[, .(N = .N, top_pct_median = round(median(top_pct, na.rm = TRUE), 3),
                  coi_median = round(median(coi, na.rm = TRUE), 3),
                  n_pop_05_median = as.numeric(median(n_pop_05, na.rm = TRUE)),
                  purelike = sum(genetic_type == "purelike", na.rm = TRUE)), by = genetic_class][order(genetic_class)]
print(anc_tab)
#        genetic_class     N top_pct_median coi_median n_pop_05_median purelike
#               <char> <int>          <num>      <num>           <num>    <int>
# 1:             mixed  3642          0.410      0.090               3      591
# 2:    pure_confirmed  3105          0.999      0.244               1     2836
# 3:       pure_likely   451          0.687      0.142               2      206
# 4:     pure_mismatch    29          0.430      0.142               4       14
# 5: pure_unverifiable   400          0.251      0.210               3      300

cat("   보고 상태 × 유전 유형:\n")
print(dcast(dt[, .N, by = .(Breed_Status, genetic_type)], Breed_Status ~ genetic_type, value.var = "N", fill = 0))
#    Breed_Status mixlike purelike
#          <char>   <int>    <int>
# 1:  Mixed breed    3051      591
# 2:     Purebred     629     3356

cat("   보고 품종과 최대 품종이 다른 예 (pure_mismatch):\n")
print(dt[genetic_class == "pure_mismatch", .N, by = .(Breed, top_pop)][order(-N)][1:min(10, .N)])
#                          Breed                   top_pop     N
#                         <char>                    <char> <int>
#  1: Staffordshire Bull Terrier american_pit_bull_terrier     3
#  2:         Labrador Retriever american_pit_bull_terrier     2
#  3:            Norwich Terrier             cairn_terrier     2
#  4:      Australian Cattle Dog             border_collie     2
#  5:                    Whippet                 greyhound     1
#  6:                    Mastiff                 chihuahua     1
#  7:         Labrador Retriever          golden_retriever     1
#  8:              Border Collie       australian_shepherd     1
#  9:              Border Collie     australian_cattle_dog     1
# 10:           Belgian Malinois       german_shepherd_dog     1

cat("\n4) CSLB\n")
cs <- fread(f_cslb)[dog_id %in% dt$IID]
cat(sprintf("   genotype 개체의 CSLB 기록 %d개 (개 %d마리), 날짜 없는 기록 %d개 제외\n",
            nrow(cs), uniqueN(cs$dog_id), cs[is.na(cslb_date), .N]))
cs <- cs[!is.na(cslb_date)]
cs <- merge(cs, dt[, .(dog_id = IID, Estimated_DOB)], by = "dog_id")
cs[, `:=`(cslb_date = as.Date(cslb_date), age_cslb = as.numeric(as.Date(cslb_date) - Estimated_DOB) / 365.25)]
setorder(cs, dog_id, cslb_date)
head(cs)
fwrite(cs, paste0(out_prefix, "_cslb_long.csv"))

slope_of <- function(score, age) if (length(score) >= 2 && diff(range(age)) >= 1) unname(coef(lm(score ~ age))[2]) else NA_real_
cs_dog <- cs[, .(n_cslb = .N, cslb_date_latest = cslb_date[.N], cslb_year_latest = cslb_year[.N],
                 cslb_score = cslb_score[.N], age_cslb = age_cslb[.N],
                 cslb_score_first = cslb_score[1], age_cslb_first = age_cslb[1],
                 cslb_score_mean = mean(cslb_score), age_cslb_mean = mean(age_cslb),
                 cslb_score_max = max(cslb_score),
                 cslb_slope = slope_of(cslb_score, age_cslb), cslb_span = diff(range(age_cslb))), by = dog_id]
dt <- merge(dt, cs_dog, by.x = "IID", by.y = "dog_id", all.x = TRUE)
print(table(n_cslb = dt$n_cslb, useNA = "ifany"))
score_min <- min(cs$cslb_score)
cat(sprintf("   CSLB 최소 점수(관측) %d, 가장 최근 점수가 최소값인 개 %.1f%%\n", score_min,
            100 * dt[!is.na(cslb_score), mean(cslb_score == score_min)]))

cat("\n5) 공변량\n")
hles <- fread(f_hles, select = c("dog_id", "st_hles_completion_date", "dd_weight_lbs", "pa_activity_level",
                                 "pa_avg_daily_active_hours", "pa_avg_daily_active_minutes", "pa_avg_activity_intensity"))
hles[, Total_Active_Minutes := pa_avg_daily_active_hours * 60 + pa_avg_daily_active_minutes]
pa_in <- hles[complete.cases(hles[, .(pa_activity_level, pa_avg_activity_intensity, Total_Active_Minutes)])]
pa_pca <- prcomp(pa_in[, .(pa_activity_level, pa_avg_activity_intensity, Total_Active_Minutes)], center = TRUE, scale. = TRUE)

pa_sign <- sign(pa_pca$rotation["pa_activity_level", 1])
pa_in[, pa_score := pa_sign * pa_pca$x[, 1]]

cat(sprintf("   활동량 PCA: HLES %d마리, PC1 설명 비율 %.1f%%, loading(방향 고정 후) %s\n", nrow(pa_in), 100 * summary(pa_pca)$importance[2, 1],
            paste(sprintf("%s %.3f", rownames(pa_pca$rotation), pa_sign * pa_pca$rotation[, 1]), collapse = ", ")))

dt <- merge(dt, pa_in[, .(IID = dog_id, pa_score)], by = "IID", all.x = TRUE)

afus <- fread(f_afus, select = c("dog_id", "afus_complete_date", "afus_dd_dog_weight_lbs"))
wt_all <- rbind(hles[, .(IID = dog_id, wt_src = "HLES", wt_date = as.Date(st_hles_completion_date), wt_kg = dd_weight_lbs * lbs2kg)],
                afus[, .(IID = dog_id, wt_src = "AFUS", wt_date = as.Date(afus_complete_date), wt_kg = afus_dd_dog_weight_lbs * lbs2kg)])
wt_all <- merge(wt_all[!is.na(wt_kg) & !is.na(wt_date) & wt_kg > 0], dt[!is.na(cslb_date_latest), .(IID, cslb_date_latest)], by = "IID")
wt_all[, wt_gap := abs(as.numeric(wt_date - cslb_date_latest)) / 365.25]
wt_sel <- wt_all[order(IID, wt_gap)][, .SD[1], by = IID][, .(IID, wt_cov = wt_kg, wt_src, wt_gap)]
dt <- merge(dt, wt_sel, by = "IID", all.x = TRUE)
print(dt[!is.na(cslb_score), .(N = .N, pa_missing = sum(is.na(pa_score)), wt_missing = sum(is.na(wt_cov)),
                               wt_gap_median = round(median(wt_gap, na.rm = TRUE), 2)), by = Breed_Status])

hs <- rbind(fread(f_hles_hs, select = c("dog_id", "hs_condition", "hs_diagnosis_year", "hs_diagnosis_month", "hs_condition_other_description"))[
              , .(IID = dog_id, code = hs_condition, yr = hs_diagnosis_year, mo = hs_diagnosis_month, other = hs_condition_other_description, src = "HLES")],
            fread(f_afus_hs, select = c("dog_id", "afus_hs_new_condition", "afus_hs_new_condition_diagnosis_year",
                                        "afus_hs_new_condition_diagnosis_month", "afus_hs_new_condition_other_description"))[
              , .(IID = dog_id, code = afus_hs_new_condition, yr = afus_hs_new_condition_diagnosis_year,
                  mo = afus_hs_new_condition_diagnosis_month, other = afus_hs_new_condition_other_description, src = "AFUS")])
recode_to <- c(dementia = dementia_codes, vestibular = vestibular_codes, seizure = seizure_codes, brain = -1L)
hs_re <- rbindlist(lapply(names(other_neuro_recode), function(k)
  hs[code == 1298 & grepl(other_neuro_recode[[k]], other, ignore.case = TRUE)][, code := recode_to[[k]]]))
cat("   기타 신경계(1298) 자유응답 재분류:", paste(sprintf("%s %d건", names(other_neuro_recode),
    sapply(names(other_neuro_recode), function(k) hs_re[code == recode_to[[k]], .N])), collapse = ", "), "\n")

hs <- rbind(hs, hs_re)
brain_codes <- c(brain_codes, -1L)
hs <- merge(hs, dt[!is.na(cslb_date_latest), .(IID, cslb_date_latest)], by = "IID")
hs[, `:=`(yr = suppressWarnings(as.integer(yr)), mo = suppressWarnings(as.integer(mo)))]
hs[, dx_date := as.Date(ifelse(is.na(yr), NA, sprintf("%d-%02d-15", yr, fifelse(is.na(mo), 6L, mo))))]
hs_ever <- copy(hs)
hs <- hs[is.na(dx_date) | dx_date <= cslb_date_latest]
dt[, `:=`(vision_cov = as.integer(IID %in% hs[code %in% vision_codes, IID]),
          vision_cov_narrow = as.integer(IID %in% hs[code %in% vision_codes_narrow, IID]),
          hearing_cov = as.integer(IID %in% hs[code %in% hearing_codes, IID]),
          seizure_cov = as.integer(IID %in% hs[code %in% seizure_codes, IID]),
          vestibular_cov = as.integer(IID %in% hs[code %in% vestibular_codes, IID]),
          dementia_dx = as.integer(IID %in% hs[code %in% dementia_codes, IID]),
          dementia_dx_ever = as.integer(IID %in% hs_ever[code %in% dementia_codes, IID]),
          oa_cov = as.integer(IID %in% hs[code %in% oa_codes, IID]),
          brain_dx = as.integer(IID %in% hs[code %in% brain_codes, IID]))]
dt[, mimic_any := as.integer(hearing_cov + seizure_cov + vestibular_cov + brain_dx > 0)]
print(dt[!is.na(cslb_score), .(vision = sum(vision_cov), vision_narrow = sum(vision_cov_narrow), hearing = sum(hearing_cov),
                               seizure = sum(seizure_cov), vestibular = sum(vestibular_cov), brain = sum(brain_dx), oa = sum(oa_cov),
                               mimic_any = sum(mimic_any), dementia = sum(dementia_dx), dementia_ever = sum(dementia_dx_ever)), by = Breed_Status])

#    Breed_Status vision vision_narrow hearing seizure vestibular brain    oa
#          <char>  <int>         <int>   <int>   <int>      <int> <int> <int>
# 1:     Purebred    136           124      81      86         20     5   207
# 2:  Mixed breed    119           109      69      82         10     1   206
#    mimic_any dementia dementia_ever
#        <int>    <int>         <int>
# 1:       180       14            14
# 2:       150       27            27

afus_pre <- merge(fread(f_afus, select = c("dog_id", "afus_complete_date"))[!is.na(afus_complete_date)],
                  dt[!is.na(cslb_date_latest), .(dog_id = IID, cslb_date_latest)], by = "dog_id")[
  as.Date(afus_complete_date) <= cslb_date_latest, .(afus_n_pre = .N), by = .(IID = dog_id)]
dt <- merge(dt, afus_pre, by = "IID", all.x = TRUE)
dt[is.na(afus_n_pre), afus_n_pre := 0L]

#                                              ac_housesoiling, ac_memory, ac_anxiety, ac_inactivity
cat("\n5-2) CSLB 연결 변수\n")
nearest <- function(x, date_col) {
  x <- merge(x[!is.na(get(date_col))], dt[!is.na(cslb_date_latest), .(IID, cslb_date_latest)], by = "IID")
  x[, gap := abs(as.numeric(as.Date(get(date_col)) - cslb_date_latest)) / 365.25]
  x[order(IID, gap)][, .SD[1], by = IID][, cslb_date_latest := NULL]
}

ac_items <- c("none", "blind", "deaf", "weightloss", "mobility_weak", "mobility_pain", "other_pain", "cleanliness",
              "confusion", "interaction_change", "sleep", "housesoiling", "anxiety", "eat_drink", "inactivity",
              "repetitive_activity", "memory", "recognition", "sundowning", "other")
ac_ccd <- c("confusion", "interaction_change", "sleep", "housesoiling", "anxiety", "inactivity",
            "repetitive_activity", "memory", "recognition", "sundowning")
ac <- fread(f_afus, select = c("dog_id", "afus_complete_date", paste0("afus_ac_aging_char_", ac_items)))
setnames(ac, c("IID", "ac_date", paste0("ac_", ac_items)))
ac <- ac[rowSums(!is.na(ac[, paste0("ac_", ac_items), with = FALSE])) > 0]
ac <- nearest(ac, "ac_date")
for (v in paste0("ac_", ac_items)) ac[, (v) := as.integer(get(v) %in% TRUE)]
ac[, ac_ccd_count := rowSums(.SD), .SDcols = paste0("ac_", ac_ccd)]
setnames(ac, "gap", "ac_gap"); ac[, ac_date := NULL]
dt <- merge(dt, ac, by = "IID", all.x = TRUE)

sl <- rbind(fread(f_hles, select = c("dog_id", "st_hles_completion_date", "de_nighttime_sleep_avg_hours", "de_daytime_sleep_avg_hours"))[
              , .(IID = dog_id, sleep_src = "HLES", sleep_date = as.Date(st_hles_completion_date),
                  night = de_nighttime_sleep_avg_hours, day = de_daytime_sleep_avg_hours)],
            fread(f_afus, select = c("dog_id", "afus_complete_date", "afus_de_nighttime_sleep_avg_hours", "afus_de_daytime_sleep_avg_hours"))[
              , .(IID = dog_id, sleep_src = "AFUS", sleep_date = as.Date(afus_complete_date),
                  night = afus_de_nighttime_sleep_avg_hours, day = afus_de_daytime_sleep_avg_hours)])
sl[night == 99, night := NA]; sl[day == 99, day := NA]
sl <- sl[!is.na(night) & !is.na(day) & night + day <= 24]
sl <- nearest(sl, "sleep_date")
sl <- sl[, .(IID, sleep_night_hr = night, sleep_day_hr = day, sleep_total_hr = night + day,
             sleep_day_ratio = fifelse(night + day > 0, day / (night + day), NA_real_), sleep_src, sleep_gap = gap)]
dt <- merge(dt, sl, by = "IID", all.x = TRUE)

# (c) 123 Treat
trt <- fread(f_trt)
acc_cols <- grep("^trt_game_accuracy_round_", names(trt), value = TRUE)
resp_cols <- grep("^trt_game_response_round_", names(trt), value = TRUE)
trt[, `:=`(trt_n_valid = rowSums(!is.na(.SD[, acc_cols, with = FALSE])),
           trt_acc = rowMeans(.SD[, acc_cols, with = FALSE], na.rm = TRUE),
           trt_n_nosearch = rowSums(.SD[, resp_cols, with = FALSE] == 4, na.rm = TRUE),
           trt_quit = as.integer(rowSums(.SD[, resp_cols, with = FALSE] == 6, na.rm = TRUE) > 0))]
trt[is.nan(trt_acc), trt_acc := NA]
trt <- nearest(trt[, .(IID = dog_id, trt_date = trt_game_date, trt_n_valid, trt_acc, trt_n_nosearch, trt_quit,
                       trt_enthusiasm = trt_game_overall_response)], "trt_date")
setnames(trt, "gap", "trt_gap"); trt[, trt_date := NULL]
dt <- merge(dt, trt, by = "IID", all.x = TRUE)
print(dt[!is.na(cslb_score), .(N = .N, AFUS_ac = sum(!is.na(ac_none)), ac_gap_median = round(median(ac_gap, na.rm = TRUE), 2),
                               sleep = sum(!is.na(sleep_night_hr)), sleep_gap_median = round(median(sleep_gap, na.rm = TRUE), 2),
                               trt = sum(!is.na(trt_acc)), trt_gap_median = round(median(trt_gap, na.rm = TRUE), 2))])

d <- dt[!is.na(cslb_score)]

for (v in paste0("ac_", ac_items)) { cat("\n", v, "\n"); print(table(d[[v]], useNA = "ifany")) }
table(ac_ccd_count = d$ac_ccd_count, useNA = "ifany")
#    0    1    2    3    4    5    6    7    8    9 <NA> 
# 4258  851  234   98   40   28   15    7    5    2 1628 

table(ac_gap_yr = cut(d$ac_gap, c(0, 0.25, 0.5, 1, 2, Inf), include.lowest = TRUE), useNA = "ifany")
  # [0,0.25] (0.25,0.5]    (0.5,1]      (1,2]    (2,Inf]       <NA> 
  #     3381        711       1111        320         15       1628 

table(sleep_src = d$sleep_src, useNA = "ifany")
table(night = d$sleep_night_hr, useNA = "ifany")
table(day = d$sleep_day_hr, useNA = "ifany")
table(total = d$sleep_total_hr, useNA = "ifany")
table(day_ratio = cut(d$sleep_day_ratio, seq(0, 1, 0.1), include.lowest = TRUE), useNA = "ifany")
table(sleep_gap_yr = cut(d$sleep_gap, c(0, 0.25, 0.5, 1, 2, Inf), include.lowest = TRUE), useNA = "ifany")

# (c) 123 Treat
table(n_valid = d$trt_n_valid, useNA = "ifany")
table(acc = round(d$trt_acc, 2), useNA = "ifany")
table(n_nosearch = d$trt_n_nosearch, useNA = "ifany")
table(quit = d$trt_quit, useNA = "ifany")
table(enthusiasm = d$trt_enthusiasm, useNA = "ifany")
table(trt_gap_yr = cut(d$trt_gap, c(0, 0.25, 0.5, 1, 2, Inf), include.lowest = TRUE), useNA = "ifany")

cat("\n5-3) 코호트·나이 정확도·보호자/가구 특성\n")
hl <- fread(f_hles, select = c("dog_id", "od_age_range_years", "od_max_education", "od_annual_income_range_usd",
                               "ss_household_dog_count", "oc_household_person_count", "de_home_area_type", "hs_general_health",
                               "ss_vet_frequency"))
hl <- hl[, .(IID = dog_id, owner_age = od_age_range_years, owner_edu = od_max_education,
             owner_income = fifelse(od_annual_income_range_usd == 98, NA_integer_, od_annual_income_range_usd),
             household_dogs = ss_household_dog_count, household_people = oc_household_person_count,
             home_area = c("urban", "suburban", "rural")[de_home_area_type], general_health = hs_general_health,
             vet_freq = ss_vet_frequency)]
dt <- merge(dt, hl, by = "IID", all.x = TRUE)
dt[, `:=`(cohort_bhs = as.integer(Cohort == "Precision - Brain Health Study"),
          dob_exact = as.integer(AgeDOB_Estimation_Method == "Estimated using owner-reported birth year and month"))]
print(dt[!is.na(cslb_score), .(N = .N, bhs = sum(cohort_bhs), dob_exact = sum(dob_exact), owner_age_NA = sum(is.na(owner_age)),
                               income_NA = sum(is.na(owner_income)), afus_pre_0 = sum(afus_n_pre == 0))])

#        N   bhs dob_exact owner_age_NA income_NA afus_pre_0
#    <int> <int>     <int>        <int>     <int>      <int>
# 1:  7166   152      5780            0       747       1914

cat("\n5-4) 혈액·소변 검사\n")
f_sr <- file.path(base_dir, "SamplesResults/v2023")
meta_s <- fread(file.path(f_sr, "DAP_2023_SamplesResults_Metadata_v1.0.csv"),
                select = c("DAP_Sample_ID", "Sample_Collection_DateTime", "Sample_Dog_Hours_Fasted"))
head(meta_s)

meta_s[, `:=`(bio_date = as.Date(substr(Sample_Collection_DateTime, 1, 10)), fasted = suppressWarnings(as.numeric(Sample_Dog_Hours_Fasted)))]
read_panel <- function(file, cols, newnames, tag) {
  x <- fread(file.path(f_sr, file))
  x <- merge(x[, c("dog_id", "DAP_Sample_ID", cols), with = FALSE], meta_s[, .(DAP_Sample_ID, bio_date, fasted)], by = "DAP_Sample_ID")
  setnames(x, c("dog_id", cols), c("IID", newnames))
  x <- nearest(x, "bio_date")
  setnames(x, c("gap", "fasted"), c(paste0("bio_", tag, "_gap"), paste0("bio_", tag, "_fasted")))
  x[, c("IID", newnames, paste0("bio_", tag, c("_gap", "_fasted"))), with = FALSE]
}
cp_cols <- grep("^krt_cp_.*_value$", names(fread(file.path(f_sr, "DAP_2023_SamplesResults_ChemistryPanel_v1.0.csv"), nrows = 0)), value = TRUE)
cp <- read_panel("DAP_2023_SamplesResults_ChemistryPanel_v1.0.csv", cp_cols, sub("^krt_cp_(.*)_value$", "bio_cp_\\1", cp_cols), "cp")
setnames(cp, "bio_cp_fasted", "bio_fasted_hr")
cbc_cols <- c("krt_cbc_total_wbcs", "krt_cbc_abs_neutrophils", "krt_cbc_abs_lymphocytes", "krt_cbc_abs_monocytes", "krt_cbc_abs_eosinophils",
              "krt_cbc_rbc", "krt_cbc_hgb", "krt_cbc_hct", "krt_cbc_mcv", "krt_cbc_mch", "krt_cbc_mchc", "krt_cbc_rdw", "krt_cbc_plt", "krt_cbc_mpv")
cbc <- read_panel("DAP_2023_SamplesResults_CBC_v1.0.csv", cbc_cols, sub("^krt_cbc_", "bio_cbc_", cbc_cols), "cbc")[, bio_cbc_fasted := NULL]
cbc[, bio_cbc_nlr := fifelse(bio_cbc_abs_lymphocytes > 0, bio_cbc_abs_neutrophils / bio_cbc_abs_lymphocytes, NA_real_)]
cbc_cols_out <- c(sub("^krt_cbc_", "bio_cbc_", cbc_cols), "bio_cbc_nlr")
ur <- fread(file.path(f_sr, "DAP_2023_SamplesResults_Urine_v1.0.csv"), select = c("dog_id", "DAP_Sample_ID", "krt_urine_sg", "krt_urine_ph", "krt_urine_protien"))
ur <- merge(ur, meta_s[, .(DAP_Sample_ID, bio_date)], by = "DAP_Sample_ID")
ur <- ur[, .(IID = dog_id, bio_date, bio_ur_sg = suppressWarnings(as.numeric(sub("^>", "", krt_urine_sg))),
             bio_ur_ph = suppressWarnings(as.numeric(krt_urine_ph)),
             bio_ur_protein = c(Negative = 0, Trace = 0.5, `1+` = 1, `2+` = 2, `3+` = 3, `4+` = 4)[krt_urine_protien])]
ur <- nearest(ur, "bio_date"); setnames(ur, "gap", "bio_ur_gap"); ur[, bio_date := NULL]
dt <- merge(dt, cp, by = "IID", all.x = TRUE)
dt <- merge(dt, cbc, by = "IID", all.x = TRUE)
dt <- merge(dt, ur, by = "IID", all.x = TRUE)
bio_vars <- c(sub("^krt_cp_(.*)_value$", "bio_cp_\\1", cp_cols), cbc_cols_out, "bio_ur_sg", "bio_ur_ph", "bio_ur_protein")
print(dt[!is.na(cslb_score), .(N = .N, chem = sum(!is.na(bio_cp_albumin)), cbc = sum(!is.na(bio_cbc_hct)), urine = sum(!is.na(bio_ur_sg)),
                               chem_age8 = sum(!is.na(bio_cp_albumin) & age_cslb >= 8),
                               chem_gap_median = round(median(bio_cp_gap, na.rm = TRUE), 2))])

cat("\n6) 나이\n")
cc <- dt[!is.na(cslb_score) & !is.na(age_cslb) & !is.na(pa_score) & !is.na(wt_cov) & !is.na(Intact_cov) &
         !is.na(owner_age) & !is.na(cohort_bhs)]

age_grp_tab <- cc[, .(N = .N, pure = sum(Breed_Status == "Purebred"), mixed = sum(Breed_Status != "Purebred"),
                      mean = round(mean(cslb_score), 1), sd = round(sd(cslb_score), 1),
                      at_floor_pct = round(100 * mean(cslb_score == score_min), 1), ge50 = sum(cslb_score >= 50)),
                  by = .(age_group = cut(age_cslb, c(0, 2, 5, 8, 10, 12, 14, Inf), right = FALSE))][order(age_group)]
print(age_grp_tab)

floor_tab <- rbindlist(lapply(age_floor_cand, function(a) cc[age_cslb >= a, .(age_min = a, N = .N,
  pure = sum(Breed_Status == "Purebred"), mixed = sum(Breed_Status != "Purebred"),
  sd = round(sd(cslb_score), 2), ge50 = sum(cslb_score >= 50))]))
cat("   나이 하한 후보별 표본 (하한 이상 개만)\n"); print(floor_tab)

hinge <- function(x, k) pmax(0, x - k)
knot_search <- rbindlist(lapply(seq(7, 13, 0.5), function(k)
  data.table(knot = k, AIC = AIC(lm(cslb_score ~ age_cslb + hinge(age_cslb, k), data = cc)))))
age_knot <- knot_search[which.min(AIC), knot]
dt[, age_hinge := hinge(age_cslb, age_knot)]
cc[, age_hinge := hinge(age_cslb, age_knot)]
cat(sprintf("   age_hinge knot = %.1f세 (AIC 기준, v2022 순종은 11세)\n", age_knot))

knot_adj <- seq(7, 13, 0.5)[which.min(sapply(seq(7, 13, 0.5), function(k)
  AIC(lm(cslb_score ~ age_cslb + hinge(age_cslb, k) + Sex + Intact_cov + pa_score + wt_cov + vision_cov + owner_age + cohort_bhs, data = cc))))]
cat(sprintf("   주 모델 공변량 보정 시 knot = %.1f세 (비보정 %.1f세와 같으면 비보정 결과 사용)\n", knot_adj, age_knot))

young <- cc[age_cslb < 5, cslb_score]
age_year_tab <- cc[, .(n = .N, mean = round(mean(cslb_score), 2), sd = round(sd(cslb_score), 2),
                       ge50_pct = round(100 * mean(cslb_score >= 50), 1), diff_vs_lt5 = round(mean(cslb_score) - mean(young), 2),
                       p_vs_lt5 = signif(t.test(cslb_score, young)$p.value, 2)), keyby = .(age = pmin(floor(age_cslb), 16))]
cat("   나이(년)별 CSLB (16 = 16세 이상):\n"); print(age_year_tab)
set.seed(1)
knot_boot <- replicate(200, { s <- cc[sample(.N, replace = TRUE)]; ks <- seq(7, 13, 0.5)
  ks[which.min(sapply(ks, function(k) AIC(lm(cslb_score ~ age_cslb + hinge(age_cslb, k), data = s))))] })
cat(sprintf("   knot 부트스트랩(200회) 95%% 구간: %.1f ~ %.1f세\n", quantile(knot_boot, .025), quantile(knot_boot, .975)))
fit <- lm(cslb_score ~ age_cslb + age_hinge + Sex + Intact_cov + pa_score + wt_cov + vision_cov, data = cc)
print(round(coef(summary(fit)), 4))

rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
core <- "age_cslb + age_hinge + Sex + Intact_cov + pa_score"
cand <- c("wt_cov", "vision_cov", "vision_cov_narrow", "hearing_cov", "seizure_cov", "vestibular_cov", "coi", "top_pct", "purelike_prob", "factor(cslb_year_latest)")
cov_tab <- rbindlist(lapply(list(all = cc, age8 = cc[age_cslb >= 8]), function(d) {
  d <- copy(d)[, y_rint := rint(cslb_score)]
  rbindlist(lapply(cand, function(v) rbindlist(lapply(c("cslb_score", "y_rint"), function(y) {
    f0 <- lm(as.formula(paste(y, "~", core)), data = d); f1 <- lm(as.formula(paste(y, "~", core, "+", v)), data = d)
    a <- anova(f0, f1)
    data.table(candidate = v, outcome = y, beta = if (grepl("factor", v)) NA_real_ else round(coef(f1)[[v]], 3),
               p = signif(a$`Pr(>F)`[2], 3), n_pos = if (all(d[[sub("factor\\((.*)\\)", "\\1", v)]] %in% 0:1)) sum(d[[v]]) else NA_integer_)
  }))))
}), idcol = "set")
cat("\n   공변량 후보 효과 (핵심 공변량 보정 후):\n")
print(dcast(cov_tab, candidate ~ set + outcome, value.var = "p"))

r0 <- resid(lm(as.formula(paste("cslb_score ~", core)), data = cc))
cat(sprintf("   CSLB 분포: 원점수 왜도 %.2f, 첨도 %.1f / 보정 후 잔차 왜도 %.2f, 첨도 %.1f (정규분포는 0, 3)\n",
            mean(scale(cc$cslb_score)^3), mean(scale(cc$cslb_score)^4), mean(scale(r0)^3), mean(scale(r0)^4)))

main_qcov <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age")
main_dcov <- c("Sex", "Intact_cov", "vision_cov", "cohort_bhs")
main_form <- paste(c(main_qcov, main_dcov), collapse = " + ")
main_tab <- rbindlist(lapply(list(all = cc, age8 = cc[age_cslb >= 8]), function(d) {
  d <- copy(d)[, y_rint := rint(cslb_score)]
  rbindlist(lapply(c("cslb_score", "y_rint"), function(y) {
    cf <- coef(summary(lm(as.formula(paste(y, "~", main_form)), data = d)))[-1, , drop = FALSE]
    data.table(outcome = y, term = rownames(cf), beta = signif(cf[, 1], 3), p = signif(cf[, 4], 3))
  }))
}), idcol = "set")
cat("\n   주 모델 공변량 효과 (p):\n")
print(dcast(main_tab, term ~ set + outcome, value.var = "p"))

floor_sens <- rbindlist(lapply(age_floor_cand, function(a) {
  d <- copy(cc[age_cslb >= a])[, y_rint := rint(cslb_score)]
  f <- main_form
  m <- lm(as.formula(paste("cslb_score ~", f)), data = d); mr <- lm(as.formula(paste("y_rint ~", f)), data = d)
  cf <- coef(m)
  data.table(age_min = a, N = nrow(d), ge50 = sum(d$cslb_score >= 50), ge50_pct = round(100 * mean(d$cslb_score >= 50), 2),
             mean = round(mean(d$cslb_score), 2), sd = round(sd(d$cslb_score), 2),
             r2_raw = round(summary(m)$r.squared, 3), r2_rint = round(summary(mr)$r.squared, 3),
             b_vision = round(cf[["vision_cov"]], 2), b_owner_age = round(cf[["owner_age"]], 3), b_pa = round(cf[["pa_score"]], 3),
             b_bhs = round(cf[["cohort_bhs"]], 2), b_age_hinge = round(cf[["age_hinge"]], 2))
}))
ref <- floor_sens[age_min == age_floor_main]
floor_sens[, `:=`(N_vs8 = round(N / ref$N, 2), ge50_vs8 = round(ge50 / ref$ge50, 2),
                  b_vision_vs8 = round(b_vision / ref$b_vision, 2), b_owner_age_vs8 = round(b_owner_age / ref$b_owner_age, 2))]
cat(sprintf("\n   나이 하한 민감도 (기준 %d세, *_vs8 = 기준 대비 비율):\n", age_floor_main)); print(floor_sens)
floor_tab <- floor_sens
for (a in setdiff(age_floor_cand, c(0, 8, 10)))
  fwrite(cc[age_cslb >= a, .(FID, IID)], sprintf("%s_keep_agefloor%d.txt", out_prefix, a), sep = " ", col.names = FALSE)

link_role <- c(setNames(rep("교란 후보(감각)", 2), paste0("ac_", c("blind", "deaf"))),
               setNames(rep("교란 후보(통증·신체)", 6), paste0("ac_", c("mobility_weak", "mobility_pain", "other_pain", "weightloss", "eat_drink", "cleanliness"))),
               setNames(rep("검증(인지)", 8), paste0("ac_", c("confusion", "repetitive_activity", "recognition", "interaction_change", "housesoiling", "memory", "anxiety", "inactivity"))),
               setNames(rep("수면·각성", 2), paste0("ac_", c("sleep", "sundowning"))),
               ac_none = "참고(증상 없음)", ac_other = "참고(기타)", ac_ccd_count = "검증(인지·수면 합계)",
               sleep_night_hr = "수면·각성", sleep_day_hr = "수면·각성", sleep_total_hr = "수면·각성", sleep_day_ratio = "수면·각성",
               trt_acc = "검증(객관적 인지 과제)", trt_n_nosearch = "검증(객관적 인지 과제)", trt_quit = "참고(과제 참여)",
               trt_enthusiasm = "참고(과제 참여)",
               hearing_cov = "교란 후보(진단 이력)", seizure_cov = "교란 후보(진단 이력)", vestibular_cov = "교란 후보(진단 이력)",
               dementia_dx = "검증(CCD 진단 이력)", dementia_dx_ever = "검증(CCD 진단 이력)",
               oa_cov = "교란 후보(진단 이력)", brain_dx = "교란 후보(진단 이력)", mimic_any = "교란 후보(진단 이력)",
               general_health = "민감도(건강 평가)", household_dogs = "민감도(가구)", household_people = "민감도(가구)",
               home_area = "민감도(가구)", owner_income = "민감도(보호자)", owner_edu = "민감도(보호자)",
               vet_freq = "참고(진단 기회)", afus_n_pre = "참고(진단 기회)", dob_exact = "참고(나이 정확도)")
link_tab <- rbindlist(lapply(list(all = cc, age8 = cc[age_cslb >= 8]), function(d) {
  d <- copy(d)[, y_rint := rint(cslb_score)]
  rbindlist(lapply(names(link_role), function(v) rbindlist(lapply(c("cslb_score", "y_rint"), function(y) {
    dd <- d[!is.na(get(v))]
    m0 <- lm(as.formula(paste(y, "~", main_form)), data = dd)
    m <- lm(as.formula(paste(y, "~", main_form, "+", v)), data = dd)
    cf <- coef(summary(m))
    is_cat <- is.character(dd[[v]])
    data.table(variable = v, role = link_role[[v]], outcome = y, n = nrow(dd),
               n_pos = if (all(dd[[v]] %in% 0:1)) sum(dd[[v]]) else NA_integer_,
               beta = if (is_cat) NA_real_ else signif(cf[v, 1], 3),
               p = signif(if (is_cat) anova(m0, m)$`Pr(>F)`[2] else cf[v, 4], 3))
  }))))
}), idcol = "set")
cat("\n   CSLB 연결 변수 (주 모델 + 변수 하나, p):\n")
print(dcast(link_tab, role + variable ~ set + outcome, value.var = "p")[order(role)])

vis_tab <- rbindlist(lapply(list(all = cc, age8 = cc[age_cslb >= 8]), function(d) {
  d <- copy(d)[, y_rint := rint(cslb_score)]
  rbindlist(lapply(c("vision_cov_narrow", "vision_cov"), function(v) rbindlist(lapply(c("cslb_score", "y_rint"), function(y) {
    m <- lm(as.formula(paste(y, "~", sub("vision_cov", v, main_form, fixed = TRUE))), data = d)
    data.table(definition = v, outcome = y, n_pos = sum(d[[v]]), beta = signif(coef(m)[[v]], 3),
               p = signif(coef(summary(m))[v, 4], 3), r2 = round(summary(m)$r.squared, 4), AIC = round(AIC(m), 1))
  }))))
}), idcol = "set")
cat("\n   시력 정의 비교 (vision_cov_narrow = 101/102/115, vision_cov = + 103/118/119):\n")
print(vis_tab)
cat(sprintf("   확장으로 추가된 개체: %d마리 (그중 8세 이상 %d)\n",
            cc[vision_cov == 1 & vision_cov_narrow == 0, .N], cc[vision_cov == 1 & vision_cov_narrow == 0 & age_cslb >= 8, .N]))

bio_tab <- rbindlist(lapply(list(all = cc, age8 = cc[age_cslb >= 8]), function(d) {
  rbindlist(lapply(bio_vars, function(v) {
    dd <- d[!is.na(get(v))]
    if (nrow(dd) < 50) return(NULL)
    dd[, `:=`(x = rint(get(v)), y_rint = rint(cslb_score))]
    r <- lapply(c("cslb_score", "y_rint"), function(y) coef(summary(lm(as.formula(paste(y, "~", main_form, "+ x")), data = dd)))["x", c(1, 4)])
    data.table(variable = v, n = nrow(dd), beta_raw_per_sd = signif(r[[1]][1], 3), p_raw = signif(r[[1]][2], 3), p_rint = signif(r[[2]][2], 3))
  }))
}), idcol = "set")
bio_tab[, `:=`(fdr_raw = signif(p.adjust(p_raw, "BH"), 3), fdr_rint = signif(p.adjust(p_rint, "BH"), 3)), by = set]
cat("\n   바이오마커 (주 모델 + 검사값, 상위 10개):\n")
print(bio_tab[order(p_raw)][, head(.SD, 10), by = set])

cat("\n7) 분석 대상\n")
dt[, set_all := !is.na(cslb_score) & !is.na(age_cslb) & !is.na(pa_score) & !is.na(wt_cov) & !is.na(Intact_cov) &
                !is.na(owner_age) & !is.na(cohort_bhs)]
dt[, `:=`(set_age8 = set_all & age_cslb >= 8, set_age10 = set_all & age_cslb >= 10,
          set_purelike = set_all & genetic_type %in% "purelike", set_mixlike = set_all & genetic_type %in% "mixlike",
          set_nomimic = set_all & mimic_any == 0, set_dob_exact = set_all & dob_exact == 1,
          set_afus = set_all & afus_n_pre >= 1, set_ac = set_all & !is.na(ac_none), set_no_bhs = set_all & cohort_bhs == 0)]
sets <- c("set_all", "set_age8", "set_age10", "set_purelike", "set_mixlike",
          "set_nomimic", "set_dob_exact", "set_afus", "set_ac", "set_no_bhs")
set_tab <- rbindlist(lapply(sets, function(s) dt[get(s) == TRUE, .(set = s, N = .N, pure = sum(Breed_Status == "Purebred"),
  mixed = sum(Breed_Status != "Purebred"), age_median = round(median(age_cslb), 1), ge50 = sum(cslb_score >= 50))]))
print(set_tab)
for (s in sets) fwrite(dt[get(s) == TRUE, .(FID, IID)], sprintf("%s_keep_%s.txt", out_prefix, sub("^set_", "", s)),
                       sep = " ", col.names = FALSE)

gcta_cov <- dt[set_all == TRUE]
fwrite(gcta_cov[, c("FID", "IID", main_qcov), with = FALSE], paste0(out_prefix, "_qcovar.txt"), sep = " ", col.names = FALSE)
fwrite(gcta_cov[, .(FID, IID, Sex, Intact = fifelse(Intact_cov == 1, "neutered", "intact"),
                    Vision = fifelse(vision_cov == 1, "vision_problem", "normal"),
                    Cohort_BHS = fifelse(cohort_bhs == 1, "brain_health_study", "other"))],
       paste0(out_prefix, "_covar.txt"), sep = " ", col.names = FALSE)
cat(sprintf("   GCTA 공변량 파일: %d마리 (qcovar: %s / covar: Sex, Intact, Vision, Cohort_BHS)\n",
            nrow(gcta_cov), paste(main_qcov, collapse = ", ")))

fwrite(gcta_cov[, .(FID, IID, hearing_cov, seizure_cov, vestibular_cov, brain_dx, oa_cov, general_health,
                    household_dogs, owner_income, owner_edu, ac_blind, ac_deaf, ac_mobility_pain, ac_mobility_weak)],
       paste0(out_prefix, "_sens_covariates.txt"), sep = " ", na = "NA")

p1 <- ggplot(cc, aes(age_cslb, cslb_score, color = Breed_Status)) +
  geom_point(alpha = 0.15, size = 0.8) + geom_smooth(method = "loess", se = TRUE) +
  geom_vline(xintercept = age_knot, linetype = 2) + geom_hline(yintercept = 50, linetype = 3, color = "red") +
  scale_color_manual(values = c("#2A4F9E", "#FF9D00")) +
  labs(title = "CSLB score by age at latest CSLB (v2023)", subtitle = sprintf("dashed: hinge knot %.1f y; red dotted: CSLB = 50", age_knot),
       x = "Age at CSLB (years)", y = "CSLB score", color = NULL) + theme_bw(base_size = 12)
ggsave(file.path(pic_dir, "00_cslb_by_age.png"), p1, width = 8, height = 5, dpi = 300)

p2 <- ggplot(dt[!is.na(top_pct)], aes(top_pct, fill = Breed_Status)) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("#2A4F9E", "#FF9D00")) +
  labs(title = "Largest single-breed ancestry proportion", x = "Top ancestry proportion", y = "Dogs", fill = NULL) + theme_bw(base_size = 12)
ggsave(file.path(pic_dir, "00_ancestry_top_pct.png"), p2, width = 8, height = 4.5, dpi = 300)

p3 <- ggplot(dt[!is.na(own_pct)], aes(own_pct)) + geom_histogram(bins = 50, fill = "#2A4F9E") +
  geom_vline(xintercept = own_cut[c("discordant", "confirmed")], linetype = 2, color = "red") +
  labs(title = "Reported-breed ancestry proportion (purebreds with breed in reference panel)",
       x = "Ancestry proportion of reported breed", y = "Dogs") + theme_bw(base_size = 12)
ggsave(file.path(pic_dir, "00_ancestry_own_pct.png"), p3, width = 8, height = 4.5, dpi = 300)

p4 <- ggplot(dt[!is.na(purelike_prob)], aes(top_pct, coi, color = purelike_prob)) +
  geom_point(size = 0.7, alpha = 0.6) + facet_wrap(~ Breed_Status) +
  scale_color_gradient(low = "#FF9D00", high = "#2A4F9E", name = "P(purelike)") +
  labs(title = "Genetic breed structure: ancestry vs inbreeding", subtitle = "colour = posterior probability of the pure-like component (2-component mixture)",
       x = "Top single-breed ancestry proportion", y = "COI (ROH-based)") + theme_bw(base_size = 12)
ggsave(file.path(pic_dir, "00_genetic_type.png"), p4, width = 10, height = 4.8, dpi = 300)

dt[, IID_chr := NULL]
fwrite(dt, paste0(out_prefix, "_all_total.csv"))

col_desc <- rbindlist(lapply(list(
  c("IID", "식별", "DAP dog_id (= genotype IID)"),
  c("FID", "식별", "genotype FID (_ed.fam, = IID)"),
  c("Sex", "식별", "성별, DAP .tsv sex_meta(X 염색체 coverage 기반 유전 추정치) 기준. .fam과 다른 개체도 제외 없이 sex_meta로 정정. 주 모델 공변량"),
  c("Sex_fam", "식별", ".fam 성별 (대조용, 주 모델에는 쓰지 않음; 일부 개체가 sex_meta·보호자 응답과 다름)"),
  c("in_v2022", "식별", "v2022 genotype release에도 있는 개체인지 (TRUE/FALSE)"),
  c("Breed_Status", "보호자 보고", "보고 품종 상태: Purebred / Mixed breed (분석에는 쓰지 않고 기술용)"),
  c("Breed", "보호자 보고", "보고 품종명 (순종만, 잡종은 Mixed)"),
  c("Estimated_DOB", "보호자 보고", "추정 생년월일 (나이 계산 기준)"),
  c("Sex_Class_at_HLES", "보호자 보고", "HLES 가입 시 성별·중성화 상태 (예: Female, spayed)"),
  c("Breed_Size_Class_at_HLES", "보호자 보고", "HLES 체중 기반 크기 등급 x AKC 순종 여부 (10단계, DAP 정의)"),
  c("Status", "보호자 보고", "현재 상태: Alive / Dead"),
  c("Cohort", "보호자 보고", "DAP 코호트: Foundation / Precision / DAP Pack / Precision - Brain Health Study"),
  c("AgeDOB_Estimation_Method", "보호자 보고", "생년월일 추정 방법 (연·월 / 연도만 / 나이로 추정)"),
  c("Intact_cov", "공변량(주 모델)", "중성화 여부 (1 = 중성화, 0 = intact), HLES 기준"),
  c("sex_owner", "보호자 보고", "보호자 응답 성별 (Sex·Sex_fam과 대조용)"),
  c("sex_meta", "유전 지표", "DAP .tsv 성별 (X 염색체 coverage 비율로 추정, F/M)"),
  c("coi", "유전 지표", "근교계수: ROH(500 kb 이상)가 덮는 유전체 비율 (DAP .tsv). 높을수록 순종형"),
  c("size_score", "유전 지표", "유전체 기반 크기 예측 점수 0(tiny)~4(giant): 키 관련 SNP 2,733개 random forest (DAP .tsv)"),
  c("top_pop", "유전 지표", "ancestry(ADMIXTURE, 참조 110집단)에서 비율이 가장 큰 집단"),
  c("top_pct", "유전 지표", "top_pop의 비율 (0~1). 클수록 한 품종에 가까움"),
  c("second_pop", "유전 지표", "두 번째로 비율이 큰 집단"),
  c("second_pct", "유전 지표", "second_pop의 비율"),
  c("n_pop_05", "유전 지표", "비율 0.05 이상인 집단 수 (잡종일수록 큼)"),
  c("own_pct", "유전 지표", "보고 품종(참조 패널에 있는 경우)의 ancestry 비율, 크기 변종 합산. 잡종·패널 밖 품종은 NA"),
  c("top_is_own", "유전 지표", "최대 비율 집단이 보고 품종인지"),
  c("genetic_class", "유전 지표", "보고 품종과 ancestry 비교(기술용): pure_confirmed(≥0.85) / pure_likely / pure_mismatch / pure_unverifiable / mixed"),
  c("purelike_prob", "유전 지표", "log(COI) 2성분 혼합분포에서 순종형 성분일 사후확률"),
  c("admix_class", "유전 지표", "보고 품종과 무관한 유전 분류 (Morrill 2022): top_pct ≥ 0.85 single_breed / 0.45~0.85 intermediate / < 0.45 highly_admixed"),
  c("genetic_type", "유전 지표", "purelike_prob ≥ 0.5 → purelike, 아니면 mixlike (효과 이질성 확인용)"),
  # CSLB
  c("n_cslb", "CSLB", "CSLB 응답 횟수 (1~4)"),
  c("cslb_date_latest", "CSLB", "가장 최근 CSLB 날짜 (다른 변수의 날짜 매칭 기준)"),
  c("cslb_year_latest", "CSLB", "가장 최근 CSLB 조사 연도"),
  c("cslb_score", "CSLB(표현형)", "가장 최근 CSLB 점수 (16~80, ≥50 = CCD). 주 표현형 (분석 시 RINT)"),
  c("age_cslb", "공변량(주 모델)", "가장 최근 CSLB 시점 나이(년)"),
  c("cslb_score_first", "CSLB", "첫 CSLB 점수"),
  c("age_cslb_first", "CSLB", "첫 CSLB 시점 나이"),
  c("cslb_score_mean", "CSLB", "전체 CSLB 평균 점수 (민감도 표현형)"),
  c("age_cslb_mean", "CSLB", "CSLB 시점 평균 나이"),
  c("cslb_score_max", "CSLB", "CSLB 최대 점수"),
  c("cslb_slope", "CSLB", "CSLB 연간 변화율(점/년), 2회 이상 & 간격 1년 이상만 (종단 표현형 후보)"),
  c("cslb_span", "CSLB", "첫~마지막 CSLB 사이 기간(년)"),
  c("pa_score", "공변량(주 모델)", "활동량: HLES 활동 수준·강도·하루 활동 시간 PCA 1번째 축 (활동 수준 loading이 양수가 되도록 부호 고정 → 클수록 활동적)"),  c("wt_cov", "공변량(주 모델)", "체중(kg): HLES·AFUS 중 가장 최근 CSLB와 가장 가까운 기록"),
  c("wt_src", "공변량", "wt_cov 출처 (HLES / AFUS)"),
  c("wt_gap", "공변량", "wt_cov 기록과 가장 최근 CSLB의 날짜 간격(년)"),
  c("vision_cov", "공변량(주 모델)", "진단 이력 시력 문제 (보호자 보고, HLES·AFUS health condition 'diagnosed with' 문항): 실명 101, 백내장 102, 어린 시기 백내장 115, 녹내장 103, PRA 118, 망막 박리 119. 가장 최근 CSLB 이전 진단 (1 = 있음)"),
  c("vision_cov_narrow", "비교용", "이전 시력 정의: 실명·백내장(101/102/115)만, 같은 기준 (00_summary_vision_def.csv에서 vision_cov와 비교)"),
  c("hearing_cov", "공변량 후보", "진단 이력 청력 문제 (보호자 보고, 난청·부분 청력 손실, 코드 201/205), 같은 기준"),
  c("seizure_cov", "공변량 후보", "진단 이력 발작·뇌전증 (보호자 보고, 코드 1215), 같은 기준"),
  c("vestibular_cov", "공변량 후보(민감도)", "진단 이력 전정 질환 (보호자 보고, 코드 1216), 같은 기준. 빙빙 돌기·균형 상실·방향 혼란이 CSLB 서성임·혼란 문항과 겹침"),
  c("dementia_dx", "검증(CCD 진단 이력)", "진단 이력 치매·노쇠 (보호자 보고, 코드 1205), 가장 최근 CSLB 이전 진단. CCD 진단 자체라 공변량으로 쓰면 안 됨"),
  c("dementia_dx_ever", "검증(CCD 진단 이력)", "진단 이력 치매·노쇠, 날짜 제한 없음 (가장 최근 CSLB 이후 진단도 포함)"),
  c("oa_cov", "공변량 후보(민감도)", "진단 이력 골관절염 (보호자 보고, 코드 1113), 가장 최근 CSLB 이전. 통증 → 쓰다듬기 회피·활동 감소 문항과 겹침"),
  c("brain_dx", "공변량 후보(민감도)", "진단 이력 뇌 질환: 소뇌 형성 부전 1201, 수두증 1202, 기타 신경계(1298) 자유응답 중 뇌경색·TIA·NCL·수막염 등"),
  c("mimic_any", "민감도(제외)", "CCD와 증상이 비슷한 질환(hearing_cov, seizure_cov, vestibular_cov, brain_dx) 중 하나라도 있음 → set_nomimic에서 제외"),
  c("afus_n_pre", "참고(진단 기회)", "가장 최근 CSLB 이전 AFUS 응답 수 (0이면 HLES 이후 새 진단이 보고될 기회가 없었음)"),
  c("owner_age", "공변량(주 모델)", "보호자 나이 등급 1(18-24) 2(25-34) 3(35-44) 4(45-54) 5(55-64) 6(65-74) 7(75+), HLES"),
  c("owner_edu", "민감도(보호자)", "보호자 최종 학력 0(학교 안 다님)~11(박사·전문 학위), HLES"),
  c("owner_income", "민감도(보호자)", "가구 연소득 등급 1(<2만 달러)~10(18만 달러 이상), 응답 거부(98) → NA"),
  c("household_dogs", "민감도(가구)", "가구 내 개 수 (HLES)"),
  c("household_people", "민감도(가구)", "가구 내 사람 수 (HLES)"),
  c("home_area", "민감도(가구)", "거주 지역 urban / suburban / rural (HLES)"),
  c("general_health", "민감도(건강 평가)", "보호자의 전반적 건강 평가 1(매우 좋음)~6(매우 나쁨), HLES. CCD 때문에 낮게 평가될 수 있어 주 모델 제외"),
  c("vet_freq", "참고(진단 기회)", "동물병원 방문 빈도 0(안 감) 1(1년에 1번 미만) 2(1년에 1번) 3(1년에 1번 넘게), HLES"),
  c("cohort_bhs", "공변량(주 모델)", "Brain Health Study 코호트 (1 = 해당). 고령·인지 상태 기준 모집으로 CSLB가 높음"),
  c("dob_exact", "민감도(나이 정확도)", "생년월일을 연·월로 앎 (1 = 예) → set_dob_exact"),
  c("ac_none", "AFUS 노화 증상", "노화 증상 없음에 체크"),
  c("ac_blind", "AFUS 노화 증상 · 교란 후보(감각)", "실명 또는 시력 저하 (보호자 관찰, 진단 아님)"),
  c("ac_deaf", "AFUS 노화 증상 · 교란 후보(감각)", "난청 또는 청력 저하 (보호자 관찰)"),
  c("ac_weightloss", "AFUS 노화 증상 · 교란 후보(신체)", "체중 감소"),
  c("ac_mobility_weak", "AFUS 노화 증상 · 교란 후보(신체)", "다리·관절 약화나 주저앉음으로 이동 능력 저하"),
  c("ac_mobility_pain", "AFUS 노화 증상 · 교란 후보(통증)", "다리·관절 통증이나 뻣뻣함으로 이동 능력 저하"),
  c("ac_other_pain", "AFUS 노화 증상 · 교란 후보(통증)", "다른 부위 통증"),
  c("ac_cleanliness", "AFUS 노화 증상 · 교란 후보(신체)", "청결 유지 능력 저하 (소변 냄새, 털 상태 나쁨)"),
  c("ac_confusion", "AFUS 노화 증상 · 검증(인지)", "방향 감각 상실·혼란 (멍하니 벽 보기, 물건 뒤에 끼임) ↔ CSLB stare, stuck"),
  c("ac_interaction_change", "AFUS 노화 증상 · 검증(인지)", "사람·동물과의 상호작용 변화 (사회적 상호작용 감소, 관심 요구 증가, 공격성) ↔ CSLB avoid"),
  c("ac_sleep", "AFUS 노화 증상 · 수면·각성", "수면 장애 (낮에 자고 밤에 안절부절). CSLB에 없는 영역"),
  c("ac_housesoiling", "AFUS 노화 증상 · 검증(인지)", "실내 배변 실수 (전에는 깨끗하던 곳) ↔ CSLB defecate_6mo"),
  c("ac_anxiety", "AFUS 노화 증상 · 검증(인지)", "사람·장소·상황에 대한 불안·두려움"),
  c("ac_eat_drink", "AFUS 노화 증상 · 교란 후보(신체)", "먹고 마시는 것에 대한 관심 감소"),
  c("ac_inactivity", "AFUS 노화 증상 · 검증(인지)", "활동 시간 감소 ↔ CSLB active_6mo (활동량 공변량과 겹침)"),
  c("ac_repetitive_activity", "AFUS 노화 증상 · 검증(인지)", "반복 행동 증가 (서성임, 원 그리며 걷기, 목적 없이 배회) ↔ CSLB pace"),
  c("ac_memory", "AFUS 노화 증상 · 검증(인지)", "학습 곤란·기억 상실 (새 재주를 못 배우거나 배운 것을 못 함). CSLB에 없는 영역"),
  c("ac_recognition", "AFUS 노화 증상 · 검증(인지)", "아는 사람·동물을 못 알아봄 ↔ CSLB recognize"),
  c("ac_sundowning", "AFUS 노화 증상 · 수면·각성", "늦은 오후·초저녁에 노화 증상이 나타나거나 심해짐 (sundowning). 일주기 관련"),
  c("ac_other", "AFUS 노화 증상", "기타 노화 증상"),
  c("ac_gap", "AFUS 노화 증상", "ac_* 를 가져온 AFUS와 가장 최근 CSLB의 날짜 간격(년). AFUS 1년차엔 문항이 없어 NA가 많음"),
  c("ac_ccd_count", "AFUS 노화 증상 · 검증(합계)", "인지·수면 영역 10문항 체크 수 (confusion, interaction_change, sleep, housesoiling, anxiety, inactivity, repetitive_activity, memory, recognition, sundowning)"),
  c("sleep_night_hr", "수면·각성", "평균 밤잠 시간(시간), HLES·AFUS 중 가장 최근 CSLB와 가장 가까운 기록 (99 = 모름 → NA)"),
  c("sleep_day_hr", "수면·각성", "평균 낮잠 시간(시간), 같은 기록"),
  c("sleep_total_hr", "수면·각성", "밤잠 + 낮잠 (24시간 초과 기록은 제외)"),
  c("sleep_day_ratio", "수면·각성", "낮잠 / (밤잠 + 낮잠): 낮에 자는 비율 (일주기 지표 후보)"),
  c("sleep_src", "수면·각성", "수면 시간 출처 (HLES / AFUS)"),
  c("sleep_gap", "수면·각성", "수면 시간 기록과 가장 최근 CSLB의 날짜 간격(년)"),
  # 123 Treat
  c("trt_n_valid", "123 Treat 인지 게임", "9회 중 상자를 찾은(정답 판정이 있는) 횟수"),
  c("trt_acc", "123 Treat 인지 게임 · 검증(객관적)", "정답률 = 맞는 상자를 찾은 횟수 / trt_n_valid (우연 수준 0.33). trt_n_valid = 0이면 NA"),
  c("trt_n_nosearch", "123 Treat 인지 게임", "아무 상자도 찾지 않은 횟수 (응답 코드 4)"),
  c("trt_quit", "123 Treat 인지 게임", "싫어해서 중단한 적 있음 (응답 코드 6, 1 = 있음)"),
  c("trt_enthusiasm", "123 Treat 인지 게임", "보호자 평가 흥미도 1(전혀 없음)~3(매우 흥분)"),
  c("trt_gap", "123 Treat 인지 게임", "게임 날짜와 가장 최근 CSLB의 날짜 간격(년)"),
  c("age_hinge", "공변량(주 모델)", "max(0, age_cslb - knot), knot = AIC로 정한 값 (11세): 고령에서 기울기 변화"),
  c("set_all", "분석 대상", "주 분석 대상 (CSLB·주 모델 공변량·유전 지표 완비, 제외 없음)"),
  c("set_age8", "분석 대상", "set_all 중 가장 최근 CSLB 시점 8세 이상"),
  c("set_age10", "분석 대상", "set_all 중 10세 이상"),
  c("set_purelike", "분석 대상", "set_all 중 genetic_type = purelike"),
  c("set_mixlike", "분석 대상", "set_all 중 genetic_type = mixlike"),
  c("set_nomimic", "분석 대상(민감도)", "set_all 중 mimic_any = 0 (청력·발작·전정·뇌 질환 진단 이력 없음)"),
  c("set_dob_exact", "분석 대상(민감도)", "set_all 중 dob_exact = 1 (생년월일 연·월 앎)"),
  c("set_afus", "분석 대상(민감도)", "set_all 중 afus_n_pre ≥ 1 (CSLB 전 AFUS 1회 이상)"),
  c("set_ac", "분석 대상(민감도)", "set_all 중 AFUS 노화 증상 문항 응답 있음 (ac_* 공변량 분석용)"),
  c("set_no_bhs", "분석 대상(민감도)", "set_all 중 Brain Health Study 코호트 제외")
), function(x) data.table(column = x[1], group = x[2], description = x[3])))

cb <- fread(file.path(base_dir, "DAP_2023_CODEBOOK_v1.1.csv"))[DataFile %like% "^SamplesResults_"]
bio_src <- c(sub("^bio_cp_(.*)$", "krt_cp_\\1_value", grep("^bio_cp_", bio_vars, value = TRUE)),
             sub("^bio_cbc_", "krt_cbc_", setdiff(grep("^bio_cbc_", bio_vars, value = TRUE), "bio_cbc_nlr")))
col_desc <- rbind(col_desc,
  data.table(column = c(grep("^bio_cp_", bio_vars, value = TRUE), setdiff(grep("^bio_cbc_", bio_vars, value = TRUE), "bio_cbc_nlr")),
             group = "바이오마커(민감도·탐색)",
             description = paste0(fifelse(grepl("^krt_cp_", bio_src), "혈액 화학 검사: ", "혈구 검사(CBC): "),
                                  cb$SurveyText[match(bio_src, cb$Variable)], " (가장 최근 CSLB와 가장 가까운 검사)")),
  data.table(column = c("bio_cbc_nlr", "bio_ur_sg", "bio_ur_ph", "bio_ur_protein", "bio_cp_gap", "bio_cbc_gap", "bio_ur_gap", "bio_fasted_hr"),
             group = "바이오마커(민감도·탐색)",
             description = c("호중구/림프구 비 (전신 염증 지표)", "요비중 (>1.060은 1.060)", "요 pH",
                             "요단백 Negative 0 / Trace 0.5 / 1+ 1 / 2+ 2 / 3+ 3", "화학 검사와 가장 최근 CSLB의 날짜 간격(년)",
                             "CBC와 가장 최근 CSLB의 날짜 간격(년)", "소변 검사와 가장 최근 CSLB의 날짜 간격(년)",
                             "화학 검사 채혈 전 금식 시간 (혈당·중성지방 해석용)")))

R <- list(
  landsberg = "Landsberg 2019 WSAVA / Canine and Feline Dementia (Springer 2017): CCD는 배제 진단, 감각 저하·신경계·내분비·대사 질환·통증을 먼저 배제 https://www.vin.com/doc/?id=9382955",
  ozawa     = "Ozawa et al. 2019 J Vet Med Sci: CCDR 웹 설문 726마리에서 시력 장애가 CCD(≥50)와 유의하게 연관 https://pmc.ncbi.nlm.nih.gov/articles/PMC6943310/",
  hopper    = "Hopper et al. 2024 PLoS One: 8세 이상에서 시각·청각 이중 장애가 인지 장애 위험 증가(OR 1.8-2.0) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11482676/",
  fefer     = "Fefer et al. 2022 J Vet Intern Med: 청력(BAER)이 나쁠수록 인지 점수(CADES) 높음 https://pubmed.ncbi.nlm.nih.gov/35932193/",
  yarb      = "Yarborough et al. 2022 Sci Rep (DAP CSLB): 신경·눈·귀 질환 이력이 CCD 위험 증가와 연관, 비활동 OR 6.47 https://www.nature.com/articles/s41598-022-15837-9",
  packer    = "Packer et al. 2018 PLoS One: 특발성 뇌전증 개에서 CCDR ≥ 50 위험 증가, 발작 빈도 높을수록 CCDR 높음 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5805257/",
  smith     = "Smith, Murrell & Mendl 2024 Anim Cogn: 골관절염 만성 통증 개에서 공간 작업기억 저하 https://pmc.ncbi.nlm.nih.gov/articles/PMC10907419/",
  mccoy     = "McCoy et al. 2023 Evol Med Public Health (DAP): 가구 소득·보호자 나이·다른 동물과의 사회적 시간 등 사회 환경이 개 건강과 연관 https://pmc.ncbi.nlm.nih.gov/articles/PMC10306367/",
  kubinyi   = "Kubinyi, Turcsán & Miklósi 2009 Behav Processes: 보호자 나이·성별·학력, 가구 내 사람·개 수가 보호자 평가 개 행동 점수와 작지만 검출 가능한 연관 https://doi.org/10.1016/j.beproc.2009.04.004",
  arons     = "Arons et al. 2013 BMC Med Res Methodol (사람): 치매 환자 대리 평가에서 보호자 나이가 평가 점수를 예측 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3844427/",
  bhs       = "DAP Brain Health Study: CCD·알츠하이머 연관 연구 위해 Precision 코호트에서 선별 초청(2023~) https://content.dogagingproject.org/canine-cognition-introducing-the-brain-health-study/",
  salvin    = "Salvin et al. 2011 Vet J: CCDR 개발(≥50 = CCD) https://pubmed.ncbi.nlm.nih.gov/20542455/",
  rc        = "Royal Canin Vet Focus (CCD 감별): 신장·간 질환, 내분비 질환, 뇌 종괴 등이 CCD와 비슷한 증상 https://prd-vetfocus-cm-01.royalcanin.com/en/scientific/canine-cognitive-dysfunction",
  mci       = "Frontiers Vet Sci 2020 (개 경도 인지 장애): 혈청 AST·ALT 약간 증가, 나트륨·염소 감소 https://www.frontiersin.org/articles/10.3389/fvets.2020.551895/text",
  aschard   = "Aschard et al. 2015 AJHG: 유전되는 공변량 보정은 GWAS 효과 추정을 왜곡할 수 있음(collider bias) https://pmc.ncbi.nlm.nih.gov/articles/PMC4320269/")
ref_of <- list(
  vision_cov = c("landsberg", "ozawa", "hopper", "yarb", "aschard"), vision_cov_narrow = c("ozawa", "hopper"),
  hearing_cov = c("landsberg", "fefer", "hopper", "yarb"), ac_blind = c("ozawa", "hopper"), ac_deaf = c("fefer", "hopper"),
  seizure_cov = c("packer", "yarb"), vestibular_cov = c("landsberg", "yarb"), brain_dx = c("landsberg", "rc"), mimic_any = c("landsberg", "rc"),
  oa_cov = c("smith", "landsberg"), ac_mobility_pain = c("smith", "landsberg"), ac_mobility_weak = c("smith"), ac_other_pain = c("landsberg"),
  dementia_dx = c("salvin"), dementia_dx_ever = c("salvin"), cohort_bhs = c("bhs"),
  owner_age = c("kubinyi", "mccoy", "arons"), owner_edu = c("kubinyi", "mccoy"), owner_income = c("mccoy"),
  household_dogs = c("mccoy", "kubinyi"), household_people = c("mccoy", "kubinyi"), home_area = c("mccoy"),
  general_health = c("yarb"), pa_score = c("yarb"), wt_cov = c("aschard"), cslb_score = c("salvin", "yarb"))
for (v in grep("^bio_", col_desc$column, value = TRUE)) ref_of[[v]] <- c("rc", "mci", "landsberg")
col_desc[, reference := vapply(column, function(v) if (is.null(ref_of[[v]])) "" else paste(unlist(R[ref_of[[v]]]), collapse = " || "), character(1))]
col_desc[column %in% c("dob_exact", "set_dob_exact"), reference := "문헌 없음: 나이를 보호자가 추정한 개(입양견 등)는 가장 강한 공변량인 나이에 측정 오차가 생김 → 오차가 결과에 영향을 주는지 확인용"]
col_desc[column %in% c("afus_n_pre", "set_afus", "vet_freq"), reference := "문헌 없음: 진단 이력은 보고 기회(AFUS 응답, 병원 방문)가 있어야 기록되므로 진단 누락 영향 확인용"]
undoc <- setdiff(names(dt), col_desc$column)
if (length(undoc)) warning("설명이 없는 컬럼: ", paste(undoc, collapse = ", "))
col_desc <- col_desc[column %in% names(dt)][match(names(dt), column)]
col_desc[, `:=`(n_nonmissing_set_all = vapply(column, function(v) sum(!is.na(dt[set_all == TRUE][[v]])), integer(1)),
                example = vapply(column, function(v) as.character(dt[set_all == TRUE & !is.na(get(v))][[v]][1]), character(1)))]
fwrite(col_desc, paste0(out_prefix, "_column_description.csv"), bom = TRUE)
fwrite(release_tab, paste0(out_prefix, "_summary_release.csv"))
fwrite(anc_tab, paste0(out_prefix, "_summary_ancestry.csv"))
fwrite(age_grp_tab, paste0(out_prefix, "_summary_age_group.csv"))
fwrite(floor_tab, paste0(out_prefix, "_summary_age_floor.csv"))
fwrite(age_year_tab, paste0(out_prefix, "_summary_age_year.csv"))
fwrite(set_tab, paste0(out_prefix, "_summary_sets.csv"))
fwrite(cov_tab, paste0(out_prefix, "_summary_covariates.csv"))
fwrite(main_tab, paste0(out_prefix, "_summary_main_model.csv"))
fwrite(vis_tab, paste0(out_prefix, "_summary_vision_def.csv"))
fwrite(link_tab, paste0(out_prefix, "_summary_linked_vars.csv"))
fwrite(bio_tab, paste0(out_prefix, "_summary_biomarker.csv"))

id_cols <- c("FID", "IID")
sa <- dt[set_all == TRUE][, cslb_rint := rint(cslb_score)]
role_cols <- list(
  main = c("cslb_score", "cslb_rint", main_qcov, main_dcov, "set_age8", "set_age10", "set_purelike", "set_mixlike"),
  sensitivity = c("vision_cov_narrow", "hearing_cov", "seizure_cov", "vestibular_cov", "brain_dx", "mimic_any", "oa_cov",
                  "ac_blind", "ac_deaf", "ac_mobility_weak", "ac_mobility_pain", "ac_other_pain", "ac_weightloss", "ac_eat_drink", "ac_cleanliness", "ac_gap",
                  "general_health", "household_dogs", "household_people", "home_area", "owner_edu", "owner_income",
                  "dob_exact", "afus_n_pre", "vet_freq", "coi", "top_pct", "purelike_prob", "size_score",
                  "set_nomimic", "set_dob_exact", "set_afus", "set_ac", "set_no_bhs"),
  validation = c("cslb_score", "dementia_dx", "dementia_dx_ever", grep("^ac_", names(dt), value = TRUE),
                 "sleep_night_hr", "sleep_day_hr", "sleep_total_hr", "sleep_day_ratio", "sleep_src", "sleep_gap",
                 grep("^trt_", names(dt), value = TRUE), "n_cslb", "cslb_score_first", "cslb_score_mean", "cslb_score_max", "cslb_slope", "cslb_span"),
  biomarker = c("cslb_score", "age_cslb", bio_vars, "bio_cp_gap", "bio_cbc_gap", "bio_ur_gap", "bio_fasted_hr"))
role_cols$validation <- setdiff(unique(role_cols$validation), c("ac_blind", "ac_deaf", "ac_mobility_weak", "ac_mobility_pain", "ac_other_pain",
                                                                "ac_weightloss", "ac_eat_drink", "ac_cleanliness", "ac_gap"))
role_cols$validation <- c(role_cols$validation, "ac_gap")
for (r in names(role_cols)) {
  out <- sa[, c(id_cols, role_cols[[r]]), with = FALSE]
  if (r == "biomarker") out <- out[rowSums(!is.na(out[, bio_vars, with = FALSE])) > 0]
  fwrite(out, sprintf("%s_df_%s.csv", out_prefix, r))
  cat(sprintf("   역할별 데이터: %-12s %5d마리 x %3d열\n", r, nrow(out), ncol(out)))
}
cat("\n저장:", out_dir, "의 00_* / 그림:", pic_dir, "\n")
