# 11-1_pheno_purify.R - Flags dogs with a sensory or neurological diagnosis before the CSLB survey. See README.md.
library(data.table)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
sv <- file.path(base_dir, "SurveyData/v2023")
out_prefix <- file.path(qc_dir, "11-1")
hles_f <- file.path(sv, "DAP_2023_HLES_health_condition_v1.0.csv"); afus_f <- file.path(sv, "DAP_2023_AFUS_health_condition_v1.0.csv")
hcan_f <- file.path(sv, "DAP_2023_HLES_cancer_conditions_v1.0.csv"); acan_f <- file.path(sv, "DAP_2023_AFUS_cancer_conditions_v1.0.csv")
cslb_f <- file.path(sv, "DAP_2023_CSLB_v1.0.csv")
pheno_f <- file.path(in_dir, "04-2_pheno.txt"); a8_keep <- file.path(in_dir, "05-1_inputs", "age8_keep.txt")
for (f in c(hles_f, afus_f, hcan_f, acan_f, cslb_f, pheno_f, a8_keep)) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
BRAIN <- -1L
codes <- list(sensory = c(101, 102, 115, 201, 205), sens_severe = c(101, 201),
              neuro = c(1215, 1216, 1210, 1298, 1201, 1202, 801, BRAIN), pain = c(1113, 1111, 1118, 1905),
              systemic = c(1311, 1301, 1903, 511, 512))

to_date <- function(y, m) {
  y <- suppressWarnings(as.integer(y)); m <- suppressWarnings(as.integer(m))
  y[!(y %in% 1990:2025)] <- NA; m[is.na(m) | m < 1 | m > 12] <- 6L
  as.IDate(ifelse(is.na(y), NA_character_, sprintf("%d-%02d-15", y, m)))
}
cc <- list(character = c("hs_diagnosis_year", "hs_diagnosis_month"))
h <- fread(hles_f, select = c("dog_id", "hs_condition", "hs_diagnosis_year", "hs_diagnosis_month"), integer64 = "double", colClasses = cc)[
  , .(dog_id = as.character(dog_id), code = as.integer(hs_condition), dx = to_date(hs_diagnosis_year, hs_diagnosis_month), src = "HLES")]
a <- fread(afus_f, select = c("dog_id", "afus_hs_new_condition", "afus_hs_new_condition_diagnosis_year", "afus_hs_new_condition_diagnosis_month"), integer64 = "double",
           colClasses = list(character = c("afus_hs_new_condition_diagnosis_year", "afus_hs_new_condition_diagnosis_month")))[
  , .(dog_id = as.character(dog_id), code = as.integer(afus_hs_new_condition), dx = to_date(afus_hs_new_condition_diagnosis_year, afus_hs_new_condition_diagnosis_month), src = "AFUS")]
hc <- fread(hcan_f, select = c("dog_id", "hs_cancer_locations_brain", "hs_initial_diagnosis_year", "hs_initial_diagnosis_month"), integer64 = "double",
            colClasses = list(character = c("hs_initial_diagnosis_year", "hs_initial_diagnosis_month")))[hs_cancer_locations_brain %in% c(TRUE, 1)][
  , .(dog_id = as.character(dog_id), code = BRAIN, dx = to_date(hs_initial_diagnosis_year, hs_initial_diagnosis_month), src = "HLES_cancer")]
ac <- fread(acan_f, select = c("dog_id", "afus_hs_new_cancer_locations_brain", "afus_hs_new_initial_diagnosis_year", "afus_hs_new_initial_diagnosis_month"), integer64 = "double",
            colClasses = list(character = c("afus_hs_new_initial_diagnosis_year", "afus_hs_new_initial_diagnosis_month")))[afus_hs_new_cancer_locations_brain %in% c(TRUE, 1)][
  , .(dog_id = as.character(dog_id), code = BRAIN, dx = to_date(afus_hs_new_initial_diagnosis_year, afus_hs_new_initial_diagnosis_month), src = "AFUS_cancer")]
cond <- unique(rbind(h, a, hc, ac)[!is.na(code)])

cs <- fread(cslb_f, select = c("dog_id", "cslb_date"), integer64 = "double")[, .(cslb_date = max(as.IDate(cslb_date))), by = .(dog_id = as.character(dog_id))]
cond <- merge(cond, cs, by = "dog_id", all.x = TRUE)
cond[, pre := is.na(dx) | is.na(cslb_date) | dx <= cslb_date]

ids <- fread(pheno_f, header = FALSE, colClasses = list(character = 1:2))[, .(FID = V1, IID = V2)]
a8 <- fread(a8_keep, header = FALSE, colClasses = "character")$V2
fl <- copy(ids)
for (k in names(codes)) {
  hit <- cond[code %in% codes[[k]], .(codes = paste(sort(unique(code)), collapse = ";"), any_pre = any(pre)), by = dog_id]
  fl[, (k) := IID %in% hit$dog_id]; fl[, (paste0(k, "_pre")) := IID %in% hit[any_pre == TRUE, dog_id]]
  fl[hit, on = .(IID = dog_id), (paste0(k, "_codes")) := i.codes]
}
fl[, `:=`(age8 = IID %in% a8,
          excl_primary = sensory_pre | neuro_pre,
          excl_primary_any = sensory | neuro,
          excl_sensory = sensory_pre,
          excl_severe = sens_severe_pre | neuro_pre,
          excl_broad = sensory_pre | neuro_pre | pain_pre | systemic_pre)]
fwrite(fl, paste0(out_prefix, "_flags.csv"))

vv <- c("sensory_pre", "neuro_pre", "pain_pre", "systemic_pre", "excl_primary", "excl_primary_any", "excl_sensory", "excl_severe", "excl_broad")
cnt <- rbindlist(lapply(c("전체", "8세 이상"), function(g) { s <- if (g == "전체") fl else fl[age8 == TRUE]
  cbind(data.table(group = g, n = nrow(s)), as.data.table(lapply(setNames(vv, vv), function(v) sum(s[[v]])))) }))
code_n <- cond[dog_id %in% ids$IID & code %in% unlist(codes), .(n_dogs = uniqueN(dog_id), n_dogs_pre = uniqueN(dog_id[pre])), by = code][order(code)]
code_n[code == BRAIN, code := NA][, label := fifelse(is.na(code), "뇌종양 (암 위치)", "")]
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("11-1 표현형 정제용 질환 표시 (가설 H8 준비; 보호자 보고 HLES · AFUS, 주 정의 = CSLB 이전 진단)\n\n■ 표시된 개 수 (_pre = 진단이 쓰인 CSLB 이전)\n"); print(cnt)
cat(sprintf("\n시점을 맞추면 주 정제에서 덜 빠지는 개: 전체 %d · 8세 이상 %d\n", fl[excl_primary_any & !excl_primary, .N], fl[age8 & excl_primary_any & !excl_primary, .N]))
cat("\n■ 코드별 유전형 개 수 (시점 무관 / CSLB 이전)\n"); print(code_n)
sink()
cat(sprintf("저장: %s_flags.csv · _summary.txt\n", out_prefix))

