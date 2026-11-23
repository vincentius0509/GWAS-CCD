# 06-7_cslb_measurement.R - CSLB measurement properties in all DAP surveys: test-retest reliability, dementia-diagnosis AUC. See README.md.
#         SurveyData/v2023/DAP_2023_HLES_health_condition_v1.0.csv · DAP_2023_AFUS_health_condition_v1.0.csv ·
#   06-7_by_age.csv · _items_by_age.csv · _consistency.csv · _retest.csv · _icc.csv · _change_by_age.csv · _validity.csv · _summary.txt

library(data.table)
library(ggplot2)
library(mgcv)
options(width = 220)
set.seed(20261112)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")    
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
out_prefix <- file.path(qc_dir, "06-7")
pic_prefix <- file.path(pic_dir, "06-7")
cslb_file <- file.path(base_dir, "SurveyData/v2023/DAP_2023_CSLB_v1.0.csv")
dog_file <- file.path(base_dir, "DogOverview/v2023/DAP_2023_DogOverview_v1.1.csv")
hles_file <- file.path(base_dir, "SurveyData/v2023/DAP_2023_HLES_health_condition_v1.0.csv")
afus_file <- file.path(base_dir, "SurveyData/v2023/DAP_2023_AFUS_health_condition_v1.0.csv")
geno_file <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023/04-2_pheno.txt")
for (f in c(cslb_file, dog_file, hles_file, afus_file, geno_file)) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
n_boot <- as.integer(Sys.getenv("N_BOOT", "1000"))

freq_it <- c("cslb_pace", "cslb_stare", "cslb_stuck", "cslb_recognize", "cslb_walk_walls", "cslb_avoid", "cslb_find_food")
chg_it <- c("cslb_pace_6mo", "cslb_stare_6mo", "cslb_defecate_6mo", "cslb_food_6mo", "cslb_recognize_6mo", "cslb_active_6mo")
items <- c(freq_it, chg_it)
w <- setNames(rep(1, 13), items); w[c("cslb_food_6mo", "cslb_recognize_6mo")] <- c(2, 3)
score_none <- sum(w[freq_it] * 1) + sum(w[chg_it] * 3)
ccd_cut <- 50

cs <- fread(cslb_file, select = c("dog_id", "cslb_date", items, "cslb_score"))
dob <- fread(dog_file, select = c("dog_id", "Estimated_DOB"))
cs <- merge(cs, dob, by = "dog_id")
cs[, `:=`(date = as.IDate(cslb_date), dob = as.IDate(Estimated_DOB))]
cs[, age := as.numeric(date - dob) / 365.25]
n_raw <- nrow(cs)
cs <- cs[complete.cases(cs[, c(items, "cslb_score", "age"), with = FALSE]) & age >= 0.5 & age < 25]
chk <- cs[, max(abs(as.matrix(.SD) %*% w - cslb_score)), .SDcols = items]
if (chk > 1e-6) warning("총점 = 가중 문항 합이 아닌 응답이 있음 (최대 차이 ", chk, ")")
setorder(cs, dog_id, date)
cs[, `:=`(visit = seq_len(.N), n_visit = .N), by = dog_id]
geno_ids <- fread(geno_file, header = FALSE)$V2
cs[, genotyped := dog_id %in% geno_ids]
cs[, age_bin := pmin(floor(age), 17)]
cs[, age_grp := cut(age, c(0, 8, 12, Inf), right = FALSE, labels = c("<8", "8-12", ">=12"))]
cs[, `:=`(none = cslb_score <= score_none, ge40 = cslb_score >= 40, ge50 = cslb_score >= ccd_cut)]
cat(sprintf("CSLB 응답 %s회 중 분석 %s회 (개 %s마리; 유전형 개 %s마리 · %s회), 반복 측정 개 %s마리\n",
            format(n_raw, big.mark = ","), format(nrow(cs), big.mark = ","), format(uniqueN(cs$dog_id), big.mark = ","),
            format(cs[genotyped == TRUE, uniqueN(dog_id)], big.mark = ","), format(cs[genotyped == TRUE, .N], big.mark = ","),
            format(cs[n_visit >= 2, uniqueN(dog_id)], big.mark = ",")))

ci_mean <- function(x) { x <- x[is.finite(x)]; n <- length(x); m <- mean(x); s <- if (n > 1) sd(x) / sqrt(n) else NA_real_
  list(n = n, mean = m, lo = m - qt(0.975, max(n - 1, 1)) * s, hi = m + qt(0.975, max(n - 1, 1)) * s) }
r_ci <- function(r, n) { z <- atanh(r); se <- 1 / sqrt(n - 3); list(lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se)) }

dog_bin <- cs[, .(score = mean(cslb_score), none = mean(none), ge40 = mean(ge40), ge50 = mean(ge50), genotyped = genotyped[1]), by = .(dog_id, age_bin)]
by_age <- dog_bin[, c(list(n_dogs = .N), ci_mean(score)[c("mean", "lo", "hi")], list(sd = sd(score), pct_none = 100 * mean(none),
                     pct_ge40 = 100 * mean(ge40), pct_ge50 = 100 * mean(ge50))), by = age_bin][order(age_bin)]
by_age_geno <- dog_bin[genotyped == TRUE, .(n_dogs_geno = .N, mean_geno = mean(score), pct_none_geno = 100 * mean(none), pct_ge50_geno = 100 * mean(ge50)), by = age_bin]
by_age <- merge(by_age, by_age_geno, by = "age_bin", all.x = TRUE)
it_long <- melt(cs[, c("dog_id", "age_bin", items), with = FALSE], id.vars = c("dog_id", "age_bin"), variable.name = "item")
it_long[, problem := fifelse(item %in% freq_it, value >= 2, value >= 4)]
items_by_age <- it_long[, .(pct_problem = 100 * mean(problem)), by = .(item, age_bin)]

first <- cs[visit == 1]
alpha_of <- function(X) { k <- ncol(X); v <- apply(X, 2, var); k / (k - 1) * (1 - sum(v) / var(rowSums(X))) }
omega_of <- function(X) {
  X <- X[, apply(X, 2, var) > 0, drop = FALSE]
  f <- tryCatch(factanal(covmat = cor(X), factors = 1, n.obs = nrow(X)), error = function(e) NULL)
  if (is.null(f)) return(NA_real_); l <- f$loadings[, 1]; s <- sum(l); s^2 / (s^2 + sum(f$uniquenesses))
}
consistency <- rbindlist(lapply(c("all", levels(cs$age_grp)), function(g) {
  d <- if (g == "all") first else first[age_grp == g]
  X <- as.matrix(d[, ..items]); Xw <- sweep(X, 2, w, `*`)
  itc <- vapply(seq_len(ncol(Xw)), function(j) cor(Xw[, j], rowSums(Xw[, -j])), 0)
  data.table(age_grp = g, n_dogs = nrow(d), alpha = alpha_of(Xw), omega = omega_of(X),
             alpha_freq = alpha_of(X[, freq_it]), alpha_change = alpha_of(X[, chg_it]),
             item_rest_min = min(itc, na.rm = TRUE), item_rest_median = median(itc, na.rm = TRUE), sd_score = sd(d$cslb_score))
}))

pr <- cs[, .(dog_id, visit, age, cslb_score, none)][, `:=`(age2 = shift(age, -1), score2 = shift(cslb_score, -1), none2 = shift(none, -1)), by = dog_id][!is.na(age2)]
pr[, `:=`(gap = age2 - age, mid = (age + age2) / 2)]
pr1 <- pr[gap >= 0.75 & gap <= 1.5]
pr1[, base_grp := cut(age, c(0, 4, 8, 10, 12, 14, Inf), right = FALSE, labels = c("<4", "4-8", "8-10", "10-12", "12-14", ">=14"))]
retest <- pr1[, { d <- .SD[, .SD[1], by = dog_id]
  r <- cor(d$cslb_score, d$score2); rs <- cor(d$cslb_score, d$score2, method = "spearman"); ci <- r_ci(r, nrow(d))
  list(n_dogs = nrow(d), pearson = r, lo95 = ci$lo, hi95 = ci$hi, spearman = rs, sd_score1 = sd(d$cslb_score),
       mean_change = mean(d$score2 - d$cslb_score), pct_none_both = 100 * mean(d$none & d$none2)) }, by = base_grp][order(base_grp)]
icc_of <- function(d) {
  d <- d[, if (.N >= 2) .SD, by = dog_id]; if (!nrow(d)) return(NULL)
  g <- d[, .(n = .N, m = mean(res), ss = sum((res - mean(res))^2)), by = dog_id]
  N <- sum(g$n); a <- nrow(g); gm <- mean(d$res)
  msb <- sum(g$n * (g$m - gm)^2) / (a - 1); msw <- sum(g$ss) / (N - a); k0 <- (N - sum(g$n^2) / N) / (a - 1)
  list(n_dogs = a, n_obs = N, k0 = k0, icc = (msb - msw) / (msb + (k0 - 1) * msw), var_between = (msb - msw) / k0, var_within = msw)
}
fit_age <- bam(cslb_score ~ s(age, k = 20), data = cs, discrete = TRUE)
cs[, res := cslb_score - predict(fit_age, cs)]
icc <- rbindlist(lapply(c("all", levels(cs$age_grp)), function(g) {
  d <- if (g == "all") cs else cs[age_grp == g]
  r <- icc_of(d[, .(dog_id, res)]); if (is.null(r)) NULL else c(list(age_grp = g), r)
}))
sb <- function(r, k) k * r / (1 + (k - 1) * r)

pr[, `:=`(rate = (score2 - cslb_score) / gap, mid_bin = pmin(floor(mid), 17), onset = none & !none2, any_none = none)]
pr_ok <- pr[gap >= 0.5 & gap <= 3.5]
change_by_age <- pr_ok[, .(rate = mean(rate), onset = mean(onset[any_none]) / mean(gap[any_none]), n_pairs = .N), by = .(dog_id, mid_bin)][
  , c(list(n_dogs = .N), setNames(ci_mean(rate)[c("mean", "lo", "hi")], c("rate", "rate_lo", "rate_hi")),
      list(onset_per_year = 100 * mean(onset, na.rm = TRUE))), by = mid_bin][order(mid_bin)]
sig_up <- change_by_age[n_dogs >= 30, .(mid_bin, up = rate_lo > 0)]
onset_age <- { u <- sig_up[order(mid_bin)]; i <- which(vapply(seq_len(nrow(u)), function(j) all(u$up[j:nrow(u)]), TRUE))[1]; if (is.na(i)) NA_integer_ else u$mid_bin[i] }
meaningful <- 1
mean_age <- { u <- change_by_age[n_dogs >= 30][order(mid_bin)]; i <- which(vapply(seq_len(nrow(u)), function(j) all(u$rate[j:nrow(u)] >= meaningful), TRUE))[1]
              if (is.na(i)) NA_integer_ else u$mid_bin[i] }
rate_young <- ci_mean(pr_ok[mid < 8, .(r = mean(rate)), by = dog_id]$r)
rate_old <- ci_mean(pr_ok[mid >= 12, .(r = mean(rate)), by = dog_id]$r)

h <- fread(hles_file, select = c("dog_id", "hs_condition", "hs_diagnosis_year"), integer64 = "double")
a <- fread(afus_file, select = c("dog_id", "afus_hs_new_condition", "afus_hs_new_condition_diagnosis_year"), integer64 = "double")
dx <- unique(rbind(h[hs_condition == 1205, .(dog_id, dx_year = hs_diagnosis_year)], a[afus_hs_new_condition == 1205, .(dog_id, dx_year = afus_hs_new_condition_diagnosis_year)]))
dx <- dx[, .(dx_year = suppressWarnings(min(as.numeric(dx_year), na.rm = TRUE))), by = dog_id]
dog_lv <- cs[, .(max_score = max(cslb_score), last_score = cslb_score[.N], last_age = age[.N], genotyped = genotyped[1]), by = dog_id]
dog_lv[, dementia := dog_id %in% dx$dog_id]
auc_of <- function(x, y) { r <- rank(c(x, y)); (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y)) }
boot_auc <- function(d, v) {
  pos <- d[dementia == TRUE, get(v)]; neg <- d[dementia == FALSE, get(v)]
  if (length(pos) < 10) return(list(auc = NA_real_, lo = NA_real_, hi = NA_real_))
  b <- replicate(n_boot, auc_of(sample(pos, replace = TRUE), sample(neg, replace = TRUE)))
  list(auc = auc_of(pos, neg), lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975)))
}
dog_lv[, last_grp := cut(last_age, c(0, 8, 12, 14, Inf), right = FALSE, labels = c("<8", "8-12", "12-14", ">=14"))]
validity <- rbindlist(lapply(c(">=8 all", "8-12", "12-14", ">=14"), function(g) {
  d <- if (g == ">=8 all") dog_lv[last_age >= 8] else dog_lv[last_grp == g]
  rbindlist(lapply(c("max_score", "last_score"), function(v) {
    b <- boot_auc(d, v)
    data.table(age_grp = g, score = v, n_dementia = d[dementia == TRUE, .N], n_other = d[dementia == FALSE, .N],
               mean_dementia = d[dementia == TRUE, mean(get(v))], mean_other = d[dementia == FALSE, mean(get(v))],
               auc = b$auc, auc_lo = b$lo, auc_hi = b$hi,
               sens50 = d[dementia == TRUE, mean(get(v) >= ccd_cut)], spec50 = d[dementia == FALSE, mean(get(v) < ccd_cut)])
  }))
}))
strat <- validity[age_grp != ">=8 all" & !is.na(auc), .(auc_age_stratified = weighted.mean(auc, n_dementia)), by = score]

geno_age <- cs[genotyped == TRUE, .(age = age[.N], n = .N), by = dog_id]
geno_age[, base_grp := cut(age, c(0, 4, 8, 10, 12, 14, Inf), right = FALSE, labels = levels(pr1$base_grp))]
mix <- merge(geno_age[, .N, by = base_grp], retest[, .(base_grp, pearson)], by = "base_grp")
rel_geno <- mix[, sum(N * pearson) / sum(N)]
k_mean <- geno_age[, mean(n)]

fwrite(by_age, paste0(out_prefix, "_by_age.csv")); fwrite(items_by_age, paste0(out_prefix, "_items_by_age.csv"))
fwrite(consistency, paste0(out_prefix, "_consistency.csv")); fwrite(retest, paste0(out_prefix, "_retest.csv")); fwrite(icc, paste0(out_prefix, "_icc.csv"))
fwrite(change_by_age, paste0(out_prefix, "_change_by_age.csv")); fwrite(validity, paste0(out_prefix, "_validity.csv"))

lab17 <- function(x) ifelse(x >= 17, "17+", as.character(x))
pa <- melt(by_age[, .(age_bin, `Mean CSLB score (95% CI)` = mean, `% dogs with score >= 50` = pct_ge50, `% dogs with no problems (<= 34)` = pct_none)],
           id.vars = "age_bin")
pa <- merge(pa, by_age[, .(age_bin, lo, hi)], by = "age_bin")
pa[variable != "Mean CSLB score (95% CI)", `:=`(lo = NA, hi = NA)]
g1 <- ggplot(pa, aes(age_bin, value)) +
  geom_line(colour = "#1f4e9c", linewidth = 0.8) + geom_point(colour = "#1f4e9c", size = 1.8) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.3, colour = "#1f4e9c", na.rm = TRUE) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = 0:17, labels = lab17) +
  labs(x = "Age at survey (years)", y = NULL, title = "CSLB by age (DAP, all dogs; one value per dog per age year)",
       subtitle = sprintf("Score 34 = no problems and no change; 50 = CCD cut-off (Salvin 2011). Dogs %s, surveys %s",
                          format(uniqueN(cs$dog_id), big.mark = ","), format(nrow(cs), big.mark = ","))) +
  theme_minimal(base_size = 11) + theme(strip.text = element_text(face = "bold", hjust = 0))
ggsave(paste0(pic_prefix, "_score_by_age.png"), g1, width = 8, height = 8, dpi = 300)

cb <- change_by_age[n_dogs >= 30]
g2 <- ggplot(cb, aes(mid_bin, rate)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = rate_lo, ymax = rate_hi, colour = rate_lo > 0), size = 0.4) +
  scale_colour_manual(values = c(`TRUE` = "#b4462b", `FALSE` = "#7f8fa6"), labels = c(`TRUE` = "95% CI above 0", `FALSE` = "CI includes 0"), name = NULL) +
  geom_hline(yintercept = meaningful, linetype = 3, colour = "grey40") +
  { if (!is.na(mean_age)) geom_vline(xintercept = mean_age - 0.5, linetype = 2, colour = "#b4462b") } +
  scale_x_continuous(breaks = 0:17, labels = lab17) +
  labs(x = "Age at mid-point of two consecutive surveys (years)", y = "Within-dog change in CSLB (points per year)",
       title = "Within-dog change in CSLB with age",
       subtitle = sprintf("Dog-level annual change, mean and 95%% CI (bins >= 30 dogs)\nCI above 0 from age %s; >= %g point/yr (dotted line) from age %s (dashed line)",
                          ifelse(is.na(onset_age), "-", onset_age), meaningful, ifelse(is.na(mean_age), "-", mean_age))) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(paste0(pic_prefix, "_change_by_age.png"), g2, width = 8, height = 5, dpi = 300)

ib <- copy(items_by_age)[, item := factor(sub("^cslb_", "", item), levels = sub("^cslb_", "", rev(items)))]
g3 <- ggplot(ib, aes(age_bin, item, fill = pct_problem)) + geom_tile(colour = "white") +
  scale_fill_gradient(low = "#f4f6fb", high = "#1f4e9c", name = "% with\nproblem") +
  scale_x_continuous(breaks = 0:17, labels = lab17) +
  labs(x = "Age at survey (years)", y = NULL, title = "Item-level problem rate by age",
       subtitle = "Problem = frequency item >= 2 (at least monthly) or change item >= 4 (worse than 6 months ago)") +
  theme_minimal(base_size = 11) + theme(panel.grid = element_blank())
ggsave(paste0(pic_prefix, "_items_by_age.png"), g3, width = 9, height = 5, dpi = 300)

g4 <- ggplot(retest, aes(base_grp, pearson)) +
  geom_pointrange(aes(ymin = lo95, ymax = hi95), colour = "#1f4e9c", size = 0.5) +
  geom_text(aes(label = sprintf("n %s", format(n_dogs, big.mark = ","))), nudge_x = 0.3, size = 3, colour = "grey30") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Age at first of two surveys (years)", y = "Correlation between surveys ~1 year apart",
       title = "Test-retest correlation of CSLB by age", subtitle = "One pair per dog (interval 0.75-1.5 years); 95% CI") +
  theme_minimal(base_size = 11)
ggsave(paste0(pic_prefix, "_retest_by_age.png"), g4, width = 7, height = 4.5, dpi = 300)

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("06-7 CSLB 측정 특성 (DAP 전체 응답; 유전형은 쓰지 않음)\n")
cat(sprintf("응답 %s회 · 개 %s마리 (반복 측정 %s마리), 유전형 개 %s마리. 점수 = 가중 문항 합 (16 ~ 80), 문제·변화 없음 = %d점, CCD 기준 %d점\n\n",
            format(nrow(cs), big.mark = ","), format(uniqueN(cs$dog_id), big.mark = ","), format(cs[n_visit >= 2, uniqueN(dog_id)], big.mark = ","),
            format(cs[genotyped == TRUE, uniqueN(dog_id)], big.mark = ","), score_none, ccd_cut))
cat("■ A. 나이별 분포 (개 단위; pct_none = 34점 이하 = 문제·변화 없음)\n"); print(f3(by_age))
cat("\n■ B. 내적 일관성 (개마다 첫 측정; alpha = 가중 총점, omega = 한 요인 factanal 근사)\n"); print(f3(consistency))
cat("\n■ C. 반복 측정 신뢰도\n  연속 두 측정 (간격 0.75 ~ 1.5년) 상관, 첫 측정 나이 구간별:\n"); print(f3(retest))
cat("  ICC(1) (나이 곡선을 뺀 잔차; var_between = 개 사이, var_within = 같은 개 안):\n"); print(f3(icc))
cat(sprintf("  유전형 개의 나이 분포로 가중한 1년 재측정 상관 ≈ %.2f, 개당 평균 측정 %.1f회 → 평균 표현형의 신뢰도 (Spearman-Brown) ≈ %.2f\n",
            rel_geno, k_mean, sb(rel_geno, k_mean)))
cat(sprintf("  → 측정 오차만으로 관측 h²는 실제의 약 %.0f%% (최근 1회) · %.0f%% (반복 평균)로 줄어든다고 기대 (06-6 h2_main vs h2_mean과 대조)\n",
            100 * rel_geno, 100 * sb(rel_geno, k_mean)))
cat("\n■ D. 나이에 따른 변화 (같은 개 안, 연속 두 측정의 연간 변화; 개 단위 평균 · 95% CI; onset = 34점 이하였던 개가 다음에 넘는 연간 비율 %)\n")
print(f3(change_by_age))
cat(sprintf("  8세 미만: %.2f점/년 (95%% CI %.2f ~ %.2f, 개 %s) / 12세 이상: %.2f점/년 (%.2f ~ %.2f, 개 %s)\n",
            rate_young$mean, rate_young$lo, rate_young$hi, format(rate_young$n, big.mark = ","), rate_old$mean, rate_old$lo, rate_old$hi, format(rate_old$n, big.mark = ",")))
cat(sprintf("  같은 개 안에서 CSLB가 유의하게 오르기 시작해 계속 오르는 나이 (95%% CI 하한 > 0이 이어지는 첫 구간, 개 30마리 이상): %s세\n",
            ifelse(is.na(onset_age), "없음", onset_age)))
cat(sprintf("  의미 있는 크기 (연 %g점 이상)로 계속 오르기 시작하는 나이: %s세  ← 개 수가 많아 아주 작은 변화도 유의해지므로 이 값을 함께 봄\n",
            meaningful, ifelse(is.na(mean_age), "없음", mean_age)))
cat("\n■ E. 진단 타당도: 보호자 보고 수의사 진단 'Dementia or senility' (HLES·AFUS 1205) vs 나머지\n"); print(f3(validity))
cat("  나이 구간 안 AUC의 진단견 수 가중 평균:\n"); print(f3(strat))
sink()
cat(sprintf("저장: %s_by_age.csv · _items_by_age.csv · _consistency.csv · _retest.csv · _icc.csv · _change_by_age.csv · _validity.csv · _summary.txt\n그림: %s_score_by_age.png · _change_by_age.png · _items_by_age.png · _retest_by_age.png\n",
            out_prefix, pic_prefix))

#      (by_age 2cbdf494 · change_by_age c9674de0 · consistency c60acbac · icc a9a1634e · items_by_age 436ef9da · retest 7d7fd52d · validity 93c91444)
