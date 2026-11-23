# 14-3_power_table.R - Sample sizes for 80% power under alternative survey designs. See README.md.
suppressMessages(library(data.table))
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
f_h2 <- file.path(d0, "06-6_h2.csv"); f_rt <- file.path(d0, "06-7_retest.csv")
for (f in c(f_h2, f_rt)) if (!file.exists(f)) stop("파일 없음: ", f)
power <- 0.80; alpha_h2 <- 0.05; alpha_snp <- 5e-8
h2_grid <- c(0.05, 0.10, 0.15, 0.20); q2_grid <- c(0.005, 0.01, 0.02)
rel_all1 <- 0.39; rel_allmean <- 0.64
sb <- function(r, k) k * r / (1 + (k - 1) * r)                # Spearman-Brown

H <- fread(f_h2)[model %in% c("h2_main", "h2_young8", "h2_age8", "h2_sex_m", "h2_sex_f")][, c_se := se * n]
c_med <- median(H$c_se); c_rng <- range(H$c_se)
R <- fread(f_rt)
if (!all(c("base_grp", "n_dogs", "pearson") %in% names(R))) stop("06-7_retest.csv 열 이름이 예상과 다름: ", paste(names(R), collapse = " "))
old <- R[as.character(base_grp) %in% c("12-14", ">=14")]
if (nrow(old) != 2) stop("06-7_retest.csv에서 12-14 · >=14 구간을 못 찾음")
rel_old1 <- old[, sum(pearson * n_dogs) / sum(n_dogs)]
designs <- data.table(design = c("All ages, single score (current)", "All ages, mean of repeated scores (2.8 per dog)",
                                 "Dogs >= 12 y, single score", "Dogs >= 12 y, mean of 3 annual scores"),
                      reliability = c(rel_all1, rel_allmean, rel_old1, sb(rel_old1, 3)))

k_h2 <- qnorm(1 - alpha_h2) + qnorm(power)
ncp_snp <- (qnorm(1 - alpha_snp / 2) + qnorm(power))^2
T1 <- CJ(design = designs$design, h2_true = h2_grid, sorted = FALSE)[designs, on = "design"]
T1[, `:=`(target = "SNP heritability > 0 (one-sided P < 0.05)", effect = h2_true, effect_obs = h2_true * reliability)]
T1[, `:=`(N_req = ceiling(k_h2 * c_med / effect_obs), N_req_lo = ceiling(k_h2 * c_rng[1] / effect_obs), N_req_hi = ceiling(k_h2 * c_rng[2] / effect_obs))]
T2 <- CJ(design = designs$design, q2_true = q2_grid, sorted = FALSE)[designs, on = "design"]
T2[, `:=`(target = "Single variant (P < 5e-8)", effect = q2_true, effect_obs = q2_true * reliability)]
T2[, N_req := ceiling(ncp_snp * (1 - effect_obs) / effect_obs)]
TT <- rbind(T1[, .(target, design, reliability, effect, effect_obs, N_req, N_req_lo, N_req_hi)],
            T2[, .(target, design, reliability, effect, effect_obs, N_req, N_req_lo = NA_real_, N_req_hi = NA_real_)])
TT[, `:=`(reliability = round(reliability, 2), effect_obs = signif(effect_obs, 3))]
fwrite(TT, file.path(d0, "14-3_power_table.csv"))

h2_hi <- fread(f_h2)[model == "h2_main", hi95]
q2_det <- ncp_snp / (7115 + ncp_snp)
rel_start <- c("8-10 y" = R[as.character(base_grp) == "8-10", pearson], "10-12 y" = R[as.character(base_grp) == "10-12", pearson], ">= 12 y" = rel_old1)
if (length(rel_start) != 3 || anyNA(rel_start)) stop("06-7_retest.csv에서 8-10 · 10-12 구간을 못 찾음")
DG <- CJ(start_age = names(rel_start), k = 1:5, sorted = FALSE)[, rel1 := rel_start[start_age]][, reliability := sb(rel1, k)]
DG[, `:=`(N_h2_5 = ceiling(k_h2 * c_med / (0.05 * reliability)), N_h2_10 = ceiling(k_h2 * c_med / (0.10 * reliability)),
          N_snp_1 = ceiling(ncp_snp * (1 - 0.01 * reliability) / (0.01 * reliability)), N_snp_2 = ceiling(ncp_snp * (1 - 0.02 * reliability) / (0.02 * reliability)))]
DG[, `:=`(surveys_total_h2_5 = N_h2_5 * k, reliability = round(reliability, 2), rel1 = round(rel1, 3))]
fwrite(DG, file.path(d0, "14-3_power_design_grid.csv"))

sink(file.path(d0, "14-3_power_summary.txt"))
cat("14-3 필요 표본 수 (80% 검정력; 방법 · 가정은 코드 머리말)\n")
cat(sprintf("GREML SE × N 상수 c: 중앙값 %.0f (범위 %.0f ~ %.0f; 06-6 모형 5개)\n", c_med, c_rng[1], c_rng[2]))
print(H[, .(model, n, se, c_se = round(c_se))], row.names = FALSE)
cat("\n신뢰도:\n"); print(designs[, .(design, reliability = round(reliability, 2))], row.names = FALSE)
cat("\n■ SNP 유전율 > 0 검출 (N_req = 중앙값 c, [lo, hi] = c 범위)\n")
print(dcast(T1, design ~ paste0("h2_", h2_true), value.var = "N_req"), row.names = FALSE)
cat("\n■ 단일 변이 (P < 5e-8) 검출\n")
print(dcast(T2, design ~ paste0("q2_", q2_true), value.var = "N_req"), row.names = FALSE)
cat("\n■ 현실적 설계 격자: 시작 나이 × 매년 설문 횟수 k (N_h2_5 · N_h2_10 = 실제 h² 5 · 10% 검출, N_snp_1 · N_snp_2 = 분산 1 · 2% 변이 검출; surveys_total = 마릿수 × k)\n")
print(DG, row.names = FALSE)
cat(sprintf("\n참고 (현재 연구): 관측 h² 95%% 상한 %.3f → 측정 신뢰도 %.2f로 되돌리면 실제 h² 상한 약 %.2f.\n", h2_hi, rel_all1, h2_hi / rel_all1))
cat(sprintf("  현재 N 7,115에서 80%% 검정력으로 잡을 수 있는 단일 변이 = 관측 q² %.2f%% → 실제 q² 약 %.1f%% (신뢰도 %.2f)\n", 100 * q2_det, 100 * q2_det / rel_all1, rel_all1))
sink()
cat(readLines(file.path(d0, "14-3_power_summary.txt")), sep = "\n")
