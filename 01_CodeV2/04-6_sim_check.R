# 04-6_sim_check.R - Checks the simulated-phenotype experiment. See README.md.
#        Vision 0.008, Sex 0.005, Cohort_BHS 0.004

library(data.table)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")   
anc_file <- file.path(base_dir, "GeneticData/v2023/DogAgingProject_GeneticData_CuratedRelease_2023.ancestry.tsv")
out_prefix <- file.path(qc_dir, "04-6")
struct_var <- 0.3; sim_seed <- 1; n_sim <- 20
min_breed_n <- 20

log_lines <- character()
say <- function(...) { x <- sprintf(...); cat(x, "\n"); log_lines <<- c(log_lines, x) }
adj_r2 <- function(f) summary(f)$adj.r.squared

fam <- fread(file.path(qc_dir, "03-3_maf0.01_hwe.fam"), header = FALSE, select = 1:2, col.names = c("FID", "IID"))
anc <- fread(anc_file, select = c("dog", "pop", "pct"))[, `:=`(dog = as.character(dog), pct = as.numeric(pct))]
A <- dcast(anc[dog %in% as.character(fam$IID)], dog ~ pop, value.var = "pct", fill = 0)
A <- A[match(as.character(fam$IID), dog)]
stopifnot(!anyNA(A$dog))
M <- as.matrix(A[, -1])
say("개체 %d마리, 참조 집단 %d개", nrow(M), ncol(M))

info <- fread(file.path(qc_dir, "00_all_total.csv"), select = c("IID", "genetic_type", "top_pop"))[match(fam$IID, IID)]
stopifnot(identical(info$IID, fam$IID), !anyNA(info$genetic_type))
ph <- fread(file.path(qc_dir, "04-2_pheno.txt"), header = FALSE, col.names = c("FID", "IID", "rint", "raw"))
qc <- fread(file.path(qc_dir, "03-5_qcovar_PC5.txt"), header = FALSE)[, 1:7]
setnames(qc, c("FID", "IID", "age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age"))
cv <- fread(file.path(qc_dir, "03-5_covar.txt"), header = FALSE, col.names = c("FID", "IID", "Sex", "Intact", "Vision", "Cohort_BHS"))
stopifnot(identical(ph$IID, fam$IID), identical(qc$IID, fam$IID), identical(cv$IID, fam$IID))
X <- data.frame(qc[, -(1:2)], lapply(cv[, -(1:2)], factor))
Mdf <- as.data.frame(M[, colSums(M) > 0, drop = FALSE]); names(Mdf) <- make.names(paste0("anc_", names(Mdf)), unique = TRUE)

say("\n[Q1] 실제 CSLB의 ancestry 설명 분산")
y <- ph$rint
f_raw <- lm(y ~ ., data = cbind(data.frame(y = y), Mdf))
f0 <- lm(y ~ ., data = cbind(data.frame(y = y), X))
f1 <- lm(y ~ ., data = cbind(data.frame(y = y), X, Mdf))
q1_p <- anova(f0, f1)$`Pr(>F)`[2]
say("CSLB RINT ~ ancestry: 수정 R2 %.4f (공변량 없이)", adj_r2(f_raw))
say("공변량만: 수정 R2 %.4f / 공변량 + ancestry: %.4f → ancestry 추가 설명 %.4f (F검정 p = %.2g)",
    adj_r2(f0), adj_r2(f1), adj_r2(f1) - adj_r2(f0), q1_p)

icc_aov <- function(v, g) {
  g <- factor(g); k <- nlevels(g); N <- length(v); ni <- table(g)
  a <- anova(lm(v ~ g)); msb <- a$`Mean Sq`[1]; msw <- a$`Mean Sq`[2]
  n0 <- (N - sum(ni^2) / N) / (k - 1)
  max(0, (msb - msw) / (msb + (n0 - 1) * msw))
}
pl <- info$genetic_type == "purelike" & !is.na(info$top_pop)
bn <- table(info$top_pop[pl]); keep_b <- names(bn)[bn >= min_breed_n]
sel <- pl & info$top_pop %in% keep_b
res0 <- resid(f0)
icc_raw <- icc_aov(y[sel], info$top_pop[sel]); icc_adj <- icc_aov(res0[sel], info$top_pop[sel])
say("순종형 %d마리 중 %d마리 이상 품종 %d개 (%d마리): 품종 간 ICC 원 RINT %.4f, 공변량 보정 잔차 %.4f",
    sum(pl), min_breed_n, length(keep_b), sum(sel), icc_raw, icc_adj)

say("\n[Q2] 결정 규칙별 PC 개수")
st <- fread(file.path(qc_dir, "04-4_sim_lambda_summary.csv"))[order(n_pc)]
rule_min <- function(cut) if (any(st$lambda_mean <= cut)) min(st[lambda_mean <= cut, n_pc]) else NA_integer_
drop_next <- c(-diff(st$lambda_mean), NA)
elbow <- st$n_pc[which(drop_next < 0.05)[1]]
rules <- data.table(rule = c("lambda <= 1.05", "lambda <= 1.10", "elbow: next drop < 0.05"),
                    n_pc = c(rule_min(1.05), rule_min(1.10), elbow))
print(st[, .(n_pc, lambda_mean, lambda_sd, drop_to_next = round(drop_next, 3))])
for (i in seq_len(nrow(rules))) say("%s → PC %d", rules$rule[i], rules$n_pc[i])
fwrite(rules, paste0(out_prefix, "_threshold_rules.csv"))

say("\n[Q3] ancestry 비율 합과 구조 효과")
s_sum <- rowSums(M)
sum_tab <- data.table(type = info$genetic_type, s = s_sum)[, .(n = .N, median = round(median(s), 3), min = round(min(s), 3),
                                                              q05 = round(quantile(s, 0.05), 3)), by = type]
print(sum_tab)
say("비율 합 중앙값(최솟값): %s", paste(sprintf("%s %.3f (%.3f)", sum_tab$type, sum_tab$median, sum_tab$min), collapse = ", "))

sim <- fread(file.path(qc_dir, "04-2_sim_pheno.txt"), header = FALSE)
stopifnot(identical(as.character(sim$V2), as.character(fam$IID)))
Mn <- M / s_sum
set.seed(sim_seed)
is_p <- info$genetic_type == "purelike"
r3 <- rbindlist(lapply(seq_len(n_sim), function(s) {
  b <- rnorm(ncol(M))
  g <- as.vector(M %*% b); g <- (g - mean(g)) / sd(g) * sqrt(struct_var)
  e <- rnorm(nrow(M), 0, sqrt(1 - struct_var))
  rep_diff <- max(abs(round(g + e, 5) - sim[[s + 2]]))
  gn <- as.vector(Mn %*% b); gn <- (gn - mean(gn)) / sd(gn) * sqrt(struct_var)
  yy <- sim[[s + 2]]
  r2g <- function(idx) summary(lm(yy[idx] ~ M[idx, colSums(M[idx, ]) > 0]))$adj.r.squared
  data.table(sim = s, repro_maxdiff = rep_diff, cor_g_norm = cor(g, gn),
             r2_all = r2g(rep(TRUE, nrow(M))), r2_purelike = r2g(is_p), r2_mixlike = r2g(!is_p))
}))
fwrite(r3, paste0(out_prefix, "_sim_r2_by_group.csv"))
if (max(r3$repro_maxdiff) > 1e-4) warning("가상 표현형 재현 실패 (최대 차이 ", max(r3$repro_maxdiff), ") → 03-4 설정(seed 등) 확인")
say("가상 표현형 재현 최대 차이 %.1e (20개)", max(r3$repro_maxdiff))
say("원래 구조 효과 vs 합 1로 맞춘 구조 효과: 상관 평균 %.4f, 최소 %.4f", mean(r3$cor_g_norm), min(r3$cor_g_norm))
say("가상 표현형 안의 실제 구조 분산(ancestry 수정 R2) 평균: 전체 %.3f, 순종형 %.3f, 혼혈형 %.3f",
    mean(r3$r2_all), mean(r3$r2_purelike), mean(r3$r2_mixlike))

say("\n[Q4] 공변량별 ancestry 설명 분산")
cov_vars <- list(age_cslb = qc$age_cslb, pa_score = qc$pa_score, wt_cov = qc$wt_cov, owner_age = qc$owner_age)
for (v in c("Sex", "Intact", "Vision", "Cohort_BHS")) {
  lv <- sort(unique(cv[[v]]))
  for (l in lv[-1]) cov_vars[[paste0(v, "=", l)]] <- as.numeric(cv[[v]] == l)
}
r4 <- rbindlist(lapply(names(cov_vars), function(v) data.table(covariate = v, adj_r2 = round(adj_r2(lm(cov_vars[[v]] ~ ., data = Mdf)), 4))))
r4 <- r4[order(-adj_r2)]
print(r4)
fwrite(r4, paste0(out_prefix, "_covar_ancestry_r2.csv"))
say("%s", paste(sprintf("%s %.3f", r4$covariate, r4$adj_r2), collapse = ", "))

say("\n[Q6] 가상 표현형 20개의 λ 정밀도")
lam <- fread(file.path(qc_dir, "04-3_sim", "lambda.tsv"))[, .(n_pc, sim, lambda)]
w <- dcast(lam, sim ~ n_pc, value.var = "lambda")
r6 <- rbindlist(lapply(setdiff(names(w), "sim"), function(k) {
  d <- w[[k]] - w[["0"]]
  data.table(n_pc = as.integer(k), mean = mean(w[[k]]), sd = sd(w[[k]]), se = sd(w[[k]]) / sqrt(nrow(w)),
             diff_vs_pc0 = mean(d), se_diff = sd(d) / sqrt(nrow(w)))
}))
print(r6[, lapply(.SD, function(x) round(x, 4))])
fwrite(r6, paste0(out_prefix, "_lambda_precision.csv"))
say("GRM + PC λ SD 범위 %.3f ~ %.3f → 평균 λ SE %.4f ~ %.4f", min(r6$sd), max(r6$sd), min(r6$se), max(r6$se))
say("PC30 - PC0 짝지은 차이 %.4f (SE %.4f)", r6[n_pc == 30, diff_vs_pc0], r6[n_pc == 30, se_diff])
pc_only <- fread(file.path(qc_dir, "03-5_pca_lambda.csv"))[label %in% c("none", "pcair")]
say("참고: PC만 넣은 회귀(03-4) λ SD 범위 %.3f ~ %.3f", min(pc_only$lambda_sd), max(pc_only$lambda_sd))

writeLines(log_lines, paste0(out_prefix, "_sim_check.txt"))
