# 06-8_age8_h2.R - Heritability in dogs >= 8 years with mean scores, and genotype-by-age GREML. See README.md.
#         SurveyData/v2023/DAP_2023_CSLB_v1.0.csv · DogOverview/v2023/DAP_2023_DogOverview_v1.1.csv

library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_keep <- as.integer(Sys.getenv("TEST_KEEP", "0"))
tag <- if (test_keep > 0) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("06-8", tag))
pic_prefix <- file.path(pic_dir, paste0("06-8", tag))
work_dir <- paste0(out_prefix, "_work")
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "20"))
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)

grm_full <- file.path(in_dir, "04-1_grm", "grm_all")
pheno_main <- file.path(in_dir, "04-2_pheno.txt")
qcov_main <- file.path(in_dir, "03-5_qcovar_PC5.txt")
cov_main <- file.path(in_dir, "03-5_covar.txt")
all_total <- file.path(in_dir, "00_all_total.csv")
h2_66 <- file.path(in_dir, "06-6_h2.csv")
cslb_file <- file.path(base_dir, "SurveyData/v2023/DAP_2023_CSLB_v1.0.csv")
dog_file <- file.path(base_dir, "DogOverview/v2023/DAP_2023_DogOverview_v1.1.csv")
for (f in c(gcta, paste0(grm_full, c(".grm.bin", ".grm.id", ".grm.N.bin")), pheno_main, qcov_main, cov_main, all_total, cslb_file, dog_file))
  if (!file.exists(f)) stop("파일 없음: ", f)
qcov_names <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5))
age_cut <- 8

sh <- function(args, log) invisible(system2(gcta, args, stdout = log, stderr = log))
read_hsq <- function(f) {
  if (!file.exists(f)) return(NULL)
  l <- readLines(f); tb <- strsplit(l, "\t")
  comp <- rbindlist(lapply(tb[lengths(tb) == 3 & !grepl("^Source", l)], function(x) data.table(term = x[1], est = as.numeric(x[2]), se = as.numeric(x[3]))))
  kv <- function(k) { i <- grep(paste0("^", k, "\t"), l); if (!length(i)) NA_character_ else sub(paste0("^", k, "\t"), "", l[i[1]]) }
  list(comp = comp, LRT = as.numeric(sub(" .*", "", kv("LRT"))), Pval = as.numeric(sub(" .*", "", kv("Pval"))), n = as.integer(kv("n")))
}
rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
wr <- function(d, f) { fwrite(d, f, sep = " ", col.names = FALSE, quote = FALSE, na = "NA"); f }

ids <- fread(paste0(grm_full, ".grm.id"), header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
if (test_keep > 0) {
  kf <- wr(ids[seq(1, .N, by = test_keep)], file.path(work_dir, "test_keep.txt"))
  grm <- file.path(work_dir, "grm_test")
  if (!file.exists(paste0(grm, ".grm.bin"))) sh(c("--grm", grm_full, "--keep", kf, "--make-grm", "--out", grm, "--threads", n_cores), file.path(work_dir, "grm_test.stdout"))
  ids <- fread(paste0(grm, ".grm.id"), header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
  cat(sprintf("TEST_KEEP=%d: %d마리만 사용 (결과 해석 금지)\n", test_keep, nrow(ids)))
} else grm <- grm_full

ph <- fread(pheno_main, header = FALSE, colClasses = list(character = 1:2)); setnames(ph, c("FID", "IID", "rint", "score"))
qc <- fread(qcov_main, header = FALSE, colClasses = list(character = 1:2)); setnames(qc, c("FID", "IID", qcov_names))
cv <- fread(cov_main, header = FALSE, colClasses = "character"); setnames(cv, c("FID", "IID", "Sex", "Intact", "Vision", "Cohort_BHS"))
at <- fread(all_total, select = c("IID", "cslb_score", "age_cslb", "age_hinge", "n_cslb"), colClasses = c(IID = "character"))
at <- at[IID %in% ph$IID][!duplicated(IID)]
knot <- at[age_hinge > 0, median(age_cslb - age_hinge)]

cs <- fread(cslb_file, select = c("dog_id", "cslb_date", "cslb_score"), colClasses = c(dog_id = "character"))
dob <- fread(dog_file, select = c("dog_id", "Estimated_DOB"), colClasses = c(dog_id = "character"))
cs <- merge(cs[dog_id %in% ph$IID], dob, by = "dog_id")
cs[, age := as.numeric(as.IDate(cslb_date) - as.IDate(Estimated_DOB)) / 365.25]
cs <- cs[!is.na(cslb_score) & !is.na(age)]
setorder(cs, dog_id, cslb_date)
per <- cs[, .(n_all = .N, last = cslb_score[.N], age_last = age[.N], n8 = sum(age >= age_cut),
              mean8 = mean(cslb_score[age >= age_cut]), age_mean8 = mean(age[age >= age_cut]),
              hinge_mean8 = mean(pmax(0, age[age >= age_cut] - knot))), by = .(IID = dog_id)]
per[, hinge_last := pmax(0, age_last - knot)]
chk <- merge(at, per, by = "IID")
if (nrow(chk) != nrow(at) || any(chk$n_all != chk$n_cslb) || any(abs(chk$last - chk$cslb_score) > 1e-6) || any(abs(chk$age_last - chk$age_cslb) > 1e-6))
  stop("CSLB 원자료에서 다시 만든 응답 수·최근 점수·나이가 00_all_total.csv와 다름")
if (any(abs(chk$hinge_last - chk$age_hinge) > 1e-6)) stop("최근 나이의 hinge가 00_all_total.csv의 age_hinge와 다름 (기준점 확인)")
if (per[age_last < age_cut, any(n8 > 0)]) stop("최근 나이 < 8세인데 8세 이후 측정이 있는 개가 있음 (age8 정의 확인)")

d <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")), list(ph[, .(FID, IID, rint, score)], qc, cv))
d <- merge(d, per, by = "IID")[IID %in% ids$IID]
d[, old := age_last >= age_cut]
a8 <- d[old == TRUE]
cat(sprintf("age8 %d마리 (8세 이후 측정 1회 %d · 2회 %d · 3회 %d · 4회 %d), 나이 기준점 %.1f세\n", nrow(a8),
            sum(a8$n8 == 1), sum(a8$n8 == 2), sum(a8$n8 == 3), sum(a8$n8 == 4), knot))
cat(sprintf("  → 8세 이후 2회 이상 측정 (m2) %d마리. 11세를 걸치는 개 (mean8에서 hinge 평균 ≠ 평균 나이의 hinge) %d마리\n",
            sum(a8$n8 >= 2), sum(abs(a8$hinge_mean8 - pmax(0, a8$age_mean8 - knot)) > 1e-6)))

pheno_info <- rbindlist(lapply(c("all", "m2"), function(s) {
  x <- if (s == "all") a8 else a8[n8 >= 2]
  data.table(subset = s, n = nrow(x), n8_mean = mean(x$n8), age_last_mean = mean(x$age_last), age_mean8_mean = mean(x$age_mean8),
             score_last_mean = mean(x$last), score_last_sd = sd(x$last), score_mean8_mean = mean(x$mean8), score_mean8_sd = sd(x$mean8),
             r_last_mean8 = cor(x$last, x$mean8), pct_ge50_last = 100 * mean(x$last >= 50), pct_ge50_mean8 = 100 * mean(x$mean8 >= 50))
}))

make_inputs <- function(name, x, y, age, hinge) {
  x <- copy(x)[, `:=`(y = rint(get(y)), age_q = get(age), hinge_q = get(hinge))]
  x[, `:=`(age_cslb = age_q, age_hinge = hinge_q)]
  keepc <- c("Sex", "Intact", "Vision", "Cohort_BHS", "n8")
  keepc <- keepc[vapply(keepc, function(j) uniqueN(x[[j]]) > 1, TRUE)]
  keepq <- qcov_names[vapply(qcov_names, function(j) isTRUE(sd(x[[j]]) > 0), TRUE)]
  list(pheno = wr(x[, .(FID, IID, y)], file.path(work_dir, paste0("pheno_", name, ".txt"))),
       qcovar = wr(x[, c("FID", "IID", keepq), with = FALSE], file.path(work_dir, paste0("qcovar_", name, ".txt"))),
       covar = wr(x[, c("FID", "IID", keepc), with = FALSE], file.path(work_dir, paste0("covar_", name, ".txt"))),
       keep = wr(x[, .(FID, IID)], file.path(work_dir, paste0("keep_", name, ".txt"))), n = nrow(x))
}
h2_models <- list(
  h2_last8    = c(make_inputs("last8", a8, "last", "age_last", "hinge_last"), desc = "age8: latest CSLB (same y as 06-6 h2_age8, + n8 covariate)"),
  h2_mean8    = c(make_inputs("mean8", a8, "mean8", "age_mean8", "hinge_mean8"), desc = "age8: mean of CSLB at >= 8 y"),
  h2_last8_m2 = c(make_inputs("last8_m2", a8[n8 >= 2], "last", "age_last", "hinge_last"), desc = "age8, >= 2 measurements at >= 8 y: latest CSLB"),
  h2_mean8_m2 = c(make_inputs("mean8_m2", a8[n8 >= 2], "mean8", "age_mean8", "hinge_mean8"), desc = "age8, >= 2 measurements at >= 8 y: mean CSLB"))

env_f <- wr(d[, .(FID, IID, fifelse(old, "age8", "young8"))], file.path(work_dir, "env_age8.txt"))
gx_q <- wr(d[, c("FID", "IID", qcov_names), with = FALSE], file.path(work_dir, "qcovar_gxe.txt"))
gx_c <- wr(d[, .(FID, IID, Sex, Intact, Vision, Cohort_BHS)], file.path(work_dir, "covar_gxe.txt"))
d[, y_within := rint(score), by = old]
gxe_models <- list(
  gxe_main   = list(pheno = wr(d[, .(FID, IID, rint)], file.path(work_dir, "pheno_gxe_main.txt")), desc = "all dogs, main RINT"),
  gxe_within = list(pheno = wr(d[, .(FID, IID, y_within)], file.path(work_dir, "pheno_gxe_within.txt")), desc = "all dogs, RINT within age group"))

run_gcta <- function(name, args) {
  out <- file.path(work_dir, name)
  if (file.exists(paste0(out, ".hsq"))) return(invisible(0))
  if (file.exists(paste0(out, ".log")) && any(grepl("not converged", readLines(paste0(out, ".log"), warn = FALSE)))) {
    cat(sprintf("  %s: 이전 실행에서 200회 안에 수렴하지 않음 → 다시 돌리지 않음\n", name)); return(invisible(1)) }
  cat("GREML:", name, "\n")
  sh(c(args, "--reml-maxit", "200", "--out", out, "--threads", n_cores), paste0(out, ".stdout"))
}
t0 <- Sys.time()
for (nm in names(h2_models)) { m <- h2_models[[nm]]
  run_gcta(nm, c("--reml", "--grm", grm, "--pheno", m$pheno, "--qcovar", m$qcovar, "--covar", m$covar, "--keep", m$keep)) }
for (nm in names(gxe_models)) { m <- gxe_models[[nm]]
  run_gcta(nm, c("--reml", "--grm", grm, "--pheno", m$pheno, "--gxe", env_f, "--qcovar", gx_q, "--covar", gx_c, "--reml-lrt", "2")) }
cat(sprintf("GREML 끝 (%.1f분)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

h2 <- rbindlist(lapply(names(h2_models), function(nm) {
  r <- read_hsq(file.path(work_dir, paste0(nm, ".hsq"))); if (is.null(r)) return(NULL)
  g <- r$comp[term == "V(G)/Vp"]
  data.table(source = "06-8", model = nm, desc = h2_models[[nm]]$desc, n = r$n, h2 = g$est, se = g$se, p_LRT = r$Pval)
}), fill = TRUE)
if (file.exists(h2_66) && test_keep == 0) {
  ref <- fread(h2_66)[model %in% c("h2_main", "h2_mean", "h2_age8", "h2_young8")]
  h2 <- rbind(ref[, .(source = "06-6", model, desc, n, h2, se, p_LRT)], h2, fill = TRUE)
}
if (nrow(h2)) h2[, `:=`(lo95 = pmax(0, h2 - 1.96 * se), hi95 = h2 + 1.96 * se)]

gxe <- rbindlist(lapply(names(gxe_models), function(nm) {
  r <- read_hsq(file.path(work_dir, paste0(nm, ".hsq"))); if (is.null(r)) return(NULL)
  g <- function(t) { x <- r$comp[term == t]; if (nrow(x)) x else data.table(est = NA_real_, se = NA_real_) }
  vg <- g("V(G)")$est; vgxe <- g("V(GxE)")$est
  data.table(model = nm, desc = gxe_models[[nm]]$desc, n = r$n,
             h2_G = g("V(G)/Vp")$est, se_G = g("V(G)/Vp")$se, h2_GxE = g("V(GxE)/Vp")$est, se_GxE = g("V(GxE)/Vp")$se,
             h2_sum = g("Sum of V(G)/Vp")$est, se_sum = g("Sum of V(G)/Vp")$se, LRT_GxE = r$LRT, p_GxE = r$Pval,
             rg_approx = if (isTRUE(vg > 1e-4 && vgxe > 1e-4)) vg / (vg + vgxe) else NA_real_)
}), fill = TRUE)

get_h2 <- function(nm, col) { x <- h2[source == "06-8" & model == nm, get(col)]; if (length(x)) x else NA_real_ }
get_gx <- function(nm, col) { x <- if (nrow(gxe)) gxe[model == nm, get(col)] else numeric(0); if (length(x)) x else NA_real_ }
D1 <- isTRUE(get_h2("h2_mean8", "p_LRT") < 0.05)
D2 <- isTRUE(get_gx("gxe_main", "p_GxE") < 0.05 && get_gx("gxe_within", "p_GxE") < 0.05)
missing_models <- c(names(h2_models), names(gxe_models))[!file.exists(file.path(work_dir, paste0(c(names(h2_models), names(gxe_models)), ".hsq")))]

fwrite(h2, paste0(out_prefix, "_h2.csv"))
fwrite(gxe, paste0(out_prefix, "_gxe.csv"))
fwrite(pheno_info, paste0(out_prefix, "_pheno_info.csv"))

pd <- rbind(h2[, .(panel = "h2 by phenotype", label = paste0(model, " (", source, ")"), est = h2, lo = lo95, hi = hi95, p = p_LRT)],
            if (nrow(gxe)) rbind(gxe[, .(panel = "G x age (all dogs)", label = paste0(model, ": V(G)/Vp"), est = h2_G, lo = pmax(0, h2_G - 1.96 * se_G), hi = h2_G + 1.96 * se_G, p = NA_real_)],
                                 gxe[, .(panel = "G x age (all dogs)", label = paste0(model, ": V(GxE)/Vp"), est = h2_GxE, lo = pmax(0, h2_GxE - 1.96 * se_GxE), hi = h2_GxE + 1.96 * se_GxE, p = p_GxE)]))
if (nrow(pd)) {
  lv <- rev(unique(pd$label)); pd[, `:=`(label = factor(label, levels = lv), panel = factor(panel, levels = c("h2 by phenotype", "G x age (all dogs)")))]
  g1 <- ggplot(pd, aes(est, label)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#1f4e9c", size = 0.4) +
    geom_text(aes(label = ifelse(is.na(p), sprintf("%.3f", est), sprintf("%.3f (p %.2g)", est, p))), nudge_y = 0.3, size = 3, colour = "grey30") +
    facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.05))) +
    labs(x = "Variance proportion with 95% CI (lower bound truncated at 0)", y = NULL, title = "CSLB SNP heritability in older dogs (GREML)",
         subtitle = "p = one-sided LRT (V(G) = 0 for h2 rows; V(GxE) = 0 for GxE rows)\nGxE rows: gxe_main and gxe_within use different y scales; compare only whether V(GxE) > 0") +
    theme_minimal(base_size = 11)
  ggsave(paste0(pic_prefix, "_h2.png"), g1, width = 8.5, height = 6, dpi = 300)
}

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("06-8 나이 든 개의 SNP 유전력 (반복 측정 평균) · 전체 표본 G×나이 (GCTA GREML, GRM = 04-1 LD pruning SNP)\n")
if (test_keep > 0) cat(sprintf("*** TEST_KEEP=%d: %d마리만 쓴 시험 결과. 해석하지 말 것 ***\n", test_keep, nrow(ids)))
cat(sprintf("나이 기준 %d세 (최근 CSLB 나이), 나이 hinge 기준점 %.1f세. 탐색 분석 (06 결과를 본 뒤 추가; 판정 규칙은 실행 전 헤더에 적음)\n\n", age_cut, knot))
cat("■ 노령견 표현형 (최근 1회 vs 8세 이후 평균)\n"); print(f3(pheno_info))
cat("\n■ h² (V(G)/Vp ± SE, LRT 한쪽 p). 06-6 줄은 비교용\n")
print(f3(h2[, .(source, model, n, h2, se, lo95, hi95, p_LRT, desc)]))
if (nrow(h2[source == "06-8"]) == 4)
  cat(sprintf("  mean8 / last8 = %.2f (전체 age8), %.2f (2회 이상 측정)\n",
              get_h2("h2_mean8", "h2") / get_h2("h2_last8", "h2"), get_h2("h2_mean8_m2", "h2") / get_h2("h2_last8_m2", "h2")))
cat("  · 06-6 h2_age8과 06-8 h2_last8은 같은 개·같은 y이고 공변량 n8 (8세 이후 측정 횟수) 하나만 다름 → 두 값의 차이는 n8 보정의 영향\n",
    "  · mean8과 last8은 같은 개라 두 추정치가 상관되어 있고 차이의 SE가 없음 → 비는 기술용이며 \"나아졌다\"는 판정에 쓰지 않음\n",
    "  · mean8 공변량: 평균 나이 + 측정별 hinge의 평균 (측정 궤적 보정). 그래도 평균이 개마다 다른 나이 구간을 섞는 한계는 남음\n", sep = "")
cat("\n■ G×나이 (전체 표본; V(G) = 나이와 상관없이 공유되는 유전 몫, V(G×E) = 나이 집단에 따라 달라지는 몫; p는 V(G×E) = 0 LRT)\n")
if (nrow(gxe)) print(f3(gxe)) else cat("  결과 없음\n")
cat("  rg_approx = V(G) / (V(G) + V(G×E)): 두 나이 집단의 유전 분산이 같다는 가정에서만 rg와 같음. 노령견 유전 분산이 더 크면\n",
    "  rg = 1이어도 V(G×E) > 0이 생겨 rg_approx가 실제보다 작게 나옴 (참고값). 두 성분이 모두 0보다 클 때만 계산\n",
    "  같은 이유로 V(G×E) > 0은 \"나이에 따라 다른 유전자\"와 \"나이에 따라 커지는 같은 유전 효과\"를 구분하지 못함\n",
    "  gxe_main (전체 RINT)과 gxe_within (집단 안 RINT)은 y 척도가 달라 V(G)/Vp 크기를 서로 비교하지 않음. 비교 대상은 V(G×E) 유무\n", sep = "")
cat("\n■ 판정 (헤더의 규칙)\n")
cat(sprintf("  D1 노령견 반복 평균으로 후속 세트 분석: %s (h2_mean8 LRT p %.3g; 참고: mean8 %.4f, last8 %.4f)\n", if (D1) "탐색적 근거 있음" else "근거 없음",
            get_h2("h2_mean8", "p_LRT"), get_h2("h2_mean8", "h2"), get_h2("h2_last8", "h2")))
cat(sprintf("  D2 나이에 따른 유전 효과 차이: %s (V(G×E) LRT p: gxe_main %.3g, gxe_within %.3g)\n", if (D2) "탐색적 근거 있음" else "근거 없음",
            get_gx("gxe_main", "p_GxE"), get_gx("gxe_within", "p_GxE")))
cat("  \"탐색적 근거\"여도 06 결과를 본 뒤 추가한 분석이므로 독립 표본 확인 전에는 \"효과가 확인됐다\"로 쓰지 않음\n")
cat("  D3 전체 표본 주 표현형 교체: 06-8로 정하지 않음 (헤더 참고)\n")
cat("  시도한 표현형·모형 수: 06-6 9개 + 06-8 h² 4개 + G×E 2개 (보고 시 함께 적음)\n")
cat(sprintf("  결과 없는 모형 (수렴 실패 등; work의 .log 확인): %s\n", if (length(missing_models)) paste(missing_models, collapse = ", ") else "없음"))
sink()
cat(sprintf("저장: %s_h2.csv · _gxe.csv · _pheno_info.csv · _summary.txt / 그림 %s_h2.png\n", out_prefix, pic_prefix))

#    h2_mean8_m2   1,609  0.014    0.049   0.39
#    gxe_main     0.013 ± 0.016    0.006 ± 0.013    0.277          0.020      0.67
#    gxe_within   0.014 ± 0.016    0.006 ± 0.013    0.282          0.020      0.69
