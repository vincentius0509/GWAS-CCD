# 06-6_h2_rg_power.R - SNP heritability (GREML), genetic correlations between strata, and power. See README.md.
#   sex: sex_m (3,657) vs sex_f (3,458)
#   06-6_h2.csv · _rg.csv · _power_h2.csv · _power_snp.csv · _power_sets.csv · _summary.txt

library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_keep <- as.integer(Sys.getenv("TEST_KEEP", "0"))
tag <- if (test_keep > 0) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("06-6", tag))
pic_prefix <- file.path(pic_dir, paste0("06-6", tag))
work_dir <- paste0(out_prefix, "_work")
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "20"))
parts <- strsplit(Sys.getenv("PARTS", "h2 rg power"), " ")[[1]]
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)

grm_full <- file.path(in_dir, "04-1_grm", "grm_all")
pheno_main <- file.path(in_dir, "04-2_pheno.txt")
qcov_main <- file.path(in_dir, "03-5_qcovar_PC5.txt")
cov_main <- file.path(in_dir, "03-5_covar.txt")
all_total <- file.path(in_dir, "00_all_total.csv")
runs_tsv <- file.path(in_dir, "05-1_runs.tsv")
s1_file <- file.path(in_dir, "06-5_stage1.csv"); s2_file <- file.path(in_dir, "06-5_stage2.csv")
for (f in c(gcta, paste0(grm_full, c(".grm.bin", ".grm.id", ".grm.N.bin")), pheno_main, qcov_main, cov_main, all_total, runs_tsv))
  if (!file.exists(f)) stop("파일 없음: ", f)
qcov_names <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5))
bk_cut <- 0.05
strata <- c("age8", "young8", "sex_m", "sex_f")
pairs <- list(age = c("young8", "age8"), sex = c("sex_m", "sex_f"))

sh <- function(args, log) invisible(system2(gcta, args, stdout = log, stderr = log))
read_hsq <- function(f) {
  if (!file.exists(f)) return(NULL)
  l <- readLines(f); tb <- strsplit(l, "\t")
  comp <- rbindlist(lapply(tb[lengths(tb) == 3 & !grepl("^Source", l)], function(x) data.table(term = x[1], est = as.numeric(x[2]), se = as.numeric(x[3]))))
  kv <- function(k) { i <- grep(paste0("^", k, "\t"), l); if (!length(i)) NA_character_ else sub(paste0("^", k, "\t"), "", l[i[1]]) }
  lrts <- grep("^(logL0|LRT|Pval)\t", l, value = TRUE)
  list(comp = comp, logL = as.numeric(kv("logL")), LRT = as.numeric(sub(" .*", "", kv("LRT"))), Pval = as.numeric(sub(" .*", "", kv("Pval"))),
       n = as.integer(kv("n")), all_p = as.numeric(sub(" .*", "", sub("^Pval\t", "", grep("^Pval\t", l, value = TRUE)))))
}

ids <- fread(paste0(grm_full, ".grm.id"), header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
if (test_keep > 0) {
  kf <- file.path(work_dir, "test_keep.txt"); fwrite(ids[seq(1, .N, by = test_keep)], kf, sep = " ", col.names = FALSE, quote = FALSE)
  grm <- file.path(work_dir, "grm_test")
  if (!file.exists(paste0(grm, ".grm.bin"))) sh(c("--grm", grm_full, "--keep", kf, "--make-grm", "--out", grm, "--threads", n_cores), file.path(work_dir, "grm_test.stdout"))
  ids <- fread(paste0(grm, ".grm.id"), header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
  cat(sprintf("TEST_KEEP=%d: %d마리만 사용 (결과 해석 금지)\n", test_keep, nrow(ids)))
} else grm <- grm_full
grm_bk <- paste0(grm, "_bK", bk_cut)
if (!file.exists(paste0(grm_bk, ".grm.bin")))
  sh(c("--grm", grm, "--make-bK", bk_cut, "--out", grm_bk, "--threads", n_cores), paste0(grm_bk, ".stdout"))

qc <- fread(qcov_main, header = FALSE, colClasses = list(character = 1:2)); setnames(qc, c("FID", "IID", qcov_names))
cv <- fread(cov_main, header = FALSE, colClasses = "character"); setnames(cv, c("FID", "IID", "Sex", "Intact", "Vision", "Cohort_BHS"))
ph <- fread(pheno_main, header = FALSE, colClasses = list(character = 1:2)); setnames(ph, c("FID", "IID", "rint", "score"))
qc_nopc <- file.path(work_dir, "qcovar_nopc.txt"); fwrite(qc[, c("FID", "IID", qcov_names[1:5]), with = FALSE], qc_nopc, sep = " ", col.names = FALSE, quote = FALSE)
rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
at <- fread(all_total, select = c("IID", "cslb_score", "cslb_score_mean", "age_cslb", "age_cslb_mean", "age_hinge", "n_cslb"), colClasses = c(IID = "character"))
at <- at[IID %in% ph$IID][!duplicated(IID)]
knot <- at[age_hinge > 0, median(age_cslb - age_hinge)]
mn <- merge(ph[, .(FID, IID)], at, by = "IID", sort = FALSE)[match(ph$IID, IID)]
if (anyNA(mn$cslb_score_mean)) stop("반복 측정 평균이 없는 개: ", sum(is.na(mn$cslb_score_mean)))
pheno_mean <- file.path(work_dir, "pheno_mean.txt")
fwrite(mn[, .(FID, IID, round(rint(cslb_score_mean), 6))], pheno_mean, sep = " ", col.names = FALSE, quote = FALSE)
qc_mean <- copy(qc)[, `:=`(age_cslb = mn$age_cslb_mean, age_hinge = pmax(0, mn$age_cslb_mean - knot))]
qcov_mean <- file.path(work_dir, "qcovar_mean.txt"); fwrite(qc_mean, qcov_mean, sep = " ", col.names = FALSE, quote = FALSE)
cat(sprintf("반복 측정: 1회 %d · 2회 %d · 3회 %d · 4회 %d마리 (나이 hinge 기준점 %.1f세)\n",
            sum(mn$n_cslb == 1), sum(mn$n_cslb == 2), sum(mn$n_cslb == 3), sum(mn$n_cslb == 4), knot))

runs <- fread(runs_tsv)
row_of <- function(rid) { r <- runs[id == rid]; if (nrow(r) != 1) stop("05-1_runs.tsv에 없음: ", rid); r }

# 2. h² (GREML)
h2_models <- list(
  h2_main  = list(grm = grm, pheno = pheno_main, mpheno = 1, qcovar = qcov_main, covar = cov_main, desc = "main: RINT, PC5"),
  h2_nopc  = list(grm = grm, pheno = pheno_main, mpheno = 1, qcovar = qc_nopc, covar = cov_main, desc = "no PCs"),
  h2_2grm  = list(mgrm = c(grm, grm_bk), pheno = pheno_main, mpheno = 1, qcovar = qcov_main, covar = cov_main, desc = "GRM + close-relative GRM (Zaitlen 2013)"),
  h2_score = list(grm = grm, pheno = pheno_main, mpheno = 2, qcovar = qcov_main, covar = cov_main, desc = "raw score"),
  h2_mean  = list(grm = grm, pheno = pheno_mean, mpheno = 1, qcovar = qcov_mean, covar = cov_main, desc = "mean of repeated CSLB (RINT)"))
clean_covar <- function(s, r) {
  k <- fread(r$covar_file, header = FALSE, colClasses = "character")
  kp <- fread(r$keep_file, header = FALSE, colClasses = "character")[[2]]
  k <- k[V2 %in% intersect(kp, ids$IID)]
  keepc <- c(TRUE, TRUE, vapply(k[, -(1:2)], function(x) uniqueN(x) > 1, TRUE))
  if (!all(keepc)) cat(sprintf("  %s: 범주 하나뿐인 공변량 %d개 뺌\n", s, sum(!keepc)))
  f <- file.path(work_dir, paste0("covar_", s, ".txt")); fwrite(k[, which(keepc), with = FALSE], f, sep = " ", col.names = FALSE, quote = FALSE); f
}
clean_qcovar <- function(s, r) {
  q <- fread(r$qcovar_file, header = FALSE, colClasses = list(character = 1:2))
  kp <- fread(r$keep_file, header = FALSE, colClasses = "character")[[2]]
  q <- q[V2 %in% intersect(kp, ids$IID)]
  keepc <- c(TRUE, TRUE, vapply(q[, -(1:2)], function(x) isTRUE(sd(x) > 0), TRUE))
  if (!all(keepc)) cat(sprintf("  %s: 상수인 양적 공변량 %d개 뺌\n", s, sum(!keepc)))
  f <- file.path(work_dir, paste0("qcovar_", s, ".txt")); fwrite(q[, which(keepc), with = FALSE], f, sep = " ", col.names = FALSE, quote = FALSE); f
}
for (s in strata) { r <- row_of(s)
  h2_models[[paste0("h2_", s)]] <- list(grm = grm, pheno = r$pheno_file, mpheno = 1, qcovar = clean_qcovar(s, r), covar = clean_covar(s, r), keep = r$keep_file, desc = paste("stratum", s)) }

run_reml <- function(name, m) {
  out <- file.path(work_dir, name)
  if (file.exists(paste0(out, ".hsq"))) return(invisible(0))
  g <- if (!is.null(m$mgrm)) { mf <- paste0(out, "_mgrm.txt"); writeLines(m$mgrm, mf); c("--mgrm", mf, "--reml-no-lrt") } else c("--grm", m$grm)
  args <- c("--reml", g, "--pheno", m$pheno, "--mpheno", m$mpheno, "--qcovar", m$qcovar, "--covar", m$covar,
            if (!is.null(m$keep)) c("--keep", m$keep), "--reml-maxit", "200", "--out", out, "--threads", n_cores)
  sh(args, paste0(out, ".stdout"))
}
if ("h2" %in% parts) {
  t0 <- Sys.time()
  for (nm in names(h2_models)) { cat("GREML:", nm, "\n"); run_reml(nm, h2_models[[nm]]) }
  cat(sprintf("GREML 끝 (%.1f분)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  miss <- names(h2_models)[!file.exists(file.path(work_dir, paste0(names(h2_models), ".hsq")))]
  if (length(miss)) warning("GREML 결과 없음 (수렴 실패 등; work의 .log 확인): ", paste(miss, collapse = ", "))
}
h2 <- rbindlist(lapply(names(h2_models), function(nm) {
  r <- read_hsq(file.path(work_dir, paste0(nm, ".hsq"))); if (is.null(r)) return(NULL)
  g <- function(t) r$comp[term == t]
  two <- nrow(g("V(G1)/Vp")) > 0
  hh <- if (two) g("V(G1)/Vp") else g("V(G)/Vp")
  data.table(model = nm, desc = h2_models[[nm]]$desc, n = r$n, h2 = hh$est, se = hh$se,
             h2_rel = if (two) g("V(G2)/Vp")$est else NA_real_, se_rel = if (two) g("V(G2)/Vp")$se else NA_real_,
             h2_sum = if (two) g("Sum of V(G)/Vp")$est else NA_real_, Vp = g("Vp")$est, LRT = r$LRT, p_LRT = if (two) pnorm(hh$est / hh$se, lower.tail = FALSE) else r$Pval,
             p_type = if (two) "Wald (one-sided)" else "LRT")
}))
if (nrow(h2)) h2[, `:=`(lo95 = pmax(0, h2 - 1.96 * se), hi95 = h2 + 1.96 * se)]

resid_of <- function(s) {
  r <- row_of(s)
  p <- fread(r$pheno_file, header = FALSE, colClasses = list(character = 1:2)); setnames(p, c("FID", "IID", "y"))
  q <- fread(r$qcovar_file, header = FALSE, colClasses = list(character = 1:2)); setnames(q, c("FID", "IID", paste0("q", seq_len(ncol(q) - 2))))
  k <- fread(r$covar_file, header = FALSE, colClasses = "character"); setnames(k, c("FID", "IID", paste0("c", seq_len(ncol(k) - 2))))
  d <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")), list(p, q, k))
  cn <- grep("^c[0-9]", names(d), value = TRUE); for (j in cn) set(d, j = j, value = factor(d[[j]]))
  cn <- cn[vapply(cn, function(j) nlevels(d[[j]]) > 1, TRUE)]
  f <- as.formula(paste("y ~", paste(c(grep("^q[0-9]", names(d), value = TRUE), cn), collapse = " + ")))
  d[, res := resid(lm(f, data = d))]
  d[, .(FID, IID, res)]
}
rg_tab <- NULL
if ("rg" %in% parts || file.exists(file.path(work_dir, "rg_age_0.hsq"))) {
  rg_tab <- rbindlist(lapply(names(pairs), function(pn) {
    a <- resid_of(pairs[[pn]][1]); b <- resid_of(pairs[[pn]][2])
    if (length(intersect(a$IID, b$IID))) stop(pn, ": 두 층에 같은 개가 있음 (겹치지 않는 층이어야 함)")
    pf <- file.path(work_dir, paste0("rg_", pn, "_pheno.txt"))
    d <- merge(ids, rbind(a[, .(FID, IID, t1 = res, t2 = NA_real_)], b[, .(FID, IID, t1 = NA_real_, t2 = res)]), by = c("FID", "IID"))
    fwrite(d, pf, sep = " ", col.names = FALSE, na = "NA", quote = FALSE)
    for (v in c("0", "1")) {
      out <- file.path(work_dir, paste0("rg_", pn, "_", v))
      failed <- file.exists(paste0(out, ".log")) && any(grepl("not converged", readLines(paste0(out, ".log"), warn = FALSE)))
      if (failed && !file.exists(paste0(out, ".hsq"))) cat(sprintf("이변량 GREML: %s rg = %s 검정은 이전 실행에서 200회 안에 수렴하지 않음 → 다시 돌리지 않음 (p = NA)\n", pn, v))
      if ("rg" %in% parts && !file.exists(paste0(out, ".hsq")) && !failed) {
        cat(sprintf("이변량 GREML: %s (%s vs %s), rg = %s 검정\n", pn, pairs[[pn]][1], pairs[[pn]][2], v))
        sh(c("--reml-bivar", "1", "2", "--reml-bivar-nocove", "--reml-bivar-lrt-rg", v, "--grm", grm, "--pheno", pf,
             "--reml-maxit", "200", "--out", out, "--threads", n_cores), paste0(out, ".stdout"))
      }
    }
    r0 <- read_hsq(file.path(work_dir, paste0("rg_", pn, "_0.hsq"))); r1 <- read_hsq(file.path(work_dir, paste0("rg_", pn, "_1.hsq")))
    if (is.null(r0)) return(NULL)
    g <- function(t) r0$comp[term == t]
    data.table(pair = pn, tr1 = pairs[[pn]][1], tr2 = pairs[[pn]][2], n = r0$n,
               h2_tr1 = g("V(G)/Vp_tr1")$est, se_tr1 = g("V(G)/Vp_tr1")$se, h2_tr2 = g("V(G)/Vp_tr2")$est, se_tr2 = g("V(G)/Vp_tr2")$se,
               rg = g("rG")$est, se_rg = g("rG")$se, p_rg0 = r0$Pval, p_rg1 = if (is.null(r1)) NA_real_ else r1$Pval)
  }))
  if (nrow(rg_tab)) {
    uh <- function(s, col) { x <- h2[model == paste0("h2_", s), get(col)]; if (length(x)) x else NA_real_ }
    rg_tab[, `:=`(h2_uni_tr1 = mapply(uh, tr1, "h2"), se_uni_tr1 = mapply(uh, tr1, "se"), h2_uni_tr2 = mapply(uh, tr2, "h2"), se_uni_tr2 = mapply(uh, tr2, "se"))]
    rg_tab[, estimable := h2_uni_tr1 > 2 * se_uni_tr1 & h2_uni_tr2 > 2 * se_uni_tr2]
    rg_tab[is.na(estimable), estimable := h2_tr1 > 2 * se_tr1 & h2_tr2 > 2 * se_tr2]
    rg_tab[, `:=`(rg_raw = rg, se_rg_raw = se_rg, p_rg0_raw = p_rg0, p_rg1_raw = p_rg1)]
    rg_tab[estimable == FALSE, `:=`(rg = NA_real_, se_rg = NA_real_, p_rg0 = NA_real_, p_rg1 = NA_real_)]
  }
}

pw_h2 <- pw_snp <- pw_sets <- NULL
if ("power" %in% parts) {
  n <- nrow(ids)
  nb <- n * (n + 1) / 2
  con <- file(paste0(grm, ".grm.bin"), "rb"); A <- readBin(con, "numeric", n = nb, size = 4); close(con)
  diag_pos <- cumsum(seq_len(n))
  offd <- A[-diag_pos]; rm(A); gc()
  var_off <- var(offd); mean_off <- mean(offd); rm(offd)
  se_theory <- sqrt(2 / (n^2 * var_off))
  se_obs <- if (nrow(h2) && "h2_main" %in% h2$model) h2[model == "h2_main", se] else NA_real_
  pw_h2 <- data.table(h2_true = c(0.02, 0.05, 0.10, 0.15, 0.20))
  pw_h2[, `:=`(power_theorySE = pnorm(h2_true / se_theory - qnorm(0.95)), power_obsSE = pnorm(h2_true / se_obs - qnorm(0.95)))]
  attr(pw_h2, "info") <- sprintf("N %d, GRM 비대각 평균 %.4f · 분산 %.2e → 이론 SE(h²) %.3f (사람 비혈연 약 316/N = %.3f), 관측 SE %.3f",
                                 n, mean_off, var_off, se_theory, 316 / n, se_obs)
  n_snp <- if (test_keep > 0) n * test_keep else n
  q2 <- c(0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01)
  alphas <- c(`5e-8` = 5e-8, `4.2e-7` = 4.2e-7, `1e-5` = 1e-5, `1e-4` = 1e-4)
  pw_snp <- CJ(q2 = q2, alpha = names(alphas))
  pw_snp[, power := pchisq(qchisq(1 - alphas[alpha], 1), 1, ncp = n_snp * q2 / (1 - q2), lower.tail = FALSE)]
  q2_80 <- vapply(alphas, function(a) uniroot(function(x) pchisq(qchisq(1 - a, 1), 1, ncp = n_snp * x / (1 - x), lower.tail = FALSE) - 0.8, c(1e-6, 0.2))$root, 0)
  attr(pw_snp, "q2_80") <- q2_80
  if (file.exists(s1_file)) {
    s1 <- fread(s1_file); s2 <- if (file.exists(s2_file)) fread(s2_file) else NULL
    z80 <- qnorm(0.8)
    pw_sets <- rbind(s1[, .(stage = "1 (main)", run = "main", set, NGENES, BETA, SE)], if (!is.null(s2)) s2[, .(stage = "2 (age)", run, set, NGENES, BETA, SE)])
    pw_sets[, m_family := fifelse(stage == "1 (main)", 14, 28)]
    pw_sets[, `:=`(mde_rank1 = (qnorm(1 - 0.05 / m_family) + z80) * SE, mde_nominal = (qnorm(0.95) + z80) * SE)]
  }
}

if (nrow(h2)) fwrite(h2, paste0(out_prefix, "_h2.csv"))
if (!is.null(rg_tab) && nrow(rg_tab)) fwrite(rg_tab, paste0(out_prefix, "_rg.csv"))
if (!is.null(pw_h2)) fwrite(pw_h2, paste0(out_prefix, "_power_h2.csv"))
if (!is.null(pw_snp)) fwrite(pw_snp, paste0(out_prefix, "_power_snp.csv"))
if (!is.null(pw_sets)) fwrite(pw_sets, paste0(out_prefix, "_power_sets.csv"))

if (nrow(h2)) {
  lv <- rev(h2$model); hp <- copy(h2)[, model := factor(model, levels = lv)]
  g1 <- ggplot(hp, aes(h2, model)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_pointrange(aes(xmin = lo95, xmax = hi95), colour = "#1f4e9c", size = 0.4) +
    geom_text(aes(label = sprintf("%.3f (p %.2g)", h2, p_LRT)), nudge_y = 0.3, size = 3, colour = "grey30") +
    labs(x = "SNP heritability (V(G)/Vp) with 95% CI (lower bound truncated at 0)", y = NULL, title = "CSLB SNP heritability (GREML)",
         subtitle = "h2_2grm shows the SNP component (V(G1)/Vp); p = one-sided test of V(G) = 0 (LRT; Wald for h2_2grm)") +
    theme_minimal(base_size = 11)
  ggsave(paste0(pic_prefix, "_h2.png"), g1, width = 8, height = 5, dpi = 300)
}
if (!is.null(pw_snp)) {
  g2 <- ggplot(pw_snp, aes(q2 * 100, power, colour = alpha)) + geom_line(linewidth = 1) + geom_point(size = 1.8) +
    geom_hline(yintercept = 0.8, linetype = 2, colour = "grey50") +
    scale_colour_manual(values = c(`5e-8` = "#1f4e9c", `4.2e-7` = "#b4462b", `1e-5` = "#6b8e23", `1e-4` = "#8a5aa8"), breaks = names(alphas), name = "alpha") +
    labs(x = "Variance explained by one SNP (%)", y = "Power", title = sprintf("Single-SNP power (N = %s)", format(if (test_keep > 0) nrow(ids) * test_keep else nrow(ids), big.mark = ","))) +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  ggsave(paste0(pic_prefix, "_power_snp.png"), g2, width = 6.5, height = 4.5, dpi = 300)
}

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("06-6 SNP 유전력 · 유전 상관 · 검정력 (GCTA GREML, GRM = 04-1 LD pruning SNP 104,441개)\n")
if (test_keep > 0) cat(sprintf("*** TEST_KEEP=%d: %d마리만 쓴 시험 결과. 해석하지 말 것 ***\n", test_keep, nrow(ids)))
cat(sprintf("가까운 친척 GRM 기준 %.2f, 반복 측정 평균의 나이 hinge 기준점 %.1f세\n\n", bk_cut, knot))
if (nrow(h2)) { cat("■ h² (V(G)/Vp ± SE, LRT p는 V(G) = 0에 대한 한쪽 검정)\n"); print(f3(h2[, .(model, desc, n, h2, se, lo95, hi95, h2_rel, se_rel, p_LRT, p_type)]))
  cat("  h2_2grm: h2 = 전체 GRM 몫 (SNP h²), h2_rel = 가까운 친척에서 추가로 보이는 몫. 개에서는 같은 품종·가정의 공유 환경,\n",
      "  희귀 변이, 비상가 효과가 섞여 있어 h2_rel를 유전 몫으로 해석하지 않는다. p_LRT는 V(G1) = 0 검정 (h2_2grm만 Wald 한쪽 검정: LRT 축소 모형이 수렴하지 않음).\n", sep = "") }
if (!is.null(rg_tab) && nrow(rg_tab)) {
  cat("\n■ 유전 상관 (이변량 GREML, 잔차 공분산 0 고정; estimable = 단변량 층별 h²가 두 층 모두 > 2 SE)\n")
  print(f3(rg_tab[, .(pair, tr1, tr2, n, h2_uni_tr1, se_uni_tr1, h2_uni_tr2, se_uni_tr2, h2_tr1, h2_tr2, estimable, rg, se_rg, p_rg0, p_rg1, rg_raw, se_rg_raw)]))
  cat("  h2_tr1·h2_tr2 = 이변량 안의 h² (층 공변량을 먼저 회귀한 잔차 기준), h2_uni = 단변량 층별 GREML. estimable = FALSE면 rg·p를 NA로 두고 원래 값은 _raw에만 남김.\n")
}
if (!is.null(pw_h2)) {
  cat("\n■ 검정력 ① h²\n  ", attr(pw_h2, "info"), "\n", sep = ""); print(f3(pw_h2))
  cat("  이론 SE는 Visscher et al. 2014의 근사 (비혈연·공변량 없음 가정, 분모 = GRM 비대각 분산). 친척·품종 구조와 PC 공변량이 있는 이 자료에서는\n",
      "  근사일 뿐이므로 판단은 관측 SE (power_obsSE)로 하고, 이론 SE는 규모 비교용으로만 본다.\n", sep = "")
  cat("\n■ 검정력 ② 단일 SNP (NCP = N q² / (1 − q²))\n"); print(dcast(f3(pw_snp), q2 ~ alpha, value.var = "power"))
  cat("  80% 검정력에 필요한 분산 설명 비율 q²: ", paste(sprintf("%s %.2f%%", names(attr(pw_snp, "q2_80")), 100 * attr(pw_snp, "q2_80")), collapse = " · "), "\n", sep = "")
}
if (!is.null(pw_sets)) {
  cat("\n■ 검정력 ③ 06-5 세트: 80% 검정력으로 검출할 수 있었던 최소 세트 β (유전자 z 평균 차이; 유전자 z의 null SD ≈ 1)\n")
  print(f3(pw_sets[, .(stage, run, set, NGENES, BETA, SE, m_family, mde_rank1, mde_nominal)]))
  cat(sprintf("  mde_rank1 = 그 단계 BH의 첫 순위 역치 (1단계 0.05/14, 2단계 0.05/28)에서 80%% 검정력. 1단계 main: %.2f ~ %.2f, 2단계: %.2f ~ %.2f (세트가 작을수록 큼)\n",
              min(pw_sets[stage == "1 (main)", mde_rank1]), max(pw_sets[stage == "1 (main)", mde_rank1]),
              min(pw_sets[stage == "2 (age)", mde_rank1]), max(pw_sets[stage == "2 (age)", mde_rank1])))
}
sink()
cat(sprintf("저장: %s_h2.csv · _rg.csv · _power_h2.csv · _power_snp.csv · _power_sets.csv · _summary.txt / 그림 %s_h2.png · _power_snp.png\n", out_prefix, pic_prefix))

#    model      n      h²       SE      95% CI          p
#    h2_age8    2,249  0.0122   0.0367  0 ~ 0.084       0.38
#    h2_young8  4,866  0.0082   0.0156  0 ~ 0.039       0.30
#    h2_sex_m   3,657  0.0000   0.0199  0 ~ 0.039       0.50
#    h2_sex_f   3,458  0.0293   0.0266  0 ~ 0.081       0.12
