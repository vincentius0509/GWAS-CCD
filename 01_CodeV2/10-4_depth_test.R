# 10-4_depth_test.R - mtDNA copy number and Y-chromosome dosage from read depth, and their association with CSLB. See README.md.
library(data.table)
library(ggplot2)
library(parallel)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_mode <- Sys.getenv("TEST", "0") == "1"; pilot <- Sys.getenv("PILOT", "0") == "1"
tag <- if (test_mode) "_test" else if (pilot) "_pilot" else ""
out_prefix <- file.path(qc_dir, paste0("10-4", tag)); pic_prefix <- file.path(pic_dir, paste0("10-4", tag)); work_dir <- paste0(out_prefix, "_work")
p103 <- file.path(qc_dir, "10-3")
cnt_dir <- if (test_mode) paste0(p103, "_test_counts") else paste0(p103, "_counts")
bins_f <- file.path(in_dir, "10-3_ref", "bins10k.tsv"); samp_f <- paste0(p103, "_samples.tsv"); rinfo_f <- paste0(p103, "_runinfo.csv")
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
grm <- file.path(in_dir, "04-1_grm", "grm_all")
a8 <- setNames(file.path(in_dir, "05-1_inputs", paste0("age8_", c("pheno", "qcovar", "covar", "keep"), ".txt")), c("pheno", "qcovar", "covar", "keep"))
main_f <- c(pheno = file.path(in_dir, "04-2_pheno.txt"), qcovar = file.path(in_dir, "03-5_qcovar_PC5.txt"), covar = file.path(in_dir, "03-5_covar.txt"))
for (f in c(bins_f, samp_f, gcta, a8, main_f)) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
par_start <- 6.65e6; n_cal <- 200

bins <- fread(bins_f)[chr %in% c(as.character(1:38), "X", "Y", "MT")]
bins[, kind := fifelse(chr == "MT", "MT", fifelse(chr == "X", "X", fifelse(chr == "Y", "Y", "A")))]
bins[, key := paste(seq, bin)]
samp <- fread(samp_f, header = FALSE, col.names = c("dog", "run", "gsex"), colClasses = "character")
samp <- samp[file.exists(file.path(cnt_dir, paste0(dog, ".stats")))]
cat(sprintf("10-3이 끝난 개 %d마리 (수컷 %d · 암컷 %d)\n", nrow(samp), sum(samp$gsex == "M"), sum(samp$gsex == "F")))
if (nrow(samp) < 4) stop("10-3 결과가 너무 적음: ", cnt_dir)
rd_cnt <- function(dog, sfx = "") {
  f <- file.path(cnt_dir, paste0(dog, sfx, ".counts.gz")); v <- numeric(nrow(bins))
  if (file.size(f) < 30) return(v)
  x <- fread(cmd = paste("gzip -dc", shQuote(f)), header = FALSE, col.names = c("seq", "bin", "n"))
  m <- match(paste(x$seq, x$bin), bins$key); v[m[!is.na(m)]] <- x$n[!is.na(m)]; v
}
rd_stats <- function(dog, sfx = "") { s <- fread(file.path(cnt_dir, paste0(dog, sfx, ".stats")), header = FALSE, col.names = c("k", "v"), colClasses = "character"); setNames(s$v, s$k) }

set.seed(1)
cal <- samp[, .SD[sample(.N, min(.N, n_cal))], by = gsex]
cm <- do.call(cbind, mclapply(cal$dog, rd_cnt, mc.cores = n_cores))
base_a <- bins$kind == "A" & bins$len == 10000 & bins$n < 0.01 & !is.na(bins$gc)
cm <- sweep(cm, 2, colMeans(cm[base_a, , drop = FALSE]), "/")
med_m <- apply(cm[, cal$gsex == "M", drop = FALSE], 1, median); med_f <- apply(cm[, cal$gsex == "F", drop = FALSE], 1, median)
bins[, `:=`(med_m = med_m, med_f = med_f)]
bins[, ok := fcase(kind == "A", base_a & med_m >= 0.5 & med_m <= 1.5 & med_f >= 0.5 & med_f <= 1.5,
                   kind == "X", bin * 1e4 >= par_start & n < 0.01 & !is.na(gc) & med_f >= 0.5 & med_f <= 1.5 & med_m / pmax(med_f, 1e-9) >= 0.3 & med_m / pmax(med_f, 1e-9) <= 0.7,
                   kind == "Y", n < 0.01 & !is.na(gc) & med_m >= 0.25 & med_m <= 0.75 & med_f < 0.05 * med_m,
                   kind == "MT", TRUE, default = FALSE)]
rm(cm); invisible(gc())
fwrite(bins[, .(seq, chr, bin, len, gc = round(gc, 4), n = round(n, 4), med_m = round(med_m, 4), med_f = round(med_f, 4), ok)], paste0(out_prefix, "_bins.csv"))
cat(sprintf("좋은 구간: 상염색체 %d / %d · X %d / %d · Y %d / %d (%.2f Mb) · chrM %d\n",
            bins[kind == "A", sum(ok)], bins[kind == "A", .N], bins[kind == "X", sum(ok)], bins[kind == "X", .N],
            bins[kind == "Y", sum(ok)], bins[kind == "Y", .N], bins[kind == "Y" & ok, sum(len)] / 1e6, bins[kind == "MT", .N]))

ia <- which(bins$kind == "A" & bins$ok); gca <- round(bins$gc[ia], 2)
tgt <- lapply(c(X = "X", Y = "Y", MT = "MT"), function(k) which(bins$kind == k & bins$ok))
gc_fit <- function(v) {
  st <- data.table(gc = gca, n = v[ia])[, .(m = mean(n), nb = .N), by = gc][nb >= 50][order(gc)]
  fit <- loess(m ~ gc, data = st, weights = nb, span = 0.3)
  function(g) { g <- pmin(pmax(g, min(st$gc)), max(st$gc)); pmax(predict(fit, data.frame(gc = g)), 1e-9) }
}
ex_sum <- function(ef, i) sum(ef(bins$gc[i]) * bins$len[i] / 1e4)
metric1 <- function(dog, gsex) {
  v <- rd_cnt(dog); s <- rd_stats(dog); ef <- gc_fit(v)
  ob <- lapply(tgt, function(i) sum(v[i])); ex <- lapply(tgt, function(i) ex_sum(ef, i))
  y_dose <- 2 * ob$Y / ex$Y; n_y <- ob$Y; top <- FALSE
  if (gsex == "M" && file.exists(file.path(cnt_dir, paste0(dog, ".top.stats")))) {
    st2 <- rd_stats(dog, ".top")
    if (as.numeric(st2["total"]) > 0) {
      v2 <- v + rd_cnt(dog, ".top"); ef2 <- gc_fit(v2)
      y_dose <- 2 * sum(v2[tgt$Y]) / ex_sum(ef2, tgt$Y); n_y <- sum(v2[tgt$Y]); top <- TRUE
    }
  }
  data.table(dog = dog, total = as.numeric(s["total"]), mapped = as.numeric(s["mapped"]), q30 = as.numeric(s["q30"]),
             auto_per10k = mean(v[ia]), gc_slope = ef(0.55) / ef(0.40),
             X_dose = ob$X / ex$X, Y_dose = y_dose, mtCN = 2 * ob$MT / ex$MT, n_Y = n_y, n_MT = ob$MT, topup = top)
}
cat("개별 지표 계산 중\n")
met <- rbindlist(mclapply(seq_len(nrow(samp)), function(j) tryCatch(metric1(samp$dog[j], samp$gsex[j]), error = function(e) NULL), mc.cores = n_cores))
met <- merge(samp, met, by = "dog")[, map_rate := mapped / total]
if (file.exists(rinfo_f)) {
  ri <- fread(rinfo_f, select = c("Run", "LoadDate", "Model"), colClasses = "character")
  if (all(c("Run", "LoadDate") %in% names(ri))) {
    ri <- ri[, .(run = Run, load_ym = substr(as.character(LoadDate), 1, 7), instrument = if ("Model" %in% names(ri)) Model else NA_character_)]
    met <- merge(met, unique(ri, by = "run"), by = "run", all.x = TRUE)
  }
}
if (!"load_ym" %in% names(met)) met[, `:=`(load_ym = NA_character_, instrument = NA_character_)]
met[, Batch := fifelse(is.na(load_ym), "unknown", load_ym)][, nb := .N, by = Batch][nb < 50, Batch := "small"][, nb := NULL]
met[, sex_ok := fifelse(gsex == "M", Y_dose > 0.5, Y_dose < 0.1)]
fwrite(met, paste0(out_prefix, "_metrics.csv"))

qn <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5)); cn <- c("Sex", "Intact", "Vision", "Cohort_BHS")
rd <- function(f, nm) { d <- fread(f, header = FALSE, colClasses = list(character = 1:2)); setnames(d, c("FID", "IID", nm)); d }
M <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")), list(rd(main_f["pheno"], c("y", "score"))[, score := NULL], rd(main_f["qcovar"], qn),
            fread(main_f["covar"], header = FALSE, colClasses = "character", col.names = c("FID", "IID", cn))))
A8 <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")), list(rd(a8["pheno"], "y"), rd(a8["qcovar"], qn),
             fread(a8["covar"], header = FALSE, colClasses = "character", col.names = c("FID", "IID", cn))))
tech_q <- c("log_total", "map_rate", "gc_slope"); tech_c <- "Batch"
prep <- function(d) {
  d <- merge(d, met, by.x = "IID", by.y = "dog")
  d[, `:=`(log_total = log(total), log_mtCN = log(mtCN))]
  d[, `:=`(mt_z = as.numeric(scale(log_mtCN)), y_z = NA_real_, x_z = NA_real_)]
  d[gsex == "M" & sex_ok == TRUE, y_z := as.numeric(scale(Y_dose))]; d[gsex == "F" & sex_ok == TRUE, x_z := as.numeric(scale(X_dose))]
  d
}
M <- prep(M); A8 <- prep(A8)

sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("10-4 가설 H7-2: read 깊이 지표 (mtDNA 복제수 · Y 양 · X 양; 설정 · 판정 규칙은 헤더에 실행 전 고정)\n")
if (test_mode) cat("*** TEST=1 시험 (4마리). 해석하지 말 것 ***\n"); if (pilot) cat("*** PILOT=1 시범 (10-3 PILOT 표본). 지표 · 점검 · A (OLS)만 ***\n")
cat(sprintf("\n개 %d마리. 조각 수 중앙값 %.2f M (상염색체 10 kb당 %.1f), 정렬률 중앙값 %.3f\n", nrow(met), median(met$total) / 1e6, median(met$auto_per10k), median(met$map_rate)))
cat(sprintf("좋은 구간: 상염색체 %d · X %d · Y %d (%.2f Mb) · chrM %d\n", length(ia), length(tgt$X), length(tgt$Y), sum(bins$len[tgt$Y]) / 1e6, length(tgt$MT)))
cat("\n■ 점검: 유전 성별별 지표 (중앙값 [2.5%, 97.5%])\n")
q3 <- function(x) sprintf("%.3f [%.3f, %.3f]", median(x, na.rm = TRUE), quantile(x, 0.025, na.rm = TRUE), quantile(x, 0.975, na.rm = TRUE))
print(met[, .(n = .N, Y_dose = q3(Y_dose), X_dose = q3(X_dose), mtCN = q3(mtCN), sex_ok = sum(sex_ok)), by = gsex])
cat(sprintf("수컷 추가 읽기 (Y_dose = 400만 조각): %d / %d마리\n", met[gsex == "M", sum(topup)], met[gsex == "M", .N]))
ym <- M[gsex == "M" & sex_ok == TRUE & age_cslb <= 3, Y_dose]
cat(sprintf("젊은 수컷 (3세 이하) %d마리 Y_dose SD %.3f → 개별 측정 잡음 상한 (Y 소실 %.0f%% 정도부터 개별 구분 가능)\n", length(ym), sd(ym), 100 * 2 * sd(ym)))
sink()

if (test_mode) { cat("시험 모드: 지표까지만\n"); quit(save = "no") }

wr <- function(d, f) { fwrite(d, f, sep = " ", col.names = FALSE, quote = FALSE); f }
greml_fix <- function(name, d, extra, qbase = c(qn, tech_q), cbase = c(cn, tech_c)) {
  out <- file.path(work_dir, name); d <- d[complete.cases(d[, c("y", extra), with = FALSE])]
  kq <- setdiff(qbase[vapply(qbase, function(j) isTRUE(sd(d[[j]]) > 0), TRUE)], extra); kc <- cbase[vapply(cbase, function(j) uniqueN(d[[j]]) > 1, TRUE)]
  qcols <- c(kq, extra)
  if (!file.exists(paste0(out, ".stdout")) || !any(grepl("Estimates of fixed effects", readLines(paste0(out, ".stdout"), warn = FALSE)))) {
    wr(d[, .(FID, IID, y)], paste0(out, ".pheno")); wr(d[, c("FID", "IID", qcols), with = FALSE], paste0(out, ".qcovar"))
    wr(d[, c("FID", "IID", kc), with = FALSE], paste0(out, ".covar")); wr(d[, .(FID, IID)], paste0(out, ".keep"))
    system2(gcta, c("--reml", "--grm", grm, "--keep", paste0(out, ".keep"), "--pheno", paste0(out, ".pheno"), "--qcovar", paste0(out, ".qcovar"),
                    "--covar", paste0(out, ".covar"), "--reml-est-fix", "--reml-no-lrt", "--reml-maxit", "200", "--out", out, "--threads", n_cores),
            stdout = paste0(out, ".stdout"), stderr = paste0(out, ".stdout"))
  }
  l <- readLines(paste0(out, ".stdout"), warn = FALSE); i <- grep("^Estimates of fixed effects", l)
  if (!length(i)) { warning(name, ": GREML 실패 (", out, ".stdout)"); return(NULL) }
  fx <- fread(text = l[(i + 2):length(l)], header = TRUE, fill = TRUE, nrows = 1 + length(qcols) + 400)
  pos <- match(extra, qcols) + 1L
  data.table(model = name, term = extra, beta = as.numeric(fx[match(paste0("X_", pos), Source), Estimate]),
             se = as.numeric(fx[match(paste0("X_", pos), Source), SE]), n = nrow(d), grm = TRUE)
}
ols_fix <- function(name, d, extra, qbase = c(qn, tech_q), cbase = c(cn, tech_c)) {
  d <- d[complete.cases(d[, c("y", extra), with = FALSE])]
  kq <- setdiff(qbase[vapply(qbase, function(j) isTRUE(sd(d[[j]]) > 0), TRUE)], extra); kc <- cbase[vapply(cbase, function(j) uniqueN(d[[j]]) > 1, TRUE)]
  s <- summary(lm(as.formula(paste("y ~", paste(c(kq, kc, extra), collapse = " + "))), data = d))$coefficients
  data.table(model = name, term = extra, beta = s[extra, 1], se = s[extra, 2], n = nrow(d), grm = FALSE)
}
fitA <- if (pilot) ols_fix else greml_fix
q_noage <- c(setdiff(qn, c("age_cslb", "age_hinge")), tech_q)
mkA <- function(d, v) { d <- copy(d); d[, y := get(v)]; d }
res <- rbindlist(list(
  fitA("A1_mtCN_age", mkA(M, "log_mtCN"), "age_cslb", qbase = q_noage),
  fitA("A2_Y_age", mkA(M[gsex == "M" & sex_ok == TRUE], "Y_dose"), "age_cslb", qbase = q_noage),
  fitA("A3_X_age", mkA(M[gsex == "F" & sex_ok == TRUE], "X_dose"), "age_cslb", qbase = q_noage)), fill = TRUE)
if (!pilot) {
  Mx <- copy(M)[, age_c := age_cslb - mean(age_cslb)][, `:=`(mt_x_age = mt_z * age_c, y_x_age = y_z * age_c)]
  res <- rbindlist(list(res,
    greml_fix("B1_age8_mtCN", A8, "mt_z"),
    greml_fix("B2_age8_Y_male", A8[gsex == "M" & sex_ok == TRUE], "y_z"),
    greml_fix("S1_all_mtCNxAge", Mx, c("mt_z", "mt_x_age")),
    greml_fix("S1_all_YxAge_male", Mx[gsex == "M" & sex_ok == TRUE], c("y_z", "y_x_age")),
    greml_fix("S2_age8_mtCN_cov1M", A8[q30 >= 1e6], "mt_z"),
    greml_fix("S2_age8_Y_male_cov1M", A8[gsex == "M" & sex_ok == TRUE & q30 >= 1e6], "y_z"),
    ols_fix("S3_age8_mtCN_noGRM", A8, "mt_z"),
    ols_fix("S3_age8_Y_male_noGRM", A8[gsex == "M" & sex_ok == TRUE], "y_z"),
    greml_fix("S4_age8_X_female", A8[gsex == "F" & sex_ok == TRUE], "x_z")), fill = TRUE)
}
res[, p := 2 * pnorm(-abs(beta / se))]
res[grepl("^B", model), p_holm := p.adjust(p, "holm")]
fwrite(res, paste0(out_prefix, "_tests.csv"))

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), append = TRUE, split = TRUE)
cat("\n■ 검정 (A: 지표 ~ 나이, 1년당; B: 8세 이상 CSLB_rint ~ 지표 1 SD당)\n"); print(f3(res))
for (k in c("A1_mtCN_age", "A2_Y_age")) { r <- res[model == k]; if (nrow(r)) cat(sprintf("■ %s: %.3g / 년 (p %.3g) → %s\n", k, r$beta, r$p, if (r$p < 0.05) "나이에 따라 변함" else "나이 관련 근거 없음")) }
for (k in c("B1_age8_mtCN", "B2_age8_Y_male")) {
  r <- res[model == k]; if (!nrow(r)) next
  cat(sprintf("■ 판정 %s: β %.3g ± %.3g, p %.3g (Holm %.3g) → %s\n", k, r$beta, r$se, r$p, r$p_holm,
              if (r$p_holm < 0.05 && r$beta < 0) "관련 있음 (탐색적; 기대 방향)" else if (r$p_holm < 0.05) "반대 방향으로 유의 (지지 안 됨, 그대로 보고)" else "지지 안 됨"))
}
sink()

pd <- rbind(M[, .(age = age_cslb, value = Y_dose, panel = "Y dose by genetic sex", grp = gsex)],
            M[, .(age = age_cslb, value = log_mtCN, panel = "log mtDNA copy number vs age", grp = gsex)],
            M[gsex == "M" & sex_ok == TRUE, .(age = age_cslb, value = Y_dose, panel = "Y dose vs age (males)", grp = gsex)])
g <- ggplot(pd, aes(age, value, colour = grp)) + geom_point(size = 0.4, alpha = 0.3) + geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  facet_wrap(~panel, scales = "free_y") + labs(x = "age at CSLB (years)", y = NULL, colour = "genetic sex") + theme_bw(base_size = 11)
ggsave(paste0(pic_prefix, "_depth.png"), g, width = 12, height = 4, dpi = 150)
cat(sprintf("저장: %s_bins.csv · _metrics.csv · _tests.csv · _summary.txt / 그림 %s_depth.png\n", out_prefix, pic_prefix))

