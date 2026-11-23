# 07-2_roh_assoc.R - Tests F_ROH against CSLB in a GREML mixed model. See README.md.
library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_mode <- Sys.getenv("TEST", "0") == "1"
tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("07-2", tag)); pic_prefix <- file.path(pic_dir, paste0("07-2", tag))
work_dir <- paste0(out_prefix, "_work")
froh_f <- if (test_mode) file.path(qc_dir, "07-1_test_froh.csv") else file.path(in_dir, "07-1_froh.csv")
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
grm <- file.path(in_dir, "04-1_grm", "grm_all")
a8 <- file.path(in_dir, "05-1_inputs", paste0("age8_", c("pheno", "qcovar", "covar", "keep"), ".txt")); names(a8) <- c("pheno", "qcovar", "covar", "keep")
main_f <- c(pheno = file.path(in_dir, "04-2_pheno.txt"), qcovar = file.path(in_dir, "03-5_qcovar_PC5.txt"), covar = file.path(in_dir, "03-5_covar.txt"))
for (f in c(froh_f, gcta, paste0(grm, c(".grm.bin", ".grm.id")), a8, main_f)) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
qn <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5)); cn <- c("Sex", "Intact", "Vision", "Cohort_BHS")
wr <- function(d, f) { fwrite(d, f, sep = " ", col.names = FALSE, quote = FALSE); f }

fr <- fread(froh_f, colClasses = c(IID = "character", FID = "character"))
rd <- function(f, nm) { d <- fread(f, header = FALSE, colClasses = list(character = 1:2)); setnames(d, c("FID", "IID", nm)); d }
A <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")), list(rd(a8["pheno"], "y"), rd(a8["qcovar"], qn), fread(a8["covar"], header = FALSE, colClasses = "character", col.names = c("FID", "IID", cn))))
M <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")), list(rd(main_f["pheno"], c("y", "score"))[, score := NULL], rd(main_f["qcovar"], qn), fread(main_f["covar"], header = FALSE, colClasses = "character", col.names = c("FID", "IID", cn))))
keep_cols <- c("IID", "F_ROH", "F_ROH_short", "F_ROH_medium", "F_ROH_long", "Fhat3", "Breed_Status")
A <- merge(A, fr[, ..keep_cols], by = "IID"); M <- merge(M, fr[, ..keep_cols], by = "IID")
cat(sprintf("8세 이상 %d마리 · 전체 %d마리 (F_ROH 있는 개)\n", nrow(A), nrow(M)))

greml_fix <- function(name, d, extra) {
  out <- file.path(work_dir, name)
  extra <- extra[vapply(extra, function(j) isTRUE(sd(d[[j]], na.rm = TRUE) > 0), TRUE)]
  if (!length(extra)) { warning(name, ": 검정할 변수가 모두 상수"); return(NULL) }
  kq <- qn[vapply(qn, function(j) isTRUE(sd(d[[j]]) > 0), TRUE)]; kc <- cn[vapply(cn, function(j) uniqueN(d[[j]]) > 1, TRUE)]
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
  fx <- fread(text = l[(i + 2):length(l)], header = TRUE, fill = TRUE, nrows = 1 + length(qcols) + 200)
  fx <- fx[grepl("^X_", Source) | Source == "mean"]
  pos <- match(extra, qcols) + 1L
  h <- fread(paste0(out, ".hsq"), fill = TRUE)
  data.table(model = name, term = extra, beta = as.numeric(fx[match(paste0("X_", pos), Source), Estimate]),
             se = as.numeric(fx[match(paste0("X_", pos), Source), SE]), n = nrow(d), h2 = as.numeric(h[Source == "V(G)/Vp", Variance]), grm = TRUE)
}
ols_fix <- function(name, d, extra) {
  kq <- qn[vapply(qn, function(j) isTRUE(sd(d[[j]]) > 0), TRUE)]; kc <- cn[vapply(cn, function(j) uniqueN(d[[j]]) > 1, TRUE)]
  m <- lm(as.formula(paste("y ~", paste(c(kq, kc, extra), collapse = " + "))), data = d); s <- summary(m)$coefficients
  data.table(model = name, term = extra, beta = s[extra, 1], se = s[extra, 2], n = nrow(d), h2 = NA_real_, grm = FALSE)
}

M[, age_c := age_cslb - mean(age_cslb)][, F_ROH_x_age := F_ROH * age_c]
res <- rbindlist(list(
  greml_fix("P_age8_FROH", A, "F_ROH"),
  greml_fix("S1_all_FROHxAge", M, c("F_ROH", "F_ROH_x_age")),
  greml_fix("S2_age8_length", A, c("F_ROH_short", "F_ROH_medium", "F_ROH_long")),
  greml_fix("S3_age8_purebred", A[Breed_Status == "Purebred"], "F_ROH"),
  greml_fix("S3_age8_mixed", A[Breed_Status == "Mixed breed"], "F_ROH"),
  ols_fix("S4_age8_noGRM", A, "F_ROH"),
  greml_fix("S5_age8_Fhat3", A, "Fhat3")), fill = TRUE)
res[, `:=`(z = beta / se)][, p := 2 * pnorm(-abs(z))]
sd_of <- function(term, model) { d <- if (grepl("^S1", model)) M else A; if (grepl("purebred", model)) d <- A[Breed_Status == "Purebred"]
  if (grepl("mixed", model)) d <- A[Breed_Status == "Mixed breed"]; sd(d[[term]], na.rm = TRUE) }
res[, `:=`(beta_per_0.1 = beta * 0.1, beta_per_sd = beta * mapply(sd_of, term, model))]
fwrite(res, paste0(out_prefix, "_results.csv"))
P <- res[model == "P_age8_FROH" & term == "F_ROH"]
H3 <- isTRUE(P$beta > 0 && P$p < 0.05)

A[, resid := resid(lm(as.formula(paste("y ~", paste(c(qn, cn), collapse = " + "))), data = A))]
A[, dec := cut(F_ROH, unique(quantile(F_ROH, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
pd <- A[, .(F_ROH = mean(F_ROH), resid = mean(resid), se = sd(resid) / sqrt(.N), n = .N), by = .(dec, Breed_Status)]
g1 <- ggplot(pd, aes(F_ROH, resid, colour = Breed_Status)) + geom_hline(yintercept = 0, colour = "grey60") + geom_pointrange(aes(ymin = resid - 1.96 * se, ymax = resid + 1.96 * se)) +
  labs(x = "F_ROH (decile mean)", y = "CSLB (RINT) residual after covariates", colour = NULL, title = "Dogs aged 8+: genomic inbreeding vs CSLB",
       subtitle = "descriptive; the test is the GREML fixed effect (GRM-adjusted)") + theme_minimal(base_size = 10)
ggsave(paste0(pic_prefix, "_froh_cslb.png"), g1, width = 7, height = 4.5, dpi = 300)

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("07-2 가설 H3: 게놈 근친도 (F_ROH)와 CSLB (2026-11-13 추가, 판정 규칙은 실행 전 헤더에 고정)\n")
if (test_mode) cat("*** TEST=1: chr38 F_ROH · 890마리. 해석하지 말 것 ***\n")
cat("\n■ 모형별 계수 (beta = F 1 단위당; beta_per_0.1 · beta_per_sd = F 0.1 · 1 SD 증가당 CSLB_rint 변화, p = 양측 Wald)\n")
print(f3(res[, .(model, term, n, beta_per_0.1, beta_per_sd, se, z, p, grm, h2)]))
cat(sprintf("\n■ 판정 (주 검정 P_age8_FROH): β %.3g (F 0.1당 %.3g), p %.3g → %s\n", P$beta, P$beta_per_0.1, P$p,
            if (H3) "H3 지지 (Holm 보정 전; 주 가설 H3 · H4 · H5 · H7 사이 보정 필요)" else "H3 지지 안 됨"))
cat("  S4 (GRM 없음)와 주 검정의 차이 = 품종 구조가 근친도 효과로 보이던 몫\n")
sink()
cat(sprintf("저장: %s_results.csv · _summary.txt / 그림 %s_froh_cslb.png\n", out_prefix, pic_prefix))

