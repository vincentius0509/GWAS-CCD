# 06-2_grm_het.R - (A) GRM over-correction diagnostics; (B) heterogeneity of SNP effects between strata. See README.md.
#                                               ancestry (purelike vs mixlike), weight (wt_t1 · wt_t2 · wt_t3)

library(data.table)
library(ggplot2)
library(parallel)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")   
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
out_prefix <- file.path(qc_dir, "06-2")
pic_prefix <- file.path(pic_dir, "06-2")
work_dir <- paste0(out_prefix, "_work")
het_pic_dir <- paste0(pic_prefix, "_het")
runs_tsv <- file.path(qc_dir, "05-1_runs.tsv")
gcta_dir <- file.path(qc_dir, "05-2_gwas")
bfile <- file.path(qc_dir, "03-3_maf0.01_hwe")
main_in <- list(pheno = file.path(qc_dir, "04-2_pheno.txt"), qcovar = file.path(qc_dir, "03-5_qcovar_PC5.txt"),
                covar = file.path(qc_dir, "03-5_covar.txt"), keep = NA_character_)
plink2 <- "/data/00_tools/plink2_linux_avx2_20231123/plink2"
for (f in c(runs_tsv, paste0(bfile, c(".bed", ".bim", ".fam")), plink2, unlist(main_in[1:3]))) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(het_pic_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- as.integer(Sys.getenv("N_CORES", "10"))
parts <- strsplit(Sys.getenv("PARTS", "diag het"), " ")[[1]]
diag_every <- as.integer(Sys.getenv("DIAG_EVERY", "10"))
n_perm <- as.integer(Sys.getenv("N_PERM", "100"))
perm_batch <- as.integer(Sys.getenv("PERM_BATCH", "25"))
r2_every <- as.integer(Sys.getenv("R2_EVERY", "1000"))
if (r2_every %% diag_every != 0) stop("R2_EVERY는 DIAG_EVERY의 배수여야 함 (R² SNP가 진단용 SNP 안에 들어가게)")
plink_threads <- as.integer(Sys.getenv("PLINK_THREADS", "4"))
grm_factor_cut <- 0.97
lead_p <- 1e-5
lead_p2 <- 1e-4
concord_p <- 1e-4
locus_kb <- 500

thr <- c("5e-8" = 5e-8, "4.2e-7" = 4.2e-7, "1e-6" = 1e-6, "1e-5" = 1e-5, "1e-4" = 1e-4)
thr_note <- c("5e-8" = "human / DAP", "4.2e-7" = "Pettersson 2026", "1e-6" = "DAP suggestive", "1e-5" = "suggestive", "1e-4" = "reference")
thr_key <- c("5e-8" = "5e8", "4.2e-7" = "4.2e7", "1e-6" = "1e6", "1e-5" = "1e5", "1e-4" = "1e4")
thr_lty <- c("5e-8" = "solid", "4.2e-7" = "dashed", "1e-6" = "longdash", "1e-5" = "dotted", "1e-4" = "dotdash")
meff <- 0.05 / 4.2e-7
null_loci <- c("5e-8" = 0, "4.2e-7" = 0.70, "1e-6" = 1.2, "1e-5" = 13.5, "1e-4" = 105)

diag_ids <- c("main", "purelike", "mixlike", "age8", "young8", "sex_m", "sex_f", "wt_t1", "wt_t2", "wt_t3")
cmps <- list(
  sex      = list(ids = c("sex_m", "sex_f"),           lab = c("Male", "Female")),
  age      = list(ids = c("young8", "age8"),           lab = c("< 8 y", ">= 8 y")),
  ancestry = list(ids = c("purelike", "mixlike"),      lab = c("purebred-like", "mixed-like")),
  weight   = list(ids = c("wt_t1", "wt_t2", "wt_t3"),  lab = c("<= 15.9 kg", "15.9-29.0 kg", "> 29.0 kg"), x = 1:3))
if (nzchar(Sys.getenv("CMP_IDS"))) cmps <- cmps[intersect(names(cmps), strsplit(Sys.getenv("CMP_IDS"), " ")[[1]])]

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), axis.title = element_text(size = 13), axis.text = element_text(size = 11))
col_chr <- c("#2A78D6", "#898781")
col_ser <- c("#2A78D6", "#EB6834", "#1BAF7A")
chi2 <- function(b, se) (b / se)^2
lam <- function(x2) median(x2[is.finite(x2)]) / qchisq(0.5, 1)
p_to_lambda <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
count_loci <- function(chr, bp, p, t) {
  h <- which(p < t); if (!length(h)) return(0L)
  o <- h[order(chr[h], bp[h])]
  sum(c(TRUE, diff(chr[o]) != 0 | diff(bp[o]) > locus_kb * 1e3))
}
clump_leads <- function(chr, bp, p, t) {
  h <- which(p < t); if (!length(h)) return(integer(0))
  h <- h[order(p[h])]; taken <- rep(FALSE, length(h)); out <- integer(0)
  for (i in seq_along(h)) if (!taken[i]) {
    out <- c(out, h[i]); taken <- taken | (chr[h] == chr[h[i]] & abs(bp[h] - bp[h[i]]) <= locus_kb * 1e3)
  }
  out
}

runs <- fread(runs_tsv, na.strings = c("", "NA", "\"\""))
need <- unique(c(diag_ids, unlist(lapply(cmps, `[[`, "ids"))))
if (!all(need %in% runs$id)) stop("실행 표에 없는 모형: ", paste(setdiff(need, runs$id), collapse = ", "))
mlma_file <- function(id) if (id == "main") runs[id == "main", ref_result] else file.path(gcta_dir, paste0(id, ".mlma.gz"))
summ_file <- function(id) if (id == "main") sub("\\.mlma\\.gz$", "_summary.txt", runs[id == "main", ref_result]) else file.path(gcta_dir, paste0(id, "_summary.txt"))
model_in <- function(id) {
  if (id == "main") return(main_in)
  i <- which(runs$id == id)[1]; r <- as.list(runs[i])
  list(pheno = r$pheno_file, qcovar = r$qcovar_file, covar = r$covar_file, keep = r$keep_file)
}
read_mlma <- function(id) {
  f <- mlma_file(id); if (!file.exists(f)) stop("결과 없음: ", f)
  fread(cmd = paste("zcat", shQuote(f)), select = c("Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p"))
}
read_vg <- function(id) {
  f <- summ_file(id); if (!file.exists(f)) return(c(vg_mean = NA, vg_min = NA, vg_max = NA, vg_zero = NA))
  l <- readLines(f); i <- grep("V\\(G\\)/Vp", l)[1]; if (is.na(i)) return(c(vg_mean = NA, vg_min = NA, vg_max = NA, vg_zero = NA))
  v <- as.numeric(sub("chr[0-9]+ ", "", regmatches(l[i + 1], gregexpr("chr[0-9]+ [0-9.]+", l[i + 1]))[[1]]))
  c(vg_mean = mean(v), vg_min = min(v), vg_max = max(v), vg_zero = sum(v == 0))
}

plink_inputs <- function(id, n_perm_k, seed) {
  inp <- model_in(id); wd <- file.path(work_dir, id); dir.create(wd, showWarnings = FALSE)
  y <- fread(inp$pheno, header = FALSE)[, 1:3]; setnames(y, c("FID", "IID", "y"))
  if (n_perm_k > 0) { set.seed(seed); for (k in seq_len(n_perm_k)) y[[paste0("perm", k)]] <- sample(y$y) }
  q <- fread(inp$qcovar, header = FALSE); setnames(q, c("FID", "IID", paste0("q", seq_len(ncol(q) - 2))))
  cv <- fread(inp$covar, header = FALSE, colClasses = "character"); setnames(cv, c("FID", "IID", paste0("c", seq_len(ncol(cv) - 2))))
  cv[, `:=`(FID = as.character(FID), IID = as.character(IID))]; q[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  y[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  cov <- merge(q, cv, by = c("FID", "IID"))
  if (!is.na(inp$keep) && nzchar(inp$keep)) { k <- fread(inp$keep, header = FALSE, colClasses = "character")[, 1:2]; setnames(k, c("FID", "IID")) }
  else k <- y[, .(FID, IID)]
  k <- k[paste(FID, IID) %in% paste(y$FID, y$IID) & paste(FID, IID) %in% paste(cov$FID, cov$IID)]
  fwrite(setnames(copy(y), c("FID", "IID"), c("#FID", "IID")), file.path(wd, "pheno.txt"), sep = "\t")
  fwrite(setnames(copy(cov), c("FID", "IID"), c("#FID", "IID")), file.path(wd, "covar.txt"), sep = "\t")
  fwrite(k, file.path(wd, "keep.txt"), sep = "\t", col.names = FALSE)
  list(wd = wd, cov = cov, keep = k, n = nrow(k))
}
run_plink <- function(wd, out, extra, pheno = "pheno.txt", cols = "cols=+a1freq") {
  args <- c("--dog", "--bfile", bfile, "--keep", file.path(wd, "keep.txt"), "--pheno", file.path(wd, pheno),
            "--covar", file.path(wd, "covar.txt"), "--covar-variance-standardize", "--glm", "hide-covar", cols,
            "--threads", plink_threads, extra, "--out", file.path(wd, out))
  st <- system2(plink2, args, stdout = file.path(wd, paste0(out, ".stdout")), stderr = file.path(wd, paste0(out, ".stdout")))
  if (st != 0) stop("plink2 실패 (", wd, "/", out, ".log)")
}
read_glm <- function(f) fread(f, select = c("ID", "A1", "TEST", "BETA", "SE", "P"))[TEST == "ADD" & !is.na(SE) & SE > 0]
align_ols <- function(g, o) {
  m <- merge(g, o[, .(SNP = ID, A1o = A1, B = BETA, SE)], by = "SNP")
  m[A1o != A1, B := -B]; m
}

bim_all <- fread(paste0(bfile, ".bim"), header = FALSE, select = c(1, 2), col.names = c("Chr", "SNP"))
diag_snps <- bim_all[seq(1, .N, by = diag_every), SNP]
diag_snp_file <- file.path(work_dir, "diag_snps.txt"); fwrite(data.table(diag_snps), diag_snp_file, col.names = FALSE)
r2_snp_file <- file.path(work_dir, "r2_snps.txt"); fwrite(bim_all[seq(1, .N, by = r2_every), .(SNP)], r2_snp_file, col.names = FALSE)
diag_lab <- sprintf("genome-wide every %dth SNP (%s SNPs)", diag_every, format(length(diag_snps), big.mark = ","))
rm(bim_all)

diag_one <- function(id) {
  tryCatch({
    pin <- plink_inputs(id, 0, 1)
    run_plink(pin$wd, "diag", c("--extract", diag_snp_file))
    g <- read_mlma(id)
    lam_all <- lam(chi2(g$b, g$se))
    gs <- g[SNP %in% diag_snps]; rm(g)
    o <- read_glm(file.path(pin$wd, "diag.y.glm.linear")); m <- align_ols(gs, o)
    yb <- fread(file.path(pin$wd, "pheno.txt")); set.seed(20261111 + match(id, diag_ids)); lam_perm <- numeric(0)
    for (b0 in seq(1, n_perm, by = perm_batch)) {
      ks <- b0:min(n_perm, b0 + perm_batch - 1); pp <- yb[, 1:2]
      for (k in ks) pp[[paste0("perm", k)]] <- sample(yb$y)
      fwrite(pp, file.path(pin$wd, "perm.txt"), sep = "\t")
      run_plink(pin$wd, "perm", c("--extract", diag_snp_file), pheno = "perm.txt", cols = "cols=p")
      for (k in ks) { f <- file.path(pin$wd, sprintf("perm.perm%d.glm.linear", k)); lam_perm <- c(lam_perm, p_to_lambda(fread(f, select = "P")$P)); file.remove(f) }
    }
    m[, ratio := chi2(b, se) / chi2(B, SE)]
    big <- m[abs(B / SE) > 2]
    st <- system2(plink2, c("--dog", "--bfile", bfile, "--keep", file.path(pin$wd, "keep.txt"), "--extract", r2_snp_file,
                            "--export", "A", "--threads", plink_threads, "--out", file.path(pin$wd, "r2")),
                  stdout = file.path(pin$wd, "r2.stdout"), stderr = file.path(pin$wd, "r2.stdout"))
    if (st != 0) stop("plink2 --export 실패")
    G <- fread(file.path(pin$wd, "r2.raw"))
    cov <- pin$cov[match(as.character(G$IID), IID)]
    X <- model.matrix(~ ., data = as.data.frame(cov[, !c("FID", "IID")]))
    if (nrow(X) != nrow(G)) stop("R²: 공변량 결측으로 개체 수가 맞지 않음")
    gm <- as.matrix(G[, -(1:6)]); colnames(gm) <- sub("_[^_]+$", "", colnames(gm))
    for (j in seq_len(ncol(gm))) { x <- gm[, j]; x[is.na(x)] <- mean(x, na.rm = TRUE); gm[, j] <- x }
    res <- qr.resid(qr(X), gm); tot <- colSums(sweep(gm, 2, colMeans(gm))^2)
    r2 <- data.table(SNP = colnames(gm), r2 = 1 - colSums(res^2) / tot)[is.finite(r2)]
    mr <- merge(m, r2, by = "SNP")
    data.table(id = id, n = pin$n, lambda_gcta_all = lam_all, n_snp_diag = nrow(m), n_perm = length(lam_perm),
               lambda_gcta = lam(chi2(m$b, m$se)), lambda_ols = lam(chi2(m$B, m$SE)),
               perm_mean = mean(lam_perm), perm_sd = sd(lam_perm), perm_min = min(lam_perm), perm_max = max(lam_perm),
               chi2_ratio = median(m$ratio, na.rm = TRUE), b_ratio = median(big$b / big$B), se_ratio = median(big$se / big$SE),
               sign_agree = mean(sign(big$b) == sign(big$B)), n_r2 = nrow(mr), one_minus_r2 = median(1 - mr$r2),
               grm_factor = median(mr$ratio / (1 - mr$r2), na.rm = TRUE), t(read_vg(id)), status = "있음")
  }, error = function(e) data.table(id = id, status = paste("오류:", conditionMessage(e))))
}

if ("diag" %in% parts) {
  cat(sprintf("(A) 진단: 모형 %d개, 진단용 SNP %s개 (%d개마다 1개), 순열 %d번 (%d개씩), R² SNP %d개마다 1개, 동시 %d개\n",
              length(diag_ids), format(length(diag_snps), big.mark = ","), diag_every, n_perm, perm_batch, r2_every, n_cores))
  dg <- rbindlist(mclapply(diag_ids, diag_one, mc.cores = n_cores, mc.preschedule = FALSE), fill = TRUE)
  dg[status == "있음", `:=`(perm_lo = perm_mean - 2 * perm_sd, perm_hi = perm_mean + 2 * perm_sd, lambda_adj = lambda_gcta / one_minus_r2)]
  dg[status == "있음", ols_above_null := lambda_ols > perm_hi]
  dg[status == "있음", verdict := fifelse(lambda_adj < perm_lo, "GRM 과보정 의심 (λ_adj < null 하한)",
                                    fifelse(grm_factor < grm_factor_cut, "GRM 보정 있음, 보정 뒤 null 범위 안 (정상)", "과보정 아님 (수축 = 1 − R² 근사)"))]
  num <- names(dg)[sapply(dg, is.numeric) & !names(dg) %in% c("n", "n_snp_diag", "n_perm", "n_r2", "vg_zero")]
  dg[, (num) := lapply(.SD, round, 4), .SDcols = num]
  fwrite(dg, paste0(out_prefix, "_grm_diag.csv"), bom = TRUE)
  sink(paste0(out_prefix, "_grm_diag.txt"))
  cat(sprintf("06-2 (A) GCTA 수축·GRM 과보정 진단 (진단용 SNP: 게놈 전체 %d개마다 1개 = %s개, 순열 %d번; 같은 표본·공변량의 plink2 일반 회귀와 비교)\n",
              diag_every, format(length(diag_snps), big.mark = ","), n_perm))
  cat("lambda_gcta / lambda_ols: 진단용 SNP의 λ | perm: OLS 표현형 순열 λ (평균 ± SD) | chi2_ratio: χ²(GCTA)/χ²(OLS) 중앙값\n")
  cat("one_minus_r2: 1 − R²(SNP ~ 공변량) 중앙값 | lambda_adj: λ_GCTA ÷ one_minus_r2 | grm_factor: (χ² 비 ÷ (1 − R²)) 중앙값, 1이면 수축은 모두 근사 몫\n\n")
  print(dg[, .(id, n, vg_mean, vg_zero, lambda_gcta_all, lambda_gcta, lambda_ols, perm = sprintf("%.3f ± %.3f", perm_mean, perm_sd),
               chi2_ratio, one_minus_r2, lambda_adj, grm_factor, b_ratio, se_ratio, ols_above_null, verdict, status)], row.names = FALSE)
  sink()
  cat(readLines(paste0(out_prefix, "_grm_diag.txt")), sep = "\n")

  pd <- dg[status == "있음"]
  pd[, id := factor(id, levels = rev(diag_ids))]
  pts <- rbind(pd[, .(id, what = "GCTA MLMA-LOCO", v = lambda_gcta)], pd[, .(id, what = "plink2 OLS (same covariates)", v = lambda_ols)],
               pd[, .(id, what = "GCTA / median(1 - R2)", v = lambda_gcta / one_minus_r2)])
  pts[, what := factor(what, levels = c("GCTA MLMA-LOCO", "plink2 OLS (same covariates)", "GCTA / median(1 - R2)"))]
  gd <- ggplot() +
    geom_rect(aes(xmin = 0.95, xmax = 1.05, ymin = -Inf, ymax = Inf), fill = "#F4F3EF") +
    geom_vline(xintercept = 1, color = "#898781", linewidth = 0.4) +
    geom_errorbarh(data = pd, aes(y = id, xmin = perm_mean - 2 * perm_sd, xmax = perm_mean + 2 * perm_sd), height = 0.35, color = "#52514E", linewidth = 0.5) +
    geom_point(data = pts, aes(v, id, color = what, shape = what), size = 2.8) +
    scale_color_manual(values = setNames(col_ser, levels(pts$what)), name = NULL) +
    scale_shape_manual(values = setNames(c(16, 17, 15), levels(pts$what)), name = NULL) +
    labs(title = "Lambda: GCTA vs OLS on the same samples and covariates",
         subtitle = sprintf("%s; bar = OLS permuted-phenotype null (mean +/- 2 SD, %d permutations); shaded = gate 0.95-1.05", diag_lab, n_perm),
         x = "lambda (median chi-square / 0.455)", y = NULL) +
    theme_qc + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                     plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
  ggsave(paste0(pic_prefix, "_grm_diag.png"), gd, width = 9, height = 0.42 * nrow(pd) + 2.4, dpi = 300)
}

het_one <- function(cn) {
  cm <- cmps[[cn]]; K <- length(cm$ids)
  ds <- lapply(cm$ids, read_mlma)
  d <- ds[[1]][, .(Chr, SNP, bp, A1, A2)]
  for (k in seq_len(K)) {
    x <- ds[[k]][, .(SNP, A1k = A1, A2k = A2, Freq, b, se, p)]
    x <- x[is.finite(b) & is.finite(se) & se > 0]
    d <- merge(d, x, by = "SNP", sort = FALSE)
    d <- d[(A1k == A1 & A2k == A2) | (A1k == A2 & A2k == A1)]
    d[A1k != A1, `:=`(b = -b, Freq = 1 - Freq)]
    setnames(d, c("Freq", "b", "se", "p"), paste0(c("freq", "b", "se", "p"), k)); d[, c("A1k", "A2k") := NULL]
  }
  ds_n <- ds[[1]][, .(SNP)]; rm(ds); setorder(d, Chr, bp)
  B <- as.matrix(d[, paste0("b", 1:K), with = FALSE]); S <- as.matrix(d[, paste0("se", 1:K), with = FALSE]); Z <- B / S
  R <- diag(K)
  for (i in 1:(K - 1)) for (j in (i + 1):K) { sel <- abs(Z[, i]) < 1.96 & abs(Z[, j]) < 1.96; R[i, j] <- R[j, i] <- cor(Z[sel, i], Z[sel, j]) }
  if (K == 2) {
    zd <- (B[, 1] - B[, 2]) / sqrt(S[, 1]^2 + S[, 2]^2 - 2 * R[1, 2] * S[, 1] * S[, 2])
    d[, `:=`(z_diff = zd, p_het = 2 * pnorm(-abs(zd)))]
  } else {
    M <- solve(R); V <- 1 / S; VM <- V %*% M
    s00 <- rowSums(VM * V); t0 <- rowSums(VM * Z); bbar <- t0 / s00
    E <- Z - bbar * V; Q <- rowSums((E %*% M) * E)
    VX <- sweep(V, 2, cm$x, "*"); VXM <- VX %*% M
    s01 <- rowSums(VM * VX); s11 <- rowSums(VXM * VX); t1 <- rowSums(VXM * Z); det <- s00 * s11 - s01^2
    slope <- (s00 * t1 - s01 * t0) / det; se_slope <- sqrt(s00 / det)
    d[, `:=`(b_mean = bbar, Q = Q, I2 = pmax(0, (Q - (K - 1)) / Q), p_het = pchisq(Q, K - 1, lower.tail = FALSE),
             trend = slope, trend_se = se_slope, p_trend = 2 * pnorm(-abs(slope / se_slope)))]
  }
  d <- d[is.finite(p_het)]
  n_drop <- nrow(ds_n) - nrow(d)
  num_cols <- names(d)[sapply(d, is.double)]
  fwrite(d[, (num_cols) := lapply(.SD, signif, 5), .SDcols = num_cols], paste0(out_prefix, "_het_", cn, ".tsv.gz"), sep = "\t")

  cnt <- rbindlist(lapply(names(thr), function(nm) data.table(thr = nm, snp = sum(d$p_het < thr[nm]), loci = count_loci(d$Chr, d$bp, d$p_het, thr[nm]),
                                                               expect = meff * thr[nm])))
  prs <- combn(K, 2)
  r_txt <- paste(apply(prs, 2, function(ij) sprintf("r(%s,%s)=%.4f", cm$ids[ij[1]], cm$ids[ij[2]], R[ij[1], ij[2]])), collapse = " ")
  zc_txt <- paste(apply(prs, 2, function(ij) sprintf("cor_z(%s,%s)=%.4f", cm$ids[ij[1]], cm$ids[ij[2]], cor(Z[, ij[1]], Z[, ij[2]]))), collapse = " ")
  ols_het <- function(leads, tag) {
    if (!nrow(leads)) return(leads)
    for (k in seq_len(K)) {
      pin <- plink_inputs(cm$ids[k], 0, 1); fwrite(leads[, .(SNP)], file.path(pin$wd, paste0(tag, cn, ".txt")), col.names = FALSE)
      run_plink(pin$wd, paste0(tag, cn), c("--extract", file.path(pin$wd, paste0(tag, cn, ".txt"))))
      o <- read_glm(file.path(pin$wd, paste0(tag, cn, ".y.glm.linear")))
      m <- align_ols(leads[, .(SNP, A1)], o)
      leads[m, on = "SNP", `:=`(bo = i.B, seo = i.SE)]; setnames(leads, c("bo", "seo"), paste0(c("b_ols", "se_ols"), k))
    }
    Bo <- as.matrix(leads[, paste0("b_ols", 1:K), with = FALSE]); So <- as.matrix(leads[, paste0("se_ols", 1:K), with = FALSE])
    if (K == 2) leads[, p_het_ols := 2 * pnorm(-abs((Bo[, 1] - Bo[, 2]) / sqrt(So[, 1]^2 + So[, 2]^2 - 2 * R[1, 2] * So[, 1] * So[, 2])))]
    else {
      M <- solve(R); Zo <- Bo / So; Vo <- 1 / So; VMo <- Vo %*% M; bb <- rowSums(VMo * Zo) / rowSums(VMo * Vo)
      Eo <- Zo - bb * Vo; leads[, p_het_ols := pchisq(rowSums((Eo %*% M) * Eo), K - 1, lower.tail = FALSE)]
    }
    leads
  }
  leads <- ols_het(d[clump_leads(d$Chr, d$bp, d$p_het, lead_p)], "lead_")
  leads[, cmp := cn]
  leads2 <- ols_het(d[clump_leads(d$Chr, d$bp, d$p_het, lead_p2)], "lead4_")
  leads2[, cmp := cn]
  sc <- rbindlist(lapply(seq_len(K), function(a) rbindlist(lapply(setdiff(seq_len(K), a), function(bk) {
    la <- clump_leads(d$Chr, d$bp, d[[paste0("p", a)]], concord_p)
    s <- sum(sign(d[[paste0("b", a)]][la]) == sign(d[[paste0("b", bk)]][la]))
    data.table(cmp = cn, from = cm$ids[a], to = cm$ids[bk], n_loci = length(la), same_sign = s, prop = s / max(1, length(la)),
               p_binom = if (length(la)) binom.test(s, length(la), 0.5, alternative = "greater")$p.value else NA_real_)
  }))))

  bim_off <- d[, .(len = max(bp)), by = Chr][order(Chr)][, off := cumsum(as.numeric(shift(len, fill = 0)))]
  x_total <- max(bim_off$off + bim_off$len); axis_df <- bim_off[, .(Chr, mid = off + len / 2)]
  pd <- d[, .(Chr, xpos = bp + bim_off$off[match(Chr, bim_off$Chr)], lp = -log10(p_het))]
  pd[, `:=`(xb = floor(xpos / (x_total / 4000)), yb = floor(lp / 0.01))]
  pd <- pd[lp >= 3 | !duplicated(pd, by = c("Chr", "xb", "yb"))]
  ps <- sort(d$p_het); n <- length(ps)
  qi <- sort(unique(c(1:min(20000, n), round(exp(seq(log(min(20001, n)), log(n), length.out = 4000))))))
  qq <- data.table(exp = -log10((qi - 0.5) / n), obs = -log10(ps[qi]))
  band <- data.table(exp = qq$exp, lo = -log10(qbeta(0.975, qi, n - qi + 1)), hi = -log10(qbeta(0.025, qi, n - qi + 1)))
  y_max <- max(-log10(5e-8), max(pd$lp), max(qq$obs)) + 0.6
  hl <- data.table(name = names(thr), y = -log10(thr), lab = sprintf("%s (%s)", names(thr), thr_note))
  lam_h <- p_to_lambda(d$p_het)
  gm <- ggplot(pd, aes(xpos, lp)) +
    geom_point(aes(color = factor(Chr %% 2)), size = 0.4, show.legend = FALSE) +
    geom_hline(data = hl, aes(yintercept = y, linetype = name), color = "#52514E", linewidth = 0.4, show.legend = FALSE) +
    geom_label(data = hl, aes(x = x_total, y = y, label = lab), inherit.aes = FALSE, hjust = 1, vjust = 0.5, size = 3.1, color = "#52514E",
               fill = "white", label.size = 0, label.padding = unit(0.12, "lines")) +
    scale_linetype_manual(values = thr_lty) +
    scale_color_manual(values = c("0" = col_chr[1], "1" = col_chr[2])) +
    scale_x_continuous(breaks = axis_df$mid, labels = ifelse(axis_df$Chr %% 2 == 1 | axis_df$Chr <= 10, axis_df$Chr, ""), expand = c(0.005, 0)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(title = "Manhattan (heterogeneity p)", subtitle = "lines: p thresholds (label = source); overlapping points merged",
         x = "Chromosome", y = expression(-log[10](p[het]))) +
    theme_qc + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                     plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
  gq <- ggplot() +
    geom_ribbon(data = band, aes(x = exp, ymin = lo, ymax = hi), fill = "#E1E0D9") +
    geom_abline(slope = 1, intercept = 0, color = "#898781", linewidth = 0.4) +
    geom_hline(data = hl, aes(yintercept = y, linetype = name), color = "#52514E", linewidth = 0.3, show.legend = FALSE) +
    scale_linetype_manual(values = thr_lty) +
    geom_point(data = qq, aes(exp, obs), color = "#2A78D6", size = 0.6) +
    annotate("text", x = max(qq$exp), y = 0.3, hjust = 1, vjust = 0, size = 3.6, label = sprintf("lambda = %.3f", lam_h)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(title = "QQ", subtitle = "grey band = 95% null interval", x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p))) +
    theme_qc + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
  ttl <- sprintf("%s: %s  |  %s  |  lambda %.3f  |  SNPs (loci) p_het < %s", cn, paste(sprintf("%s (%s)", cm$ids, cm$lab), collapse = " vs "),
                 if (K == 2) "z_diff" else "Cochran Q (df 2)", lam_h, paste(sprintf("%s: %d (%d)", cnt$thr, cnt$snp, cnt$loci), collapse = " / "))
  gg <- gridExtra::arrangeGrob(gm, gq, ncol = 2, widths = c(3.3, 1), top = grid::textGrob(ttl, gp = grid::gpar(fontsize = 13)))
  ggsave(file.path(het_pic_dir, paste0(cn, "_manhattan.png")), gg, width = 18, height = 5.5, dpi = 300)

  if (K == 2) {
    ui <- unique(c(clump_leads(d$Chr, d$bp, d$p1, lead_p), clump_leads(d$Chr, d$bp, d$p2, lead_p)))
    sd2 <- d[ui][, from := fifelse(p1 < lead_p & p2 < lead_p, "both", fifelse(p1 < lead_p, cm$lab[1], cm$lab[2]))]
    sd2[, from := factor(from, levels = c(cm$lab, "both"))]
    lim <- max(abs(c(sd2$b1 + 1.96 * sd2$se1, sd2$b1 - 1.96 * sd2$se1, sd2$b2 + 1.96 * sd2$se2, sd2$b2 - 1.96 * sd2$se2)), 0.1)
    ge <- ggplot(sd2, aes(b1, b2)) +
      geom_hline(yintercept = 0, color = "#C3C2B7", linewidth = 0.3) + geom_vline(xintercept = 0, color = "#C3C2B7", linewidth = 0.3) +
      geom_abline(slope = 1, intercept = 0, color = "#898781", linetype = "dashed", linewidth = 0.4) +
      geom_errorbar(aes(ymin = b2 - 1.96 * se2, ymax = b2 + 1.96 * se2, color = from), width = 0, alpha = 0.5) +
      geom_errorbarh(aes(xmin = b1 - 1.96 * se1, xmax = b1 + 1.96 * se1, color = from), height = 0, alpha = 0.5) +
      geom_point(aes(color = from, shape = p_het < 1e-3), size = 2.6, stroke = 0.8) +
      scale_color_manual(values = setNames(col_ser, c(cm$lab, "both")), name = "p < 1e-5 in", drop = FALSE) +
      scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16), labels = c("FALSE" = "p_het >= 1e-3", "TRUE" = "p_het < 1e-3"), name = NULL) +
      coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
      labs(title = sprintf("Lead SNP effects: %s vs %s", cm$lab[1], cm$lab[2]),
           subtitle = "loci with p < 1e-5 in either stratum; bars = 95% CI; dashed = equal effect (b in within-stratum SD units)",
           x = sprintf("b (%s)", cm$lab[1]), y = sprintf("b (%s)", cm$lab[2])) +
      theme_qc + theme(panel.grid.minor = element_blank(), legend.position = "right", plot.subtitle = element_text(hjust = 0.5, size = 9, color = "#52514E"))
    ggsave(file.path(het_pic_dir, paste0(cn, "_effects.png")), ge, width = 8.5, height = 7, dpi = 300)
  } else {
    top <- head(d[clump_leads(d$Chr, d$bp, d$p_het, 1e-3)], 15)
    if (nrow(top)) {
      top[, locus := factor(sprintf("%s  (p_het %.1e)", SNP, p_het), levels = rev(sprintf("%s  (p_het %.1e)", SNP, p_het)))]
      lt <- rbindlist(lapply(seq_len(K), function(k) top[, .(locus, stratum = cm$lab[k], b = get(paste0("b", k)), se = get(paste0("se", k)))]))
      lt[, stratum := factor(stratum, levels = cm$lab)]
      ge <- ggplot(lt, aes(b, locus, color = stratum)) +
        geom_vline(xintercept = 0, color = "#898781", linewidth = 0.4) +
        geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0, position = position_dodge(width = 0.6)) +
        geom_point(size = 2.2, position = position_dodge(width = 0.6)) +
        scale_color_manual(values = setNames(col_ser, cm$lab), name = NULL) +
        labs(title = "Top heterogeneity loci: effect by weight tertile", subtitle = "lead SNP per 500 kb, ordered by Cochran Q p; bars = 95% CI",
             x = "b (within-stratum SD units)", y = NULL) +
        theme_qc + theme(panel.grid.minor = element_blank(), legend.position = "bottom", axis.text.y = element_text(size = 9),
                         plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
      ggsave(file.path(het_pic_dir, paste0(cn, "_effects.png")), ge, width = 9, height = 0.38 * nrow(top) + 2.4, dpi = 300)
    }
  }
  list(sum = data.table(cmp = cn, ids = paste(cm$ids, collapse = " vs "), test = if (K == 2) "z_diff" else "Cochran Q (df 2)", n_snp = nrow(d),
                        n_snp_dropped = n_drop, r_null = r_txt, cor_z = zc_txt, lambda_het = round(lam_h, 4),
                        dcast(cnt[, .(x = 1, key = paste0("snp_", thr_key[thr]), v = snp)], x ~ key, value.var = "v")[, -1],
                        dcast(cnt[, .(x = 1, key = paste0("loci_", thr_key[thr]), v = loci)], x ~ key, value.var = "v")[, -1],
                        n_lead = nrow(leads), n_lead_ols_1e4 = if (nrow(leads)) sum(leads$p_het_ols < 1e-4, na.rm = TRUE) else 0L),
       cnt =cnt[, cmp := cn], leads = leads, leads2 = leads2, sc = sc)
}

if ("het" %in% parts) {
  cat(sprintf("(B) 이질성: 비교 %s, 동시 %d개\n", paste(names(cmps), collapse = ", "), min(n_cores, length(cmps))))
  hr <- mclapply(names(cmps), function(cn) tryCatch(het_one(cn), error = function(e) list(err = paste(cn, conditionMessage(e)))),
                 mc.cores = min(n_cores, length(cmps)), mc.preschedule = FALSE)
  bad <- sapply(hr, function(x) inherits(x, "try-error") || !is.null(x$err))
  if (any(bad)) cat("오류:", sapply(hr[bad], function(x) if (inherits(x, "try-error")) as.character(x) else x$err), sep = "\n")
  hr <- hr[!bad]
  hs <- rbindlist(lapply(hr, `[[`, "sum"), fill = TRUE)
  cnt <- rbindlist(lapply(hr, `[[`, "cnt"))
  ld <- rbindlist(lapply(hr, `[[`, "leads"), fill = TRUE)
  sc <- rbindlist(lapply(hr, `[[`, "sc"))
  fwrite(hs, paste0(out_prefix, "_het_summary.csv"), bom = TRUE)
  fwrite(ld, paste0(out_prefix, "_het_leads.csv"), bom = TRUE)
  ld2 <- rbindlist(lapply(hr, `[[`, "leads2"), fill = TRUE)
  fwrite(ld2, paste0(out_prefix, "_het_leads_1e4.csv"), bom = TRUE)
  fwrite(sc, paste0(out_prefix, "_sign_concordance.csv"), bom = TRUE)
  sink(paste0(out_prefix, "_het_summary.txt"))
  cat("06-2 (B) 층 사이 이질성 검정 (b는 층 안 표준편차 단위; r = 두 층 모두 |z| < 1.96인 SNP의 z 상관)\n")
  cat(sprintf("참고: 유효 검정 수 %s 가정 기대값 (0.05 / 4.2e-7에서 역산) → %s\n", format(round(meff), big.mark = ","),
              paste(sprintf("%s %.2g", names(thr), meff * thr), collapse = " / ")))
  cat(sprintf("우연 기대 locus 수 (실측: 주 분석 표본 표현형 순열 GWAS 20회 평균, 2026-11-11): %s\n\n",
              paste(sprintf("%s %g", names(null_loci), null_loci), collapse = " / ")))
  for (cn in hs$cmp) {
    h <- hs[cmp == cn]; cat(sprintf("== %s: %s (%s) ==\n", cn, h$ids, h$test))
    cat(sprintf("SNP %s (층 안 단형성 등으로 뺀 SNP %d), λ_het %.4f, %s, %s\n", format(h$n_snp, big.mark = ","), h$n_snp_dropped, h$lambda_het, h$r_null, h$cor_z))
    cat("p_het 기준별 SNP (locus):", paste(cnt[cmp == cn, sprintf("%s %d (%d)", thr, snp, loci)], collapse = " / "), "\n")
    l <- ld[cmp == cn]
    if (nrow(l)) {
      cat(sprintf("대표 SNP (p_het < %g, 500 kb마다 1개) %d개 — p_het_ols = plink2 일반 회귀로 다시 계산:\n", lead_p, nrow(l)))
      cols <- c("Chr", "SNP", "bp", grep(sprintf("^(freq|b|p)[1-%d]$", length(strsplit(h$ids, " vs ")[[1]])), names(l), value = TRUE), "p_het", "p_het_ols")
      print(l[order(p_het), ..cols][, lapply(.SD, function(x) if (is.numeric(x) && !is.integer(x)) signif(x, 3) else x)], row.names = FALSE)
    } else cat(sprintf("대표 SNP (p_het < %g) 없음\n", lead_p))
    cat("효과 방향 일치 (층별 p < 1e-4 locus 대표 SNP):\n"); print(sc[cmp == cn, .(from, to, n_loci, same_sign, prop = round(prop, 3), p_binom = signif(p_binom, 3))], row.names = FALSE)
    cat("\n")
  }
  cat(sprintf("== (2026-11-14 추가) 대표 SNP p_het < %g (500 kb마다 1개; 전체 표 _het_leads_1e4.csv; 우연 기대 locus 약 %g개 = 주 분석 순열 null, 참고) ==\n",
              lead_p2, null_loci["1e-4"]))
  for (cn in hs$cmp) {
    l2 <- ld2[cmp == cn]
    cat(sprintf("  %-9s 대표 SNP %3d개 · 일반 회귀로도 p_het_ols < %g: %3d개\n", cn, nrow(l2), lead_p2,
                if (nrow(l2) && "p_het_ols" %in% names(l2)) sum(l2$p_het_ols < lead_p2, na.rm = TRUE) else 0L))
  }
  sink()
  cat(readLines(paste0(out_prefix, "_het_summary.txt")), sep = "\n")
}
cat(sprintf("저장: %s_grm_diag.*, %s_het_* (1e-4 대표 SNP는 _het_leads_1e4.csv), %s_sign_concordance.csv, 그림 %s_grm_diag.png, %s/\n", out_prefix, out_prefix, out_prefix, pic_prefix, het_pic_dir))

#    sex       4,842,200*  0.932  0.023   1         1       8       78      chr31:23,436,068  1.6e-7 (8.1e-8)
#    age       4,842,736   0.942  0.019   0         1       8       75      chr5:39,898,938   9.3e-7 (5.3e-7)
#    ancestry  4,842,736   0.948  0.013   0         1       12      77      chr17:28,639,756  8.1e-7 (4.7e-7)
#    weight    4,842,736   0.889  ≈ 0     0         0       6       70      chr28:40,923,038  5.3e-6 (8.0e-6)
#      sex       78                        8                   75                              1.6e-7 (chr31)
#      age       75                        8                   75                              9.3e-7
#      ancestry  77                        12                  76                              8.1e-7
#      weight    70                        6                   67                              5.3e-6
