# 13-3_finemap.R - Conditional analysis (COJO) and fine-mapping (SuSiE) at lead loci. See README.md.
if (Sys.getenv("OPENBLAS_NUM_THREADS") != "1") {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  if (length(f) == 1) { Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"); quit(save = "no", status = system2(file.path(R.home("bin"), "Rscript"), shQuote(f))) }
}
suppressMessages({ library(data.table); library(ggplot2); library(parallel); library(susieR) })
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
pic_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/05_GWAS_2023"
test_mode <- Sys.getenv("TEST", "0") == "1"; tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("13-3", tag)); fig_dir <- file.path(pic_dir, paste0("13-3", tag, "_finemap"))
in13 <- file.path(qc_dir, "13-1"); in132 <- file.path(qc_dir, "13-2")
n_cores <- as.integer(Sys.getenv("N_CORES", "8"))
plink1 <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
bfile <- file.path(d0, "03-3_maf0.01_hwe")
gl <- fread("/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/06-1_gene_loc.txt", header = FALSE,
            col.names = c("GENE", "chr", "start", "end", "strand", "sym"))
cojo_p <- 1e-4
p_keep <- 1e-4; p_gone <- 0.01
ld_hi <- 0.5; ld_lo <- 0.1
susie_L <- 10; cs_cov <- 0.95; cs_purity <- 0.5; cs_narrow <- 10; s_flag <- 0.1
fam_order <- c("main", "wt_t1", "wt_t2", "wt_t3", "sex_m", "sex_f", "purelike", "mixlike", "young8", "age8", "int_age",
               "agefloor6", "agefloor7", "agefloor9", "age10", "agefloor11", "agefloor12", "bin_all", "bin_age8")
approx_models <- c("int_age", "bin_all", "bin_age8")
gfile <- function(id) {
  if (id == "main") return("04-5_gwas/rint_maf0.01_PC5.mlma.gz")
  if (id %in% c("int_age", "bin_age8", "bin_all")) return(sprintf("05-3_regenie/%s.regenie.gz", id))
  sprintf("05-2_gwas/%s.mlma.gz", id)
}
ld_grp <- function(m) fifelse(m %in% c("main", "int_age", "bin_all"), "pooled", fifelse(m == "bin_age8", "age8", m))
age_series <- c("young8", "main", "agefloor6", "agefloor7", "age8", "agefloor9", "age10", "agefloor11", "agefloor12")
wt_series <- c("wt_t1", "wt_t2", "wt_t3")
series_order <- c(age_series, wt_series, "sex_m", "sex_f", "purelike", "mixlike", "int_age", "bin_all", "bin_age8")
series_of <- function(m) fifelse(m %in% age_series, "나이", fifelse(m %in% wt_series, "체중", "기타"))
keep_of <- function(g) if (g == "pooled") file.path(d0, "04-1_grm", "grm_all.grm.id") else file.path(d0, "05-1_inputs", paste0(g, "_keep.txt"))
cov_of <- function(g) if (g == "pooled") file.path(d0, c("03-5_covar.txt", "03-5_qcovar_PC5.txt")) else file.path(d0, "05-1_inputs", paste0(g, c("_covar.txt", "_qcovar.txt")))
cov_design <- function(g) {
  cf <- fread(cov_of(g)[1], header = FALSE, colClasses = "character"); qf <- fread(cov_of(g)[2], header = FALSE)
  setnames(cf, c("FID", "IID", paste0("c", seq_len(ncol(cf) - 2)))); setnames(qf, c("FID", "IID", paste0("q", seq_len(ncol(qf) - 2))))
  qf[, IID := as.character(IID)]; d <- merge(qf[, -"FID"], cf[, -"FID"], by = "IID")
  keepc <- names(d)[-1][sapply(d[, -1], function(v) uniqueN(v) > 1)]
  D <- model.matrix(~ ., data = as.data.frame(d[, ..keepc])); rownames(D) <- d$IID; D
}
for (f in c(paste0(in13, c("_compare.csv", "_region_leads.csv")), paste0(in132, "_regions.csv"), paste0(bfile, ".bed"), plink1, gcta,
            file.path(d0, sapply(fam_order, gfile)), sapply(unique(ld_grp(fam_order)), keep_of), unlist(lapply(unique(ld_grp(fam_order)), cov_of))))
  if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
wd <- tempfile("r13_3_"); dir.create(wd)

R <- fread(paste0(in132, "_regions.csv"))[, .(region_id, chr, index, index_model, index_bp, lo, hi, panels, multi_signal, in94)]
if (test_mode) R <- R[region_id %in% c("chr4_r8", "chr7_r8", "chr5_r25", "chr2_r33")]
if (nzchar(Sys.getenv("REGIONS"))) { R <- R[region_id %in% strsplit(Sys.getenv("REGIONS"), ",")[[1]]]; tag <- "_test"; test_mode <- TRUE
  out_prefix <- file.path(qc_dir, paste0("13-3", tag)); fig_dir <- file.path(pic_dir, paste0("13-3", tag, "_finemap")); dir.create(fig_dir, showWarnings = FALSE) }
setorder(R, chr, index_bp); regs <- R$region_id
cmp <- fread(paste0(in13, "_compare.csv"))[region_id %in% regs]
leads13 <- unique(cmp[, .(region_id, model, SNP, bp, p)], by = c("region_id", "model", "SNP"))
lp <- fread(paste0(in13, "_region_leads.csv"))[region_id %in% regs]
J <- unique(cmp[, .(region_id, model)])[R[, .(region_id, chr, lo, hi)], on = "region_id", nomatch = NULL]
J[, grp := ld_grp(model)]; J[, approx := model %in% approx_models]
J[, has_1e5 := paste(region_id, model) %in% cmp[thr == 1e-5, paste(region_id, model)]]
R[, panels := sapply(region_id, function(r) paste(intersect(fam_order, J[region_id == r, model]), collapse = ","))]
n_dogs <- sapply(setNames(unique(ld_grp(fam_order)), unique(ld_grp(fam_order))), function(g) length(readLines(keep_of(g))))
cat(sprintf("13-3 위치 %d개 · 작업 (위치 × 모형) %d개 · 모형 %d개 · LD 표본 %d개 · 대표 SNP 짝 %d개\n", nrow(R), nrow(J), uniqueN(J$model),
            uniqueN(J$grp), nrow(lp)))

cat(sprintf("[%s] 위치 bfile %d개\n", format(Sys.time(), "%T"), nrow(R)))
rb <- mclapply(seq_len(nrow(R)), function(j) {
  out <- file.path(wd, paste0("reg_", R$region_id[j]))
  system2(plink1, c("--bfile", bfile, "--dog", "--threads", "1", "--keep-allele-order", "--chr", R$chr[j], "--from-bp", R$lo[j], "--to-bp", R$hi[j], "--make-bed", "--out", out),
          stdout = FALSE, stderr = FALSE)
  file.exists(paste0(out, ".bed"))
}, mc.cores = min(n_cores, nrow(R)))
if (!all(unlist(rb) %in% TRUE)) stop("위치 bfile 실패: ", paste(R$region_id[!unlist(rb) %in% TRUE], collapse = ", "))
BIM <- rbindlist(lapply(regs, function(r) fread(file.path(wd, paste0("reg_", r, ".bim")), header = FALSE,
                                                col.names = c("chr", "SNP", "cm", "bp", "bA1", "bA2"))[, region_id := r]))

win_f <- file.path(wd, "win.txt"); fwrite(R[, .(chr, lo, hi)], win_f, sep = "\t", col.names = FALSE)
models <- intersect(fam_order, unique(c(J$model, "young8", "age8", "sex_m", "sex_f", "purelike", "mixlike", "wt_t1", "wt_t2", "wt_t3")))
model_n <- sapply(setNames(models, models), function(m) n_dogs[[ld_grp(m)]])
read_win <- function(id) {
  f <- file.path(d0, gfile(id)); rg <- grepl("regenie", f); test <- if (id == "int_age") "ADD-INT_SNPxage_cslb" else "ADD"
  aw <- sprintf("awk 'NR==FNR {n[$1]++; s[$1,n[$1]]=$2; e[$1,n[$1]]=$3; next} FNR==1 {print; next} %s {c=$1; p=$%d; for (k=1; k<=n[c]; k++) if (p>=s[c,k] && p<=e[c,k]) {print; break}}' %s -",
                if (rg) sprintf("$8==\"%s\"", test) else "", if (rg) 2 else 3, win_f)
  x <- fread(cmd = sprintf("gzip -dc %s | %s", f, aw))
  if (rg) x[, .(SNP = ID, A1 = ALLELE1, A2 = ALLELE0, freq = A1FREQ, b = BETA, se = SE, p = 10^(-LOG10P), N = as.numeric(N))]
  else x[, .(SNP, A1, A2, freq = Freq, b, se, p, N = as.numeric(model_n[[id]]))]
}
cat(sprintf("[%s] GWAS %d개 읽는 중\n", format(Sys.time(), "%T"), length(models)))
G <- mclapply(setNames(models, models), read_win, mc.cores = min(n_cores, length(models)))
bad <- names(G)[sapply(G, function(x) inherits(x, "try-error") || !is.data.table(x) || !nrow(x))]
if (length(bad)) stop("GWAS 읽기 실패 또는 범위 안 SNP 0개: ", paste(bad, collapse = ", "))
for (m in models) G[[m]] <- G[[m]][is.finite(b) & is.finite(se) & se > 0]

run_gcta <- function(args, out) {
  system2(gcta, c(args, "--thread-num", "1", "--out", out), stdout = FALSE, stderr = FALSE)
  invisible(out)
}
write_ma <- function(r, m) {
  f <- file.path(wd, sprintf("ma_%s_%s.txt", r, m))
  if (!file.exists(f)) fwrite(G[[m]][SNP %in% BIM[region_id == r, SNP], .(SNP, A1, A2, freq, b, se, p, N)], f, sep = " ")
  f
}
for (k in seq_len(nrow(J))) write_ma(J$region_id[k], J$model[k])
JS <- rbind(J[, .(region_id, model, grp, pthr = cojo_p)], J[has_1e5 == TRUE, .(region_id, model, grp, pthr = 1e-5)])
cat(sprintf("[%s] COJO-slct %d개 (1e-4 %d · 1e-5 %d)\n", format(Sys.time(), "%T"), nrow(JS), nrow(J), sum(J$has_1e5)))
SL <- rbindlist(mclapply(seq_len(nrow(JS)), function(k) {
  r <- JS$region_id[k]; m <- JS$model[k]; pt <- JS$pthr[k]; out <- file.path(wd, sprintf("slct_%s_%s_%g", r, m, pt))
  run_gcta(c("--bfile", file.path(wd, paste0("reg_", r)), "--keep", keep_of(JS$grp[k]), "--cojo-file", write_ma(r, m), "--cojo-slct", "--cojo-p", pt), out)
  f <- paste0(out, ".jma.cojo")
  if (!file.exists(f)) {
    ok <- file.exists(paste0(out, ".log")) && any(grepl("No SNP|no SNP|not.*significant|0 SNPs", readLines(paste0(out, ".log"))))
    return(data.table(region_id = r, model = m, pthr = pt, SNP = NA_character_, cojo_note = if (ok) "기준 넘는 SNP 없음" else "GCTA 출력 없음"))
  }
  x <- fread(f); x[, .(region_id = r, model = m, pthr = pt, SNP, bp, refA, freq, b, se, p, bJ, bJ_se, pJ, LD_r, cojo_note = "")]
}, mc.cores = min(n_cores, nrow(JS))), fill = TRUE)
bad <- SL[cojo_note == "GCTA 출력 없음"]
if (nrow(bad)) { print(bad); stop("COJO-slct 실패 (로그는 작업 폴더 ", wd, ")") }

dirs <- unique(rbind(lp[, .(region_id, model = model_b, target = lead_b, cond = lead_a, cond_model = model_a)],
                     lp[, .(region_id, model = model_a, target = lead_a, cond = lead_b, cond_model = model_b)]))
dirs <- dirs[target != cond]
cj <- unique(dirs[, .(region_id, model, cond)]); cj[, grp := ld_grp(model)]
cat(sprintf("[%s] COJO-cond %d개 (짝 방향 %d개)\n", format(Sys.time(), "%T"), nrow(cj), nrow(dirs)))
CD <- rbindlist(mclapply(seq_len(nrow(cj)), function(k) {
  r <- cj$region_id[k]; m <- cj$model[k]; s <- cj$cond[k]; out <- file.path(wd, sprintf("cond_%s_%s_%d", r, m, k))
  writeLines(s, paste0(out, ".snp"))
  run_gcta(c("--bfile", file.path(wd, paste0("reg_", r)), "--keep", keep_of(cj$grp[k]), "--cojo-file", write_ma(r, m), "--cojo-cond", paste0(out, ".snp")), out)
  f <- paste0(out, ".cma.cojo")
  if (!file.exists(f)) return(data.table(region_id = r, model = m, cond = s, SNP = NA_character_, cond_note = "GCTA 출력 없음"))
  x <- fread(f); x[, .(region_id = r, model = m, cond = s, SNP, pC, bC, bC_se, cond_note = "")]
}, mc.cores = min(n_cores, nrow(cj))), fill = TRUE)

SJ <- unique(J[, .(region_id, grp)])
lead_snps <- unique(c(leads13$SNP, R$index))
cat(sprintf("[%s] SuSiE: 위치 × 표본 %d개 (모형 %d개)\n", format(Sys.time(), "%T"), nrow(SJ), nrow(J)))
SU <- mclapply(seq_len(nrow(SJ)), function(k) {
  r <- SJ$region_id[k]; g <- SJ$grp[k]; out <- file.path(wd, sprintf("geno_%s_%s", r, g))
  system2(plink1, c("--bfile", file.path(wd, paste0("reg_", r)), "--dog", "--threads", "1", "--keep-allele-order", "--keep", keep_of(g), "--recode", "A", "--out", out),
          stdout = FALSE, stderr = FALSE)
  raw <- fread(paste0(out, ".raw")); unlink(paste0(out, ".raw")); iid <- as.character(raw$IID)
  X <- as.matrix(raw[, -(1:6)]); rm(raw)
  bim <- BIM[region_id == r]
  cnt <- sub("^.*_", "", colnames(X)); if (!identical(cnt, bim$bA1)) stop("recode A 대립유전자가 bim A1과 다름: ", r, " ", g)
  colnames(X) <- bim$SNP
  mu <- colMeans(X, na.rm = TRUE); na <- which(is.na(X), arr.ind = TRUE); if (nrow(na)) X[na] <- mu[na[, 2]]
  sdv <- apply(X, 2, sd); X <- X[, sdv > 0, drop = FALSE]
  D <- cov_design(g); if (!all(iid %in% rownames(D))) stop("공변량 없는 개체: ", g)
  Q <- qr.Q(qr(D[iid, , drop = FALSE]))
  Xr <- X - Q %*% crossprod(Q, X); sr <- sqrt(colSums(Xr^2) / (nrow(Xr) - 1))
  Xr <- sweep(Xr[, sr > 1e-8, drop = FALSE], 2, sr[sr > 1e-8], "/")
  X <- scale(X)
  res <- list(); pip <- list(); cs <- list(); lr <- NULL; last_key <- ""
  for (m in J[region_id == r & grp == g, model]) {
    gw <- merge(G[[m]], bim[, .(SNP, bA1, bA2)], by = "SNP")[SNP %in% colnames(Xr)]
    gw[, flip := fifelse(A1 == bA1 & A2 == bA2, 1, fifelse(A1 == bA2 & A2 == bA1, -1, NA_real_))]
    gw <- gw[!is.na(flip)]; gw[, z := flip * b / se]
    idx <- match(gw$SNP, colnames(Xr)); key <- paste(idx, collapse = ",")
    if (!identical(key, last_key)) {
      Rm <- crossprod(Xr[, idx, drop = FALSE]) / (nrow(Xr) - 1); Rm <- (Rm + t(Rm)) / 2; attr(Rm, "eigen") <- eigen(Rm, symmetric = TRUE); last_key <- key
    }
    n <- median(gw$N)
    fit <- try(suppressMessages(susie_rss(z = gw$z, R = Rm, n = n, L = susie_L, coverage = cs_cov, min_abs_corr = cs_purity, max_iter = 500)), silent = TRUE)
    s_est <- try(estimate_s_rss(gw$z, Rm, n), silent = TRUE)
    kr <- try(kriging_rss(gw$z, Rm, n, s = if (is.numeric(s_est)) s_est else 0)$conditional_dist, silent = TRUE)
    n_out <- if (is.data.frame(kr)) sum(kr$logLR > 2 & abs(kr$z) > 2, na.rm = TRUE) else NA_integer_
    if (inherits(fit, "try-error")) {
      res[[m]] <- data.table(region_id = r, model = m, grp = g, n_snp = nrow(gw), n_susie = n, susie_ok = FALSE, susie_note = as.character(fit)[1])
      next
    }
    pp <- data.table(region_id = r, model = m, SNP = gw$SNP, bp = bim$bp[match(gw$SNP, bim$SNP)], z = gw$z, pip = fit$pip)
    css <- fit$sets$cs; pur <- fit$sets$purity
    if (length(css)) {
      cs[[m]] <- rbindlist(lapply(seq_along(css), function(i) {
        ii <- css[[i]]; data.table(region_id = r, model = m, cs = i, SNP = gw$SNP[ii], bp = pp$bp[ii], pip = fit$pip[ii], cs_purity = pur$min.abs.corr[i])
      }))
      pp[, cs := 0L]; for (i in seq_along(css)) pp[css[[i]], cs := i]
    } else pp[, cs := 0L]
    pip[[m]] <- pp
    res[[m]] <- data.table(region_id = r, model = m, grp = g, n_snp = nrow(gw), n_susie = n, susie_ok = TRUE, susie_converged = isTRUE(fit$converged),
                           n_cs = length(css), max_pip = max(fit$pip), top_pip_snp = gw$SNP[which.max(fit$pip)], n_pip_ge_0.1 = sum(fit$pip >= 0.1),
                           s_rss = if (is.numeric(s_est)) s_est else NA_real_, n_kriging_outlier = n_out, susie_note = "")
  }
  ls_ <- intersect(lead_snps, colnames(X))
  if (length(ls_) >= 2) { rr <- crossprod(X[, ls_, drop = FALSE]) / (nrow(X) - 1); lr <- data.table(region_id = r, grp = g, a = rep(ls_, length(ls_)), b = rep(ls_, each = length(ls_)), r = as.vector(rr)) }
  list(res = rbindlist(res, fill = TRUE), pip = rbindlist(pip, fill = TRUE), cs = rbindlist(cs, fill = TRUE), lr = lr)
}, mc.cores = min(n_cores, nrow(SJ)))
bad <- which(sapply(SU, function(x) inherits(x, "try-error") || !is.list(x)))
if (length(bad)) { print(SJ[bad]); print(SU[[bad[1]]]); stop("SuSiE 작업 실패") }
RS <- rbindlist(lapply(SU, `[[`, "res"), fill = TRUE); PIP <- rbindlist(lapply(SU, `[[`, "pip"), fill = TRUE)
CSN <- rbindlist(lapply(SU, `[[`, "cs"), fill = TRUE); LR <- rbindlist(lapply(SU, `[[`, "lr"), fill = TRUE); rm(SU)
if (any(!RS$susie_ok)) { cat("SuSiE 실패 작업:\n"); print(RS[susie_ok == FALSE]) }

cs_sum <- if (nrow(CSN)) CSN[order(-pip), .(cs_size = .N, cs_purity = cs_purity[1], top_snp = SNP[1], top_bp = bp[1], top_pip = pip[1],
                                           has_index = any(SNP %in% R[region_id == .BY$region_id, index]),
                                           has_own_lead = any(SNP %in% leads13[region_id == .BY$region_id & model == .BY$model, SNP])),
                                       by = .(region_id, model, cs)] else data.table()
if (nrow(cs_sum)) cs_sum[, narrow := cs_size <= cs_narrow]
RUN <- merge(J[, .(region_id, model, grp, approx, has_1e5)], SL[pthr == cojo_p, .(n_cojo = sum(!is.na(SNP)), cojo_snps = paste(na.omit(SNP), collapse = ";")),
                                                                    by = .(region_id, model)], by = c("region_id", "model"), all.x = TRUE)
RUN <- merge(RUN, SL[pthr == 1e-5, .(n_cojo_1e5 = sum(!is.na(SNP))), by = .(region_id, model)], by = c("region_id", "model"), all.x = TRUE)
RUN <- merge(RUN, RS[, !"grp"], by = c("region_id", "model"), all.x = TRUE)
RUN[, s_caution := !is.na(s_rss) & s_rss >= s_flag]
RUN[, model := factor(model, fam_order)]; setorder(RUN, region_id, model); RUN[, model := as.character(model)]

dirs[, grp := ld_grp(model)]
dirs[CD, on = .(region_id, model, cond, target = SNP), pC := i.pC]
r2_of <- function(reg, g, ta, cb) { x <- LR[region_id == reg & grp == g & a == ta & b == cb, r]; if (length(x)) x[1]^2 else NA_real_ }
dirs[, r2_model := mapply(r2_of, region_id, grp, target, cond)]
dirs[, p_target := { v <- G[[model[1]]][SNP == target[1], p]; if (length(v)) v[1] else NA_real_ }, by = .(model, target)]
dirs[, cond_status := fifelse(is.na(pC) & !is.na(r2_model) & r2_model > 0.9, "사라짐 (공선)",
                       fifelse(is.na(pC), "계산 안 됨", fifelse(pC < p_keep, "남음", fifelse(pC >= p_gone, "사라짐", "약해짐"))))]
PC <- copy(lp[, .(region_id, model_a, lead_a, model_b, lead_b, dist_kb, r2_13_1 = r2, ld_signal)])
PC[dirs, on = .(region_id, model_b = model, lead_b = target, lead_a = cond), `:=`(pC_b_given_a = i.pC, status_b = i.cond_status, p_b = i.p_target)]
PC[dirs, on = .(region_id, model_a = model, lead_a = target, lead_b = cond), `:=`(pC_a_given_b = i.pC, status_a = i.cond_status, p_a = i.p_target)]
PC[lead_a == lead_b, `:=`(status_a = "같은 SNP", status_b = "같은 SNP")]
gone <- function(s) s %in% c("사라짐", "사라짐 (공선)", "같은 SNP")
PC[, pair_call := fifelse(gone(status_a) & gone(status_b), "같은 신호",
                  fifelse(status_a == "남음" & status_b == "남음", "다른 신호",
                  fifelse((status_a == "남음" & gone(status_b)) | (status_b == "남음" & gone(status_a)), "한쪽 신호 포함", "불분명")))]

CP <- data.table(region_id = character(), model_a = character(), cs_a = integer(), top_a = character(), model_b = character(), cs_b = integer(),
                 top_b = character(), n_shared = integer(), r2_top = numeric(), cs_call = character())
if (nrow(cs_sum)) {
  tops <- unique(cs_sum$top_snp); plp <- file.path(wd, "cs_tops")
  writeLines(tops, paste0(plp, ".txt"))
  system2(plink1, c("--bfile", bfile, "--dog", "--threads", min(n_cores, 8), "--keep", keep_of("pooled"), "--extract", paste0(plp, ".txt"), "--r2",
                    "--ld-window-kb", "10000", "--ld-window", "999999", "--ld-window-r2", "0", "--out", plp), stdout = FALSE, stderr = FALSE)
  tl <- if (file.exists(paste0(plp, ".ld"))) fread(paste0(plp, ".ld")) else data.table(SNP_A = character(), SNP_B = character(), R2 = numeric())
  tl <- rbind(tl[, .(a = SNP_A, b = SNP_B, r2 = R2)], tl[, .(a = SNP_B, b = SNP_A, r2 = R2)])
  CP <- rbindlist(c(list(CP), lapply(unique(cs_sum$region_id), function(reg) {
    x <- cs_sum[region_id == reg]; if (uniqueN(x$model) < 2) return(NULL)
    pr <- CJ(i = seq_len(nrow(x)), j = seq_len(nrow(x)))[i < j][x$model[i] != x$model[j]]
    if (!nrow(pr)) return(NULL)
    rbindlist(lapply(seq_len(nrow(pr)), function(q) {
      xa <- x[pr$i[q]]; xb <- x[pr$j[q]]
      sa <- CSN[region_id == reg & model == xa$model & cs == xa$cs, SNP]; sb <- CSN[region_id == reg & model == xb$model & cs == xb$cs, SNP]
      ta <- xa$top_snp; tb <- xb$top_snp
      r2t <- if (ta == tb) 1 else tl[a == ta & b == tb, r2][1]
      data.table(region_id = reg, model_a = xa$model, cs_a = xa$cs, top_a = ta, model_b = xb$model, cs_b = xb$cs, top_b = tb,
                 n_shared = length(intersect(sa, sb)), r2_top = r2t)
    }))
  })), fill = TRUE)
  if (nrow(CP)) CP[, cs_call := fifelse(n_shared > 0, "CS 겹침", fifelse(is.na(r2_top), "LD 없음", fifelse(r2_top >= ld_hi, "같은 신호 (r² 높음)",
                                    fifelse(r2_top < ld_lo, "다른 신호", "불분명"))))]
}

RG <- RUN[, .(n_models = .N, max_cojo = max(n_cojo, na.rm = TRUE), max_cs = suppressWarnings(max(n_cs, na.rm = TRUE)),
              models_cojo2 = paste(model[n_cojo >= 2], collapse = ","), models_cs2 = paste(model[!is.na(n_cs) & n_cs >= 2], collapse = ","),
              n_s_caution = sum(s_caution)), by = region_id]
RG <- merge(RG, PC[, .(n_pairs = .N, n_same = sum(pair_call == "같은 신호"), n_diff = sum(pair_call == "다른 신호"),
                       n_nested = sum(pair_call == "한쪽 신호 포함"), n_unclear = sum(pair_call == "불분명")), by = region_id], by = "region_id", all.x = TRUE)
for (v in c("n_pairs", "n_same", "n_diff", "n_nested", "n_unclear")) set(RG, which(is.na(RG[[v]])), v, 0L)
RG <- merge(R[, .(region_id, chr, index, index_bp, multi_signal_13_1 = multi_signal)], RG, by = "region_id")
RG[, region_call := fifelse(max_cojo >= 2 | max_cs >= 2 | n_diff > 0, "여러 신호", fifelse(max_cojo == 0 & max_cs <= 0, "신호 약함", "한 신호"))]
RG[R, on = "region_id", in94 := i.in94]
setorder(RG, chr, index_bp)

LC <- cmp[, .(model, SNP, thr, region_id, chr, bp, p, tier)]
LC[, in94 := region_id %in% R[in94 == TRUE, region_id]]
near_lead <- function(reg, m, bps) {
  L <- unique(LC[region_id == reg & model == m, .(SNP, bp)]); if (nrow(L) == 1) return(rep(L$SNP, length(bps)))
  L$SNP[sapply(bps, function(b) which.min(abs(L$bp - b)))]
}
SLa <- SL[!is.na(SNP)][, lead := near_lead(region_id[1], model[1], bp), by = .(region_id, model)]
LC[, n_cojo := SLa[pthr == thr[1]][.SD, on = .(region_id, model, lead = SNP), .N, by = .EACHI]$N, by = thr]
if (nrow(cs_sum)) cs_sum[, lead := near_lead(region_id[1], model[1], top_bp), by = .(region_id, model)]
csl <- if (nrow(cs_sum)) cs_sum[, .(n_cs = .N, cs_size = min(cs_size), cs_top_pip = max(top_pip)), by = .(region_id, model, lead)] else
  data.table(region_id = character(), model = character(), lead = character(), n_cs = integer(), cs_size = integer(), cs_top_pip = numeric())
LC[csl, on = .(region_id, model, SNP = lead), `:=`(n_cs = i.n_cs, cs_size = i.cs_size, cs_top_pip = i.cs_top_pip)]
LC[is.na(n_cs) & paste(region_id, model) %in% RUN[susie_ok %in% TRUE, paste(region_id, model)], n_cs := 0L]
LC[PIP, on = .(region_id, model, SNP), lead_pip := i.pip]
LC[, max_pip_near := mapply(function(reg, m, b) { x <- PIP[region_id == reg & model == m & abs(bp - b) <= 2.5e5, pip]; if (length(x)) max(x) else NA_real_ },
                            region_id, model, bp)]
pc_long <- rbind(PC[, .(region_id, model = model_a, SNP = lead_a, other = model_b, pair_call)], PC[, .(region_id, model = model_b, SNP = lead_b, other = model_a, pair_call)])
pcs <- pc_long[, .(n_partner = uniqueN(other), n_same = sum(pair_call == "같은 신호"), n_diff = sum(pair_call == "다른 신호"),
                   n_unclear = sum(pair_call %in% c("불분명", "한쪽 신호 포함")), same_with_main = any(pair_call == "같은 신호" & other == "main")),
               by = .(region_id, model, SNP)]
LC[pcs, on = .(region_id, model, SNP), `:=`(n_partner = i.n_partner, n_same = i.n_same, n_diff = i.n_diff, n_unclear = i.n_unclear, same_with_main = i.same_with_main)]
for (v in c("n_partner", "n_same", "n_diff", "n_unclear")) set(LC, which(is.na(LC[[v]])), v, 0L)
LC[is.na(same_with_main), same_with_main := FALSE]
LC[, partner_call := fifelse(n_partner == 0, "혼자", fifelse(n_same > 0, "같은 신호 공유", fifelse(n_diff > 0, "다른 신호만", "불분명만")))]
LC[, `:=`(series = series_of(model), approx = model %in% approx_models)]
LC[, model := factor(model, series_order)]; setorder(LC, model, -thr, chr, bp); LC[, model := as.character(model)]
LCM <- LC[, .(n_loci = .N, pct_cojo_ge2 = 100 * mean(n_cojo >= 2), pct_cs = 100 * mean(n_cs >= 1, na.rm = TRUE),
              med_cs_size = if (any(n_cs >= 1, na.rm = TRUE)) as.numeric(median(cs_size, na.rm = TRUE)) else NA_real_,
              med_lead_pip = median(lead_pip, na.rm = TRUE), med_max_pip = median(max_pip_near, na.rm = TRUE),
              pct_alone = 100 * mean(partner_call == "혼자"), pct_share_same = 100 * mean(partner_call == "같은 신호 공유"),
              pct_diff_only = 100 * mean(partner_call == "다른 신호만"), pct_same_main = 100 * mean(same_with_main),
              series = series[1], approx = approx[1]), by = .(model, thr)]
LCM[, thr := factor(format(thr, scientific = TRUE), c("1e-05", "1e-04"))]; LCM[, model := factor(model, series_order)]; setorder(LCM, thr, model)

LC[cmp, on = .(model, SNP, thr), `:=`(het_age_p = i.het_age_p, het_sex_p = i.het_sex_p, het_ancestry_p = i.het_ancestry_p, het_weight_p = i.het_weight_p,
                                    het_weight_trend_p = i.het_weight_trend_p, int_age_p = i.int_age_p, int_age_b = i.int_age_b, het_p_wc = i.het_p_wc,
                                    b_disc = i.b, family = i.family)]
for (m in unique(LC$model)) LC[model == m, a1_disc := G[[m]][.SD, on = "SNP", x.A1]]
LC[, s_disc := sign(b_disc)]
z_str <- c("young8", "age8", "sex_m", "sex_f", "purelike", "mixlike", "wt_t1", "wt_t2", "wt_t3")
for (sm in z_str) {
  tmp <- G[[sm]][, .(SNP, A1s = A1, zs = b / se)]
  LC[tmp, on = "SNP", `:=`(zz = i.zs, aa = i.A1s)]
  LC[, (paste0("z_", sm)) := s_disc * fifelse(aa == a1_disc, zz, -zz)]; LC[, c("zz", "aa") := NULL]
}
wc_fam <- function(m, f) m %in% c(age = list(c("young8", "age8", "agefloor9", "age10", "agefloor11", "agefloor12")), sex = list(c("sex_m", "sex_f")),
                                  genetic = list(c("purelike", "mixlike")), weight = list(c("wt_t1", "wt_t2", "wt_t3")))[[f]]
LC[, `:=`(p_age = fifelse(wc_fam(model, "age"), het_p_wc, het_age_p), p_sex = fifelse(wc_fam(model, "sex"), het_p_wc, het_sex_p),
          p_anc = fifelse(wc_fam(model, "genetic"), het_p_wc, het_ancestry_p), p_wt = fifelse(wc_fam(model, "weight"), het_p_wc, het_weight_p),
          p_int = fifelse(model %in% c("int_age", age_series[age_series != "main"]), NA_real_, int_age_p))]
LC[, sel_bias := model %in% c("agefloor6", "agefloor7", "int_age")]
LC[, side_age := fifelse(z_age8 >= z_young8, "8세+ 쪽", "8세 미만 쪽")]
LC[, side_sex := fifelse(z_sex_f >= z_sex_m, "암컷 쪽", "수컷 쪽")]
LC[, side_anc := fifelse(z_purelike >= z_mixlike, "순종형 쪽", "잡종형 쪽")]
LC[, side_wt := c("wt_t1 쪽", "wt_t2 쪽", "wt_t3 쪽")[max.col(cbind(z_wt_t1, z_wt_t2, z_wt_t3), "first")]]
LC[, side_int := fifelse(sign(int_age_b) == s_disc, "나이 들수록 커짐", "나이 들수록 작아짐")]
lab_ax <- function(p, side, nm) fifelse(!is.na(p) & p < 0.05, paste0(nm, ": ", side), NA_character_)
LC[, spec_axes := apply(cbind(lab_ax(p_age, side_age, "나이"), lab_ax(p_wt, side_wt, "체중"), lab_ax(p_sex, side_sex, "성별"),
                              lab_ax(p_anc, side_anc, "유전 배경"), lab_ax(p_int, side_int, "나이 상호작용")), 1, function(v) paste(na.omit(v), collapse = " · "))]
LC[, n_axes := (!is.na(p_age) & p_age < 0.05) + (!is.na(p_wt) & p_wt < 0.05) + (!is.na(p_sex) & p_sex < 0.05) + (!is.na(p_anc) & p_anc < 0.05) +
                (!is.na(p_int) & p_int < 0.05)]
LC[, specass := fifelse(n_axes == 0, "층 차이 없음 (공통)", fifelse(n_axes >= 2, "여러 축",
                    fifelse(!is.na(p_age) & p_age < 0.05, paste0("나이 (", side_age, ")"), fifelse(!is.na(p_wt) & p_wt < 0.05, "체중",
                    fifelse(!is.na(p_sex) & p_sex < 0.05, "성별", fifelse(!is.na(p_anc) & p_anc < 0.05, "유전 배경", "나이 상호작용"))))))]
LC[, shared := fifelse(partner_call == "같은 신호 공유", "공유", fifelse(partner_call == "불분명만", "불분명", "공유 안 됨"))]
spec_lev <- c("층 차이 없음 (공통)", "나이 (8세+ 쪽)", "나이 (8세 미만 쪽)", "나이 상호작용", "체중", "성별", "유전 배경", "여러 축")
SPM <- LC[, c(list(n_loci = .N), setNames(lapply(spec_lev, function(k) sum(specass == k)), spec_lev),
              list(ref_chance_per_axis = round(0.05 * .N, 1), n_sel_bias = sum(sel_bias))), by = .(model, thr, shared)]
SPM[, `:=`(series = series_of(model), thr = factor(format(thr, scientific = TRUE), c("1e-05", "1e-04")))]
SPM[, model := factor(model, series_order)]; setorder(SPM, thr, model, shared)

gene_rows <- function(cc, lo, hi) {
  g <- gl[chr == cc & end >= lo & start <= hi][order(start)]; if (!nrow(g)) return(g)
  pad <- (hi - lo) * 0.06; row_end <- numeric(0); r <- integer(nrow(g))
  for (k in seq_len(nrow(g))) { j <- which(row_end < g$start[k] - pad)
    if (length(j)) { r[k] <- j[1]; row_end[j[1]] <- g$end[k] + pad } else { row_end <- c(row_end, g$end[k] + pad); r[k] <- length(row_end) } }
  g[, row := r][row <= 6]
}
cs_cols <- c("0" = "grey75", "1" = "#B2182B", "2" = "#2166AC", "3" = "#E17C05", "4" = "#1B7837", "5" = "#762A83", "6" = "#8C510A", "7" = "#01665E",
             "8" = "#C51B7D", "9" = "#4D4D4D", "10" = "#35978F")
plot_region <- function(reg) {
  x <- R[region_id == reg]; pans <- intersect(strsplit(x$panels, ",")[[1]], RUN[region_id == reg & susie_ok %in% TRUE, model]); cc <- x$chr
  pp <- PIP[region_id == reg]; if (!nrow(pp)) return(data.table(region_id = reg, plotted = FALSE))
  rn <- RUN[region_id == reg]
  strip <- rn[, .(model, lbl = sprintf("%s\nCOJO %d · CS %s%s", model, n_cojo, ifelse(is.na(n_cs), "NA", n_cs),
                                       ifelse(s_caution, sprintf("\ns=%.2f (caution)", s_rss), "")))]
  lev <- c(strip$lbl[match(pans, strip$model)], "genes")
  pp[strip, on = "model", fac := i.lbl]; pp[, fac := factor(fac, lev)]; pp[, csf := factor(as.character(cs), names(cs_cols))]
  co <- SL[region_id == reg & !is.na(SNP), .(model, SNP)][pp, on = .(model, SNP), nomatch = NULL]
  gr <- gene_rows(cc, x$lo, x$hi); if (nrow(gr)) gr[, fac := factor("genes", lev)]
  p <- ggplot() + geom_vline(xintercept = x$index_bp, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_point(data = pp[cs == 0], aes(bp, pip), colour = "grey75", size = 0.6) +
    geom_point(data = pp[cs > 0], aes(bp, pip, colour = csf), size = 1.4) +
    geom_point(data = co, aes(bp, pip), shape = 24, size = 2.4, colour = "black", fill = NA) +
    scale_colour_manual(values = cs_cols, name = "SuSiE CS", drop = TRUE) +
    facet_grid(fac ~ ., scales = "free_y", switch = "y") +
    scale_x_continuous(labels = function(b) sprintf("%.2f", b / 1e6), limits = c(x$lo, x$hi), expand = expansion(mult = 0.01)) +
    labs(title = sprintf("%s  chr%d:%.2f Mb  index %s  [13-3: %s]", reg, cc, x$index_bp / 1e6, x$index, RG[region_id == reg, region_call]),
         subtitle = "y = SuSiE PIP (LD = each model's own GWAS sample); colour = credible set; open triangle = COJO-selected SNP (p < 1e-4); dashed = index",
         x = sprintf("chr%d position (Mb)", cc), y = "PIP") +
    theme_bw(base_size = 9) + theme(strip.text.y.left = element_text(angle = 0, size = 7), strip.placement = "outside", panel.grid.minor = element_blank(),
                                    plot.subtitle = element_text(size = 7), plot.title = element_text(size = 10))
  if (nrow(gr)) p <- p + geom_segment(data = gr, aes(x = pmax(start, x$lo), xend = pmin(end, x$hi), y = -row, yend = -row), linewidth = 1.2, colour = "navy") +
    geom_text(data = gr, aes(x = pmin(pmax((pmax(start, x$lo) + pmin(end, x$hi)) / 2, x$lo + 0.04 * (x$hi - x$lo)), x$hi - 0.04 * (x$hi - x$lo)), y = -row + 0.45, label = sym),
              size = 1.9, colour = "navy")
  ggsave(file.path(fig_dir, paste0(reg, ".png")), p, width = 10, height = 1.2 * (length(pans) + 1) + 1.2, dpi = 200, limitsize = FALSE)
  data.table(region_id = reg, plotted = TRUE)
}
cat(sprintf("[%s] 그림 %d장\n", format(Sys.time(), "%T"), nrow(R)))
PS <- rbindlist(mclapply(regs, function(reg) tryCatch(plot_region(reg), error = function(e) data.table(region_id = reg, plotted = FALSE, error = conditionMessage(e))),
                         mc.cores = min(n_cores, length(regs))), fill = TRUE)
if ("error" %in% names(PS) && any(!is.na(PS$error))) { print(PS[!is.na(error)]); stop("그림 실패") }
xt <- PC[, .N, by = .(ld_signal, pair_call)]
xt[, ld_signal := factor(ld_signal, c("같은 신호", "불분명", "다른 신호 후보"), c("same (r2 >= 0.5)", "unclear", "distinct candidate (r2 < 0.1)"))]
xt[, pair_call := factor(pair_call, c("같은 신호", "한쪽 신호 포함", "불분명", "다른 신호"), c("same", "one nested", "unclear", "distinct"))]
g2 <- ggplot(xt, aes(pair_call, ld_signal, fill = N)) + geom_tile(colour = "white") + geom_text(aes(label = N), size = 3.5) +
  scale_fill_gradient(low = "#F0F0F0", high = "#46ACC8") +
  labs(x = "13-3 conditional analysis (COJO-cond, both directions)", y = "13-1 lead-pair LD", title = "Lead pairs in the same region: LD vs conditional analysis",
       subtitle = sprintf("%d pairs; remains: pC < %g, gone: pC >= %g or collinear", nrow(PC), p_keep, p_gone)) + theme_bw(base_size = 10)
ggsave(file.path(pic_dir, paste0("13-3", tag, "_pair_cond.png")), g2, width = 7.5, height = 3.8, dpi = 200)
lp4 <- melt(LCM[, .(model, thr, `% loci with >= 2 COJO signals` = pct_cojo_ge2, `% loci with SuSiE CS` = pct_cs, `median lead PIP` = med_lead_pip,
                    `% alone in region` = pct_alone, `% sharing a signal with another model` = pct_share_same, `% only distinct signals nearby` = pct_diff_only)],
            id.vars = c("model", "thr"), variable.factor = TRUE)
g4 <- ggplot(lp4, aes(model, value, colour = thr, group = thr)) + geom_point(size = 1.8) + geom_line(data = lp4[model %in% age_series], linewidth = 0.4) +
  geom_line(data = lp4[model %in% wt_series], linewidth = 0.4) + facet_wrap(~variable, ncol = 2, scales = "free_y") +
  geom_vline(xintercept = c(length(age_series), length(age_series) + length(wt_series)) + 0.5, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c(`1e-05` = "#B2182B", `1e-04` = "#46ACC8"), name = "locus threshold") +
  labs(x = "age series (young8 < 8 y, main = all, minimum age 6 ~ 12 y) | weight tertiles | others", y = NULL,
       title = "13-1 loci per model: signal structure (COJO, SuSiE) and sharing with other models (descriptive)") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(pic_dir, paste0("13-3", tag, "_loci_by_model.png")), g4, width = 11, height = 9, dpi = 200)
spec_en <- c("no stratum difference (common)", "age: stronger in 8+ y", "age: stronger in < 8 y", "age interaction (continuous)", "weight", "sex",
             "breed background", ">= 2 axes")
sp5 <- LC[, .(model = factor(model, series_order), thr = factor(format(thr, scientific = TRUE), c("1e-05", "1e-04")),
              grp = factor(fifelse(shared == "공유 안 됨", "not shared (alone / distinct only)", "shared or unclear"),
                           c("not shared (alone / distinct only)", "shared or unclear")),
              cls = factor(spec_en[match(specass, spec_lev)], spec_en))]
g5 <- ggplot(sp5, aes(model, fill = cls)) + geom_bar(position = "fill", width = 0.8) + facet_grid(thr ~ grp) +
  geom_vline(xintercept = c(length(age_series), length(age_series) + length(wt_series)) + 0.5, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("grey80", "#B2182B", "#F4A582", "#E17C05", "#46ACC8", "#762A83", "#1B7837", "grey35"), name = NULL, drop = FALSE) +
  scale_y_continuous(labels = function(v) paste0(100 * v, "%")) +
  labs(x = "age series | weight tertiles | others", y = "share of loci",
       title = "Where is each 13-1 locus 'specific'? same-SNP stratum difference p < 0.05 (nominal; ~5% per axis expected by chance, see limitation)") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), legend.position = "bottom", plot.title = element_text(size = 9))
ggsave(file.path(pic_dir, paste0("13-3", tag, "_spec_by_model.png")), g5, width = 13, height = 8, dpi = 200)

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
fwrite(f3(SL), paste0(out_prefix, "_cojo_slct.csv")); fwrite(f3(RUN), paste0(out_prefix, "_runs.csv")); fwrite(f3(PC), paste0(out_prefix, "_pair_cond.csv"))
fwrite(f3(cs_sum), paste0(out_prefix, "_susie_cs.csv")); fwrite(f3(CSN), paste0(out_prefix, "_susie_cs_snps.csv"))
fwrite(f3(PIP), paste0(out_prefix, "_susie_pip.tsv.gz"), sep = "\t"); fwrite(f3(CP), paste0(out_prefix, "_cs_pairs.csv")); fwrite(f3(RG), paste0(out_prefix, "_regions.csv"))
fwrite(f3(LC), paste0(out_prefix, "_loci.csv")); fwrite(f3(LCM), paste0(out_prefix, "_loci_by_model.csv")); fwrite(SPM, paste0(out_prefix, "_spec_by_model.csv"))
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("13-3 조건부 분석 (COJO) · 미세 지도 (SuSiE) — 서술 (우연 기대 · 다중 검정 보정 없음)\n"); if (test_mode) cat("*** TEST=1 시험 (위치 4개). 해석하지 말 것 ***\n")
cat(sprintf("위치 %d개 (94곳 %d) · 작업 (위치 × 모형) %d개 · 모형 %d개 (근사: %s) · LD 표본 = 각 모형 GWAS 표본\n", nrow(R), sum(R$in94), nrow(J), uniqueN(J$model),
            paste(intersect(approx_models, J$model), collapse = ", ")))
cat(sprintf("\n■ COJO-slct (p < %g): 모형별 독립 신호 수 분포\n", cojo_p)); print(RUN[, .N, by = n_cojo][order(n_cojo)])
cat("모형별 (작업 수 · 신호 ≥ 2인 작업 수 · CS ≥ 2인 작업 수 · CS 0인 작업 수):\n")
print(RUN[, .(jobs = .N, cojo_ge2 = sum(n_cojo >= 2), cs_ge2 = sum(n_cs >= 2, na.rm = TRUE), cs0 = sum(n_cs == 0, na.rm = TRUE), approx = approx[1]),
          by = model][order(match(model, fam_order))], nrows = 50)
cat(sprintf("\n■ SuSiE (L = %d, CS %g%%, 순도 ≥ %.1f): CS 수 분포\n", susie_L, 100 * cs_cov, cs_purity)); print(RUN[, .N, by = n_cs][order(n_cs)])
cat(sprintf("수렴 안 됨 %d · 실패 %d · s ≥ %.1f (z와 LD 불일치 주의) %d · kriging 이상치 있는 작업 %d / %d\n", sum(RUN$susie_converged %in% FALSE),
            sum(RUN$susie_ok %in% FALSE), s_flag, sum(RUN$s_caution), sum(RUN$n_kriging_outlier > 0, na.rm = TRUE), nrow(RUN)))
cat("s 분포 (작업별):\n"); print(summary(RUN$s_rss))
cat(sprintf("CS 없는 작업 %d개의 최고 PIP: 중앙값 %.2f · 범위 %.2f ~ %.2f (CS가 없어도 신호는 있을 수 있음 — PIP가 LD로 묶인 여러 SNP에 나뉨)\n",
            sum(RUN$n_cs == 0, na.rm = TRUE), median(RUN[n_cs == 0, max_pip]), min(RUN[n_cs == 0, max_pip]), max(RUN[n_cs == 0, max_pip])))
if (nrow(cs_sum)) {
  cat(sprintf("CS %d개: 크기 중앙값 %g (≤ %d 좁음 %d개) · 순도 중앙값 %.2f · 최고 PIP 중앙값 %.2f · index 포함 %d · 그 모형 13-1 대표 SNP 포함 %d\n", nrow(cs_sum),
              median(cs_sum$cs_size), cs_narrow, sum(cs_sum$narrow), median(cs_sum$cs_purity), median(cs_sum$top_pip), sum(cs_sum$has_index), sum(cs_sum$has_own_lead)))
}
cat(sprintf("COJO 신호 수와 SuSiE CS 수가 같은 작업 %d / %d\n", sum(RUN$n_cojo == RUN$n_cs, na.rm = TRUE), sum(!is.na(RUN$n_cs))))
cat(sprintf("\n■ 모형 간 조건부 (대표 SNP 짝 %d개; 남음 pC < %g · 사라짐 pC ≥ %g 또는 공선)\n", nrow(PC), p_keep, p_gone))
print(dcast(PC[, .N, by = .(ld_signal, pair_call)], ld_signal ~ pair_call, value.var = "N", fill = 0L))
cat("방향별 상태 (짝 × 2):\n"); print(table(c(PC$status_a, PC$status_b), useNA = "ifany"))
if (nrow(CP)) { cat(sprintf("\n■ 모형 간 SuSiE CS 비교 (%d쌍)\n", nrow(CP))); print(CP[, .N, by = cs_call]) }
cat("\n■ 위치 판정 (여러 신호 = COJO ≥ 2 · CS ≥ 2 · 다른 신호 짝 중 하나; 한 신호; 신호 약함 = COJO · CS 모두 0)\n")
print(RG[, .N, by = region_call]); cat("13-1 다른 신호 후보 위치 (multi_signal) × 13-3 판정:\n"); print(table(RG$multi_signal_13_1, RG$region_call))
cat("94곳 (13-1 1 ~ 2순위):\n"); print(RG[in94 == TRUE, .N, by = region_call])
cat("\n94곳 위치별 (chr 순서; 전체 830곳은 _regions.csv):\n")
print(RG[in94 == TRUE, .(region_id, index, multi_13_1 = multi_signal_13_1, region_call, n_models, max_cojo, max_cs, models_cojo2, models_cs2, n_pairs, n_same, n_diff,
                         n_nested, n_unclear, n_s_caution)], nrows = 200)
cat(sprintf("\n■ locus별 · 모형별 경향 (⑧; 13-1 locus %d개 = 1e-5 %d · 1e-4 %d) — 서술\n", nrow(LC), sum(LC$thr == 1e-5), sum(LC$thr == 1e-4)))
cat("  pct_cojo_ge2 = 그 기준 COJO 독립 신호 ≥ 2인 locus % · pct_cs = SuSiE CS 있는 % · med_cs_size = CS 크기 중앙값 (SNP) · med_lead_pip = 대표 SNP PIP 중앙값 ·\n")
cat("  med_max_pip = 대표 SNP ± 250 kb 최고 PIP 중앙값 · pct_alone = 같은 위치에 다른 모형 대표 SNP 없음 · pct_share_same = 조건부로 같은 신호인 다른 모형 있음 ·\n")
cat("  pct_diff_only = 다른 모형이 있지만 모두 다른 신호 · pct_same_main = main과 같은 신호\n")
for (sr in c("나이", "체중", "기타")) {
  cat(sprintf("[%s 계열]%s\n", sr, switch(sr, 나이 = " young8 (8세 미만) → main (전체) → 6 · 7 · 8 · 9 · 10 · 11 · 12세 이상 순", 체중 = " wt_t1 (가벼움) → wt_t3 (무거움) 순", "")))
  print(f3(LCM[series == sr, !"series"]), nrows = 100)
}
cat("\n■ locus 특이성 프로필 (⑨; 같은 SNP의 층 간 차이 p < 0.05, 명목) — 서술. ref_chance_per_axis = 우연이면 축마다 표시될 locus 수 (참고, 판정에 안 씀)\n")
cat("  공유 = 같은 위치 다른 모형과 조건부로 같은 신호 · 공유 안 됨 = 혼자 또는 다른 신호만 (나이 하한이 높은 모형에서 늘어난 쪽) · n_sel_bias = agefloor6 · 7 · int_age (나이 축 선택 편향 가능)\n")
cat("  나이 계열 모형 · int_age의 locus에는 나이 상호작용 축을 쓰지 않음 (하위 표본에서 고른 SNP는 상호작용이 저절로 커짐)\n")
for (sr in c("나이", "체중", "기타")) {
  cat(sprintf("[%s 계열]\n", sr)); print(SPM[series == sr, !"series"], nrows = 200)
}
cat("기준 × 공유 여부 전체 (축별 표시 수; 우연 기대 = 0.05 × n, 참고):\n")
print(LC[, .(n = .N, ref_chance = round(0.05 * .N, 1), age_8plus = sum(p_age < 0.05 & side_age == "8세+ 쪽", na.rm = TRUE),
             age_under8 = sum(p_age < 0.05 & side_age == "8세 미만 쪽", na.rm = TRUE), int_age = sum(p_int < 0.05, na.rm = TRUE), weight = sum(p_wt < 0.05, na.rm = TRUE),
             sex = sum(p_sex < 0.05, na.rm = TRUE), ancestry = sum(p_anc < 0.05, na.rm = TRUE), none = sum(n_axes == 0)), by = .(thr, shared)][order(-thr, shared)])
cat("나이 계열 1e-4, 공유 안 된 locus의 나이 축 방향 (8세+ 쪽 / 8세 미만 쪽):\n")
print(LC[thr == 1e-4 & series == "나이" & shared == "공유 안 됨", .(n = .N, age_8plus = sum(p_age < 0.05 & side_age == "8세+ 쪽", na.rm = TRUE),
                                                                 age_under8 = sum(p_age < 0.05 & side_age == "8세 미만 쪽", na.rm = TRUE),
                                                                 int_up = sum(p_int < 0.05 & side_int == "나이 들수록 커짐", na.rm = TRUE),
                                                                 int_down = sum(p_int < 0.05 & side_int == "나이 들수록 작아짐", na.rm = TRUE)), by = model][order(match(model, series_order))])
cat("기준별 전체:\n")
print(f3(LC[, .(n_loci = .N, pct_cojo_ge2 = 100 * mean(n_cojo >= 2), pct_cs = 100 * mean(n_cs >= 1, na.rm = TRUE), med_lead_pip = median(lead_pip, na.rm = TRUE),
                pct_alone = 100 * mean(partner_call == "혼자"), pct_share_same = 100 * mean(partner_call == "같은 신호 공유"),
                pct_diff_only = 100 * mean(partner_call == "다른 신호만")), by = thr]))
sink()
unlink(wd, recursive = TRUE)
cat(sprintf("저장: %s_* (cojo_slct · runs · pair_cond · susie_cs · susie_cs_snps · susie_pip · cs_pairs · regions · loci · loci_by_model · spec_by_model · summary) / 그림 %s/ · 13-3%s_pair_cond.png · _loci_by_model.png · _spec_by_model.png\n",
            out_prefix, fig_dir, tag))

