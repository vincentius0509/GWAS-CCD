# 04-7_model_compare.R - Compares the 0-30 PC models on the real phenotype (lambda, loci, QQ tail) to choose five PCs. See README.md.
library(data.table)
library(ggplot2)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
gwas_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023/04-5_gwas")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
out_prefix <- file.path(qc_dir, "04-7_model")
pic_prefix <- file.path(pic_dir, "04-7_model")
tmp_dir <- file.path(qc_dir, "04-7_tmp")
bfile <- file.path(qc_dir, "03-3_maf0.01_hwe")
plink2 <- "/data/00_tools/plink2_linux_avx2_20231123/plink2"
plink19 <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
n_threads <- 32
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)

pc_sets <- c(0, 5, 10, 15, 20, 25, 30)
run_name <- function(k) sprintf("rint_maf0.01_PC%d", k)
lambda_range <- c(0.95, 1.05)
p_lead <- 1e-5
p_clump2 <- 1e-2
p_support <- 1e-3
clump_r2 <- 0.1
clump_kb <- 1000
locus_kb <- 500
n_region <- 8
region_kb <- 500

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))
pc_cols <- setNames(c("#85B6EE", "#6DA0E6", "#518CDE", "#2A78D6", "#1E5CA8", "#14427D", "#0A2A54"), paste0("PC", pc_sets))
col_lone <- "#EB6834"; col_chr <- c("#2A78D6", "#898781")
seq_ramp <- c("#86B6EF", "#5598E7", "#2A78D6", "#1C5CAB", "#104281", "#0D366B")

bim <- fread(paste0(bfile, ".bim"), header = FALSE, select = c(1, 2, 4), col.names = c("Chr", "SNP", "bp"))
snp <- NULL; P <- list(); B <- list(); SE <- list()
for (k in pc_sets) {
  f <- file.path(gwas_dir, paste0(run_name(k), ".mlma.gz"))
  if (!file.exists(f)) stop("파일 없음: ", f, " (04-5 PC_LIST 실행 확인)")
  d <- fread(cmd = paste("zcat", shQuote(f)), header = TRUE)
  if (!identical(names(d)[1:9], c("Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p"))) stop("헤더가 예상과 다름: ", f)
  if (nrow(d) != nrow(bim)) warning("PC", k, ": 결과 SNP ", nrow(d), " ≠ bim ", nrow(bim))
  if (is.null(snp)) {
    snp <- d[, .(Chr, SNP, bp, A1, A2, Freq)]
  } else if (!identical(d$SNP, snp$SNP) || !identical(d$A1, snp$A1)) {
    stop("PC", k, ": SNP 순서나 A1이 PC", pc_sets[1], "과 다름 → 모델 간 비교 전에 맞춰야 함")
  }
  bad <- d[, sum(!is.finite(p) | !is.finite(b) | !is.finite(se) | se <= 0)]
  if (bad > 0) warning("PC", k, ": p·b·se가 숫자가 아니거나 se ≤ 0인 SNP ", bad, "개 (λ·상관에서 제외)")
  P[[paste0("PC", k)]] <- d$p; B[[paste0("PC", k)]] <- d$b; SE[[paste0("PC", k)]] <- d$se
  cat("읽음:", basename(f), nrow(d), "SNP\n")
}
rm(d); invisible(gc())
models <- names(P)
snp[, maf := pmin(Freq, 1 - Freq)]
snp[, idx := .I]
if (any(diff(snp$Chr) < 0) || snp[, any(diff(bp) < 0), by = Chr][, any(V1)]) stop("결과가 염색체·위치 순으로 정렬되어 있지 않음")

tail_q <- c(0.1, 0.01, 1e-3, 1e-4)
maf_bins <- c(0.01, 0.05, 0.2, 0.5)
snp[, maf_bin := cut(maf, maf_bins, include.lowest = TRUE, right = FALSE,
                     labels = c("MAF 0.01-0.05", "MAF 0.05-0.2", "MAF 0.2-0.5"))]

lam_rows <- list()
for (m in models) {
  chi <- (B[[m]] / SE[[m]])^2; ok <- is.finite(chi)
  r <- data.table(model = m, n_pc = as.integer(sub("PC", "", m)), what = "all SNP (median)",
                  lambda = median(chi[ok]) / qchisq(0.5, 1), n_snp = sum(ok))
  for (q in tail_q) r <- rbind(r, data.table(model = m, n_pc = r$n_pc[1], what = sprintf("tail top %g", q),
                                              lambda = quantile(chi[ok], 1 - q, names = FALSE) / qchisq(1 - q, 1), n_snp = sum(ok)))
  for (bn in levels(snp$maf_bin)) {
    s <- ok & snp$maf_bin == bn
    r <- rbind(r, data.table(model = m, n_pc = r$n_pc[1], what = bn, lambda = median(chi[s]) / qchisq(0.5, 1), n_snp = sum(s)))
  }
  lam_rows[[m]] <- r
}
lam_detail <- rbindlist(lam_rows)
fwrite(lam_detail, paste0(out_prefix, "_lambda_detail.csv"))

sum_info <- rbindlist(lapply(pc_sets, function(k) {
  f <- file.path(gwas_dir, paste0(run_name(k), "_summary.txt"))
  if (!file.exists(f)) return(data.table(model = paste0("PC", k), lambda_045 = NA_real_, vgvp_mean = NA_real_, vgvp_min = NA_real_, vgvp_max = NA_real_))
  x <- readLines(f, encoding = "UTF-8")
  l <- as.numeric(sub(".*: ([0-9.]+) \\(SNP.*", "\\1", grep("^λ \\(median", x, value = TRUE)[1]))
  i <- grep("^염색체별 V\\(G\\)/Vp", x)
  v <- if (length(i)) suppressWarnings(as.numeric(sub("^chr[0-9]+ ", "", regmatches(x[i + 1], gregexpr("chr[0-9]+ [^ ]+", x[i + 1]))[[1]]))) else NA_real_
  if (length(i) && length(v) != 38) warning("PC", k, ": summary의 염색체별 V(G)/Vp가 38개가 아님 (", length(v), ")")
  data.table(model = paste0("PC", k), lambda_045 = l, vgvp_mean = mean(v, na.rm = TRUE), vgvp_min = min(v, na.rm = TRUE), vgvp_max = max(v, na.rm = TRUE))
}))
chk <- merge(lam_detail[what == "all SNP (median)", .(model, lambda)], sum_info, by = "model")
if (chk[!is.na(lambda_045), any(abs(lambda - lambda_045) > 1e-3)]) { print(chk); warning("04-5 summary λ와 다시 계산한 λ가 0.001 넘게 다름") }

run_cmd <- function(args, log_name, exe = plink2) {
  out <- suppressWarnings(system2(exe, args, stdout = TRUE, stderr = TRUE))
  writeLines(out, file.path(tmp_dir, log_name))
  st <- attr(out, "status")
  if (!is.null(st) && st != 0) stop("plink2 실패 (로그: ", file.path(tmp_dir, log_name), ")")
  invisible(out)
}
clump_rows <- list()
for (m in models) {
  p <- P[[m]]
  if (!any(p < p_lead, na.rm = TRUE)) { cat(m, ": p <", p_lead, "SNP 없음 → clump 없음\n"); next }
  inp <- file.path(tmp_dir, paste0(m, "_p.tsv"))
  fwrite(data.table(ID = snp$SNP, P = p)[P < p_clump2], inp, sep = "\t")
  outp <- file.path(tmp_dir, m)
  unlink(paste0(outp, ".clumps"))
  run_cmd(c("--bfile", bfile, "--dog", "--clump", inp, "--clump-id-field", "ID", "--clump-p-field", "P",
            "--clump-p1", p_lead, "--clump-p2", p_clump2, "--clump-r2", clump_r2, "--clump-kb", clump_kb,
            "--clump-unphased", "--threads", n_threads, "--out", outp), paste0(m, "_clump.log"))
  cf <- paste0(outp, ".clumps")
  if (!file.exists(cf)) stop(m, ": clump 결과 파일이 없음 (로그 확인)")
  cl <- fread(cf)
  setnames(cl, sub("^#", "", names(cl)))
  bin_cols <- grep("^S[0-9.e-]+$", names(cl), value = TRUE)
  bin_cut <- as.numeric(sub("^S", "", bin_cols))
  if (!length(bin_cols)) stop(m, ": clump 결과에 p값 구간 열(S0.001 등)이 없음")
  cl[, n_support := rowSums(.SD), .SDcols = bin_cols[bin_cut <= p_support]]
  sp_ids <- strsplit(gsub("\\([^)]*\\)", "", cl$SP2), ",")
  rng <- t(vapply(seq_len(nrow(cl)), function(i) {
    ids <- setdiff(sp_ids[[i]], "."); pos <- c(cl$POS[i], snp$bp[match(ids, snp$SNP)])
    c(min(pos, na.rm = TRUE), max(pos, na.rm = TRUE))
  }, numeric(2)))
  li <- match(cl$ID, snp$SNP)
  clump_rows[[m]] <- data.table(model = m, n_pc = as.integer(sub("PC", "", m)), Chr = cl$CHROM, lead = cl$ID, bp = cl$POS,
                                p = cl$P, b = B[[m]][li], se = SE[[m]][li], maf = snp$maf[li],
                                n_clumped = cl$TOTAL, n_support = cl$n_support,
                                start = rng[, 1], end = rng[, 2], span_kb = round((rng[, 2] - rng[, 1]) / 1e3, 1),
                                lone = cl$n_support == 0)
}
clumps <- rbindlist(clump_rows)
setorder(clumps, n_pc, p)
fwrite(clumps, paste0(out_prefix, "_clumps.csv"))

head(clumps)

loci <- data.table()
if (nrow(clumps)) {
  cand <- clumps[order(p)]
  cand[, locus := NA_integer_]
  centers <- data.table(locus = integer(), Chr = integer(), bp = numeric(), lead = character(), lead_model = character(), lead_p = numeric())
  for (i in seq_len(nrow(cand))) {
    hit <- centers[Chr == cand$Chr[i] & abs(bp - cand$bp[i]) <= locus_kb * 1e3]
    if (nrow(hit)) { cand$locus[i] <- hit$locus[which.min(abs(hit$bp - cand$bp[i]))]; next }
    id <- nrow(centers) + 1L
    centers <- rbind(centers, data.table(locus = id, Chr = cand$Chr[i], bp = cand$bp[i], lead = cand$lead[i], lead_model = cand$model[i], lead_p = cand$p[i]))
    cand$locus[i] <- id
  }
  centers[, locus_id := sprintf("L%02d", locus)]
  centers[, label := sprintf("%s chr%d:%.2f Mb", locus_id, Chr, bp / 1e6)]
  writeLines(centers$lead, file.path(tmp_dir, "locus_leads.txt"))
  unlink(file.path(tmp_dir, "locus.ld"))
  run_cmd(c("--bfile", bfile, "--dog", "--r2", "--ld-snp-list", file.path(tmp_dir, "locus_leads.txt"),
            "--ld-window-kb", max(locus_kb, region_kb), "--ld-window", 999999, "--ld-window-r2", 0, "--threads", n_threads,
            "--out", file.path(tmp_dir, "locus")), "locus_ld.log", exe = plink19)
  ld <- fread(file.path(tmp_dir, "locus.ld"))[, .(ID_A = SNP_A, ID_B = SNP_B, R2)]
  if (!all(centers$lead %in% ld$ID_A)) warning("r²가 계산되지 않은 lead SNP: ", paste(setdiff(centers$lead, ld$ID_A), collapse = ", "))
  rows <- list()
  for (j in seq_len(nrow(centers))) {
    w <- snp[Chr == centers$Chr[j] & abs(bp - centers$bp[j]) <= locus_kb * 1e3, idx]
    li <- match(centers$lead[j], snp$SNP)
    r2l <- ld[ID_A == centers$lead[j]]
    r2v <- r2l$R2[match(snp$SNP[w], r2l$ID_B)]
    hi <- w[!is.na(r2v) & r2v >= 0.5 & w != li]
    for (m in models) {
      pm <- P[[m]][w]; jm <- which.min(pm)
      own <- cand[locus == centers$locus[j] & model == m]
      lead_lp <- -log10(P[[m]][li]); ld_lp <- if (length(hi)) max(-log10(P[[m]][hi])) else NA_real_
      rows[[length(rows) + 1]] <- data.table(locus_id = centers$locus_id[j], label = centers$label[j], model = m,
                                             n_pc = as.integer(sub("PC", "", m)),
                                             min_p = pm[jm], min_p_snp = snp$SNP[w[jm]], min_p_bp = snp$bp[w[jm]],
                                             lead_p = P[[m]][li], lead_b = B[[m]][li], lead_se = SE[[m]][li],
                                             has_lead = nrow(own) > 0, lead_lone = if (nrow(own)) own$lone[1] else NA,
                                             n_ld05 = length(hi), ld05_max_lp = ld_lp, gap_lp = lead_lp - ld_lp)
    }
  }
  loci <- rbindlist(rows)
  fwrite(loci, paste0(out_prefix, "_loci.csv"))
  wide <- dcast(loci, locus_id + label ~ model, value.var = "min_p")
  setcolorder(wide, c("locus_id", "label", models))
  rep_n <- loci[, .(n_model_1e5 = sum(min_p < 1e-5), n_model_1e4 = sum(min_p < 1e-4), n_lead = sum(has_lead),
                    n_lone = sum(lead_lone, na.rm = TRUE), best_p = min(min_p),
                    sign_agree = all(sign(lead_b) == sign(lead_b[model == models[1]])),
                    n_ld05 = n_ld05[1], gap_min = min(gap_lp), gap_max = max(gap_lp)), by = locus_id]
  wide <- merge(wide, rep_n, by = "locus_id")
  wide <- merge(wide, centers[, .(locus_id, lead, lead_model, lead_maf = snp$maf[match(lead, snp$SNP)])], by = "locus_id")
  setorder(wide, best_p)
  fwrite(wide, paste0(out_prefix, "_loci_wide.csv"))
}

Z <- sapply(models, function(m) B[[m]] / SE[[m]])
zc <- cor(Z, use = "pairwise.complete.obs")
fwrite(data.table(model = models, round(zc, 4)), paste0(out_prefix, "_z_cor.csv"))
b_cor <- NULL
if (nrow(loci) >= 3 * length(models)) {
  ref_b <- loci[model == models[1], setNames(lead_b, locus_id)]
  b_cor <- loci[, .(b_cor_vs_first = cor(lead_b, ref_b[locus_id])), by = model]
}

ms <- lam_detail[, .(lambda = lambda[what == "all SNP (median)"],
                     lambda_tail_1e3 = lambda[what == "tail top 0.001"],
                     lambda_maf_low = lambda[what == "MAF 0.01-0.05"],
                     lambda_maf_high = lambda[what == "MAF 0.2-0.5"]), by = .(model, n_pc)]
ms[, lambda_pass := lambda >= lambda_range[1] & lambda <= lambda_range[2]]
ms[, `:=`(n_p5e8 = sapply(model, function(m) sum(P[[m]] < 5e-8, na.rm = TRUE)),
          n_p1e6 = sapply(model, function(m) sum(P[[m]] < 1e-6, na.rm = TRUE)),
          n_p1e5 = sapply(model, function(m) sum(P[[m]] < 1e-5, na.rm = TRUE)),
          n_p1e4 = sapply(model, function(m) sum(P[[m]] < 1e-4, na.rm = TRUE)))]
if (nrow(clumps)) {
  cs <- clumps[, .(n_locus = .N, n_lone = sum(lone), median_support = as.numeric(median(n_support))), by = model]
  ms <- merge(ms, cs, by = "model", all.x = TRUE)
  rp <- loci[has_lead == TRUE, .(locus_id, model)][, .(n_locus_rep_all = sum(sapply(locus_id, function(l)
    all(loci[locus_id == l & model != .BY$model, min_p] < 1e-4)))), by = model]
  ms <- merge(ms, rp, by = "model", all.x = TRUE)
}
ms <- merge(ms, sum_info[, .(model, vgvp_mean, vgvp_min, vgvp_max)], by = "model", all.x = TRUE)
if (!is.null(b_cor)) ms <- merge(ms, b_cor, by = "model", all.x = TRUE)
ms[, z_cor_vs_first := zc[model, models[1]]]
for (cc in c("n_locus", "n_lone", "n_locus_rep_all")) if (cc %in% names(ms)) set(ms, which(is.na(ms[[cc]])), cc, 0L)
setorder(ms, n_pc)
fwrite(ms, paste0(out_prefix, "_model_summary.csv"))

pc_lab <- setNames(sprintf("%s (lambda %.3f)", ms$model, ms$lambda), ms$model)

n_all <- length(P[[1]])
keep_i <- sort(unique(c(1:20000, round(exp(seq(log(20001), log(n_all), length.out = 4000))))))
keep_i <- keep_i[keep_i <= n_all]
qq <- rbindlist(lapply(models, function(m) {
  ps <- sort(P[[m]][is.finite(P[[m]])]); n <- length(ps); i <- keep_i[keep_i <= n]
  data.table(model = m, exp = -log10((i - 0.5) / n), obs = -log10(ps[i]))
}))
band <- data.table(i = keep_i)[, `:=`(exp = -log10((i - 0.5) / n_all),
                                       lo = -log10(qbeta(0.975, i, n_all - i + 1)), hi = -log10(qbeta(0.025, i, n_all - i + 1)))]
qq2 <- rbind(qq[, panel := "All SNPs"], qq[exp > 3][, panel := "Tail (expected -log10 p > 3)"])
band2 <- rbind(band[, panel := "All SNPs"], band[exp > 3][, panel := "Tail (expected -log10 p > 3)"])
qq2[, model := factor(model, levels = models)]
g <- ggplot() +
  geom_ribbon(data = band2, aes(x = exp, ymin = lo, ymax = hi), fill = "#E1E0D9") +
  geom_abline(slope = 1, intercept = 0, color = "#898781", linewidth = 0.4) +
  geom_line(data = qq2, aes(exp, obs, color = model), linewidth = 0.6) +
  facet_wrap(~panel, scales = "free") +
  scale_color_manual(values = pc_cols, labels = pc_lab, name = "Model") +
  labs(title = "QQ plot: CSLB RINT, GCTA MLMA-LOCO, 7 PC settings",
       x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p))) +
  theme_qc + theme(legend.position = "right")
ggsave(paste0(pic_prefix, "_qq.png"), g, width = 13, height = 6, dpi = 300)

n_xbin <- 4000; y_bin <- 0.01
chr_len <- snp[, .(len = max(bp)), by = Chr][order(Chr)][, off := cumsum(as.numeric(shift(len, fill = 0)))]
snp[chr_len, x := bp + i.off, on = "Chr"]
x_bin <- floor(snp$x / (max(snp$x) / n_xbin))
mh <- rbindlist(lapply(models, function(m) {
  lp <- -log10(P[[m]])
  d <- data.table(i = seq_along(lp), Chr = snp$Chr, xb = x_bin, yb = floor(lp / y_bin), lp = lp)
  keep <- d[lp >= 3 | !duplicated(d, by = c("Chr", "xb", "yb")), i]
  data.table(model = m, x = snp$x[keep], Chr = snp$Chr[keep], lp = lp[keep])
}))
mh[, model := factor(model, levels = models, labels = pc_lab[models])]
ld_pts <- if (nrow(clumps)) clumps[, .(model = factor(model, levels = models, labels = pc_lab[models]),
                                        x = snp$x[match(lead, snp$SNP)], lp = -log10(p), lone)] else NULL
axis_df <- chr_len[, .(Chr, mid = off + len / 2)]
g <- ggplot(mh, aes(x, lp)) +
  geom_point(aes(color = factor(Chr %% 2)), size = 0.35, show.legend = FALSE) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "#52514E", linewidth = 0.4) +
  geom_hline(yintercept = -log10(1e-5), linetype = "dotted", color = "#52514E", linewidth = 0.4) +
  scale_color_manual(values = c("0" = col_chr[1], "1" = col_chr[2])) +
  scale_x_continuous(breaks = axis_df$mid, labels = ifelse(axis_df$Chr %% 2 == 1 | axis_df$Chr <= 10, axis_df$Chr, ""), expand = c(0.005, 0)) +
  facet_wrap(~model, ncol = 1, strip.position = "right") +
  labs(title = "Manhattan: 7 PC settings (orange triangle = lone lead, no p < 1e-3 SNP in its LD clump)",
       subtitle = "all 4.84M SNPs drawn; points overlapping on screen merged (4,000 x-bins, 0.01 -log10 p bins), all p < 1e-3 kept",
       x = "Chromosome", y = expression(-log[10](p))) +
  theme_qc + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = 10))
if (!is.null(ld_pts) && any(ld_pts$lone)) g <- g + geom_point(data = ld_pts[lone == TRUE], aes(x, lp), shape = 24, size = 2.2, fill = col_lone, color = "#52514E", stroke = 0.3)
ggsave(paste0(pic_prefix, "_manhattan.png"), g, width = 14, height = 17, dpi = 300)

if (nrow(loci)) {
  hm <- copy(loci)
  ord <- wide$label
  hm[, label := factor(label, levels = rev(ord))]
  hm[, model := factor(model, levels = models)]
  hm[, lp := -log10(min_p)]
  hm[, txt := sprintf("%.1f%s", lp, ifelse(!is.na(lead_lone) & lead_lone, "*", ""))]
  lp_mid <- min(hm$lp) + 0.5 * diff(range(hm$lp))
  g <- ggplot(hm, aes(model, label, fill = lp)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = txt, color = lp > lp_mid), size = 3.2, show.legend = FALSE) +
    scale_color_manual(values = c("FALSE" = "#0B0B0B", "TRUE" = "#FFFFFF")) +
    scale_fill_gradientn(colors = seq_ramp, name = expression(min~-log[10](p))) +
    labs(title = "Locus replication across PC settings", subtitle = "cell = smallest p within +/-500 kb of locus lead; * = that model's lead is a lone point",
         x = NULL, y = NULL) +
    theme_qc + theme(panel.grid = element_blank())
  ggsave(paste0(pic_prefix, "_locus_heatmap.png"), g, width = 10, height = max(4, 0.32 * nrow(wide) + 2), dpi = 300)
}

if (nrow(loci)) {
  top <- wide[seq_len(min(n_region, nrow(wide)))]
  rg <- rbindlist(lapply(seq_len(nrow(top)), function(j) {
    li <- match(top$lead[j], snp$SNP)
    w <- snp[Chr == snp$Chr[li] & abs(bp - snp$bp[li]) <= region_kb * 1e3, idx]
    r2 <- ld[ID_A == top$lead[j]]
    r2v <- r2$R2[match(snp$SNP[w], r2$ID_B)]
    r2v[w == li] <- 1
    rbindlist(lapply(models, function(m) data.table(label = top$label[j], model = m, pos = (snp$bp[w] - snp$bp[li]) / 1e3,
                                                    lp = -log10(P[[m]][w]), r2 = r2v, is_lead = w == li)))
  }))
  rg[, label := factor(label, levels = top$label)]
  rg[, model := factor(model, levels = models)]
  setorder(rg, r2, na.last = FALSE)
  g <- ggplot(rg, aes(pos, lp)) +
    geom_point(aes(color = r2), size = 0.6) +
    geom_point(data = rg[is_lead == TRUE], shape = 23, size = 2.2, fill = col_lone, color = "#0B0B0B", stroke = 0.3) +
    geom_hline(yintercept = -log10(1e-5), linetype = "dotted", color = "#52514E", linewidth = 0.3) +
    scale_color_gradientn(colors = seq_ramp, limits = c(0, 1), na.value = "#C9C8C1", name = expression(r^2~with~lead)) +
    facet_grid(label ~ model, scales = "free_y") +
    labs(title = sprintf("Top %d loci: regional association (+/-%d kb), diamond = locus lead SNP", nrow(top), region_kb),
         x = "Position relative to locus lead (kb)", y = expression(-log[10](p))) +
    theme_qc + theme(strip.text.y = element_text(size = 8, angle = 0), strip.text.x = element_text(size = 10),
                     axis.text = element_text(size = 7), panel.grid.minor = element_blank())
  ggsave(paste0(pic_prefix, "_regional.png"), g, width = 22, height = 2.3 * nrow(top) + 1.5, dpi = 250, limitsize = FALSE)
}

y_max <- max(-log10(5e-8), max(mh$lp), max(qq$obs)) + 0.6
x_max_qq <- max(qq$exp) * 1.02
lab_loc <- if (exists("cand") && nrow(clumps)) merge(cand[, .(model, lead, p, locus)], centers[, .(locus, locus_id)], by = "locus") else NULL
for (m in models) {
  lab <- pc_lab[[m]]
  d_mh <- mh[model == lab]
  gm <- ggplot(d_mh, aes(x, lp)) +
    geom_point(aes(color = factor(Chr %% 2)), size = 0.4, show.legend = FALSE) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "#52514E", linewidth = 0.4) +
    geom_hline(yintercept = -log10(1e-5), linetype = "dotted", color = "#52514E", linewidth = 0.4) +
    scale_color_manual(values = c("0" = col_chr[1], "1" = col_chr[2])) +
    scale_x_continuous(breaks = axis_df$mid, labels = ifelse(axis_df$Chr %% 2 == 1 | axis_df$Chr <= 10, axis_df$Chr, ""), expand = c(0.005, 0)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(title = "Manhattan", subtitle = "dashed p = 5e-8, dotted p = 1e-5; overlapping points merged",
         x = "Chromosome", y = expression(-log[10](p))) +
    theme_qc + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                     plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
  if (!is.null(lab_loc) && nrow(lab_loc[model == m])) {
    lm <- lab_loc[model == m][order(p)][!duplicated(locus_id)]
    lm[, `:=`(x = snp$x[match(lead, snp$SNP)], lp = -log10(p))]
    gm <- gm + geom_text(data = lm, aes(x, lp, label = locus_id), nudge_y = 0.3, size = 3, color = "#0B0B0B")
  }
  if (!is.null(ld_pts) && nrow(ld_pts[model == lab & lone == TRUE]))
    gm <- gm + geom_point(data = ld_pts[model == lab & lone == TRUE], aes(x, lp), shape = 24, size = 2.4, fill = col_lone, color = "#52514E", stroke = 0.3)
  d_qq <- qq[model == m]
  lt <- ms[model == m]
  gq <- ggplot() +
    geom_ribbon(data = band, aes(x = exp, ymin = lo, ymax = hi), fill = "#E1E0D9") +
    geom_abline(slope = 1, intercept = 0, color = "#898781", linewidth = 0.4) +
    geom_point(data = d_qq, aes(exp, obs), color = "#2A78D6", size = 0.6) +
    annotate("text", x = 0.2, y = y_max - 0.4, hjust = 0, vjust = 1, size = 3.6, color = "#0B0B0B",
             label = sprintf("lambda = %.3f\ntail lambda (top 0.1%%) = %.3f", lt$lambda, lt$lambda_tail_1e3)) +
    scale_x_continuous(limits = c(0, x_max_qq), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(title = "QQ", subtitle = "grey band = 95% null interval",
         x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p))) +
    theme_qc + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
  ttl <- sprintf("CSLB RINT, GCTA MLMA-LOCO, %s PCs  |  lambda %.3f (%s)  |  SNPs p < 5e-8 / 1e-5 / 1e-4: %d / %d / %d  |  loci (p < 1e-5): %d  |  lone leads: %d",
                 sub("PC", "", m), lt$lambda, ifelse(lt$lambda_pass, "within 0.95-1.05", "outside 0.95-1.05"),
                 lt$n_p5e8, lt$n_p1e5, lt$n_p1e4, if ("n_locus" %in% names(lt)) lt$n_locus else 0L, if ("n_lone" %in% names(lt)) lt$n_lone else 0L)
  gg <- gridExtra::arrangeGrob(gm, gq, ncol = 2, widths = c(3.3, 1),
                               top = grid::textGrob(ttl, gp = grid::gpar(fontsize = 15)))
  ggsave(sprintf("%s_model_PC%s.png", pic_prefix, sub("PC", "", m)), gg, width = 18, height = 5.5, dpi = 300)
}

sink(paste0(out_prefix, "_summary.txt"))
cat("04-7 모델 비교: CSLB RINT, GCTA MLMA-LOCO, PC", paste(pc_sets, collapse = "·"), "\n")
cat(sprintf("locus 정의: plink2 clump (index p < %g, r² ≥ %g, ±%d kb, 묶음 p < %g), 모델 간 합치기 ±%d kb\n",
            p_lead, clump_r2, clump_kb, p_clump2, locus_kb))
cat(sprintf("단독 점: lead p < %g인데 같은 clump 안에 p < %g SNP가 없음\n\n", p_lead, p_support))
cat("== 모델별 요약 (기준 1: λ, 기준 2: locus 재현, 기준 4: 단독 점) ==\n")
pr <- ms[, .(model, lambda = round(lambda, 4), pass = ifelse(lambda_pass, "통과", "불통과"),
             tail_1e3 = round(lambda_tail_1e3, 3), maf_low = round(lambda_maf_low, 3), maf_high = round(lambda_maf_high, 3),
             p5e8 = n_p5e8, p1e5 = n_p1e5, p1e4 = n_p1e4)]
if ("n_locus" %in% names(ms)) pr <- cbind(pr, ms[, .(locus = n_locus, lone = n_lone, rep_all = n_locus_rep_all)])
pr <- cbind(pr, ms[, .(vgvp = round(vgvp_mean, 4), z_cor_PC0 = round(z_cor_vs_first, 4))])
print(pr, row.names = FALSE)
cat("  p5e8·p1e5·p1e4 = 기준보다 p가 작은 SNP 수 (전처리가 엄격해 1e-4도 함께 봄, 연관 없을 때 1e-4 기대치 ≈ SNP 수 × 1e-4), tail_1e3 = 상위 0.1% 꼬리 λ, maf_low/high = MAF 0.01-0.05 / 0.2-0.5 구간 λ,",
    "rep_all = 이 모델 locus 중 나머지 6개 모델 모두에서 p < 1e-4인 수, vgvp = 염색체별 V(G)/Vp 평균\n\n")
if (nrow(loci)) {
  cat("== 합친 locus (best p 순, 칸 = locus ±500 kb 안 최소 p) ==\n")
  pw <- copy(wide)
  for (m in models) set(pw, j = m, value = sprintf("%.1e", pw[[m]]))
  print(pw[, c("label", models), with = FALSE], row.names = FALSE)
  cat("\n== 합친 locus: 재현 · LD 모양 ==\n")
  print(pw[, .(label, n_model_1e5, n_model_1e4, n_lone, sign_agree, lead_maf = round(lead_maf, 3), lead_model,
               n_ld05, gap_min = round(gap_min, 2), gap_max = round(gap_max, 2))], row.names = FALSE)
  cat("  n_model_1e5·1e4 = locus ±500 kb 안 최소 p가 기준보다 작은 모델 수, n_ld05 = lead와 r² ≥ 0.5인 SNP 수,\n",
      " gap = lead -log10 p − r² ≥ 0.5 SNP 중 최대 -log10 p (7개 모델 중 최소·최대; 크면 lead 하나만 튄 모양, NA = r² ≥ 0.5 SNP 없음)\n\n")
}
cat("== 모델 간 게놈 전체 z 상관 ==\n"); print(round(zc, 4))
sink()
cat(readLines(paste0(out_prefix, "_summary.txt")), sep = "\n")

