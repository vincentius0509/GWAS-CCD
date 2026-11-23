# 05-4_sens_compare.R - Compares every model with the primary model (lambda, loci, effect sizes). See README.md.
library(data.table)
library(ggplot2)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")    
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
out_prefix <- file.path(qc_dir, "05-4")
pic_prefix <- file.path(pic_dir, "05-4")
run_pic_dir <- paste0(pic_prefix, "_runs")
gcta_dir <- file.path(qc_dir, "05-2_gwas")
regenie_dir <- file.path(qc_dir, "05-3_regenie")
runs_tsv <- file.path(qc_dir, "05-1_runs.tsv")
ref_cand <- file.path(qc_dir, c("04-7_loci_wide.csv", "04-7_loci_wide.csv", "04-7_model_loci_wide.csv"))
dir.create(run_pic_dir, showWarnings = FALSE, recursive = TRUE)
locus_kb <- 500
make_run_plots <- TRUE
if (make_run_plots && !requireNamespace("gridExtra", quietly = TRUE)) {
  message("gridExtra 패키지 없음 → 모형별 그림(5-3)은 건너뜀 (install.packages(\"gridExtra\") 후 다시 실행)"); make_run_plots <- FALSE
}

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), axis.title = element_text(size = 13), axis.text = element_text(size = 11))
col_chr <- c("#2A78D6", "#898781")
seq_ramp <- c("#86B6EF", "#5598E7", "#2A78D6", "#1C5CAB", "#104281", "#0D366B")
group_order <- c("main", "qc", "phenotype", "age", "genetic", "strata", "drop_cov", "add_cov", "subset")
p_to_lambda <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }

runs <- fread(runs_tsv)[engine != "check"]
runs[, group := factor(group, levels = group_order)]
setorder(runs, group)
main_file <- runs[id == "main", ref_result]
if (!length(main_file) || !file.exists(main_file)) stop("주 분석 결과 없음: ", main_file)
read_mlma <- function(f) fread(cmd = paste("zcat", shQuote(f)), select = c("Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p"))
mn <- read_mlma(main_file)
if (any(diff(mn$Chr) < 0) || mn[, any(diff(bp) < 0), by = Chr][, any(V1)]) stop("주 분석 결과가 염색체·위치 순이 아님")
N <- nrow(mn)
N_main <- runs[id == "main", n]
cat("주 분석:", basename(main_file), N, "SNP\n")

ref_file <- ref_cand[file.exists(ref_cand)][1]
if (is.na(ref_file)) stop("04-7 기준 locus 파일 없음: ", paste(basename(ref_cand), collapse = " / "))
ref <- fread(ref_file)[, .(locus_id, label, lead)]
ref[, li := match(lead, mn$SNP)]
if (anyNA(ref$li)) stop("기준 locus lead SNP가 주 분석 결과에 없음")
ref[, `:=`(Chr = mn$Chr[li], bp = mn$bp[li])]
win <- lapply(seq_len(nrow(ref)), function(j) which(mn$Chr == ref$Chr[j] & abs(mn$bp - ref$bp[j]) <= locus_kb * 1e3))
in_ref <- logical(N); for (w in win) in_ref[w] <- TRUE
cat("기준 locus:", nrow(ref), "개 (", basename(ref_file), ")\n")

align_gcta <- function(d) {
  i <- match(mn$SNP, d$SNP)
  flip <- ifelse(d$A1[i] == mn$A1, 1, ifelse(d$A1[i] == mn$A2, -1, NA))
  list(p = d$p[i], z = flip * d$b[i] / d$se[i], b = flip * d$b[i])
}
read_regenie <- function(f, test_pat) {
  d <- fread(cmd = paste("zcat", shQuote(f)))
  tests <- unique(d$TEST)
  tsel <- grep(test_pat, tests, value = TRUE)[1]
  if (is.na(tsel)) stop(basename(f), ": TEST에 ", test_pat, " 없음 (있는 것: ", paste(tests, collapse = ", "), ")")
  d <- d[TEST == tsel & !is.na(LOG10P)]
  i <- match(mn$SNP, d$ID)
  flip <- ifelse(d$ALLELE1[i] == mn$A1, 1, ifelse(d$ALLELE1[i] == mn$A2, -1, NA))
  list(p = 10^(-d$LOG10P[i]), z = flip * d$BETA[i] / d$SE[i], b = flip * d$BETA[i], test = tsel)
}
res <- list(); status <- character()
for (k in seq_len(nrow(runs))) {
  r <- runs[k]
  out <- tryCatch({
    if (r$engine == "existing") {
      if (!file.exists(r$ref_result)) NULL else if (r$id == "main") list(p = mn$p, z = mn$b / mn$se, b = mn$b) else align_gcta(read_mlma(r$ref_result))
    } else if (r$engine == "gcta") {
      f <- file.path(gcta_dir, paste0(r$id, ".mlma.gz")); if (!file.exists(f)) NULL else align_gcta(read_mlma(f))
    } else if (r$engine == "regenie_bt") {
      f <- file.path(regenie_dir, paste0(r$id, ".regenie.gz")); if (!file.exists(f)) NULL else read_regenie(f, "^ADD$")
    } else if (r$engine == "regenie_int") {
      f <- file.path(regenie_dir, paste0(r$id, ".regenie.gz")); if (!file.exists(f)) NULL else read_regenie(f, "^ADD-INT_SNPx")
    } else if (r$engine == "maf_filter") {
      if (!file.exists(r$ref_result)) NULL else {
        keep <- mn$SNP %in% fread(r$ref_result, header = FALSE, select = 2)$V2
        list(p = ifelse(keep, mn$p, NA), z = ifelse(keep, mn$b / mn$se, NA), b = ifelse(keep, mn$b, NA))
      }
    }
  }, error = function(e) { message(r$id, ": 읽기 실패 - ", conditionMessage(e)); NULL })
  status[r$id] <- if (is.null(out)) "결과 없음" else "있음"
  if (!is.null(out)) res[[r$id]] <- out
  cat(sprintf("%-15s %-12s %s\n", r$id, r$engine, status[r$id]))
}
runs[, status := status[id]]
ok_ids <- runs[status == "있음", id]

z_main <- mn$b / mn$se
lead_main_z <- z_main[ref$li]
loc_rows <- list(); new_rows <- list(); sum_rows <- list()
for (rid in ok_ids) {
  r <- as.list(runs[which(runs$id == rid)]); x <- res[[rid]]
  is_int <- r$engine == "regenie_int"
  same_scale <- r$engine %in% c("gcta", "existing", "maf_filter") && r$y == "rint"
  lr <- rbindlist(lapply(seq_len(nrow(ref)), function(j) {
    w <- win[[j]]; pw <- x$p[w]
    jm <- if (all(is.na(pw))) NA_integer_ else which.min(pw)
    data.table(id = rid, locus_id = ref$locus_id[j], label = ref$label[j],
               min_p = if (is.na(jm)) NA_real_ else pw[jm], min_p_snp = if (is.na(jm)) NA_character_ else mn$SNP[w[jm]],
               lead_p = x$p[ref$li[j]], lead_z = x$z[ref$li[j]], main_lead_z = lead_main_z[j],
               lead_b = x$b[ref$li[j]], main_lead_b = mn$b[ref$li[j]])
  }))
  if (is_int) lr[, `:=`(lead_z = NA_real_, lead_b = NA_real_)]
  lr[, same_sign := sign(lead_z) == sign(main_lead_z)]
  loc_rows[[rid]] <- lr
  hit <- which(x$p < 1e-5 & !in_ref)
  nl <- data.table()
  if (length(hit)) {
    h <- data.table(i = hit, Chr = mn$Chr[hit], bp = mn$bp[hit], p = x$p[hit])
    h[, grp := cumsum(c(TRUE, diff(Chr) != 0 | diff(bp) > locus_kb * 1e3))]
    nl <- h[order(p), .(lead = mn$SNP[i[1]], Chr = Chr[1], bp = bp[1], p = p[1], n_snp_p1e5 = .N,
                        span_kb = round((max(bp) - min(bp)) / 1e3), maf = pmin(mn$Freq[i[1]], 1 - mn$Freq[i[1]])), by = grp][, id := rid]
    new_rows[[rid]] <- nl
  }
  bs <- if (same_scale) lr[!is.na(lead_b), sum(lead_b * main_lead_b) / sum(main_lead_b^2)] else NA_real_
  sum_rows[[rid]] <- data.table(
    id = rid, group = as.character(r$group), engine = r$engine, test = if (!is.null(x$test)) x$test else "", n = r$n,
    n_snp = sum(is.finite(x$p)), lambda = p_to_lambda(x$p),
    n_p5e8 = sum(x$p < 5e-8, na.rm = TRUE), n_p1e5 = sum(x$p < 1e-5, na.rm = TRUE), n_p1e4 = sum(x$p < 1e-4, na.rm = TRUE),
    z_cor_main = if (is_int) NA_real_ else cor(x$z, z_main, use = "pairwise.complete.obs"),
    z_cor_expect = if (r$engine %in% c("gcta", "existing", "maf_filter")) sqrt(min(r$n, N_main) / N_main) else NA_real_,
    ref_1e5 = lr[, sum(min_p < 1e-5, na.rm = TRUE)], ref_1e4 = lr[, sum(min_p < 1e-4, na.rm = TRUE)],
    ref_sign_agree = if (is_int) NA_real_ else lr[, mean(same_sign, na.rm = TRUE)],
    ref_z_ratio = if (is_int) NA_real_ else lr[, median(lead_z / main_lead_z, na.rm = TRUE)],
    ref_b_slope = bs, n_new_loci = nrow(nl),
    best_new = if (nrow(nl)) sprintf("chr%d:%.2fMb p=%.1e", nl$Chr[1], nl$bp[1] / 1e6, nl$p[1]) else "",
    desc = r$desc)
}
summ <- rbindlist(sum_rows)
summ <- rbind(summ, runs[status != "있음", .(id, group = as.character(group), engine, n, desc)], fill = TRUE)
summ[, group := factor(group, levels = group_order)]
setorder(summ, group)
loc <- rbindlist(loc_rows)
newl <- rbindlist(new_rows, fill = TRUE)
fwrite(summ, paste0(out_prefix, "_run_summary.csv"))
fwrite(loc, paste0(out_prefix, "_ref_loci.csv"))
if (nrow(loc)) {
  wide <- dcast(loc, locus_id + label ~ id, value.var = "min_p")
  setcolorder(wide, c("locus_id", "label", intersect(runs$id, names(wide))))
  fwrite(wide, paste0(out_prefix, "_ref_loci_minp.csv"))
}
fwrite(newl, paste0(out_prefix, "_new_loci.csv"))

sink(paste0(out_prefix, "_summary.txt"))
cat("05-4 민감도 분석 비교 (주 분석:", basename(main_file), "/ 기준 locus", nrow(ref), "개:", basename(ref_file), ")\n")
cat("결과 있음", length(ok_ids), "/ 전체", nrow(runs), "모형\n\n")
options(width = 250)
print(summ[, .(id, group, n, lambda = round(lambda, 3), p5e8 = n_p5e8, p1e5 = n_p1e5, p1e4 = n_p1e4, z_cor = round(z_cor_main, 3), z_exp = round(z_cor_expect, 2),
               ref_1e5, ref_1e4, sign = round(ref_sign_agree, 2), z_ratio = round(ref_z_ratio, 2), b_slope = round(ref_b_slope, 2),
               new = n_new_loci, best_new)], row.names = FALSE, nrows = 500)
cat("\n  ref_1e5·1e4 = 기준 locus 중 ±500 kb 안 최소 p가 기준보다 작은 수, sign = lead 효과 방향 일치 비율,",
    "z_ratio = lead z / 주 분석 lead z (중앙값), b_slope = lead b 기울기 (같은 RINT 척도 GCTA 모형만), new = 기준 locus 밖 p < 1e-5 신호 수\n",
    " z_exp = 표본 겹침만으로 기대되는 z 상관 sqrt(n / 7,115): 하위집합은 z_cor가 이 값 근처면 정상(표본이 작아서), 많이 낮을 때만 신호 패턴 차이\n",
    " 이분형(bin_*): z = log OR / SE (REGENIE Firth). 연속형 RINT와 척도가 달라 z_ratio·b_slope는 직접 비교하지 않고, 방향(sign)은\n",
    "  'CSLB가 높아지는 대립유전자 = CCD 양성 위험이 커지는 대립유전자'가 기대되나 50점 기준 이분화로 일부 locus는 방향이 다를 수 있음 (참고 지표)\n",
    " 상호작용(int_age): SNP × age_cslb 항의 p만 봄 (z 상관·방향 없음)\n")
if (nrow(newl)) {
  cat("\n== 새 locus (기준 locus 밖 p < 1e-5) ==\n")
  print(newl[order(p), .(id, lead, Chr, Mb = round(bp / 1e6, 2), p = signif(p, 3), n_snp_p1e5, maf = round(maf, 3))], row.names = FALSE, nrows = 500)
}
sink()
cat(readLines(paste0(out_prefix, "_summary.txt")), sep = "\n")

if (nrow(loc)) {
  hm <- merge(loc, summ[, .(id, group)], by = "id")
  hm[, id := factor(id, levels = summ$id)]
  hm[, label := factor(label, levels = rev(ref$label))]
  hm[, lp := -log10(min_p)]
  g <- ggplot(hm, aes(id, label, fill = lp)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_point(data = hm[min_p < 1e-5], shape = 16, size = 1.2, color = "#FFFFFF") +
    geom_text(data = hm[same_sign == FALSE], label = "x", size = 3, color = "#EB6834") +
    scale_fill_gradientn(colors = seq_ramp, name = expression(min~-log[10](p)), na.value = "#E1E0D9") +
    facet_grid(~group, scales = "free_x", space = "free_x") +
    labs(title = "Reference loci (04-7) across sensitivity analyses",
         subtitle = "cell = smallest p within +/-500 kb of locus lead; white dot = p < 1e-5; orange x = lead effect opposite to main; grey = no result",
         x = NULL, y = NULL) +
    theme_qc + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9), panel.grid = element_blank(),
                     strip.text = element_text(size = 9), panel.spacing = unit(2, "pt"))
  ggsave(paste0(pic_prefix, "_locus_heatmap.png"), g, width = max(10, 0.32 * nrow(summ) + 5), height = 0.32 * nrow(ref) + 3, dpi = 300, limitsize = FALSE)
}
if (nrow(summ[!is.na(lambda)])) {
  mt <- melt(summ[!is.na(lambda), .(id, group, lambda, z_cor_main, ref_1e5 = as.numeric(ref_1e5))], id.vars = c("id", "group"),
             variable.name = "metric")
  mt[, metric := factor(metric, levels = c("lambda", "z_cor_main", "ref_1e5"),
                        labels = c("lambda", "z correlation with main", sprintf("reference loci with p < 1e-5 (of %d)", nrow(ref))))]
  mt[, id := factor(id, levels = rev(summ$id))]
  ref_lines <- data.table(metric = factor(levels(mt$metric), levels = levels(mt$metric)), v = c(1, 1, summ[id == "main", ref_1e5]))
  g <- ggplot(mt[!is.na(value)], aes(value, id)) +
    geom_vline(data = ref_lines, aes(xintercept = v), color = "#898781", linetype = "dashed", linewidth = 0.4) +
    geom_point(color = "#2A78D6", size = 2.4) +
    facet_wrap(~metric, scales = "free_x", nrow = 1) +
    labs(title = "Sensitivity analyses vs main model", subtitle = "dashed = 1 (lambda, correlation) or main model value", x = NULL, y = NULL) +
    theme_qc + theme(axis.text.y = element_text(size = 9), panel.grid.minor = element_blank())
  ggsave(paste0(pic_prefix, "_run_metrics.png"), g, width = 14, height = 0.22 * nrow(summ[!is.na(lambda)]) + 2, dpi = 300, limitsize = FALSE)
}
if (make_run_plots) {
  chr_len <- mn[, .(len = max(bp)), by = Chr][order(Chr)][, off := cumsum(as.numeric(shift(len, fill = 0)))]
  xpos <- mn$bp + chr_len$off[match(mn$Chr, chr_len$Chr)]
  x_bin <- floor(xpos / (max(xpos) / 4000))
  axis_df <- chr_len[, .(Chr, mid = off + len / 2)]
  for (rid in setdiff(ok_ids, c("main", "maf0.05", "maf0.1", grep("^pc[0-9]+$", ok_ids, value = TRUE)))) {
    x <- res[[rid]]; lp <- -log10(x$p); ok <- is.finite(lp)
    d <- data.table(i = which(ok), Chr = mn$Chr[ok], xb = x_bin[ok], yb = floor(lp[ok] / 0.01), lp = lp[ok])
    d <- d[lp >= 3 | !duplicated(d, by = c("Chr", "xb", "yb"))]
    ps <- sort(x$p[ok]); n <- length(ps)
    qi <- sort(unique(c(1:min(20000, n), round(exp(seq(log(min(20001, n)), log(n), length.out = 4000))))))
    qq <- data.table(exp = -log10((qi - 0.5) / n), obs = -log10(ps[qi]))
    band <- data.table(exp = qq$exp, lo = -log10(qbeta(0.975, qi, n - qi + 1)), hi = -log10(qbeta(0.025, qi, n - qi + 1)))
    y_max <- max(-log10(5e-8), max(d$lp), max(qq$obs)) + 0.6
    s <- as.list(summ[which(summ$id == rid)])
    gm <- ggplot(d, aes(xpos[i], lp)) +
      geom_point(aes(color = factor(Chr %% 2)), size = 0.4, show.legend = FALSE) +
      geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "#52514E", linewidth = 0.4) +
      geom_hline(yintercept = -log10(1e-5), linetype = "dotted", color = "#52514E", linewidth = 0.4) +
      scale_color_manual(values = c("0" = col_chr[1], "1" = col_chr[2])) +
      scale_x_continuous(breaks = axis_df$mid, labels = ifelse(axis_df$Chr %% 2 == 1 | axis_df$Chr <= 10, axis_df$Chr, ""), expand = c(0.005, 0)) +
      scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
      labs(title = "Manhattan", subtitle = "dashed p = 5e-8, dotted p = 1e-5; overlapping points merged", x = "Chromosome", y = expression(-log[10](p))) +
      theme_qc + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                       plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
    gq <- ggplot() +
      geom_ribbon(data = band, aes(x = exp, ymin = lo, ymax = hi), fill = "#E1E0D9") +
      geom_abline(slope = 1, intercept = 0, color = "#898781", linewidth = 0.4) +
      geom_point(data = qq, aes(exp, obs), color = "#2A78D6", size = 0.6) +
      annotate("text", x = 0.2, y = y_max - 0.4, hjust = 0, vjust = 1, size = 3.6, label = sprintf("lambda = %.3f", s$lambda)) +
      scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
      labs(title = "QQ", subtitle = "grey band = 95% null interval", x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p))) +
      theme_qc + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
    ttl <- sprintf("%s (%s, N = %s)  |  lambda %.3f  |  SNPs p < 5e-8 / 1e-5 / 1e-4: %d / %d / %d  |  ref loci p < 1e-5: %d of %d  |  new loci: %d",
                   rid, if (nzchar(s$test)) s$test else s$engine, format(s$n, big.mark = ","), s$lambda, s$n_p5e8, s$n_p1e5, s$n_p1e4,
                   s$ref_1e5, nrow(ref), s$n_new_loci)
    gg <- gridExtra::arrangeGrob(gm, gq, ncol = 2, widths = c(3.3, 1), top = grid::textGrob(ttl, gp = grid::gpar(fontsize = 14)))
    ggsave(file.path(run_pic_dir, paste0(rid, ".png")), gg, width = 18, height = 5.5, dpi = 300)
  }
}
