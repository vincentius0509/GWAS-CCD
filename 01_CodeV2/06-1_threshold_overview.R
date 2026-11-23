# 06-1_threshold_overview.R - Counts SNPs below each significance threshold and draws Manhattan/QQ plots for every model. See README.md.
#     p < 5e-8 0 / 4.2e-7 0.70 / 1e-6 1.2 / 1e-5 13.5 / 1e-4 105

library(data.table)
library(ggplot2)
library(parallel)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")    
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
out_prefix <- file.path(qc_dir, "06-1")
pic_prefix <- file.path(pic_dir, "06-1")
run_pic_dir <- paste0(pic_prefix, "_runs")
runs_tsv <- file.path(qc_dir, "05-1_runs.tsv")
gcta_dir <- file.path(qc_dir, "05-2_gwas")
regenie_dir <- file.path(qc_dir, "05-3_regenie")
bim_file <- file.path(qc_dir, "03-3_maf0.01_hwe.bim")
catalog_csv <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023/05-5_model_catalog.csv")
dir.create(run_pic_dir, showWarnings = FALSE, recursive = TRUE)
for (f in c(runs_tsv, bim_file)) if (!file.exists(f)) stop("파일 없음: ", f)
n_cores <- as.integer(Sys.getenv("N_CORES", "24"))

thr <- c("5e-8" = 5e-8, "4.2e-7" = 4.2e-7, "1e-6" = 1e-6, "1e-5" = 1e-5, "1e-4" = 1e-4)
thr_note <- c("5e-8" = "human / DAP", "4.2e-7" = "Pettersson 2026", "1e-6" = "DAP suggestive", "1e-5" = "suggestive", "1e-4" = "reference")
thr_key <- c("5e-8" = "5e8", "4.2e-7" = "4.2e7", "1e-6" = "1e6", "1e-5" = "1e5", "1e-4" = "1e4")
snp_cols <- paste0("snp_", thr_key[names(thr)]); loci_cols <- paste0("loci_", thr_key[names(thr)])
meff <- 0.05 / 4.2e-7
null_loci <- c("5e-8" = 0, "4.2e-7" = 0.70, "1e-6" = 1.2, "1e-5" = 13.5, "1e-4" = 105)
if (!setequal(names(null_loci), names(thr))) stop("thr를 바꾸면 null_loci도 같은 기준으로 맞춰야 함")
locus_kb <- 500
use06 <- c("main", "purelike", "mixlike", "int_age", "age8", "bin_age8",
           "young8", "sex_m", "sex_f", "wt_t1", "wt_t2", "wt_t3")
group_order <- c("main", "qc", "phenotype", "age", "genetic", "strata", "drop_cov", "add_cov", "subset")

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), axis.title = element_text(size = 13), axis.text = element_text(size = 11))
col_chr <- c("#2A78D6", "#898781")
thr_lty <- c("5e-8" = "solid", "4.2e-7" = "dashed", "1e-6" = "longdash", "1e-5" = "dotted", "1e-4" = "dotdash")
p_to_lambda <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
count_loci <- function(chr, bp, p, t) {
  h <- which(p < t); if (!length(h)) return(0L)
  o <- h[order(chr[h], bp[h])]
  sum(c(TRUE, diff(chr[o]) != 0 | diff(bp[o]) > locus_kb * 1e3))
}

runs <- fread(runs_tsv, na.strings = c("", "NA"))[engine != "check"]
if (!all(use06 %in% runs$id)) message("실행 표에 없는 06 모형: ", paste(setdiff(use06, runs$id), collapse = ", "), " (05-1 다시 실행 필요)")
if (nzchar(Sys.getenv("RUN_IDS"))) runs <- runs[id %in% strsplit(Sys.getenv("RUN_IDS"), " ")[[1]]]
runs[, use_06 := id %in% use06]
runs[, group := factor(group, levels = group_order)]
runs[, ord := fifelse(use_06, match(id, use06), 100L + as.integer(group))]
setorder(runs, ord, group)
if (file.exists(catalog_csv)) {
  ct <- fread(catalog_csv, select = c("id", "role", "gate"))
  runs <- merge(runs, ct, by = "id", all.x = TRUE, sort = FALSE)
} else runs[, `:=`(role = NA_character_, gate = NA)]
setorder(runs, ord, group)

bim <- fread(bim_file, header = FALSE, select = c(1, 2, 4), col.names = c("Chr", "SNP", "bp"))
chr_len <- bim[, .(len = max(bp)), by = Chr][order(Chr)][, off := cumsum(as.numeric(shift(len, fill = 0)))]
axis_df <- chr_len[, .(Chr, mid = off + len / 2)]
x_total <- max(chr_len$off + chr_len$len)
main_file <- runs[id == "main", ref_result]
cat(sprintf("모형 %d개 (06 사용 %d, 미사용 %d), 동시 %d개\n", nrow(runs), sum(runs$use_06), sum(!runs$use_06), n_cores))

runs[runs$use_06==TRUE][,1]

read_res <- function(r) {
  if (r$engine == "existing") {
    if (!file.exists(r$ref_result)) return(NULL)
    d <- fread(cmd = paste("zcat", shQuote(r$ref_result)), select = c("Chr", "SNP", "bp", "p")); d[, test := ""]
  } else if (r$engine == "gcta") {
    f <- file.path(gcta_dir, paste0(r$id, ".mlma.gz")); if (!file.exists(f)) return(NULL)
    d <- fread(cmd = paste("zcat", shQuote(f)), select = c("Chr", "SNP", "bp", "p")); d[, test := ""]
  } else if (r$engine %in% c("regenie_bt", "regenie_int")) {
    f <- file.path(regenie_dir, paste0(r$id, ".regenie.gz")); if (!file.exists(f)) return(NULL)
    d <- fread(cmd = paste("zcat", shQuote(f)), select = c("CHROM", "GENPOS", "ID", "TEST", "LOG10P"))
    tsel <- grep(if (r$engine == "regenie_bt") "^ADD$" else "^ADD-INT_SNPx", unique(d$TEST), value = TRUE)[1]
    if (is.na(tsel)) stop(r$id, ": REGENIE TEST 없음")
    d <- d[TEST == tsel & !is.na(LOG10P), .(Chr = CHROM, SNP = ID, bp = GENPOS, p = 10^(-LOG10P), test = tsel)]
  } else if (r$engine == "maf_filter") {
    if (!file.exists(r$ref_result) || !file.exists(main_file)) return(NULL)
    keep <- fread(r$ref_result, header = FALSE, select = 2)$V2
    d <- fread(cmd = paste("zcat", shQuote(main_file)), select = c("Chr", "SNP", "bp", "p"))[SNP %in% keep]; d[, test := ""]
  } else return(NULL)
  d[is.finite(p)]
}

plot_run <- function(r, d, s) {
  d[, xpos := bp + chr_len$off[match(Chr, chr_len$Chr)]]
  d[, lp := -log10(p)]
  pd <- d[, .(Chr, xpos, lp, xb = floor(xpos / (x_total / 4000)), yb = floor(lp / 0.01))]
  pd <- pd[lp >= 3 | !duplicated(pd, by = c("Chr", "xb", "yb"))]
  ps <- sort(d$p); n <- length(ps)
  qi <- sort(unique(c(1:min(20000, n), round(exp(seq(log(min(20001, n)), log(n), length.out = 4000))))))
  qq <- data.table(exp = -log10((qi - 0.5) / n), obs = -log10(ps[qi]))
  band <- data.table(exp = qq$exp, lo = -log10(qbeta(0.975, qi, n - qi + 1)), hi = -log10(qbeta(0.025, qi, n - qi + 1)))
  y_max <- max(-log10(5e-8), max(pd$lp), max(qq$obs)) + 0.6
  hl <- data.table(name = names(thr), y = -log10(thr), lab = sprintf("%s (%s)", names(thr), thr_note))
  gm <- ggplot(pd, aes(xpos, lp)) +
    geom_point(aes(color = factor(Chr %% 2)), size = 0.4, show.legend = FALSE) +
    geom_hline(data = hl, aes(yintercept = y, linetype = name), color = "#52514E", linewidth = 0.4, show.legend = FALSE) +
    geom_label(data = hl, aes(x = x_total, y = y, label = lab), inherit.aes = FALSE, hjust = 1, vjust = 0.5, size = 3.1, color = "#52514E",
               fill = "white", label.size = 0, label.padding = unit(0.12, "lines")) +
    scale_linetype_manual(values = thr_lty) +
    scale_color_manual(values = c("0" = col_chr[1], "1" = col_chr[2])) +
    scale_x_continuous(breaks = axis_df$mid, labels = ifelse(axis_df$Chr %% 2 == 1 | axis_df$Chr <= 10, axis_df$Chr, ""), expand = c(0.005, 0)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(title = "Manhattan", subtitle = "lines: p thresholds (label = source); overlapping points merged", x = "Chromosome", y = expression(-log[10](p))) +
    theme_qc + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                     plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
  gq <- ggplot() +
    geom_ribbon(data = band, aes(x = exp, ymin = lo, ymax = hi), fill = "#E1E0D9") +
    geom_abline(slope = 1, intercept = 0, color = "#898781", linewidth = 0.4) +
    geom_hline(data = hl, aes(yintercept = y, linetype = name), color = "#52514E", linewidth = 0.3, show.legend = FALSE) +
    scale_linetype_manual(values = thr_lty) +
    geom_point(data = qq, aes(exp, obs), color = "#2A78D6", size = 0.6) +
    annotate("text", x = max(qq$exp), y = 0.3, hjust = 1, vjust = 0, size = 3.6, label = sprintf("lambda = %.3f", s$lambda)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(title = "QQ", subtitle = "grey band = 95% null interval", x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p))) +
    theme_qc + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
  cnt_txt <- paste(sprintf("%s: %d (%d)", names(thr), unlist(s[snp_cols]), unlist(s[loci_cols])), collapse = " / ")
  ttl <- sprintf("%s [%s]  %s, N = %s  |  lambda %.3f  |  SNPs (loci) p < %s",
                 r$id, if (r$use_06) "used in 06" else "not used in 06", if (nzchar(s$test)) s$test else if (r$engine == "maf_filter") "GCTA MLMA-LOCO, SNP set filter" else "GCTA MLMA-LOCO", format(r$n, big.mark = ","), s$lambda,
                 cnt_txt)
  gg <- gridExtra::arrangeGrob(gm, gq, ncol = 2, widths = c(3.3, 1), top = grid::textGrob(ttl, gp = grid::gpar(fontsize = 13)))
  ggsave(file.path(run_pic_dir, paste0(r$id, ".png")), gg, width = 18, height = 5.5, dpi = 300)
}

one_run <- function(k) {
  r <- as.list(runs[k])
  out <- tryCatch({
    d <- read_res(r)
    if (is.null(d)) return(data.table(id = r$id, status = "결과 없음"))
    s <- list(id = r$id, status = "있음", test = d$test[1], n_snp = nrow(d), lambda = p_to_lambda(d$p), min_p = min(d$p))
    for (nm in names(thr)) {
      s[[paste0("snp_", thr_key[nm])]] <- sum(d$p < thr[nm])
      s[[paste0("loci_", thr_key[nm])]] <- count_loci(d$Chr, d$bp, d$p, thr[nm])
    }
    plot_run(r, d, s)
    as.data.table(s)
  }, error = function(e) data.table(id = r$id, status = paste("오류:", conditionMessage(e))))
  out
}

res <- mclapply(seq_len(nrow(runs)), one_run, mc.cores = n_cores, mc.preschedule = FALSE)
res <- lapply(seq_along(res), function(i) if (inherits(res[[i]], "try-error") || is.null(res[[i]])) data.table(id = runs$id[i], status = "오류 (프로세스 종료)") else res[[i]])
summ <- rbindlist(res, fill = TRUE)
summ <- merge(runs[, .(id, use_06, role, gate, group, engine, n, ord, desc)], summ, by = "id", all.x = TRUE, sort = FALSE)
setorder(summ, ord, group)
cat(summ[status != "있음", sprintf("%s: %s\n", id, status)], sep = "")

exp_loci <- null_loci[names(thr)]
cnt_cols <- c(rbind(snp_cols, loci_cols))
out <- summ[, c("id", "use_06", "role", "gate", "group", "engine", "test", "n", "status", "n_snp", "lambda", "min_p", cnt_cols, "desc"), with = FALSE]
out[, `:=`(lambda = round(lambda, 4), min_p = signif(min_p, 3))]
fwrite(out, paste0(out_prefix, "_threshold_counts.csv"), bom = TRUE)

sink(paste0(out_prefix, "_threshold_counts.txt"))
cat("06-1 유의 기준별 SNP·locus 수 (locus = 이웃 간격 > 500 kb이면 나눈 거리 묶음)\n")
cat(sprintf("기준: %s\n", paste(sprintf("%s (%s)", names(thr), thr_note), collapse = " / ")))
cat(sprintf("우연 기대 locus 수 (실측 null: 주 분석 표본 표현형 순열 GWAS 20회 평균): %s\n",
            paste(sprintf("%s → %g", names(thr), exp_loci), collapse = " / ")))
cat(sprintf("참고: 유효 검정 수 %s 가정 기대값 (0.05 / 4.2e-7에서 역산): %s\n", format(round(meff), big.mark = ","),
            paste(sprintf("%s → %.2g", names(thr), meff * thr), collapse = " / ")))
cat("표 칸: SNP 수 (locus 수)\n")
fmt <- function(x) {
  y <- x[, .(id, role, gate, n, lambda, min_p)]
  for (i in seq_along(thr)) y[[paste0("p<", names(thr)[i])]] <- sprintf("%d (%d)", x[[snp_cols[i]]], x[[loci_cols[i]]])
  y
}
cat("\n== 06에 가져가는 모형 ==\n"); print(fmt(out[use_06 == TRUE & status == "있음"]), row.names = FALSE)
cat("\n== 06에 가져가지 않는 모형 ==\n"); print(fmt(out[use_06 == FALSE & status == "있음"]), row.names = FALSE)
if (nrow(out[status != "있음"])) { cat("\n== 결과 없음 / 오류 ==\n"); print(out[status != "있음", .(id, use_06, status)], row.names = FALSE) }
sink()
cat(readLines(paste0(out_prefix, "_threshold_counts.txt")), sep = "\n")

ov <- melt(out[status == "있음", c("id", "use_06", loci_cols), with = FALSE][, use := fifelse(use_06, "used in 06", "not used in 06")][, use_06 := NULL],
           id.vars = c("id", "use"), measure.vars = loci_cols, variable.name = "thr", value.name = "loci")
thr_lab <- setNames(paste("p <", names(thr)), loci_cols)
ov[, thr := factor(thr_lab[as.character(thr)], levels = thr_lab)]
ov[, id := factor(id, levels = rev(out[status == "있음", id]))]
ov[, use := factor(use, levels = c("used in 06", "not used in 06"), labels = c("used\nin 06", "not used\nin 06"))]
ov[, fill_v := if (max(loci) > 0) loci / max(loci) else 0, by = thr]
g <- ggplot(ov, aes(thr, id)) +
  geom_tile(aes(fill = fill_v), color = "white", linewidth = 0.6) +
  geom_text(aes(label = loci), size = 3.2, color = ifelse(ov$fill_v > 0.6, "white", "#1F1E1C")) +
  scale_fill_gradient(low = "#F4F3EF", high = "#1C5CAB", limits = c(0, 1), guide = "none") +
  facet_grid(use ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Independent loci (500 kb distance clumps) by p threshold",
       subtitle = sprintf("number = loci; shading scaled within each threshold; permutation-null mean loci ~ %s",
                          paste(sprintf("%g", exp_loci), collapse = " / ")), x = NULL, y = NULL) +
  theme_qc + theme(panel.grid = element_blank(), axis.text.y = element_text(size = 9),
                   plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"), strip.text.y = element_text(size = 10, angle = 0))
ggsave(paste0(pic_prefix, "_threshold_overview.png"), g, width = 8, height = 0.24 * nrow(out[status == "있음"]) + 2.2, dpi = 300, limitsize = FALSE)
cat(sprintf("저장: %s_threshold_counts.csv / .txt, %s_threshold_overview.png, %s/{id}.png\n", out_prefix, pic_prefix, run_pic_dir))
