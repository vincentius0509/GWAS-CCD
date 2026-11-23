# 13-2_regional_ld.R - Regional plots across models and LD comparison between purebred-like and mixed-like dogs. See README.md.
#   03-3_maf0.01_hwe (bfile) · 04-1_grm/grm_all.grm.id (pooled) · 05-1_inputs/purelike_keep.txt · mixlike_keep.txt · 03_GWAS_Pure/06-1_gene_loc.txt

suppressMessages({ library(data.table); library(ggplot2); library(parallel) })
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
pic_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/05_GWAS_2023"
test_mode <- Sys.getenv("TEST", "0") == "1"; tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("13-2", tag)); fig_dir <- file.path(pic_dir, paste0("13-2", tag, "_regional"))
in13 <- file.path(qc_dir, "13-1")
n_cores <- as.integer(Sys.getenv("N_CORES", "8"))
plink1 <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
bfile <- file.path(d0, "03-3_maf0.01_hwe")
keeps <- c(pooled = file.path(d0, "04-1_grm", "grm_all.grm.id"), purelike = file.path(d0, "05-1_inputs", "purelike_keep.txt"),
           mixlike = file.path(d0, "05-1_inputs", "mixlike_keep.txt"))
strata3 <- names(keeps)
ld_grp <- function(m) fifelse(m %in% c("main", "int_age", "bin_all"), "pooled", fifelse(m == "bin_age8", "age8", m))
plink_thr <- 8
gl <- fread("/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/06-1_gene_loc.txt", header = FALSE,
            col.names = c("GENE", "chr", "start", "end", "strand", "sym"))
flank <- 2.5e5; ld_kb <- 1000
r2_strong <- 0.8; r2_hit <- 0.2
ld_hi <- 0.5; ld_lo <- 0.1; cor_flag <- 0.8
r2_cols <- c("#D0D0D0", "#46ACC8", "#E17C05", "#B2182B")
fam_order <- c("main", "wt_t1", "wt_t2", "wt_t3", "sex_m", "sex_f", "purelike", "mixlike", "young8", "age8", "int_age",
               "agefloor6", "agefloor7", "agefloor9", "age10", "agefloor11", "agefloor12", "bin_all", "bin_age8")
age_series <- c("young8", "main", "agefloor6", "agefloor7", "age8", "agefloor9", "age10", "agefloor11", "agefloor12")
wt_series <- c("wt_t1", "wt_t2", "wt_t3")
series_order <- c(age_series, wt_series, "sex_m", "sex_f", "purelike", "mixlike", "int_age", "bin_all", "bin_age8")
series_of <- function(m) fifelse(m %in% age_series, "나이", fifelse(m %in% wt_series, "체중", "기타"))
gfile <- function(id) {
  if (id == "main") return("04-5_gwas/rint_maf0.01_PC5.mlma.gz")
  if (id %in% c("int_age", "bin_age8", "bin_all")) return(sprintf("05-3_regenie/%s.regenie.gz", id))
  sprintf("05-2_gwas/%s.mlma.gz", id)
}
for (f in c(paste0(in13, c("_priority.csv", "_compare.csv", "_region_leads.csv")), paste0(bfile, ".bed"), keeps, plink1, file.path(d0, sapply(fam_order, gfile))))
  if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
wd <- tempfile("r13_2_"); dir.create(wd)

pri <- fread(paste0(in13, "_priority.csv")); cmp <- fread(paste0(in13, "_compare.csv"))
regs <- unique(cmp$region_id); regs94 <- unique(pri[tier <= 2, region_id])
if (test_mode) regs <- intersect(c("chr4_r8", "chr7_r8", "chr5_r25", "chr2_r33"), unique(cmp$region_id))[1:4]
regs <- regs[!is.na(regs)]
ld_leads <- unique(cmp[region_id %in% regs, .(region_id, model, SNP, chr, bp, p)], by = c("region_id", "model", "SNP"))
R <- ld_leads[, {
  mi <- which(model == "main"); k <- if (length(mi)) mi[which.min(p[mi])] else which.min(p)
  .(chr = chr[1], index = SNP[k], index_model = model[k], index_bp = bp[k], index_p = p[k], lo = max(1, min(bp) - flank), hi = max(bp) + flank,
    panels = paste(intersect(fam_order, c("main", unique(model))), collapse = ","))
}, by = region_id]
R[, ord := match(region_id, regs)]; setorder(R, chr, index_bp); R[, in94 := region_id %in% regs94]
R[cmp[, .(region_tier = min(tier, na.rm = TRUE), has_1e5 = any(thr == 1e-5)), by = region_id], on = "region_id", `:=`(region_tier = i.region_tier, has_1e5 = i.has_1e5)]
cat(sprintf("13-2 대상 위치 %d개 (94곳 %d · 1e-5 locus 있는 곳 %d) · 패널 모형 %d개\n", nrow(R), sum(R$in94), sum(R$has_1e5),
            length(unique(unlist(strsplit(R$panels, ","))))))
sub_grp <- setdiff(intersect(fam_order, ld_grp(unique(unlist(strsplit(R$panels, ","))))), names(keeps))
keeps <- c(keeps, setNames(file.path(d0, "05-1_inputs", paste0(sub_grp, "_keep.txt")), sub_grp))
for (f in keeps) if (!file.exists(f)) stop("파일 없음: ", f)
n_dogs <- sapply(keeps, function(f) length(readLines(f)))

win_f <- file.path(wd, "win.txt"); fwrite(R[, .(chr, lo, hi)], win_f, sep = "\t", col.names = FALSE)
models <- intersect(fam_order, unique(unlist(strsplit(R$panels, ","))))
read_win <- function(id) {
  f <- file.path(d0, gfile(id)); rg <- grepl("regenie", f); test <- if (id == "int_age") "ADD-INT_SNPxage_cslb" else "ADD"
  aw <- sprintf("awk 'NR==FNR {n[$1]++; s[$1,n[$1]]=$2; e[$1,n[$1]]=$3; next} FNR==1 {print; next} %s {c=$1; p=$%d; for (k=1; k<=n[c]; k++) if (p>=s[c,k] && p<=e[c,k]) {print; break}}' %s -",
                if (rg) sprintf("$8==\"%s\"", test) else "", if (rg) 2 else 3, win_f)
  x <- fread(cmd = sprintf("gzip -dc %s | %s", f, aw))
  if (rg) x[, .(Chr = CHROM, SNP = ID, bp = GENPOS, p = 10^(-LOG10P))] else x[, .(Chr, SNP, bp, p)]
}
cat(sprintf("[%s] GWAS %d개 읽는 중 (코어 %d)\n", format(Sys.time(), "%T"), length(models), n_cores))
G <- mclapply(setNames(models, models), read_win, mc.cores = min(n_cores, length(models)))
bad <- names(G)[sapply(G, function(x) inherits(x, "try-error") || !is.data.table(x) || !nrow(x))]
if (length(bad)) stop("GWAS 읽기 실패 또는 범위 안 SNP 0개: ", paste(bad, collapse = ", "))

writeLines(unique(R$index), file.path(wd, "idx.txt"))
lp <- fread(paste0(in13, "_region_leads.csv"))[region_id %in% regs]
writeLines(unique(c(lp$lead_a, lp$lead_b, R$index)), file.path(wd, "leads.txt"))
run_plink <- function(args, out) {
  system2(plink1, c("--bfile", bfile, "--dog", "--threads", plink_thr, args, "--out", out), stdout = FALSE, stderr = FALSE)
  if (!file.exists(paste0(out, ".ld"))) stop("plink 실패: ", out, ".log")
  fread(paste0(out, ".ld"))
}
R[, multi_signal := region_id %in% lp[ld_signal == "다른 신호 후보", region_id]]
n_job <- min(n_cores, length(keeps)); plink_thr <- max(1L, n_cores %/% n_job)
cat(sprintf("[%s] LD 표본 %d개 (동시 %d개 · 각 %d스레드): %s\n", format(Sys.time(), "%T"), length(keeps), n_job, plink_thr,
            paste(sprintf("%s %d", names(keeps), n_dogs), collapse = " · ")))
loc <- cmp[region_id %in% regs, .(model, SNP, thr, region_id, chr, bp, p, tier)]; loc[, grp := ld_grp(model)]   # 13-1 locus (⑥)
LDres <- mclapply(setNames(names(keeps), names(keeps)), function(g) {
  sl <- unique(c(R$index, loc[grp == g, SNP])); writeLines(sl, file.path(wd, paste0("idx_", g, ".txt")))
  x <- run_plink(c("--keep", keeps[[g]], "--r2", "--ld-snp-list", file.path(wd, paste0("idx_", g, ".txt")), "--ld-window-kb", ld_kb, "--ld-window", "999999",
                   "--ld-window-r2", "0"), file.path(wd, paste0("idx_", g)))
  ld <- unique(rbind(x[, .(index = SNP_A, SNP = SNP_B, bp = BP_B, r2 = R2)], data.table(index = sl, SNP = sl, bp = NA_integer_, r2 = 1)),
               by = c("index", "SNP"))
  setkey(ld, index)
  pp <- NULL
  if (g %in% strata3) {
    y <- run_plink(c("--keep", keeps[[g]], "--extract", file.path(wd, "leads.txt"), "--r2", "--ld-window-kb", "10000", "--ld-window", "999999", "--ld-window-r2", "0"),
                   file.path(wd, paste0("leads_", g)))
    pp <- rbind(y[, .(a = SNP_A, b = SNP_B, r2 = R2)], y[, .(a = SNP_B, b = SNP_A, r2 = R2)])
  }
  list(ld = ld, lp = pp)
}, mc.cores = n_job)
bad <- names(LDres)[sapply(LDres, function(x) inherits(x, "try-error") || !is.list(x))]; if (length(bad)) stop("LD 계산 실패: ", paste(bad, collapse = ", "))
LD <- lapply(LDres, `[[`, "ld"); LP <- lapply(LDres[strata3], `[[`, "lp"); rm(LDres)

LS <- rbindlist(lapply(seq_len(nrow(R)), function(j) {
  ix <- R$index[j]; base <- LD$pooled[.(ix), nomatch = NULL][SNP != ix, .(SNP, bp, r2p = r2)]; pg <- ld_grp(strsplit(R$panels[j], ",")[[1]])
  rbindlist(lapply(names(keeps), function(g) {
    x <- LD[[g]][.(ix), nomatch = NULL][SNP != ix]; m <- merge(base, x[, .(SNP, r2)], by = "SNP")
    dd <- sort(abs(x[r2 >= r2_strong & is.finite(bp), bp] - R$index_bp[j]), decreasing = TRUE) / 1e3
    data.table(region_id = R$region_id[j], index = ix, group = g, n_dogs = n_dogs[[g]], in_panels = g %in% pg,
               n_window = nrow(x), n_r2_ge_0.8 = sum(x$r2 >= r2_strong, na.rm = TRUE),
               n_r2_ge_0.2 = sum(x$r2 >= r2_hit, na.rm = TRUE), mean_r2 = mean(x$r2, na.rm = TRUE), max_r2 = if (nrow(x)) max(x$r2, na.rm = TRUE) else NA_real_,
               strong_far_kb = if (length(dd)) dd[1] else 0, strong_far2_kb = if (length(dd) >= 2) dd[2] else 0,
               cor_with_pooled = if (g == "pooled") 1 else { mh <- m[pmax(r2p, r2, na.rm = TRUE) >= r2_hit]; if (nrow(mh) >= 10) cor(mh$r2p, mh$r2, use = "complete.obs") else NA_real_ })
  }))
}))
LS[, ld_diff_flag := group != "pooled" & !is.na(cor_with_pooled) & cor_with_pooled < cor_flag]
LS[, in94 := region_id %in% regs94]

LL <- rbindlist(lapply(seq_len(nrow(loc)), function(k) {
  l <- loc[k]; x <- LD[[l$grp]][.(l$SNP), nomatch = NULL][SNP != l$SNP]
  st <- x[r2 >= r2_strong]; dd <- sort(abs(st$bp - l$bp), decreasing = TRUE) / 1e3
  pa <- G[[l$model]][SNP %in% st$SNP, p]
  data.table(n_window = nrow(x), n_r2_ge_0.8 = nrow(st), n_r2_ge_0.2 = sum(x$r2 >= r2_hit), strong_far_kb = if (length(dd)) dd[1] else 0,
             strong_far2_kb = if (length(dd) >= 2) dd[2] else 0, n_strong_in_gwas = length(pa), n_strong_p1e3 = sum(pa < 1e-3))
}))
LL <- cbind(loc, LL); LL[, frac_strong_p1e3 := fifelse(n_strong_in_gwas > 0, n_strong_p1e3 / n_strong_in_gwas, NA_real_)]
LL[, in94 := region_id %in% regs94]; LL[, model := factor(model, fam_order)]; setorder(LL, model, -thr, chr, bp); LL[, model := as.character(model)]
LM <- LL[, .(n_loci = .N, med_n_r2_ge_0.8 = as.numeric(median(n_r2_ge_0.8)), pct_no_strong = 100 * mean(n_r2_ge_0.8 == 0),
             med_strong_far2_kb = as.numeric(median(strong_far2_kb)), med_frac_strong_p1e3 = median(frac_strong_p1e3, na.rm = TRUE),
             pct_ld_support = 100 * mean(n_strong_p1e3 >= 1)), by = .(model, thr)]
LM[, thr := factor(format(thr, scientific = TRUE), c("1e-05", "1e-04"))]; LM[, series := series_of(model)]
LM[, model := factor(model, series_order)]; setorder(LM, thr, model)
cls <- function(r) fifelse(is.na(r), "LD 없음", fifelse(r >= ld_hi, "같은 신호", fifelse(r < ld_lo, "다른 신호 후보", "불분명")))
for (g in strata3) {
  lp[LP[[g]], on = .(lead_a = a, lead_b = b), (paste0("r2_", g)) := i.r2]
  lp[lead_a == lead_b, (paste0("r2_", g)) := 1]
}
lp[, `:=`(class_pooled = cls(r2_pooled), class_purelike = cls(r2_purelike), class_mixlike = cls(r2_mixlike))]
lp[, ld_sensitive := !(class_pooled == class_purelike & class_pooled == class_mixlike)]
lp[, r2_check := abs(r2_pooled - r2) < 1e-3]

gene_rows <- function(cc, lo, hi) {
  g <- gl[chr == cc & end >= lo & start <= hi][order(start)]; if (!nrow(g)) return(g)
  pad <- (hi - lo) * 0.06; row_end <- numeric(0); r <- integer(nrow(g))
  for (k in seq_len(nrow(g))) { j <- which(row_end < g$start[k] - pad)
    if (length(j)) { r[k] <- j[1]; row_end[j[1]] <- g$end[k] + pad } else { row_end <- c(row_end, g$end[k] + pad); r[k] <- length(row_end) } }
  g[, row := r][row <= 6]
}
plot_region <- function(j) {
  x <- R[j]; pans <- strsplit(x$panels, ",")[[1]]; cc <- x$chr
  pts <- rbindlist(lapply(pans, function(m) {
    g <- G[[m]][Chr == cc & bp >= x$lo & bp <= x$hi]; grp <- ld_grp(m)
    if (!nrow(g)) stop("패널 SNP 0개: ", m)
    g <- merge(g, LD[[grp]][.(x$index), nomatch = NULL][, .(SNP, r2)], by = "SNP", all.x = TRUE)
    g[, `:=`(panel = m, logp = -log10(p))]
  }))
  lab <- ld_leads[region_id == x$region_id][pts, on = .(model = panel, SNP), nomatch = NULL][, .(panel = model, SNP, bp = i.bp, logp, r2)]
  strip <- pts[, .(minp = min(p)), by = panel][, lbl := sprintf("%s\nmin p %.1e", panel, minp)]
  lev <- c(strip$lbl[match(pans, strip$panel)], "genes")
  pts[strip, on = "panel", fac := i.lbl]; lab[strip, on = "panel", fac := i.lbl]
  gr <- gene_rows(cc, x$lo, x$hi); if (nrow(gr)) gr[, fac := "genes"]
  pts[, fac := factor(fac, lev)]; lab[, fac := factor(fac, lev)]; if (nrow(gr)) gr[, fac := factor(fac, lev)]
  thr_l <- data.table(fac = factor(setdiff(lev, "genes"), lev))
  p <- ggplot() +
    geom_hline(data = thr_l, aes(yintercept = 4), linetype = "dotted", colour = "grey55", linewidth = 0.3) +
    geom_hline(data = thr_l, aes(yintercept = 5), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = x$index_bp, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_point(data = pts[is.na(r2)], aes(bp, logp), colour = "grey85", size = 0.5) +
    geom_point(data = pts[!is.na(r2)][order(r2)], aes(bp, logp, colour = r2), size = 0.9) +
    geom_point(data = pts[SNP == x$index], aes(bp, logp), shape = 23, size = 2.6, fill = "red", colour = "black") +
    geom_point(data = lab, aes(bp, logp, fill = r2), shape = 24, size = 2.4, colour = "black") +
    geom_label(data = lab, aes(bp, logp, label = sprintf("%s  r2=%.2f", SNP, r2)), size = 2.2, vjust = -0.5, label.size = 0.15,
               label.padding = unit(0.08, "lines"), fill = "white", alpha = 0.85) +
    scale_colour_gradientn(colours = r2_cols, limits = c(0, 1), name = expression(r^2 ~ "with index")) +
    scale_fill_gradientn(colours = r2_cols, limits = c(0, 1), guide = "none") +
    facet_grid(fac ~ ., scales = "free_y", switch = "y") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
    scale_x_continuous(labels = function(b) sprintf("%.2f", b / 1e6), limits = c(x$lo, x$hi), expand = expansion(mult = 0.01)) +
    labs(title = sprintf("%s  chr%d:%.2f Mb  index %s (%s, p=%.1e)%s", x$region_id, cc, x$index_bp / 1e6, x$index, x$index_model, x$index_p,
                         if (x$multi_signal) "\n[13-1: distinct-signal candidate pair here -> peaks not in LD with index look grey; conditional analysis in 13-3]" else ""),
         subtitle = "colour = r2 with index, LD in each model's own GWAS sample (main / int_age / bin_all: pooled 7,115 dogs); triangle = each model's lead; dotted p=1e-4, dashed p=1e-5",
         x = sprintf("chr%d position (Mb)", cc), y = expression(-log[10](p))) +
    theme_bw(base_size = 9) + theme(strip.text.y.left = element_text(angle = 0, size = 7), strip.placement = "outside", legend.position = "right",
                                    panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7), plot.title = element_text(size = 10))
  if (nrow(gr)) p <- p + geom_segment(data = gr, aes(x = pmax(start, x$lo), xend = pmin(end, x$hi), y = -row, yend = -row), linewidth = 1.2, colour = "navy") +
    geom_text(data = gr, aes(x = pmin(pmax((pmax(start, x$lo) + pmin(end, x$hi)) / 2, x$lo + 0.04 * (x$hi - x$lo)), x$hi - 0.04 * (x$hi - x$lo)), y = -row + 0.45, label = sym),
              size = 1.9, colour = "navy")
  ggsave(file.path(fig_dir, paste0(x$region_id, ".png")), p, width = 10, height = 1.3 * (length(pans) + 1) + 1.2, dpi = 200, limitsize = FALSE)
  data.table(region_id = x$region_id, n_snp_plotted = nrow(pts), n_r2 = sum(!is.na(pts$r2)))
}
cat(sprintf("[%s] 그림 %d장\n", format(Sys.time(), "%T"), nrow(R)))
PS <- rbindlist(mclapply(seq_len(nrow(R)), function(j) tryCatch(plot_region(j), error = function(e) data.table(region_id = R$region_id[j], error = conditionMessage(e))),
                         mc.cores = min(n_cores, nrow(R))), fill = TRUE)
if ("error" %in% names(PS) && any(!is.na(PS$error))) { print(PS[!is.na(error)]); stop("그림 실패") }

LSp <- merge(LS[group %in% strata3 & in94 == TRUE], R[, .(region_id, chr, index_bp)], by = "region_id")
LSp[, region := factor(region_id, unique(LSp[order(chr, index_bp), region_id]))]
LSl <- melt(LSp, id.vars = c("region", "group"), measure.vars = c("mean_r2", "n_r2_ge_0.8", "strong_far2_kb"), variable.factor = FALSE)
LSl[, variable := factor(variable, c("mean_r2", "n_r2_ge_0.8", "strong_far2_kb"),
                         c("mean r2 with index (1 Mb window)", "SNPs with r2 >= 0.8", "spread of r2 >= 0.8 SNPs (2nd farthest, kb from index)"))]
g2 <- ggplot(LSl, aes(region, value, colour = group)) + geom_point(size = 1.2, alpha = 0.85) + facet_wrap(~variable, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(pooled = "grey30", purelike = "#B2182B", mixlike = "#2166AC")) + labs(x = NULL, y = NULL, colour = "LD sample") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), panel.grid.minor = element_blank())
ggsave(file.path(pic_dir, paste0("13-2", tag, "_ld_strata.png")), g2, width = max(8, 0.12 * uniqueN(LSp$region) + 3), height = 8, dpi = 200, limitsize = FALSE)
LSs <- LS[group != "pooled" & !is.na(cor_with_pooled)]
lab_lev <- sprintf("%s\n(n=%d)", names(keeps), n_dogs)
LSs[, lab := factor(sprintf("%s\n(n=%d)", group, n_dogs), lab_lev)]
g3 <- ggplot(LSs, aes(lab, cor_with_pooled)) + geom_hline(yintercept = cor_flag, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.6, colour = "grey40") +
  geom_jitter(aes(colour = in_panels), width = 0.15, height = 0, size = 0.8, alpha = 0.8) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey65"), labels = c(`TRUE` = "model has a lead here", `FALSE` = "no lead here"), name = NULL) +
  labs(x = "LD sample (GWAS sample of each model)", y = "r2 correlation with pooled LD\n(SNPs with r2 >= 0.2 in either, 1 Mb around index)",
       title = sprintf("LD around index SNPs of %d regions: each model's sample vs pooled (7,115 dogs)", nrow(R)), subtitle = sprintf("dashed = %.1f (flag 'LD differs', descriptive)", cor_flag)) +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(pic_dir, paste0("13-2", tag, "_ld_samples.png")), g3, width = max(8, 0.65 * length(keeps) + 2), height = 5, dpi = 200)
LLp <- melt(LL[, .(model = factor(model, series_order), thr = factor(format(thr, scientific = TRUE), c("1e-05", "1e-04")), n_r2_ge_0.8, strong_far2_kb, frac_strong_p1e3)],
            id.vars = c("model", "thr"), variable.factor = FALSE)
LLp[, variable := factor(variable, c("n_r2_ge_0.8", "strong_far2_kb", "frac_strong_p1e3"),
                         c("SNPs with r2 >= 0.8 to lead", "spread of r2 >= 0.8 SNPs (2nd farthest, kb)", "LD support: share of r2 >= 0.8 SNPs with GWAS p < 1e-3"))]
g4 <- ggplot(LLp[!is.na(value)], aes(model, value, fill = thr)) + geom_boxplot(outlier.size = 0.4, width = 0.7, position = position_dodge(0.8)) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") + scale_y_continuous(trans = "sqrt") +
  scale_fill_manual(values = c(`1e-05` = "#B2182B", `1e-04` = "#46ACC8"), name = "locus threshold") +
  geom_vline(xintercept = c(length(age_series), length(age_series) + length(wt_series)) + 0.5, linetype = "dashed", colour = "grey50") +
  labs(x = "age series (young8 < 8 y, main = all, then minimum age 6 ~ 12 y) | weight tertiles | others", y = "value (sqrt scale)",
       title = "LD around each 13-1 locus lead, in that model's GWAS sample (descriptive)") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(pic_dir, paste0("13-2", tag, "_locus_ld.png")), g4, width = 11, height = 9, dpi = 200)

R <- merge(R, PS, by = "region_id"); setorder(R, chr, index_bp); R[, ord := NULL]
fwrite(R, paste0(out_prefix, "_regions.csv"))
f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
fwrite(f3(LS), paste0(out_prefix, "_ld_strata.csv"))
fwrite(f3(lp), paste0(out_prefix, "_lead_ld_strata.csv"))
fwrite(f3(LL), paste0(out_prefix, "_locus_ld.csv")); fwrite(f3(LM), paste0(out_prefix, "_locus_ld_by_model.csv"))
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("13-2 여러 모형 regional plot · 층별 LD (서술; 검정 · 우연 기대 없음)\n"); if (test_mode) cat("*** TEST=1 시험 (위치 4개). 해석하지 말 것 ***\n")
cat(sprintf("위치 %d개 (13-1 전체; 1 ~ 2순위 94곳 %d · 1e-5 locus 있는 곳 %d) · 패널 수 중앙값 %g (최대 %d) · 그림 폴더 %s\n", nrow(R), sum(R$in94), sum(R$has_1e5),
            median(lengths(strsplit(R$panels, ","))),
            max(lengths(strsplit(R$panels, ","))), fig_dir))
cat(sprintf("index가 main 대표 SNP인 위치 %d개 · 아닌 위치 %d개\n", sum(R$index_model == "main"), sum(R$index_model != "main")))
cat(sprintf("13-1 다른 신호 후보 짝이 있는 위치 (그림 제목 표시): %d개\n", sum(R$multi_signal)))
pr_ls <- function(D, lab) {
cat(sprintf("\n■ LD 표본별 (index 기준, 창 %d kb; 이전 연구 기준 r² ≥ %.1f 강한 LD · ≥ %.1f hit) — %s 위치 %d개 중앙값\n", ld_kb, r2_strong, r2_hit, lab, uniqueN(D$region_id)))
cat("  각 모형 패널은 그 모형 GWAS 표본 LD로 색칠 (main · int_age · bin_all = pooled, bin_age8 = age8). n_lead_here = 그 표본 모형 대표 SNP가 있는 위치 수\n")
cat("  strong_far = 가장 먼 r² ≥ 0.8 SNP까지 거리, strong_far2 = 2번째로 먼 것 (먼 SNP 1개에 흔들리지 않는 값; 끊기지 않는 덩어리 길이 아님)\n")
print(f3(D[, .(n_dogs = n_dogs[1], n_lead_here = sum(in_panels), n_window = median(n_window), n_r2_ge_0.8 = median(n_r2_ge_0.8), n_r2_ge_0.2 = median(n_r2_ge_0.2),
               mean_r2 = median(mean_r2), strong_far_kb = median(strong_far_kb), strong_far2_kb = median(strong_far2_kb),
               cor_med = median(cor_with_pooled, na.rm = TRUE), cor_min = suppressWarnings(min(cor_with_pooled, na.rm = TRUE)),
               cor_na = sum(is.na(cor_with_pooled)), n_flag = sum(ld_diff_flag), n_flag_lead_here = sum(ld_diff_flag & in_panels)), by = group]), nrows = 100)
cat(sprintf("LD 다름 표시 (pooled와 r² 상관 < %.1f, 서술): %d개 (위치 × 표본) · 그 표본 모형 대표 SNP가 있는 위치에서 %d개\n", cor_flag,
            sum(D$ld_diff_flag), sum(D$ld_diff_flag & D$in_panels)))
}
pr_ls(LS[in94 == TRUE], "13-1 1 ~ 2순위 (이전 대상)")
if (any(LS[in94 == TRUE]$ld_diff_flag)) print(f3(LS[in94 == TRUE & ld_diff_flag == TRUE][order(region_id, group)]), nrows = 300)
pr_ls(LS, "전체")
cat("전체에서 그 표본 모형 대표 SNP가 있는 위치의 LD 다름 (표본별):\n"); print(LS[ld_diff_flag & in_panels, .N, by = group][order(-N)])
fs <- dcast(LS[group %in% strata3 & in94 == TRUE], region_id ~ group, value.var = c("strong_far_kb", "strong_far2_kb"))
fs <- fs[pmax(abs(strong_far_kb_purelike - strong_far_kb_pooled), abs(strong_far_kb_mixlike - strong_far_kb_pooled)) >= 30]
cat(sprintf("\n(94곳) 세 층 중 가장 먼 강한 LD SNP 거리가 pooled와 30 kb 이상 다른 위치 %d개 (2번째 거리도 같이; 1번째만 다르면 먼 SNP 1개 때문)\n", nrow(fs)))
if (nrow(fs)) print(f3(fs), nrows = 100)
cat(sprintf("\n■ 대표 SNP 짝 (13-1 _region_leads, 이 위치들) %d개: 13-1 r²와 일치 %d / %d\n", nrow(lp), sum(lp$r2_check, na.rm = TRUE), sum(!is.na(lp$r2_check))))
print(dcast(melt(lp[, .(class_pooled, class_purelike, class_mixlike)], measure.vars = 1:3, variable.factor = FALSE)[, .N, by = .(variable, value)],
            value ~ variable, value.var = "N", fill = 0L))
cat(sprintf("층에 따라 구분이 바뀌는 짝 (LD 기준에 민감): %d개\n", sum(lp$ld_sensitive)))
if (any(lp$ld_sensitive)) print(f3(lp[ld_sensitive == TRUE, .(region_id, model_a, lead_a, model_b, lead_b, dist_kb, r2_pooled, r2_purelike, r2_mixlike)]), nrows = 300)
cat(sprintf("\n■ locus별 LD (⑥; 모든 13-1 locus %d개, 그 모형 GWAS 표본, 대표 SNP ± %d kb) — 모형 × 기준 요약 (서술)\n", nrow(LL), ld_kb))
cat("  med_n_r2_ge_0.8 = 강한 LD 짝 수 중앙값 · pct_no_strong = 강한 LD 짝 없는 locus % · med_strong_far2_kb = 퍼진 거리 (2번째) 중앙값 ·\n")
cat("  med_frac_strong_p1e3 = 강한 LD 짝 중 GWAS p < 1e-3 비율 중앙값 · pct_ld_support = 강한 LD 짝 중 p < 1e-3이 1개 이상인 locus %\n")
for (sr in c("나이", "체중", "기타")) {
  cat(sprintf("[%s 계열]%s\n", sr, switch(sr, 나이 = " young8 (8세 미만) → main (전체) → 6 · 7 · 8 · 9 · 10 · 11 · 12세 이상 순", 체중 = " wt_t1 (가벼움) → wt_t3 (무거움) 순", "")))
  print(f3(LM[series == sr, !"series"]), nrows = 100)
}
cat("기준별 전체:\n")
print(f3(LL[, .(n_loci = .N, med_n_r2_ge_0.8 = as.numeric(median(n_r2_ge_0.8)), pct_no_strong = 100 * mean(n_r2_ge_0.8 == 0),
                med_strong_far2_kb = as.numeric(median(strong_far2_kb)), med_frac_strong_p1e3 = median(frac_strong_p1e3, na.rm = TRUE),
                pct_ld_support = 100 * mean(n_strong_p1e3 >= 1)), by = .(thr)]))
sink()
unlink(wd, recursive = TRUE)
cat(sprintf("저장: %s_regions.csv · _ld_strata.csv · _lead_ld_strata.csv · _locus_ld.csv · _locus_ld_by_model.csv · _summary.txt / 그림 %s/*.png · %s · _ld_samples.png · _locus_ld.png\n", out_prefix, fig_dir,
            file.path(pic_dir, paste0("13-2", tag, "_ld_strata.png"))))

