# 13-8_pathway.R - Pathway trends across age and body-weight models. See README.md.
suppressMessages({ library(data.table); library(ggplot2); library(parallel) })
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
pic_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/05_GWAS_2023"
test_mode <- Sys.getenv("TEST", "0") == "1"; tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("13-8", tag))
pre <- function(x) file.path(qc_dir, x)
n_cores <- as.integer(Sys.getenv("N_CORES", "8"))
plink1 <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
bfile <- file.path(d0, "03-3_maf0.01_hwe")
g22 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure"
f_in <- c(cmp = pre("13-1_compare.csv"), l33 = pre("13-3_loci.csv"), cs = pre("13-3_susie_cs.csv"), css = pre("13-3_susie_cs_snps.csv"),
          l34 = pre("13-4_loci_breed.csv"), l35 = pre("13-5_loci_cov.csv"), sp = pre("06-4_sets_primary.txt"), sa = pre("06-4_sets_aux.txt"),
          si = pre("06-4_set_info.csv"), gl = file.path(g22, "06-1_gene_loc.txt"), gi = file.path(g22, "06-1_gene_info.csv"))
for (f in c(f_in, paste0(bfile, ".bed"), plink1)) if (!file.exists(f)) stop("파일 없음: ", f)
gsa_f <- Sys.getenv("GSA_FILE", pre("13-7_gsa_all.csv.gz")); has_a2 <- file.exists(gsa_f)
null_dir <- file.path(d0, "06-5_gsa")
p_cut <- 0.05; r2_min <- 0.5; win_kb <- 10; ld_kb <- 1000; min_hit <- 3; aux_top <- 10
age_series <- c("young8", "main", "agefloor6", "agefloor7", "age8", "agefloor9", "age10", "agefloor11", "agefloor12")
wt_series <- c("wt_t1", "wt_t2", "wt_t3")
run_order <- c(age_series, wt_series, "sex_m", "sex_f", "purelike", "mixlike", "int_age", "int_wt", "bin_all", "bin_age8")
series_of <- function(m) fifelse(m %in% age_series, "나이", fifelse(m %in% wt_series, "체중", "기타"))
ld_grp <- function(m) fifelse(m %in% c("main", "int_age", "bin_all"), "pooled", fifelse(m == "bin_age8", "age8", m))
keep_of <- function(g) if (g == "pooled") file.path(d0, "04-1_grm", "grm_all.grm.id") else file.path(d0, "05-1_inputs", paste0(g, "_keep.txt"))
wd <- tempfile("r13_8_"); dir.create(wd)
f3 <- function(x) { x <- copy(x); for (j in names(x)) if (is.numeric(x[[j]]) && !is.integer(x[[j]])) set(x, j = j, value = signif(x[[j]], 3)); x }

read_sets <- function(f) { s <- strsplit(readLines(f), "[\t ]+"); setNames(lapply(s, function(x) as.integer(x[-1])), sapply(s, `[`, 1)) }
SP <- read_sets(f_in[["sp"]]); SA <- read_sets(f_in[["sa"]])
si <- fread(f_in[["si"]]); fam_of <- setNames(si$family, si$set)
gl <- fread(f_in[["gl"]], header = FALSE, col.names = c("GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"))

cmp <- fread(f_in[["cmp"]])[, .(model, SNP, thr, region_id, chr, bp)]
L33 <- fread(f_in[["l33"]])[, .(model, SNP, thr, specass, shared)]
L35 <- fread(f_in[["l35"]])[, .(model, SNP, thr, own_kind, wt_source)]
L34 <- unique(fread(f_in[["l34"]])[, .(SNP, breed_conc, top_breed_far2_kb)], by = "SNP")
LC <- Reduce(function(a, b) merge(a, b, by = c("model", "SNP", "thr"), all.x = TRUE), list(cmp, L33, L35))
LC <- merge(LC, L34, by = "SNP", all.x = TRUE)
LC[, thr_set := fifelse(thr < 5e-5, "1e-5", "1e-4")]
if (test_mode) LC <- LC[SNP %in% unique(LC$SNP)[round(seq(1, uniqueN(LC$SNP), length.out = 100))]]
cs <- fread(f_in[["cs"]]); css <- fread(f_in[["css"]])
csl <- merge(css, cs[, .(region_id, model, cs, lead)], by = c("region_id", "model", "cs"))
cs_has <- csl[, .(has_lead = any(SNP == lead[1])), by = .(region_id, model, cs, lead)][has_lead == TRUE]
cs_snp <- merge(csl, cs_has[, .(region_id, model, cs)], by = c("region_id", "model", "cs"))[, .(model, SNP = lead, snp_in = SNP)]
bimp <- fread(paste0(bfile, ".bim"), header = FALSE, select = c(1, 2, 4), col.names = c("chr", "SNP", "bp")); setkey(bimp, SNP)
LC[, grp := ld_grp(model)]
need <- unique(LC[!paste(model, SNP) %in% cs_snp[, paste(model, SNP)], .(grp, SNP)])
cat(sprintf("[%s] A1: 대표 SNP %d개 · locus %d개 (CS 사용 %d · r² 사용 %d), LD 표본 %d개\n", format(Sys.time(), "%T"), uniqueN(LC$SNP), nrow(LC),
            LC[paste(model, SNP) %in% cs_snp[, paste(model, SNP)], .N], LC[!paste(model, SNP) %in% cs_snp[, paste(model, SNP)], .N], uniqueN(need$grp)))
R2 <- rbindlist(mclapply(unique(need$grp), function(g) {
  writeLines(need[grp == g, unique(SNP)], file.path(wd, paste0("idx_", g, ".txt")))
  out <- file.path(wd, paste0("ld_", g))
  system2(plink1, c("--bfile", bfile, "--dog", "--threads", "1", "--keep", keep_of(g), "--r2", "--ld-snp-list", file.path(wd, paste0("idx_", g, ".txt")),
                    "--ld-window-kb", ld_kb, "--ld-window", "999999", "--ld-window-r2", r2_min, "--out", out), stdout = FALSE, stderr = FALSE)
  if (!file.exists(paste0(out, ".ld"))) stop("plink LD 실패: ", g)
  x <- fread(paste0(out, ".ld"))[, .(grp = g, SNP = SNP_A, snp_in = SNP_B)]
  rbind(x, data.table(grp = g, SNP = need[grp == g, unique(SNP)], snp_in = need[grp == g, unique(SNP)]))
}, mc.cores = min(n_cores, uniqueN(need$grp))))
MS <- unique(LC[, .(model, SNP, grp)])
rng <- rbind(merge(MS, R2, by = c("grp", "SNP"), allow.cartesian = TRUE)[!paste(model, SNP) %in% cs_snp[, paste(model, SNP)], .(SNP, snp_in)],
             merge(MS[, .(model, SNP)], cs_snp, by = c("model", "SNP"), allow.cartesian = TRUE)[, .(SNP, snp_in)])
rng <- unique(rng); rng[, `:=`(c_in = bimp[.(snp_in), chr], b_in = bimp[.(snp_in), bp])]
lg <- rbindlist(lapply(split(rng, by = "SNP"), function(x) {
  c0 <- x$c_in[1]; gg <- gl[CHR == c0]; lo <- min(x$b_in) - win_kb * 1e3; hi <- max(x$b_in) + win_kb * 1e3
  gg <- gg[STOP >= lo & START <= hi]
  if (nrow(gg)) gg <- gg[vapply(seq_len(.N), function(i) any(x$b_in >= gg$START[i] - win_kb * 1e3 & x$b_in <= gg$STOP[i] + win_kb * 1e3), TRUE)]
  data.table(SNP = x$SNP[1], n_snp_range = nrow(x), span_kb = (max(x$b_in) - min(x$b_in)) / 1e3, GENE = if (nrow(gg)) gg$GENE else NA_integer_,
             SYMBOL = if (nrow(gg)) gg$SYMBOL else NA_character_)
}))
SN <- lg[, .(n_snp_range = n_snp_range[1], span_kb = span_kb[1], n_genes = sum(!is.na(GENE)), genes = paste(na.omit(SYMBOL), collapse = ";")), by = SNP]
hit_mat <- function(S) { gs <- split(lg[!is.na(GENE), GENE], lg[!is.na(GENE), SNP])
  M <- sapply(S, function(s) vapply(SN$SNP, function(x) any(gs[[x]] %in% s), TRUE)); if (is.null(dim(M))) M <- matrix(M, nrow = nrow(SN)); rownames(M) <- SN$SNP; M }
HP <- hit_mat(SP)
SN <- cbind(SN, as.data.table(HP)[, lapply(.SD, as.integer)])
SN[, prim_sets := apply(HP, 1, function(v) paste(colnames(HP)[v], collapse = ";"))]
cat(sprintf("[%s] 유전자 연결: 대표 SNP %d개 중 유전자 있음 %d개 (유전자 수 중앙값 %g), 주 가설 세트에 하나라도 걸림 %d개\n", format(Sys.time(), "%T"),
            nrow(SN), sum(SN$n_genes > 0), median(SN$n_genes), sum(rowSums(HP) > 0)))
LC[, series := series_of(model)]
grp_def <- list(
  `계열 × 공유 안 됨` = function(x) fifelse(x$shared %in% "공유 안 됨" & x$model != "main", paste0(x$series, " 계열 공유 안 됨"), NA_character_),
  `13-3 특이성 분류` = function(x) fifelse(x$specass %in% c("층 차이 없음 (공통)", "") | is.na(x$specass), NA_character_, x$specass),
  `13-5 자기 축 판정` = function(x) fifelse(x$own_kind %in% "다른 축으로 옮겨감", x$own_kind, NA_character_),
  `13-5 체중 해석` = function(x) fifelse(x$wt_source %in% c("품종 체격 쪽", "개체 체중 쪽", "둘 다"), x$wt_source, NA_character_))
GM <- rbindlist(lapply(names(grp_def), function(k) unique(LC[, .(SNP, thr_set, gdef = k, level = grp_def[[k]](.SD))][!is.na(level)])))
thr_sets <- c("1e-4", "1e-5")
loc_test <- function(Hm, sn_ids, lev_snps, n_genes) {
  y_in <- sn_ids %in% lev_snps; lg1 <- log1p(n_genes)
  rbindlist(lapply(colnames(Hm), function(s) {
    h <- Hm[, s]; a <- sum(h & y_in); if (a < min_hit || sum(h) == a) return(data.table(set = s, n_in = sum(y_in), hit_in = a, hit_out = sum(h & !y_in), or = NA_real_, p = NA_real_))
    f <- tryCatch(suppressWarnings(glm(h ~ y_in + lg1, family = binomial)), error = function(e) NULL)
    cf <- if (is.null(f)) c(NA, NA) else summary(f)$coefficients["y_inTRUE", c(1, 4)]
    data.table(set = s, n_in = sum(y_in), hit_in = a, hit_out = sum(h & !y_in), or = exp(cf[1]), p = cf[2])
  }))
}
bc_snp <- L34[breed_conc %in% TRUE, SNP]
A1 <- rbindlist(lapply(thr_sets, function(ts) {
  U <- SN$SNP %in% LC[thr_set == ts, SNP]
  rbindlist(lapply(split(unique(GM[thr_set == ts, .(gdef, level)]), by = c("gdef", "level")), function(gv) {
    ls <- GM[thr_set == ts & gdef == gv$gdef & level == gv$level, SNP]
    r_all <- loc_test(HP[U, , drop = FALSE], SN$SNP[U], ls, SN$n_genes[U])[, `:=`(thr_set = ts, gdef = gv$gdef, level = gv$level)]
    keep <- U & !SN$SNP %in% bc_snp
    r_nb <- loc_test(HP[keep, , drop = FALSE], SN$SNP[keep], ls, SN$n_genes[keep])[, .(set, or_nobreed = or, p_nobreed = p)]
    merge(r_all, r_nb, by = "set", all.x = TRUE)
  }))
}))
A1[, `:=`(family = fam_of[set], trend = !is.na(p) & p < p_cut & or > 1, holds_nobreed = !is.na(p_nobreed) & p_nobreed < p_cut & or_nobreed > 1)]
A1[, q_ref := p.adjust(p, "BH"), by = .(thr_set, gdef, level)]
setcolorder(A1, c("thr_set", "gdef", "level", "set", "family")); setorder(A1, thr_set, gdef, level, p)
HA <- hit_mat(SA[lengths(SA) >= 10 & lengths(SA) <= 500])
A1aux <- rbindlist(lapply(thr_sets, function(ts) { U <- SN$SNP %in% LC[thr_set == ts, SNP]; H <- HA[U, , drop = FALSE]
 rbindlist(lapply(split(unique(GM[thr_set == ts, .(gdef, level)]), by = c("gdef", "level")), function(gv) {
  ls <- GM[thr_set == ts & gdef == gv$gdef & level == gv$level, SNP]; y_in <- SN$SNP[U] %in% ls; if (sum(y_in) < 5) return(NULL)
  hi <- colSums(H[y_in, , drop = FALSE]); ho <- colSums(H[!y_in, , drop = FALSE]); ok <- hi >= min_hit
  if (!any(ok)) return(NULL)
  ft <- vapply(which(ok), function(j) fisher.test(matrix(c(hi[j], sum(y_in) - hi[j], ho[j], sum(!y_in) - ho[j]), 2), alternative = "greater")$p.value, 0)
  data.table(thr_set = ts, gdef = gv$gdef, level = gv$level, set = names(ft), hit_in = hi[ok], n_in = sum(y_in), hit_out = ho[ok], p = ft,
             q_ref = p.adjust(ft, "BH", n = ncol(H)))[order(p)][seq_len(min(aux_top, .N))]
 }))
}))

if (has_a2) {
  GA <- fread(cmd = paste("zcat", shQuote(gsa_f)))
  P <- GA[set %in% names(SP)]
  wide <- function(v, w) dcast(P[variant == v & window == w], set ~ model, value.var = c("BETA", "SE", "P"))
  B10 <- P[variant == "primary" & window == 10L, .(set, model, BETA, SE, P)]
  B50 <- P[variant == "primary" & window == 50L, .(set, model, P50 = P, BETA50 = BETA)]
  BNB <- P[variant == "nobreed" & window == 10L, .(set, model, Pnb = P, BETAnb = BETA)]
  SBR <- Reduce(function(a, b) merge(a, b, by = c("set", "model"), all.x = TRUE), list(B10, B50, BNB))
  SBR[, `:=`(family = fam_of[set], z = BETA / SE, model = factor(model, run_order))]; setorder(SBR, set, model); SBR[, model := as.character(model)]
  getb <- function(s, m, col = "BETA") SBR[set == s & model == m][[col]]
  diff_z <- function(s, a, b, sfx = "") { ba <- getb(s, a, paste0("BETA", sfx)); bb <- getb(s, b, paste0("BETA", sfx)); sa <- getb(s, a, "SE"); sb <- getb(s, b, "SE")
    z <- (ba - bb) / sqrt(sa^2 + sb^2); 2 * pnorm(-abs(z)) }
  trend <- function(s, sfx = "") { b <- sapply(wt_series, function(m) getb(s, m, paste0("BETA", sfx))); se <- sapply(wt_series, function(m) getb(s, m, "SE"))
    w <- 1 / se^2; x <- 1:3; xb <- sum(w * x) / sum(w); sl <- sum(w * (x - xb) * b) / sum(w * (x - xb)^2); sls <- sqrt(1 / sum(w * (x - xb)^2))
    mu <- sum(w * b) / sum(w); Q <- sum(w * (b - mu)^2)
    c(slope = sl, slope_p = 2 * pnorm(-abs(sl / sls)), Q_p = pchisq(Q, 2, lower.tail = FALSE)) }
  ST <- rbindlist(lapply(names(SP), function(s) {
    tr <- trend(s); tr50 <- trend(s, "50"); trnb <- trend(s, "nb")
    rbind(
      data.table(set = s, test = "int_wt (SNP × 체중)", stat = getb(s, "int_wt"), p = getb(s, "int_wt", "P"), p50 = getb(s, "int_wt", "P50"), pnb = getb(s, "int_wt", "Pnb")),
      data.table(set = s, test = "체중 3분위 추세", stat = tr[["slope"]], p = tr[["slope_p"]], p50 = tr50[["slope_p"]], pnb = trnb[["slope_p"]]),
      data.table(set = s, test = "체중 3분위 Q (참고)", stat = NA_real_, p = tr[["Q_p"]], p50 = tr50[["Q_p"]], pnb = trnb[["Q_p"]]),
      data.table(set = s, test = "int_age (SNP × 나이)", stat = getb(s, "int_age"), p = getb(s, "int_age", "P"), p50 = getb(s, "int_age", "P50"), pnb = getb(s, "int_age", "Pnb")),
      data.table(set = s, test = "age8 − young8", stat = getb(s, "age8") - getb(s, "young8"), p = diff_z(s, "age8", "young8"), p50 = diff_z(s, "age8", "young8", "50"), pnb = diff_z(s, "age8", "young8", "nb")),
      data.table(set = s, test = "sex_f − sex_m", stat = getb(s, "sex_f") - getb(s, "sex_m"), p = diff_z(s, "sex_f", "sex_m"), p50 = diff_z(s, "sex_f", "sex_m", "50"), pnb = diff_z(s, "sex_f", "sex_m", "nb")),
      data.table(set = s, test = "purelike − mixlike", stat = getb(s, "purelike") - getb(s, "mixlike"), p = diff_z(s, "purelike", "mixlike"), p50 = diff_z(s, "purelike", "mixlike", "50"), pnb = diff_z(s, "purelike", "mixlike", "nb")),
      data.table(set = s, test = "main", stat = getb(s, "main"), p = getb(s, "main", "P"), p50 = getb(s, "main", "P50"), pnb = getb(s, "main", "Pnb")))
  }))
  ST[, `:=`(family = fam_of[set], trend = !grepl("참고", test) & p < p_cut & !is.na(p50) & p50 < p_cut, holds_nobreed = !is.na(pnb) & pnb < p_cut)]
  ST[, q_ref := p.adjust(p, "BH"), by = test]
  nl <- list.files(null_dir, pattern = "^null[0-9]+__primary\\.gsa\\.out$", full.names = TRUE)
  if (length(nl)) {
    NL <- rbindlist(lapply(nl, function(f) { d <- fread(f, skip = "VARIABLE"); if ("FULL_NAME" %in% names(d)) d[, VARIABLE := FULL_NAME]; d[, .(set = VARIABLE, P)] }))
    ST <- merge(ST, NL[, .(null_p05_ref = mean(P < p_cut)), by = set], by = "set", all.x = TRUE)
  }
  setorder(ST, test, p)
  AX <- GA[variant == "aux" & window == 10L]
  ax_tr <- AX[model %in% wt_series, { if (.N == 3) { o <- match(wt_series, model); b <- BETA[o]; se <- SE[o]; w <- 1 / se^2; x <- 1:3; xb <- sum(w * x) / sum(w)
              sl <- sum(w * (x - xb) * b) / sum(w * (x - xb)^2); .(stat = sl, p = 2 * pnorm(-abs(sl * sqrt(sum(w * (x - xb)^2))))) } else .(stat = NA_real_, p = NA_real_) }, by = set]
  AXT <- rbind(AX[model == "int_wt", .(test = "int_wt", set, NGENES, stat = BETA, p = P)],
               AX[model == "int_age", .(test = "int_age", set, NGENES, stat = BETA, p = P)],
               merge(ax_tr, unique(AX[, .(set, NGENES)], by = "set"), by = "set")[, .(test = "체중 3분위 추세", set, NGENES, stat, p)])
  AXT[, `:=`(q_ref = p.adjust(p, "BH"), bonf_ref = p < p_cut / .N), by = test]
  AXT <- AXT[order(test, p)][, head(.SD, 15), by = test]
}

if (has_a2) {
  hm <- copy(SBR); hm[, `:=`(model = factor(model, run_order), set = factor(set, rev(names(SP))))]
  g1 <- ggplot(hm, aes(model, set, fill = pmax(-4, pmin(4, z)))) + geom_tile(colour = "white") +
    geom_text(data = hm[P < p_cut], aes(label = "*"), size = 4, vjust = 0.75) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "set z\n(beta/SE)") +
    geom_vline(xintercept = c(length(age_series), length(age_series) + length(wt_series)) + 0.5, linetype = "dashed", colour = "grey40") +
    labs(x = "age series | weight tertiles | others", y = NULL, title = "MAGMA gene-set effects of the 14 pre-specified sets across models (* nominal p < 0.05)") +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(file.path(pic_dir, paste0("13-8", tag, "_set_heatmap.png")), g1, width = 10, height = 6, dpi = 200)
  wtp <- SBR[model %in% wt_series][, model := factor(model, wt_series)]
  iw <- ST[test == "int_wt (SNP × 체중)", .(set, lab = sprintf("int_wt p %.2g", p))]
  g2 <- ggplot(wtp, aes(model, BETA)) + geom_hline(yintercept = 0, colour = "grey50") + geom_pointrange(aes(ymin = BETA - 1.96 * SE, ymax = BETA + 1.96 * SE), size = 0.25) +
    geom_line(aes(group = 1), linewidth = 0.3) + facet_wrap(~set, ncol = 5, scales = "free_y") +
    geom_text(data = iw, aes(x = 2, y = Inf, label = lab), vjust = 1.3, size = 2.5, inherit.aes = FALSE) +
    labs(x = "weight tertile (light -> heavy)", y = "MAGMA set beta (95% CI)", title = "Pathway effect by body-weight tertile, with the genome-wide SNP x weight interaction (int_wt) set p") +
    theme_bw(base_size = 8)
  ggsave(file.path(pic_dir, paste0("13-8", tag, "_weight_trend.png")), g2, width = 11, height = 6.5, dpi = 200)
}
a1p <- A1[!is.na(or)]
if (nrow(a1p)) {
  a1p[, lab := paste0(thr_set, " · ", gdef, ": ", level)]
  g3 <- ggplot(a1p, aes(log2(or), set, colour = p < p_cut)) + geom_vline(xintercept = 0, colour = "grey50") + geom_point() + facet_wrap(~lab, ncol = 3) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60"), name = "nominal p < 0.05") +
    labs(x = "log2 OR (locus group vs other 13-1 loci; adjusted for log genes in range)", y = NULL,
         title = "Pre-specified pathway sets among 13-1 loci of each specificity group (locus -> genes -> set)") + theme_bw(base_size = 7) + theme(legend.position = "bottom")
  ggsave(file.path(pic_dir, paste0("13-8", tag, "_locus_sets.png")), g3, width = 11, height = max(6, 1.6 * ceiling(uniqueN(a1p$lab) / 3) + 2), dpi = 200, limitsize = FALSE)
}

fwrite(f3(SN), paste0(out_prefix, "_locus_genes.csv")); fwrite(f3(A1), paste0(out_prefix, "_locus_sets.csv")); fwrite(f3(A1aux), paste0(out_prefix, "_locus_aux_top.csv"))
if (has_a2) { fwrite(f3(SBR), paste0(out_prefix, "_set_by_run.csv")); fwrite(f3(ST), paste0(out_prefix, "_set_tests.csv")); fwrite(f3(AXT), paste0(out_prefix, "_aux_top.csv")) }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("13-8 생물학적 경로 경향 (판정: 명목 p < 0.05 + 사전 조건; BH q · null 비율 = 참고 열)\n"); if (test_mode) cat("*** TEST=1 시험. 해석하지 말 것 ***\n")
cat(sprintf("\n■ A1 locus → 유전자 → 세트: 대표 SNP %d개, 유전자 있음 %d개, SNP 범위 넓이 중앙값 %.1f kb, 유전자 수 중앙값 %g\n",
            nrow(SN), sum(SN$n_genes > 0), median(SN$span_kb), median(SN$n_genes)))
cat("주 가설 세트별 대표 SNP 수 (세트 유전자를 하나라도 가짐):\n"); print(sort(colSums(HP), decreasing = TRUE))
cat("집단 크기 (대표 SNP 수; 기준 집합별 — 1e-4 집합은 1e-5 locus 포함):\n"); print(dcast(GM[, .N, by = .(thr_set, gdef, level)], gdef + level ~ thr_set, value.var = "N", fill = 0L))
cat("\n경향 (집단 × 주 가설 세트, OR > 1 · 명목 p < 0.05; holds_nobreed = 품종 집중 SNP를 빼도 p < 0.05; q_ref = 참고):\n")
print(f3(A1[trend == TRUE, .(thr_set, gdef, level, set, family, n_in, hit_in, hit_out, or, p, or_nobreed, p_nobreed, holds_nobreed, q_ref)]), nrows = 200)
cat(sprintf("(참고: 집단 × 세트 검정 %d개 중 hit ≥ %d로 계산된 것 %d개 — 우연이면 그 중 약 %.0f개가 p < 0.05)\n", nrow(A1), min_hit, sum(!is.na(A1$p)), p_cut * sum(!is.na(A1$p))))
cat("\n보조 GO/KEGG (설명용; 집단마다 Fisher 한쪽 p 상위 — 판정 아님):\n"); print(f3(A1aux[, head(.SD, 5), by = .(thr_set, gdef, level)]), nrows = 200)
if (has_a2) {
  cat("\n■ A2 게놈 전체 세트 효과 (MAGMA 13-7): 검정별 결과 (trend = p < 0.05 & ±50 kb p < 0.05; holds_nobreed = 하플로타입 유전자 뺀 판 p < 0.05)\n")
  for (t in unique(ST$test)) { cat(sprintf("[%s]\n", t)); print(f3(ST[test == t, .(set, family, stat, p, p50, pnb, trend, holds_nobreed, q_ref, null_p05_ref = if ("null_p05_ref" %in% names(ST)) null_p05_ref else NA)]), nrows = 20) }
  cat("\n경향으로 표시된 것 (검정 × 세트):\n"); print(f3(ST[trend == TRUE, .(test, set, family, stat, p, p50, pnb, holds_nobreed, q_ref)]), nrows = 100)
  cat("\n나이 · 체중 계열 세트 z (서술; 나이 하한 계열은 표본이 겹침):\n")
  print(dcast(SBR[model %in% c(age_series, wt_series)], set ~ factor(model, c(age_series, wt_series)), value.var = "z", fun.aggregate = function(v) if (length(v)) round(v[1], 2) else NA_real_))
  cat("\n보조 GO/KEGG 상위 (설명용):\n"); print(f3(AXT), nrows = 60)
} else cat("\n■ A2: 13-7 결과 (13-7_gsa_all.csv.gz) 없음 → 건너뜀\n")
sink()
unlink(wd, recursive = TRUE)
cat(sprintf("저장: %s_* / 그림 13-8%s_*.png\n", out_prefix, tag))

