# 13-1_suggestive_loci.R - Lead loci (p < 1e-5 and 1e-4) per model, cross-stratum consistency, robustness across covariate models. See README.md.
suppressMessages(library(data.table))
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
out_prefix <- file.path(qc_dir, "13-1")
thrs <- c(`1e-5` = 1e-5, `1e-4` = 1e-4)
null_loci <- c(`1e-5` = 13.5, `1e-4` = 105)
gap <- 5e5; flank <- 2.5e5
p_rep <- 0.025; p_spec <- 0.05; q_cut <- 0.05; rob_ratio <- 0.5
gl <- fread("/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/06-1_gene_loc.txt", header = FALSE,
            col.names = c("GENE", "chr", "start", "end", "strand", "sym"))

fam <- c(main = "main", wt_t1 = "weight", wt_t2 = "weight", wt_t3 = "weight", sex_m = "sex", sex_f = "sex", purelike = "genetic", mixlike = "genetic",
         young8 = "age", age8 = "age", int_age = "age", agefloor6 = "age", agefloor7 = "age", agefloor9 = "age", age10 = "age", agefloor11 = "age",
         agefloor12 = "age", bin_all = "scale", bin_age8 = "scale")
disc <- data.table(id = names(fam), family = unname(fam))
disc[, in06 := id %in% c("main", "purelike", "mixlike", "int_age", "age8", "bin_age8", "young8", "sex_m", "sex_f", "wt_t1", "wt_t2", "wt_t3")]
rep_partners <- list(wt_t1 = c("wt_t2", "wt_t3"), wt_t2 = c("wt_t1", "wt_t3"), wt_t3 = c("wt_t1", "wt_t2"), sex_m = "sex_f", sex_f = "sex_m",
                     purelike = "mixlike", mixlike = "purelike", young8 = "age8", age8 = "young8",
                     agefloor9 = "young8", age10 = "young8", agefloor11 = "young8", agefloor12 = "young8")
trend_pairs <- rbind(CJ(a = c("wt_t1", "wt_t2", "wt_t3"), b = c("wt_t1", "wt_t2", "wt_t3"))[a != b],
                     data.table(a = c("sex_m", "sex_f", "purelike", "mixlike", "young8", "age8"), b = c("sex_f", "sex_m", "mixlike", "purelike", "age8", "young8")))
het_cmp <- c("sex", "age", "ancestry", "weight")
rob_cov <- c(paste0("PC", c(0, 10, 15, 20, 25, 30)), "pc_all", "pc_qmat", "add_ancestry", "add_mimic_dx", "add_oa", "add_health", "add_household", "vision_narrow")
rob_sub <- c("nomimic", "dob_exact", "afus", "no_bhs", "add_ac", "add_ses")
rob_diag <- c("raw", paste0("drop_", c("age", "pa", "wt", "owner_age", "sex", "intact", "vision", "bhs")))
floors <- c(young8 = 4, agefloor6 = 6, agefloor7 = 7, age8 = 8, agefloor9 = 9, age10 = 10, agefloor11 = 11, agefloor12 = 12)
gfile <- function(id) {
  if (id == "main") return("04-5_gwas/rint_maf0.01_PC5.mlma.gz")
  if (grepl("^PC[0-9]+$", id)) return(sprintf("04-5_gwas/rint_maf0.01_%s.mlma.gz", id))
  if (id %in% c("int_age", "bin_age8", "bin_all")) return(sprintf("05-3_regenie/%s.regenie.gz", id))
  sprintf("05-2_gwas/%s.mlma.gz", id)
}
all_ids <- c(disc$id, rob_cov, rob_sub, rob_diag)
het_files <- file.path(d0, sprintf("06-2_het_%s.tsv.gz", het_cmp))
for (f in c(file.path(d0, sapply(all_ids, gfile)), het_files, file.path(d0, "06-2_het_summary.csv"))) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
tmp_snp <- tempfile(fileext = ".txt")

hs <- fread(file.path(d0, "06-2_het_summary.csv"))
rr <- rbindlist(lapply(hs$r_null, function(s) {
  m <- regmatches(s, gregexpr("r\\([^,]+,[^)]+\\)=-?[0-9.]+", s))[[1]]
  rbindlist(lapply(m, function(x) { y <- regmatches(x, regexec("r\\(([^,]+),([^)]+)\\)=(-?[0-9.]+)", x))[[1]]; data.table(a = y[2], b = y[3], r = as.numeric(y[4])) }))
}))
r_of <- function(x, y) {
  if (y == "young8" && x %in% c("agefloor9", "age10", "agefloor11", "agefloor12")) x <- "age8"
  v <- rr$r[(rr$a == x & rr$b == y) | (rr$a == y & rr$b == x)]
  if (!length(v)) stop("06-2 r_null에 없음: ", x, " · ", y)
  v[1]
}

read_gwas <- function(id, mode) {
  f <- file.path(d0, gfile(id)); rg <- grepl("regenie", f)
  test <- if (id == "int_age") "ADD-INT_SNPxage_cslb" else "ADD"
  cond <- if (mode == "top") { if (rg) sprintf("$8==\"%s\" && $12 > %g", test, -log10(max(thrs))) else sprintf("$9 < %g", max(thrs)) } else
    { if (rg) sprintf("$8==\"%s\" && ($3 in s)", test) else "($2 in s)" }
  cmd <- if (mode == "top") sprintf("gzip -dc %s | awk 'NR==1 || (%s)'", f, cond) else
    sprintf("gzip -dc %s | awk 'NR==FNR {s[$1]; next} FNR==1 || (%s)' %s -", f, cond, tmp_snp)
  x <- fread(cmd = cmd)
  if (rg) x[, .(Chr = CHROM, SNP = ID, bp = GENPOS, ea = ALLELE1, maf = pmin(A1FREQ, 1 - A1FREQ), b = BETA, se = SE, p = 10^(-LOG10P))]
  else x[, .(Chr, SNP, bp, ea = A1, maf = pmin(Freq, 1 - Freq), b, se, p)]
}
near_genes <- function(cc, pos, k = 3) {
  g <- gl[chr == cc]; if (!nrow(g)) return("")
  g[, dist := fifelse(pos >= start & pos <= end, 0, pmin(abs(pos - start), abs(pos - end)))]
  g <- g[order(dist)][seq_len(min(k, .N))]
  paste(sprintf("%s(%s)", g$sym, ifelse(g$dist == 0, "in", paste0(round(g$dist / 1e3), "kb"))), collapse = ", ")
}
make_loci <- function(x, thr) {
  s <- x[p < thr][order(Chr, bp)]; if (!nrow(s)) return(NULL)
  s[, locus := cumsum(c(1, diff(bp) > gap | diff(Chr) != 0))]
  s[, c(list(n_sig = .N, start = min(bp), end = max(bp)), .SD[which.min(p)]), by = locus][, locus := NULL][]
}
wc_z <- function(z, z_sel) {
  f <- function(m) -(dnorm(z - m, log = TRUE) - log(pnorm(m - z_sel) + pnorm(-m - z_sel)))
  optimize(f, c(-abs(z) - 3, abs(z) + 3))$minimum
}
het <- function(b, se) { ok <- is.finite(b) & is.finite(se); if (sum(ok) < 2) return(NA_real_); b <- b[ok]; w <- 1 / se[ok]^2; m <- sum(w * b) / sum(w)
  pchisq(sum(w * (b - m)^2), length(b) - 1, lower.tail = FALSE) }

cat("발견 모형 19개 읽는 중\n")
top <- lapply(setNames(disc$id, disc$id), read_gwas, mode = "top")
loci <- rbindlist(lapply(names(thrs), function(tn) rbindlist(lapply(disc$id, function(id) {
  l <- make_loci(top[[id]], thrs[[tn]]); if (is.null(l)) return(NULL)
  tp <- top[[id]]
  l[, `:=`(thr = tn, model = id, n_p1e4_flank = sapply(seq_len(.N), function(j) tp[Chr == l$Chr[j] & abs(bp - l$bp[j]) <= flank, .N]))]
}))))
loci <- merge(loci, disc, by.x = "model", by.y = "id", sort = FALSE)
loci[, z := sign(b) * qnorm(p / 2, lower.tail = FALSE)]
loci[, b_wc := se * mapply(wc_z, z, qnorm(thrs[thr] / 2, lower.tail = FALSE))]
loci[, has_1e5 := p < 1e-5]
loci[, genes := mapply(near_genes, Chr, bp)]
loci[, n_model := .N, by = .(thr, model)][, null_exp := null_loci[thr]][, fdr_count := pmin(1, null_exp / n_model)]
setorder(loci, Chr, bp); loci[, `:=`(lo = bp - flank, hi = bp + flank)]
loci[, region_id := { r <- integer(.N); k <- 0L; cur <- -Inf
  for (j in seq_len(.N)) { if (lo[j] > cur) { k <- k + 1L; cur <- hi[j] } else cur <- max(cur, hi[j]); r[j] <- k }; r }, by = Chr]
loci[, region_id := sprintf("chr%d_r%d", Chr, region_id)][, c("lo", "hi") := NULL]
loci[, row := .I]

writeLines(unique(loci$SNP), tmp_snp)
cat(sprintf("lead %d개를 모형 %d개 + 06-2 이질성 표 4개에서 찾는 중\n", length(unique(loci$SNP)), length(all_ids)))
look <- rbindlist(lapply(all_ids, function(id) read_gwas(id, "snp")[, model := id]))
look <- unique(look, by = c("SNP", "model"))
wide <- function(v) { m <- dcast(look, SNP ~ model, value.var = v); x <- as.matrix(m[, -1, with = FALSE]); rownames(x) <- m$SNP; x }
Mb <- wide("b"); Mse <- wide("se"); Mp <- wide("p"); Mea <- wide("ea")
gv <- function(M, s, m) M[cbind(match(s, rownames(M)), match(m, colnames(M)))]
ab <- function(s, m, ea) { b <- gv(Mb, s, m); e <- gv(Mea, s, m); ifelse(is.na(e), NA_real_, ifelse(e == ea, b, -b)) }
hl <- rbindlist(lapply(seq_along(het_cmp), function(i) {
  x <- fread(cmd = sprintf("gzip -dc %s | awk 'NR==FNR {s[$1]; next} FNR==1 || ($1 in s)' %s -", het_files[i], tmp_snp))
  x[, .(SNP, cmp = het_cmp[i], p_het, p_wtrend = if ("p_trend" %in% names(x)) p_trend else NA_real_)]
}))
H <- dcast(hl, SNP ~ cmp, value.var = "p_het"); setnames(H, het_cmp, paste0("het_", het_cmp, "_p"))
H <- merge(H, hl[cmp == "weight", .(SNP, het_weight_trend_p = p_wtrend)], by = "SNP", all.x = TRUE)
ref_ea <- gv(Mea, look$SNP, "main")
fwrite(look[, .(SNP, model, b_main_allele = signif(ifelse(is.na(ref_ea) | ea == ref_ea, b, -b), 4), se = signif(se, 4), p = signif(p, 3))],
       paste0(out_prefix, "_lookup.csv"))

cmp <- rbindlist(lapply(seq_len(nrow(loci)), function(j) {
  x <- loci[j]; o <- data.table(row = x$row)
  if (x$model %in% names(rep_partners)) {
    ms <- rep_partners[[x$model]]; s <- rep(x$SNP, length(ms))
    bB <- ab(s, ms, x$ea); sB <- gv(Mse, s, ms); pB <- gv(Mp, s, ms); rv <- vapply(ms, function(m) r_of(x$model, m), 0)
    zal <- sign(x$z) * bB / sB; u <- (zal - rv * abs(x$z)) / sqrt(1 - rv^2); pm <- pnorm(-u)
    rp <- if (all(is.na(pm))) NA_real_ else min(1, sum(!is.na(pm)) * min(pm, na.rm = TRUE))
    o[, `:=`(rep_partners = paste(sprintf("%s b=%.3f p=%.2g", ms, bB, pB), collapse = "; "), rep_p = rp,
             het_p_naive = het(c(x$b, bB), c(x$se, sB)), het_p_wc = het(c(x$b_wc, bB), c(x$se, sB)))]
  }
  if (x$model == "main" || (x$family == "age" && !x$model %in% c("int_age"))) {
    tb <- ab(rep(x$SNP, length(floors)), names(floors), x$ea); ok <- is.finite(tb) & names(floors) != "young8"
    o[, `:=`(age_traj = paste(sprintf("%s:%.3f", c("<8", 6:12), tb), collapse = " "),
             age_slope = if (sum(ok) >= 3) unname(coef(lm(tb[ok] ~ floors[ok]))[2]) else NA_real_,
             int_age_b = ab(x$SNP, "int_age", x$ea), int_age_p = gv(Mp, x$SNP, "int_age"))]
  }
  if (x$model %in% c("main", "bin_all", "bin_age8", "age8")) {
    pr <- c(main = "bin_all", bin_all = "main", bin_age8 = "age8", age8 = "bin_age8")[[x$model]]
    o[, scale_partner := sprintf("%s b=%.3f p=%.2g", pr, ab(x$SNP, pr, x$ea), gv(Mp, x$SNP, pr))]
  }
  if (x$model == "main") {
    ratio <- ab(rep(x$SNP, length(rob_cov)), rob_cov, x$ea) / x$b
    o[, `:=`(rob_min_ratio = min(ratio, na.rm = TRUE), rob_worst = rob_cov[which.min(ratio)], rob_stable = all(is.finite(ratio)) && all(ratio >= rob_ratio))]
  }
  o
}), fill = TRUE)
for (v in c("rep_partners", "rep_p", "het_p_naive", "het_p_wc", "age_traj", "age_slope", "int_age_b", "int_age_p", "scale_partner", "rob_min_ratio", "rob_worst", "rob_stable"))
  if (!v %in% names(cmp)) cmp[, (v) := NA]
cm <- merge(loci, cmp, by = "row")
cm <- merge(cm, H, by = "SNP", all.x = TRUE, sort = FALSE)
cm[, sexlink := is.na(het_sex_p) & !is.na(het_age_p)]
cm[!is.na(rep_p), q_rep := p.adjust(rep_p, "BH"), by = .(thr, model)]
cm[, rep_status := fifelse(is.na(rep_p), NA_character_, fifelse(rep_p < p_rep, "명목 재현", "재현 안 됨"))]
for (h in het_cmp) cm[model == "main", (paste0("het_", h, "_q")) := p.adjust(get(paste0("het_", h, "_p")), "BH"), by = thr]
cm[model == "main", int_age_q := p.adjust(int_age_p, "BH"), by = thr]
spec_lab <- c(het_cmp, "int_age"); spec_pv <- c(paste0("het_", het_cmp, "_p"), "int_age_p")
hitm <- matrix(sapply(spec_pv, function(v) !is.na(cm[[v]]) & cm[[v]] < p_spec), nrow = nrow(cm))
spec_s <- apply(hitm, 1, function(v) paste(spec_lab[v], collapse = "+"))
cm[, spec_main := fifelse(model == "main" & spec_s != "", spec_s, NA_character_)]
cm[!is.na(rep_p), spec_desc := rep_p >= p_rep & !is.na(het_p_wc) & het_p_wc < 0.05]
cm[, tier := fifelse(!is.na(rep_p) & rep_p < p_rep, 1L, fifelse(model == "main" & rob_stable %in% TRUE, 2L, 3L))]
cm[, flags := gsub("^; |; $", "", paste(fifelse(sexlink, "sexlink", ""), fifelse(maf < 0.05, "maf<0.05", ""), sep = "; "))]

plink1 <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
ld_bfile <- file.path(d0, "03-3_maf0.01_hwe"); ld_keep <- file.path(d0, "04-1_grm", "grm_all.grm.id")
ld_hi <- 0.5; ld_lo <- 0.1
rl <- unique(cm[, .(region_id, model, SNP, Chr, bp, ea, b, p)], by = c("region_id", "model", "SNP"))
rl <- rl[, if (uniqueN(model) >= 2) .SD, by = region_id]
rpair <- rl[rl, on = .(region_id), allow.cartesian = TRUE, nomatch = NULL][model != i.model]
setnames(rpair, c("model", "SNP", "bp", "ea", "b", "p", "i.model", "i.SNP", "i.bp", "i.ea", "i.b", "i.p"),
         c("model_a", "lead_a", "bp_a", "ea_a", "b_a", "p_a", "model_b", "lead_b", "bp_b", "ea_b", "b_b", "p_b"))
rpair[, i.Chr := NULL]
ld <- data.table(SNP_A = character(), SNP_B = character(), R2 = numeric())
ld_snps <- unique(c(rpair$lead_a, rpair$lead_b))
if (length(ld_snps) >= 2) {
  ld_dir <- tempfile("ld13_"); dir.create(ld_dir); writeLines(ld_snps, file.path(ld_dir, "snps.txt"))
  system2(plink1, c("--bfile", ld_bfile, "--keep", ld_keep, "--extract", file.path(ld_dir, "snps.txt"), "--dog", "--r2",
                    "--ld-window", "99999", "--ld-window-kb", "10000", "--ld-window-r2", "0", "--threads", "4", "--out", file.path(ld_dir, "ld")),
          stdout = FALSE, stderr = FALSE)
  if (!file.exists(file.path(ld_dir, "ld.ld"))) stop("plink LD 실패: ", ld_dir)
  ld <- fread(file.path(ld_dir, "ld.ld"))[, .(SNP_A, SNP_B, R2)]
  unlink(ld_dir, recursive = TRUE)
}
ld <- rbind(ld, ld[, .(SNP_A = SNP_B, SNP_B = SNP_A, R2)])
rpair[ld, on = .(lead_a = SNP_A, lead_b = SNP_B), r2 := i.R2]
rpair[lead_a == lead_b, r2 := 1]
rpair[, `:=`(dist_kb = round(abs(bp_a - bp_b) / 1e3, 1),
             p_b_at_a = gv(Mp, lead_a, model_b), p_a_at_b = gv(Mp, lead_b, model_a),
             same_dir_at_a = sign(ab(lead_a, model_b, ea_a)) == sign(b_a))]
rpair[, ld_signal := fifelse(is.na(r2), "LD 없음", fifelse(r2 >= ld_hi, "같은 신호", fifelse(r2 < ld_lo, "다른 신호 후보", "불분명")))]
rs_row <- rpair[, .(region_n_model = uniqueN(model_b) + 1L, region_r2_min = min(r2, na.rm = TRUE), region_r2_max = max(r2, na.rm = TRUE),
                    region_other = paste(sprintf("%s %s r2=%s p@this=%.2g%s", model_b, fifelse(lead_b == lead_a, "(같은 SNP)", lead_b),
                                                 fifelse(is.na(r2), "NA", sprintf("%.2f", r2)), p_b_at_a, fifelse(same_dir_at_a %in% FALSE, " 반대방향", "")),
                                         collapse = "; ")), by = .(model = model_a, SNP = lead_a)]
rs_row[!is.finite(region_r2_min), region_r2_min := NA][!is.finite(region_r2_max), region_r2_max := NA]
cm <- merge(cm, rs_row, by = c("model", "SNP"), all.x = TRUE, sort = FALSE)
cm[is.na(region_n_model), region_n_model := 1L]
rpair_out <- rpair[model_a < model_b][, .(region_id, chr = Chr, model_a, lead_a, bp_a, p_a = signif(p_a, 3), p_b_at_a = signif(p_b_at_a, 3),
                   model_b, lead_b, bp_b, p_b = signif(p_b, 3), p_a_at_b = signif(p_a_at_b, 3), dist_kb, r2 = round(r2, 3), ld_signal, same_dir_at_a)]
setorder(rpair_out, chr, region_id, model_a, model_b)
fwrite(rpair_out, paste0(out_prefix, "_region_leads.csv"))
setorder(cm, thr, model, Chr, bp)
fwrite(cm[, .(thr, model, family, in06, region_id, chr = Chr, bp, lead = SNP, ea, maf = round(maf, 3), p = signif(p, 3), b = signif(b, 3), b_wc = signif(b_wc, 3),
              se = signif(se, 3), has_1e5, n_sig, start, end, n_p1e4_flank, genes, n_model, null_exp, fdr_count = signif(fdr_count, 3))],
       paste0(out_prefix, "_loci.csv"))
fwrite(cm[, !c("row", "z", "Chr"), with = FALSE][, chr := cm$Chr], paste0(out_prefix, "_compare.csv"))

pri <- cm[order(tier, -has_1e5, thr != "1e-5", p), .(tier, thr, model, region_id, chr = Chr, bp, lead = SNP, p = signif(p, 3), b = signif(b, 3), b_wc = signif(b_wc, 3),
          maf = round(maf, 3), rep_status, rep_p = signif(rep_p, 3), q_rep = signif(q_rep, 3), rep_partners, rob_stable, spec_main, spec_desc,
          int_age_p = signif(int_age_p, 3), region_n_model, region_r2_min = round(region_r2_min, 2), flags, genes)]
fwrite(pri, paste0(out_prefix, "_priority.csv"))
rob_snps <- unique(cm[model == "main", SNP])
rob_tab <- dcast(look[SNP %in% rob_snps & model %in% c("main", rob_cov, rob_sub, rob_diag)], SNP ~ model, value.var = "b")
rob_tab <- merge(cm[model == "main", .(thr_min = fifelse(any(thr == "1e-5"), "1e-5", "1e-4")), by = SNP], rob_tab, by = "SNP")
fwrite(rob_tab, paste0(out_prefix, "_robust_main.csv"))

rs <- cm[!is.na(rep_p), .(n_lead = .N, n_nominal = sum(rep_p < p_rep), expected = p_rep * .N,
                          p_excess = binom.test(sum(rep_p < p_rep), .N, p_rep, alternative = "greater")$p.value, n_q05 = sum(q_rep < q_cut)), by = .(thr, model)]
rs <- rbind(rs, rs[, .(model = "합계", n_lead = sum(n_lead), n_nominal = sum(n_nominal), expected = sum(expected),
                       p_excess = binom.test(sum(n_nominal), sum(n_lead), p_rep, alternative = "greater")$p.value, n_q05 = sum(n_q05)), by = thr])
rs[, fdr_upper := pmin(1, expected / pmax(n_nominal, 1))]
fwrite(rs, paste0(out_prefix, "_replication_summary.csv"))

ss <- rbindlist(lapply(names(thrs), function(tn) {
  m <- cm[model == "main" & thr == tn]
  rbindlist(lapply(seq_along(spec_lab), function(j) {
    pv <- m[[spec_pv[j]]]; n <- sum(!is.na(pv)); k <- sum(pv < p_spec, na.rm = TRUE)
    data.table(thr = tn, comparison = spec_lab[j], n_lead = n, n_p05 = k, expected = p_spec * n,
               p_excess = if (n) binom.test(k, n, p_spec, alternative = "greater")$p.value else NA_real_,
               n_q05 = sum(m[[sub("_p$", "_q", spec_pv[j])]] < q_cut, na.rm = TRUE))
  }))
}))
fwrite(ss, paste0(out_prefix, "_specificity_summary.csv"))

trend <- rbindlist(lapply(names(thrs), function(tn) {
  tr <- rbindlist(lapply(seq_len(nrow(trend_pairs)), function(j) {
    a <- trend_pairs$a[j]; b <- trend_pairs$b[j]; A <- loci[thr == tn & model == a]; if (!nrow(A)) return(NULL)
    bB <- ab(A$SNP, b, A$ea); sB <- gv(Mse, A$SNP, b); ok <- is.finite(bB) & sB > 0; A <- A[ok]; bB <- bB[ok]; sB <- sB[ok]; if (!nrow(A)) return(NULL)
    r <- r_of(a, b); zal <- sign(A$z) * bB / sB; u <- (zal - r * abs(A$z)) / sqrt(1 - r^2)
    e_same <- mean(pnorm(r * abs(A$z) / sqrt(1 - r^2)))
    al <- sign(A$b) * bB; w <- 1 / sB^2; ivw <- sum(w * al) / sum(w)
    data.table(thr = tn, a = a, b = b, r = r, n_lead = nrow(A), same_sign = mean(zal > 0), exp_same = e_same,
               p_sign = binom.test(sum(zal > 0), nrow(A), e_same)$p.value, z_trend = sum(u) / sqrt(nrow(A)), ivw_z = ivw * sqrt(sum(w)))
  }))
  tr[, p_trend := 2 * pnorm(-abs(z_trend))][, p_holm := p.adjust(p_trend, "holm")]
}))
fwrite(trend, paste0(out_prefix, "_trend.csv"))

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("13-1 suggestive locus (1e-5 · 1e-4): 독립 층 재현 · 특이성 · 강건성 · 우선순위 (판정 규칙은 헤더; 2026-11-15 개정 = 우연 기대 · BH q를 참고 값으로, 결과를 본 뒤 결정)\n")
cat(sprintf("판정: 명목 재현 = rep_p < %.3f (같은 방향 · p < 0.05) · 층 차이 / 나이 상호작용 후보 = 명목 p < %.2f. 우연 기대 수 · BH q는 참고 값 (걸러내지 않음)\n",
            p_rep, p_spec))
cat(sprintf("참고: 우연 기대 locus 수 (주 분석 표본 순열) 1e-5 %.1f · 1e-4 %.0f\n\n", null_loci[1], null_loci[2]))
cat("■ ① 모형별 locus 수 (참고: fdr_ = 순열 기대 / 관측)\n")
nl <- dcast(cm[, .(n = .N), by = .(model, thr)], model ~ thr, value.var = "n", fill = 0L)
nl <- merge(disc[, .(model = id, family)], nl, by = "model", all.x = TRUE, sort = FALSE)
for (tn in names(thrs)) { if (!tn %in% names(nl)) nl[, (tn) := 0L]; nl[is.na(get(tn)), (tn) := 0L]; nl[, (paste0("fdr_", tn)) := round(pmin(1, null_loci[tn] / pmax(get(tn), 1)), 2)] }
print(nl)
cat("\n■ ② 독립 층 재현 (모형 × 기준; 참고 열: expected = 우연 기대 수 · p_excess = 한쪽 이항 · n_q05 = BH q < 0.05 수 · fdr_upper = 기대 / 관측)\n")
print(f3(rs))
cat("\n■ ③ 방향 일치 경향 (참고: 우연 대비 검정, 모형 특이성 판정에 안 씀. A의 lead 전부를 B에서; z_trend = r 보정 u 합 / √n, 기준마다 Holm; ivw_z = 첫 판 방식, r 보정 없음)\n")
print(f3(trend))
for (tn in names(thrs)) {
  cat(sprintf("\n■ ④ 주 분석 lead (%s) 특이성 · 강건성: 06-2 p_het · int_age 명목 p < %.2f면 표시, 안정 = 공변량 모형 14개 모두 같은 방향 · |b| ≥ %.1f배\n", tn, p_spec, rob_ratio))
  m <- cm[model == "main" & thr == tn]
  cat(sprintf("  lead %d개 · 안정 %d개 · 층 차이 / 나이 상호작용 후보 (명목, 하나 이상) %d개 · 성 연관 의심 %d개\n", nrow(m), sum(m$rob_stable %in% TRUE), sum(!is.na(m$spec_main)), sum(m$sexlink)))
  cat("  비교별 명목 표시 수 · 참고 (우연 기대 수 · 한쪽 이항 p · BH q < 0.05 수):\n"); print(f3(ss[thr == tn, !"thr"]))
  if (tn == "1e-5" || any(!is.na(m$spec_main)))
    print(f3(m[if (tn == "1e-5") TRUE else !is.na(spec_main), .(region_id, lead = SNP, p, b, b_wc, maf, rob_stable, rob_min_ratio, het_sex_p, het_age_p, het_ancestry_p, het_weight_p, int_age_p, spec_main, flags, genes)]), nrows = 200)
}
cat("\n■ ⑤ 층 lead의 특이 후보 (서술; 명목 재현이 아니고 het_p_wc < 0.05) 수\n")
print(cm[!is.na(rep_p), .(n_lead = .N, spec_desc = sum(spec_desc, na.rm = TRUE)), by = .(thr, family)])
cat("\n■ ⑥ 우선순위 (1 명목 재현 · 2 주 분석 안정 · 3 나머지)\n")
print(dcast(cm[, .N, by = .(thr, tier)], thr ~ tier, value.var = "N", fill = 0L))
cat(sprintf("  13-2 대상 (1 ~ 2순위) region %d개\n", uniqueN(cm[tier <= 2, region_id])))
cat(sprintf("\n■ ⑦ 같은 위치의 다른 모형 대표 SNP (서술; r² ≥ %.1f 같은 신호 · < %.1f 다른 신호 후보; LD = 주 분석 7,115마리)\n", ld_hi, ld_lo))
cat(sprintf("  대표 SNP가 둘 이상 모형에 있는 위치 %d개 / 전체 %d개 · 모형이 다른 대표 SNP 짝 %d개\n", uniqueN(rl$region_id), uniqueN(cm$region_id), nrow(rpair_out)))
print(rpair_out[, .N, by = ld_signal][order(-N)])
rsum <- rpair_out[, .(n_model = uniqueN(c(model_a, model_b)), n_pair = .N, r2_min = min(r2, na.rm = TRUE), n_diff = sum(ld_signal == "다른 신호 후보"),
                      n_opp = sum(same_dir_at_a %in% FALSE)), by = .(region_id, chr)]
cat("  짝별 요약 (p_b_at_a = B 모형에서 본 A의 대표 SNP p; 표본이 겹치는 모형끼리는 독립 근거 아님):\n")
print(f3(rpair_out[, .(N = .N, med_r2 = median(r2, na.rm = TRUE), med_dist_kb = median(dist_kb), med_p_b_at_a = median(p_b_at_a, na.rm = TRUE),
                       p_b_at_a_lt05 = mean(p_b_at_a < 0.05, na.rm = TRUE), same_dir = mean(same_dir_at_a, na.rm = TRUE)), by = ld_signal]))
cat("  13-2 대상 (1 ~ 2순위) 위치 중 대표 SNP가 둘 이상 모형에 있는 곳:\n")
print(f3(rsum[region_id %in% cm[tier <= 2, region_id]][order(chr, region_id)]), nrows = 300)
cat("\n1순위 locus (q_rep는 참고):\n"); print(f3(pri[tier == 1, .(tier, thr, model, lead, p, b_wc, maf, rep_p, q_rep, rep_partners, spec_desc, flags, genes)]), nrows = 200)
sink()
cat(sprintf("저장: %s_loci.csv · _compare.csv · _priority.csv · _replication_summary.csv · _specificity_summary.csv · _region_leads.csv · _trend.csv · _robust_main.csv · _lookup.csv · _summary.txt\n", out_prefix))

