# 13-4_breed_ld.R - LD around lead SNPs within breeds. See README.md.
#   04-1_grm/grm_all.grm.id, 03-3_maf0.01_hwe (bfile), 03_GWAS_Pure/00_Breed_Count.csv

suppressMessages({ library(data.table); library(ggplot2); library(parallel) })
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
pic_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/05_GWAS_2023"
test_mode <- Sys.getenv("TEST", "0") == "1"; tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("13-4", tag))
in13 <- file.path(qc_dir, "13-1"); in132 <- file.path(qc_dir, "13-2")
n_cores <- as.integer(Sys.getenv("N_CORES", "8"))
plink1 <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
bfile <- file.path(d0, "03-3_maf0.01_hwe")
pooled_keep <- file.path(d0, "04-1_grm", "grm_all.grm.id")
anc_f <- file.path(d0, "08-1_ancestry_wide.csv"); size_f <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/00_Breed_Count.csv"
pure_cut <- 0.85; min_n <- 20
ld_kb <- 1000; r2_strong <- 0.8; r2_hit <- 0.2; cor_flag <- 0.8; ld_hi <- 0.5; ld_lo <- 0.1
rand_n <- c(20, 30, 50, 100, 200, 400); rand_rep <- 3; seed <- 20261115
in133 <- file.path(qc_dir, paste0("13-3", tag))
conc_share <- 0.5; conc_enrich <- 3; conc_min_mac <- 10
cs_maf <- 0.05
for (f in c(paste0(in132, c("_regions.csv", "_ld_strata.csv")), paste0(in13, "_region_leads.csv"), paste0(in133, c("_loci.csv", "_susie_cs.csv", "_susie_cs_snps.csv")), anc_f, pooled_keep, size_f, paste0(bfile, ".bed"), plink1))
  if (!file.exists(f)) stop("파일 없음: ", f)
wd <- tempfile("r13_4_"); dir.create(wd)

ids <- fread(pooled_keep, header = FALSE, col.names = c("FID", "IID"))
A <- fread(anc_f)[IID %in% ids$IID]; M <- as.matrix(A[, -(1:2)])
A[, `:=`(top = apply(M, 1, max), breed = colnames(M)[max.col(M, "first")])]
cat(sprintf("분석 표본 %d마리 중 ancestry 있음 %d · 유전체 순종형 (최대 품종 비율 ≥ %.2f) %d\n", nrow(ids), nrow(A), pure_cut, sum(A$top >= pure_cut)))
sz <- fread(size_f)[, .(breed = gsub("[^a-z]+", "_", tolower(Breed)), size_grp = get(paste0("size_", "class")))]
B <- A[top >= pure_cut, .N, by = breed][N >= min_n][order(-N)]
B <- merge(B, unique(sz, by = "breed"), by = "breed", all.x = TRUE)[order(-N)]
if (test_mode) B <- B[c(1:3, (nrow(B) - 2):nrow(B))]
cat(sprintf("품종 %d개 (≥ %d마리) · %d마리\n", nrow(B), min_n, sum(B$N)))
keeps <- list(pooled = pooled_keep)
for (b in B$breed) { f <- file.path(wd, paste0("keep_", b, ".txt")); fwrite(A[top >= pure_cut & breed == b, .(FID, IID)], f, sep = "\t", col.names = FALSE); keeps[[b]] <- f }
set.seed(seed)
for (n in rand_n) for (k in seq_len(rand_rep)) {
  nm <- sprintf("random_n%d_r%d", n, k); f <- file.path(wd, paste0("keep_", nm, ".txt"))
  fwrite(ids[sample(.N, n)], f, sep = "\t", col.names = FALSE); keeps[[nm]] <- f
}
n_of <- sapply(keeps, function(f) length(readLines(f)))

R <- fread(paste0(in132, "_regions.csv"))[, .(region_id, chr, index, index_bp, in94)]
if (test_mode) R <- R[region_id %in% c("chr4_r8", "chr7_r8", "chr5_r25", "chr2_r33")]
setorder(R, chr, index_bp)
cmp <- fread(paste0(in13, "_compare.csv"))[region_id %in% R$region_id]
U <- unique(rbind(R[, .(SNP = index, chr, bp = index_bp, region_id)], unique(cmp[, .(SNP, chr, bp, region_id)])), by = "SNP")
U[, is_index := SNP %in% R$index]; setorder(U, chr, bp)
lp <- fread(paste0(in13, "_region_leads.csv"))[region_id %in% R$region_id]
writeLines(U$SNP, file.path(wd, "idx.txt")); writeLines(unique(c(lp$lead_a, lp$lead_b, R$index)), file.path(wd, "leads.txt"))
run_plink <- function(args, out) {
  system2(plink1, c("--bfile", bfile, "--dog", "--threads", "1", args, "--out", out), stdout = FALSE, stderr = FALSE)
  if (file.exists(paste0(out, ".ld"))) fread(paste0(out, ".ld")) else NULL
}
cat(sprintf("[%s] LD 표본 %d개 (pooled 1 · 품종 %d · 무작위 %d) · 기준 SNP %d개 (위치 index %d)\n", format(Sys.time(), "%T"), length(keeps), nrow(B),
            length(rand_n) * rand_rep, nrow(U), sum(U$is_index)))
LDres <- mclapply(names(keeps), function(g) {
  x <- run_plink(c("--keep", keeps[[g]], "--r2", "with-freqs", "--ld-snp-list", file.path(wd, "idx.txt"), "--ld-window-kb", ld_kb, "--ld-window", "999999",
                   "--ld-window-r2", "0"), file.path(wd, paste0("idx_", g)))
  y <- run_plink(c("--keep", keeps[[g]], "--extract", file.path(wd, "leads.txt"), "--r2", "--ld-window-kb", "10000", "--ld-window", "999999", "--ld-window-r2", "0"),
                 file.path(wd, paste0("leads_", g)))
  list(ld = if (is.null(x)) NULL else x[, .(index = SNP_A, maf_index = MAF_A, SNP = SNP_B, bp = BP_B, r2 = R2)],
       lp = if (is.null(y)) NULL else rbind(y[, .(a = SNP_A, b = SNP_B, r2 = R2)], y[, .(a = SNP_B, b = SNP_A, r2 = R2)]))
}, mc.cores = min(n_cores, length(keeps)))
names(LDres) <- names(keeps)
bad <- names(LDres)[sapply(LDres, function(x) inherits(x, "try-error") || !is.list(x) || is.null(x$ld))]
if (length(bad)) stop("LD 계산 실패: ", paste(bad, collapse = ", "))

base_all <- LDres$pooled$ld[SNP != index & !is.na(r2), .(index, SNP, r2p = r2)]; setkey(base_all, index)
for (g in names(LDres)) setkey(LDres[[g]]$ld, index)
LS <- rbindlist(mclapply(names(keeps), function(g) {
  L <- LDres[[g]]$ld
  rbindlist(lapply(seq_len(nrow(U)), function(j) {
    ix <- U$SNP[j]; x <- L[.(ix), nomatch = NULL][SNP != ix]
    maf <- if (nrow(x)) x$maf_index[1] else NA_real_
    mono <- !nrow(x) || is.na(maf) || maf == 0 || all(is.na(x$r2))
    x <- x[!is.na(r2)]
    m <- merge(base_all[.(ix), nomatch = NULL][, .(SNP, r2p)], x[, .(SNP, r2)], by = "SNP"); mh <- m[pmax(r2p, r2) >= r2_hit]
    dd <- sort(abs(x[r2 >= r2_strong, bp] - U$bp[j]), decreasing = TRUE) / 1e3
    data.table(region_id = U$region_id[j], index = ix, is_index = U$is_index[j], group = g, n_dogs = n_of[[g]], maf_index = maf, index_mono = mono,
               n_window = nrow(x), n_r2_ge_0.8 = sum(x$r2 >= r2_strong), n_r2_ge_0.2 = sum(x$r2 >= r2_hit), mean_r2 = if (nrow(x)) mean(x$r2) else NA_real_,
               max_r2 = if (nrow(x)) max(x$r2) else NA_real_, strong_far_kb = if (length(dd)) dd[1] else 0, strong_far2_kb = if (length(dd) >= 2) dd[2] else 0,
               n_cor_snp = nrow(mh), cor_with_pooled = if (g == "pooled") 1 else if (!mono && nrow(mh) >= 10) cor(mh$r2p, mh$r2) else NA_real_)
  }))
}, mc.cores = min(n_cores, length(keeps))))
LS[, kind := fifelse(group == "pooled", "pooled", fifelse(grepl("^random_n", group), "random", "breed"))]
LS[, ld_diff_flag := kind == "breed" & !is.na(cor_with_pooled) & cor_with_pooled < cor_flag]
RN <- LS[kind == "random", .(ref_cor_med = median(cor_with_pooled, na.rm = TRUE), ref_cor_q10 = quantile(cor_with_pooled, 0.1, na.rm = TRUE),
                              ref_n_r2_ge_0.2 = as.numeric(median(n_r2_ge_0.2)), ref_mean_r2 = median(mean_r2, na.rm = TRUE)), by = .(n_dogs)]
near_n <- function(n) RN$n_dogs[which.min(abs(log(RN$n_dogs) - log(n)))]
LB <- LS[kind == "breed"]; LB[, ref_n := sapply(n_dogs, near_n)]
LB <- merge(LB, RN[, .(ref_n = n_dogs, ref_cor_med)], by = "ref_n", all.x = TRUE)
LB <- merge(LB, B[, .(group = breed, size_grp)], by = "group", all.x = TRUE)
LB[, cor_minus_ref := cor_with_pooled - ref_cor_med]
nbin <- function(n) cut(n, c(0, 49, 99, 199, Inf), labels = c("20-49", "50-99", "100-199", ">=200"))

cls <- function(r) fifelse(is.na(r), "LD 없음", fifelse(r >= ld_hi, "같은 신호", fifelse(r < ld_lo, "다른 신호 후보", "불분명")))
LPB <- rbindlist(lapply(c("pooled", B$breed, grep("^random_n", names(keeps), value = TRUE)), function(g) {
  y <- LDres[[g]]$lp; x <- copy(lp[, .(region_id, model_a, lead_a, model_b, lead_b, dist_kb, r2_13_1 = r2, ld_signal)])
  if (!is.null(y)) x[y, on = .(lead_a = a, lead_b = b), r2_g := i.r2]
  if (!"r2_g" %in% names(x)) x[, r2_g := NA_real_]
  x[lead_a == lead_b, r2_g := 1]; x[, `:=`(group = g, n_dogs = n_of[[g]], signal_g = cls(r2_g))]
}), fill = TRUE)
LPB[, changed := group != "pooled" & signal_g != "LD 없음" & signal_g != ld_signal]
LPB[, self_pair := lead_a == lead_b]
LPU <- unique(LPB[self_pair == FALSE][, `:=`(s1 = pmin(lead_a, lead_b), s2 = pmax(lead_a, lead_b))], by = c("group", "s1", "s2"))

age_series <- c("young8", "main", "agefloor6", "agefloor7", "age8", "agefloor9", "age10", "agefloor11", "agefloor12")
wt_series <- c("wt_t1", "wt_t2", "wt_t3")
series_order <- c(age_series, wt_series, "sex_m", "sex_f", "purelike", "mixlike", "int_age", "bin_all", "bin_age8")
series_of <- function(m) fifelse(m %in% age_series, "나이", fifelse(m %in% wt_series, "체중", "기타"))
SB <- LB[, .(n_breed_poly = sum(!index_mono), pct_breed_mono = 100 * mean(index_mono), breed_cor_med = median(cor_with_pooled, na.rm = TRUE),
             breed_n08_med = as.numeric(median(n_r2_ge_0.8[!index_mono])), breed_far2_med = as.numeric(median(strong_far2_kb[!index_mono]))), by = .(SNP = index)]
SB <- merge(SB, LS[kind == "pooled", .(SNP = index, maf_pooled = maf_index, pooled_n08 = n_r2_ge_0.8, pooled_far2 = strong_far2_kb)], by = "SNP", all.x = TRUE)
LCB <- merge(cmp[, .(model, SNP, thr, region_id, chr, bp, p)], SB, by = "SNP", all.x = TRUE)
LCB[, series := series_of(model)]; LCB[, model := factor(model, series_order)]; setorder(LCB, model, -thr, chr, bp); LCB[, model := as.character(model)]
LBM <- LCB[, .(n_loci = .N, med_maf_pooled = median(maf_pooled, na.rm = TRUE), med_pct_breed_mono = median(pct_breed_mono, na.rm = TRUE),
               med_breed_cor = median(breed_cor_med, na.rm = TRUE), med_pooled_far2 = as.numeric(median(pooled_far2, na.rm = TRUE)),
               med_breed_far2 = median(breed_far2_med, na.rm = TRUE), med_breed_n08 = median(breed_n08_med, na.rm = TRUE), series = series[1]), by = .(model, thr)]
LBM[, thr := factor(format(thr, scientific = TRUE), c("1e-05", "1e-04"))]; LBM[, model := factor(model, series_order)]; setorder(LBM, thr, model)

clst <- rbindlist(lapply(B$breed, function(b) fread(keeps[[b]], header = FALSE)[, .(FID = V1, IID = V2, CLST = b)]))
fwrite(clst, file.path(wd, "clst.txt"), sep = "\t", col.names = FALSE); fwrite(clst[, .(FID, IID)], file.path(wd, "pure_keep.txt"), sep = "\t", col.names = FALSE)
CSI <- fread(paste0(in133, "_susie_cs.csv"))[region_id %in% R$region_id]; CSN <- fread(paste0(in133, "_susie_cs_snps.csv"))[region_id %in% R$region_id]
writeLines(unique(c(cmp$SNP, CSN$SNP)), file.path(wd, "loci.txt"))
system2(plink1, c("--bfile", bfile, "--dog", "--threads", min(n_cores, 8), "--keep-allele-order", "--keep", file.path(wd, "pure_keep.txt"), "--extract", file.path(wd, "loci.txt"),
                  "--within", file.path(wd, "clst.txt"), "--freq", "--out", file.path(wd, "conc")), stdout = FALSE, stderr = FALSE)
if (!file.exists(file.path(wd, "conc.frq.strat"))) stop("품종별 빈도 실패: ", file.path(wd, "conc.log"))
FS <- fread(file.path(wd, "conc.frq.strat"))                                  # CHR SNP CLST A1 A2 MAC NCHROBS (A1 = bim A1, --keep-allele-order)
FS[, c1 := MAC]; FS[, c2 := NCHROBS - MAC]
tot <- FS[, .(t1 = sum(c1), tN = sum(NCHROBS)), by = SNP]
FS <- merge(FS, tot, by = "SNP"); FS[, minor_is_a1 := t1 <= tN / 2]
FS[, cm := fifelse(minor_is_a1, c1, c2)]
CC <- FS[, {
  tm <- sum(cm); k <- which.max(cm)
  .(pure_mac = tm, pure_maf = tm / sum(NCHROBS), n_breed_carry = sum(cm > 0), top_breed = CLST[k], top_breed_freq = cm[k] / NCHROBS[k],
    top_share = if (tm > 0) cm[k] / tm else NA_real_, expected = NCHROBS[k] / sum(NCHROBS))
}, by = SNP]
CC[, enrich := top_share / expected]; CC[, breed_conc := !is.na(top_share) & pure_mac >= conc_min_mac & top_share >= conc_share & enrich >= conc_enrich]
LCB <- merge(LCB, CC, by = "SNP", all.x = TRUE)
TB <- LB[, .(SNP = index, top_breed = group, top_breed_mono = index_mono, top_breed_n08 = n_r2_ge_0.8, top_breed_far2_kb = strong_far2_kb, top_breed_cor = cor_with_pooled)]
LCB <- merge(LCB, TB, by = c("SNP", "top_breed"), all.x = TRUE)
FB <- FS[SNP %in% CSN$SNP, .(SNP, breed = CLST, maf_b = pmin(c1, c2) / NCHROBS)]
CSB <- merge(CSN[, .(region_id, model, cs, SNP, pip)], FB, by = "SNP", allow.cartesian = TRUE)[, .(
  n_cs_snp = .N, frac_mono = mean(maf_b == 0), n_snp_maf05 = sum(maf_b >= cs_maf), pip_w_maf = sum(pip * maf_b) / sum(pip), maf_top_pip = maf_b[which.max(pip)]),
  by = .(region_id, model, cs, breed)]
CSB[, usable := n_snp_maf05 > 0]; CSB <- merge(CSB, B[, .(breed, n_dogs = N, size_grp)], by = "breed", all.x = TRUE)
CSS <- CSB[, .(n_breeds = .N, n_usable = sum(usable), med_frac_mono = median(frac_mono), med_pip_w_maf = median(pip_w_maf),
               usable_breeds = paste(head(breed[usable][order(-pip_w_maf[usable])], 5), collapse = ",")), by = .(region_id, model, cs)]
CSS <- merge(CSI[, .(region_id, model, cs, cs_size, narrow, top_snp, top_pip)], CSS, by = c("region_id", "model", "cs"), all.x = TRUE)[order(-narrow, cs_size)]
L33 <- fread(paste0(in133, "_loci.csv"))[, .(model, SNP, thr, specass, spec_axes, shared, partner_call)]
LCB <- merge(LCB, L33, by = c("model", "SNP", "thr"), all.x = TRUE)
LCB[, model := factor(model, series_order)]; setorder(LCB, model, -thr, chr, bp); LCB[, model := as.character(model)]
BCM <- LCB[, .(n_loci = .N, n_breed_conc = sum(breed_conc, na.rm = TRUE), pct_breed_conc = 100 * mean(breed_conc, na.rm = TRUE),
               med_top_share = median(top_share, na.rm = TRUE), med_n_breed_carry = as.numeric(median(n_breed_carry, na.rm = TRUE)),
               top_breeds = paste(head(names(sort(table(top_breed[breed_conc %in% TRUE]), decreasing = TRUE)), 3), collapse = ","),
               series = series_of(model[1])), by = .(model, thr, shared)]
BCM[, thr := factor(format(thr, scientific = TRUE), c("1e-05", "1e-04"))]; BCM[, model := factor(model, series_order)]; setorder(BCM, thr, model, shared)

s132 <- fread(paste0(in132, "_ld_strata.csv"))
wt_regs <- unique(s132[group %in% c("wt_t1", "wt_t2", "wt_t3") & ld_diff_flag == TRUE, region_id])

lev_b <- B[order(-N), sprintf("%s (n=%d)", breed, N)]
LB[, blab := factor(sprintf("%s (n=%d)", group, n_dogs), lev_b)]
rib <- RN[, .(n_dogs, ref_cor_med, ref_cor_q10)]
LB[, size_grp := fifelse(is.na(size_grp), "unknown", size_grp)]
g1 <- ggplot(LB[!is.na(cor_with_pooled)], aes(blab, cor_with_pooled)) +
  geom_hline(yintercept = cor_flag, linetype = "dashed", colour = "grey50") +
  geom_boxplot(aes(fill = size_grp), outlier.shape = NA, width = 0.65, alpha = 0.6) +
  geom_point(data = unique(LB[, .(blab, ref_cor_med)]), aes(blab, ref_cor_med), shape = 4, size = 2, colour = "black") +
  scale_fill_manual(values = c(Giant = "#B2182B", Large = "#E17C05", Standard = "#46ACC8", Medium = "#1B7837", `Toy and Small` = "#762A83", unknown = "grey70"),
                    name = "size class") +
  labs(x = NULL, y = "r2 correlation with pooled LD (SNPs with r2 >= 0.2 in either)", title = "Breed-specific LD around 13-2 index SNPs vs pooled (7,115 dogs)",
       subtitle = sprintf("genomic purebreds (ancestry >= %.2f), breeds >= %d dogs; x = median for random samples of similar size (reference: sampling noise only)",
                          pure_cut, min_n)) +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(pic_dir, paste0("13-4", tag, "_breed_cor.png")), g1, width = max(8, 0.28 * nrow(B) + 3), height = 6, dpi = 200)
R[, rlab := factor(region_id, rev(R$region_id))]
hm <- merge(LB[is_index == TRUE, .(region_id, blab, cor_with_pooled, index_mono)], R[in94 == TRUE, .(region_id, rlab)], by = "region_id")
hm[, rlab := factor(rlab, rev(R[in94 == TRUE, region_id]))]
g2 <- ggplot(hm, aes(blab, rlab, fill = cor_with_pooled)) + geom_tile(colour = "white", linewidth = 0.2) +
  geom_point(data = hm[index_mono == TRUE], aes(blab, rlab), shape = 4, size = 0.8, colour = "grey30", inherit.aes = FALSE) +
  scale_fill_gradientn(colours = c("#B2182B", "#E17C05", "#F7F7F7", "#46ACC8"), values = c(0, 0.5, 0.8, 1), limits = c(0, 1), na.value = "grey85",
                       name = "r2 cor\nwith pooled") +
  labs(x = NULL, y = NULL, title = "Breed LD vs pooled LD, 94 priority regions (grey = not computed: <10 SNPs or index monomorphic (x))") +
  theme_minimal(base_size = 8) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6), axis.text.y = element_text(size = 5))
ggsave(file.path(pic_dir, paste0("13-4", tag, "_breed_heatmap.png")), g2, width = max(8, 0.22 * nrow(B) + 3), height = max(5, 0.11 * uniqueN(hm$rlab) + 2), dpi = 200,
       limitsize = FALSE)
lb4 <- melt(LBM[, .(model, thr, `median breed r2-cor with pooled` = med_breed_cor, `% breeds where lead is monomorphic` = med_pct_breed_mono,
                    `pooled MAF of lead` = med_maf_pooled, `spread of r2 >= 0.8 SNPs within breeds (2nd, kb)` = med_breed_far2)],
            id.vars = c("model", "thr"), variable.factor = TRUE)
g3 <- ggplot(lb4, aes(model, value, colour = thr, group = thr)) + geom_point(size = 1.8) + geom_line(data = lb4[model %in% age_series], linewidth = 0.4) +
  geom_line(data = lb4[model %in% wt_series], linewidth = 0.4) + facet_wrap(~variable, ncol = 2, scales = "free_y") +
  geom_vline(xintercept = c(length(age_series), length(age_series) + length(wt_series)) + 0.5, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c(`1e-05` = "#B2182B", `1e-04` = "#46ACC8"), name = "locus threshold") +
  labs(x = "age series (young8 < 8 y, main = all, minimum age 6 ~ 12 y) | weight tertiles | others", y = "median over loci",
       title = "13-1 loci per model: breed-specific LD around the lead (descriptive)") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(pic_dir, paste0("13-4", tag, "_breed_by_model.png")), g3, width = 11, height = 7, dpi = 200)
bc4 <- LCB[!is.na(breed_conc), .(pct = 100 * mean(breed_conc), n = .N), by = .(model, thr, grp = fifelse(shared == "공유 안 됨", "not shared", "shared / unclear"))]
bc4[, `:=`(model = factor(model, series_order), thr = factor(format(thr, scientific = TRUE), c("1e-05", "1e-04")))]
g4 <- ggplot(bc4, aes(model, pct, fill = grp)) + geom_col(position = position_dodge(0.8), width = 0.75) + facet_wrap(~thr, ncol = 1) +
  geom_vline(xintercept = c(length(age_series), length(age_series) + length(wt_series)) + 0.5, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c(`not shared` = "#B2182B", `shared / unclear` = "#46ACC8"), name = NULL) +
  labs(x = "age series | weight tertiles | others", y = "% loci whose minor allele is concentrated in one breed",
       title = sprintf("Breed concentration of 13-1 lead alleles (top breed holds >= %.0f%% of copies among 45 breeds and >= %gx its share)", 100 * conc_share, conc_enrich)) +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom", plot.title = element_text(size = 9))
ggsave(file.path(pic_dir, paste0("13-4", tag, "_breed_conc.png")), g4, width = 11, height = 7, dpi = 200)

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
fwrite(B, paste0(out_prefix, "_breeds.csv")); fwrite(f3(rbind(LS[kind == "pooled"], LB[, names(LS), with = FALSE])), paste0(out_prefix, "_ld_breed.csv"))
fwrite(f3(LS[kind == "random"]), paste0(out_prefix, "_ld_random.csv")); fwrite(f3(LPB[!grepl("^random_n", group)]), paste0(out_prefix, "_lead_ld_breed.csv"))
fwrite(f3(LCB), paste0(out_prefix, "_loci_breed.csv")); fwrite(f3(LBM), paste0(out_prefix, "_loci_breed_by_model.csv")); fwrite(f3(BCM), paste0(out_prefix, "_breed_conc_by_model.csv"))
fwrite(f3(CSB), paste0(out_prefix, "_cs_breed.csv")); fwrite(f3(CSS), paste0(out_prefix, "_cs_breed_summary.csv"))
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("13-4 품종별 LD (서술; 우연 기대는 참고 열만)\n"); if (test_mode) cat("*** TEST=1 시험 (위치 4개 · 품종 6개). 해석하지 말 것 ***\n")
cat(sprintf("기준 SNP %d개 (13-1 locus 대표 SNP + 위치 %d곳 index, 창 %d kb) · 품종 %d개 (유전체 순종형 ≥ %.2f, ≥ %d마리; %d마리) · 무작위 참고 표본 %s × %d번\n",
            nrow(U), nrow(R), ld_kb, nrow(B),
            pure_cut, min_n, sum(B$N), paste(rand_n, collapse = "/"), rand_rep))
cat("\n■ 무작위 표본 참고 (분석 표본에서 품종 무관하게 뽑음 — 표본 크기만으로 생기는 차이; 기준 SNP 중앙값)\n"); print(f3(RN[order(n_dogs)]))
cat("\n■ 품종별 (기준 SNP 전체 중앙값; 기준 SNP 단형성은 뺌. ref_cor_med = 비슷한 크기 무작위 표본의 상관 중앙값, 참고; n_regions = 계산된 기준 SNP 수)\n")
PB <- LB[index_mono == FALSE, .(n_dogs = n_dogs[1], size_grp = size_grp[1], n_regions = .N, n_mono = nrow(LB[group == .BY$group & index_mono == TRUE]),
                                n_r2_ge_0.8 = as.numeric(median(n_r2_ge_0.8)), n_r2_ge_0.2 = as.numeric(median(n_r2_ge_0.2)), mean_r2 = median(mean_r2, na.rm = TRUE),
                                strong_far2_kb = as.numeric(median(strong_far2_kb)), cor_med = median(cor_with_pooled, na.rm = TRUE), ref_cor_med = ref_cor_med[1],
                                cor_minus_ref_med = median(cor_minus_ref, na.rm = TRUE),
                                n_flag = sum(ld_diff_flag)), by = group][order(-n_dogs)]
print(f3(PB), nrows = 100)
cat(sprintf("pooled: r² ≥ 0.8 수 중앙값 %g · r² ≥ 0.2 수 %g · 평균 r² %.4f · 퍼진 거리 (2번째) %.1f kb\n", LS[kind == "pooled", median(n_r2_ge_0.8)],
            LS[kind == "pooled", median(n_r2_ge_0.2)], LS[kind == "pooled", median(mean_r2)], LS[kind == "pooled", median(strong_far2_kb)]))
cat("\n체격 분류별 (품종 단위 중앙값의 중앙값):\n")
print(f3(PB[, .(n_breeds = .N, dogs = sum(n_dogs), n_r2_ge_0.8 = median(n_r2_ge_0.8), strong_far2_kb = median(strong_far2_kb), cor_med = median(cor_med, na.rm = TRUE),
                ref_cor_med = median(ref_cor_med)), by = size_grp][order(size_grp)]))
cat(sprintf("\n■ LD 다름 표시 (pooled와 상관 < %.1f, 서술): %d / %d (기준 SNP × 품종, 계산된 것) — 무작위 표본에서는 %d / %d (참고)\n", cor_flag, sum(LB$ld_diff_flag),
            sum(!is.na(LB$cor_with_pooled)), LS[kind == "random", sum(!is.na(cor_with_pooled) & cor_with_pooled < cor_flag)], LS[kind == "random", sum(!is.na(cor_with_pooled))]))
cat("  품종 크기 구간별 (⑨ a): 표시 비율 · 상관 중앙값 · cor_minus_ref 중앙값 (품종 − 비슷한 크기 무작위) | ref_pct_flag = 같은 구간 무작위 표본의 표시 비율 (참고)\n")
fb <- LB[!is.na(cor_with_pooled), .(n_breeds = uniqueN(group), n = .N, n_flag = sum(ld_diff_flag), pct_flag = 100 * mean(ld_diff_flag), cor_med = median(cor_with_pooled),
                                    cor_minus_ref_med = median(cor_minus_ref, na.rm = TRUE)), by = .(n_bin = nbin(n_dogs))]
fr <- LS[kind == "random" & !is.na(cor_with_pooled), .(ref_pct_flag = 100 * mean(cor_with_pooled < cor_flag)), by = .(n_bin = nbin(n_dogs))]
print(f3(merge(fb, fr, by = "n_bin", all.x = TRUE)[order(n_bin)]))
cat("94곳 위치 (index) 품종 상관 분포 (품종 상관 중앙값이 낮은 순 10 · 높은 순 5; 품종 목록은 _ld_breed.csv):\n")
PR <- LB[is_index == TRUE & region_id %in% R[in94 == TRUE, region_id], .(n_breeds_computed = sum(!is.na(cor_with_pooled)), n_flag = sum(ld_diff_flag), cor_med = median(cor_with_pooled, na.rm = TRUE),
             cor_min = suppressWarnings(min(cor_with_pooled, na.rm = TRUE)), cor_max = suppressWarnings(max(cor_with_pooled, na.rm = TRUE))), by = region_id][order(cor_med)]
print(f3(rbind(head(PR, 10), tail(PR, 5))))
cat(sprintf("\n■ 13-2에서 체중 층 LD가 pooled와 달랐던 위치 %d곳 (%s): 품종별 상관 (낮은 순 8개)\n", length(wt_regs), paste(wt_regs, collapse = ", ")))
for (r in intersect(wt_regs, R$region_id)) {
  x <- LB[is_index == TRUE & region_id == r & !is.na(cor_with_pooled)][order(cor_with_pooled)][1:min(8, .N)]
  cat(sprintf("  %s: %s\n", r, paste(sprintf("%s (%s, n=%d) %.2f", x$group, x$size_grp, x$n_dogs, x$cor_with_pooled), collapse = " · ")))
}
cat(sprintf("\n■ 대표 SNP 짝 %d개의 13-1 구분이 품종에서 바뀌는 비율 (LD 없음 제외; 작은 표본일수록 r²가 흔들림)\n", nrow(lp)))
print(f3(LPB[group %in% B$breed, .(n_dogs = n_dogs[1], n_pairs = sum(signal_g != "LD 없음"), n_changed = sum(changed), pct = 100 * mean(changed[signal_g != "LD 없음"])),
             by = group][order(-n_dogs)]), nrows = 100)
cat("참고: 같은 크기 무작위 표본 (크기마다 3번 합침)\n")
print(f3(LPB[grepl("^random_n", group), .(n_pairs = sum(signal_g != "LD 없음"), n_changed = sum(changed), pct = 100 * mean(changed[signal_g != "LD 없음"])),
             by = n_dogs][order(n_dogs)]))
cat(sprintf("SNP 짝 단위 (⑨ c; 같은 SNP 짝이 여러 모형 짝에 나오면 한 번, 같은 SNP끼리인 자기 짝 %d개 제외): %d짝\n",
            nrow(lp[lead_a == lead_b]), uniqueN(LPU[group == "pooled", .(s1, s2)])))
print(f3(LPU[group %in% B$breed, .(n_dogs = n_dogs[1], n_pairs = sum(signal_g != "LD 없음"), n_changed = sum(changed), pct = 100 * mean(changed[signal_g != "LD 없음"])),
             by = group][order(-n_dogs)]), nrows = 100)
cat("참고: 같은 크기 무작위 표본 (SNP 짝 단위)\n")
print(f3(LPU[grepl("^random_n", group), .(n_pairs = sum(signal_g != "LD 없음"), n_changed = sum(changed), pct = 100 * mean(changed[signal_g != "LD 없음"])),
             by = n_dogs][order(n_dogs)]))
cat(sprintf("\n■ 모형 × 기준 경향 (⑦; 13-1 locus %d개) — locus마다 품종 값의 중앙값을 다시 모형 × 기준 중앙값으로 (서술)\n", nrow(LCB)))
cat("  med_maf_pooled = 대표 SNP 전체 표본 MAF · med_pct_breed_mono = 대표 SNP가 단형성인 품종 % · med_breed_cor = 품종 r² 상관 (pooled와) ·\n")
cat("  med_pooled_far2 / med_breed_far2 = 강한 LD 퍼진 거리 (2번째, kb) 전체 표본 / 품종 안 · med_breed_n08 = 품종 안 r² ≥ 0.8 SNP 수\n")
for (sr in c("나이", "체중", "기타")) {
  cat(sprintf("[%s 계열]%s\n", sr, switch(sr, 나이 = " young8 (8세 미만) → main (전체) → 6 · 7 · 8 · 9 · 10 · 11 · 12세 이상 순", 체중 = " wt_t1 (가벼움) → wt_t3 (무거움) 순", "")))
  print(f3(LBM[series == sr, !"series"]), nrows = 100)
}
cat(sprintf("\n■ 품종 집중 (⑧; 순종형 45개 품종; 가장 많이 가진 품종 비율 ≥ %.1f · 기대 대비 ≥ %g배) — 서술\n", conc_share, conc_enrich))
cat("  공유 여부와 특이성 분류는 13-3 ⑨ (공유 안 됨 = 혼자 또는 다른 신호만). top_breeds = 품종 집중 locus에서 가장 많이 나온 품종 (최대 3개)\n")
for (sr in c("나이", "체중", "기타")) { cat(sprintf("[%s 계열]\n", sr)); print(f3(BCM[series == sr, !"series"]), nrows = 200) }
cat("기준 × 공유 여부 전체:\n")
print(f3(LCB[, .(n = .N, pct_breed_conc = 100 * mean(breed_conc, na.rm = TRUE), med_top_share = median(top_share, na.rm = TRUE)), by = .(thr, shared)][order(-thr, shared)]))
cat("13-3 특이성 분류별 품종 집중 비율 (전체 locus):\n")
print(f3(LCB[, .(n = .N, pct_breed_conc = 100 * mean(breed_conc, na.rm = TRUE)), by = .(specass)][order(-n)]))
cat("품종 집중 locus가 가장 많이 몰린 품종 (상위 10):\n"); print(head(LCB[breed_conc %in% TRUE, .N, by = top_breed][order(-N)], 10))
uc <- unique(CC, by = "SNP")[SNP %in% cmp$SNP]
cat(sprintf("순종형 소수 대립유전자 총 개수 (pure_mac): 최소 %d · 중앙값 %g · 하한 %d 미만이라 품종 집중에서 빠진 대표 SNP %d개 (⑨ b)\n",
            min(uc$pure_mac), median(uc$pure_mac), conc_min_mac, sum(uc$pure_mac < conc_min_mac & uc$top_share >= conc_share & uc$enrich >= conc_enrich, na.rm = TRUE)))
cat("\n■ 품종 집중 locus의 그 품종 안 LD (⑨ d; 대표 SNP 단위, 서술): top_breed_far2_kb = 가장 많이 가진 품종 안 강한 LD 퍼진 거리 (2번째, kb)\n")
ul <- unique(LCB, by = "SNP")
print(f3(ul[, .(n = .N, med_top_breed_far2_kb = as.numeric(median(top_breed_far2_kb[top_breed_mono %in% FALSE], na.rm = TRUE)),
                med_breed_far2_all = median(breed_far2_med, na.rm = TRUE), med_top_breed_n08 = as.numeric(median(top_breed_n08[top_breed_mono %in% FALSE], na.rm = TRUE)),
                med_top_breed_cor = median(top_breed_cor, na.rm = TRUE), pct_far2_ge_500kb = 100 * mean(top_breed_far2_kb[top_breed_mono %in% FALSE] >= 500, na.rm = TRUE)),
            by = .(breed_conc)]))
cat("  품종 집중이면서 그 품종 안 퍼진 거리 ≥ 500 kb인 대표 SNP (품종 병목 하플로타입 덩어리 후보):\n")
print(f3(ul[breed_conc %in% TRUE & top_breed_far2_kb >= 500, .(SNP, region_id, top_breed, top_share, enrich, top_breed_far2_kb, top_breed_n08,
                                                              models = sapply(SNP, function(x) paste(unique(LCB[SNP == x, model]), collapse = ",")))][order(-top_breed_far2_kb)]), nrows = 100)
cat("\n■ 나이 · 체중 계열 locus: 13-3 특이성 분류 × 품종 집중 (⑨ d; 기준 × 계열 × 공유 여부, 서술)\n")
print(f3(LCB[series %in% c("나이", "체중"), .(n = .N, n_breed_conc = sum(breed_conc, na.rm = TRUE), pct_breed_conc = 100 * mean(breed_conc, na.rm = TRUE)),
             by = .(thr, series, shared, specass)][order(-thr, series, shared, -n)]), nrows = 200)
cat(sprintf("\n■ 13-3 SuSiE 신뢰 집합 (CS) %d개의 품종별 확인 가능성 (⑨ d; CS SNP 중 하나라도 품종 안 MAF ≥ %.2f = 확인 가능, 서술)\n", nrow(CSS), cs_maf))
cat("  n_usable = 확인 가능한 품종 수 (45개 중) · med_frac_mono = 품종별 CS SNP 단형성 비율의 중앙값 · usable_breeds = PIP 가중 MAF가 큰 품종 (최대 5)\n")
print(f3(CSS), nrows = 200)
cat("품종 크기 구간별 확인 가능 비율 (CS × 품종):\n")
print(f3(CSB[, .(n = .N, pct_usable = 100 * mean(usable), med_frac_mono = median(frac_mono)), by = .(n_bin = nbin(n_dogs))][order(n_bin)]))
sink()
unlink(wd, recursive = TRUE)
cat(sprintf("저장: %s_* (breeds · ld_breed · ld_random · lead_ld_breed · loci_breed · loci_breed_by_model · breed_conc_by_model · cs_breed(_summary) · summary) / 그림 13-4%s_breed_cor.png · _breed_heatmap.png · _breed_by_model.png\n", out_prefix, tag))

