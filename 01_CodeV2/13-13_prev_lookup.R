# 13-13_prev_lookup.R - Looks up the loci of an earlier analysis in the present data. See README.md.
#   03-3_maf0.01_hwe.fam · 04-2_pheno.txt · 03-5_qcovar_PC5.txt · 03-5_covar.txt · 04-1_grm/ · 05-1_inputs/purelike_* · 04-5_gwas/rint_maf0.01_PC5.mlma.gz ·

suppressMessages({ library(data.table); library(parallel) })
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
out_prefix <- file.path(qc_dir, "13-13")
pre <- function(x) file.path(qc_dir, x)
inp <- function(x) file.path(d0, x)
n_cores <- as.integer(Sys.getenv("N_CORES", "12"))
plink2 <- "/data/00_tools/plink2_linux_avx2_20231123/plink2"
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
g22 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure"
rel <- "/data/02_ByProject/DogAgingProject/GeneticData/v2023/DogAgingProject_GeneticData_CuratedRelease_2023_ed"
f_in <- c(bed = paste0(rel, ".bed"), vmiss = inp("01-1_raw.vmiss"), afreq = inp("01-1_raw.afreq"), dup = inp("01-1_dup_pos_exclude.txt"), g05 = inp("01-2_geno0.05.pvar"),
          fam = inp("03-3_maf0.01_hwe.fam"), pheno = inp("04-2_pheno.txt"), qcovar = inp("03-5_qcovar_PC5.txt"), covar = inp("03-5_covar.txt"),
          grm = inp("04-1_grm/grm_all.grm.bin"), pl_keep = inp("05-1_inputs/purelike_keep.txt"), pl_pheno = inp("05-1_inputs/purelike_pheno.txt"),
          pl_qcovar = inp("05-1_inputs/purelike_qcovar.txt"), pl_covar = inp("05-1_inputs/purelike_covar.txt"),
          main = inp("04-5_gwas/rint_maf0.01_PC5.mlma.gz"), pl_gwas = inp("05-2_gwas/purelike.mlma.gz"), genes = inp("06-3_magma/main.genes.out"),
          gl = file.path(g22, "06-1_gene_loc.txt"), orth = pre("06-4_ortholog_map.csv"), prev_fam = file.path(g22, "00_purebred_data.fam"))
for (f in c(f_in, plink2, gcta)) if (!file.exists(f)) stop("파일 없음: ", f)
grm_dir <- inp("04-1_grm")
wd <- tempfile("r13_13_"); dir.create(wd)
p_cut <- 0.05; win_kb <- 50; gene_kb <- 10
f3 <- function(x) { x <- copy(x); for (j in names(x)) if (is.numeric(x[[j]]) && !is.integer(x[[j]])) set(x, j = j, value = signif(x[[j]], 3)); x }
t0 <- Sys.time()

PREV <- data.table(
  SNP = c("1:112125573:C:A", "1:112125589:A:G", "7:79899038:G:A", "11:38208882:G:A", "11:38211005:G:A", "22:57833607:G:A",
          "5:25529184:C:T", "24:12787767:A:G", "24:12787779:T:C", "26:13787793:A:G", "27:44543118:A:T"),
  prev_table = c(rep("Table 1 (p < 5e-8)", 6), rep("Table 2 (p < 1e-4)", 5)),
  prev_gene = c("ENSCAFG00000004882", "ENSCAFG00000004882", "ENSCAFG00000019076", "ENSCAFG00000047955", "ENSCAFG00000047955", "-",
                "GUCY1A2", "PLCB4", "PLCB4", "NOS1", "CACNA1C"),
  prev_A1 = c("A", "G", "A", "A", "A", "A", "T", "A", "T", "G", "T"),
  prev_beta = c(-3.6554, -3.6400, 6.1008, 4.7163, 4.6566, -5.8488, -1.4793, -0.9766, -0.9437, -3.6821, -1.5759),
  prev_p = c(NA, NA, NA, NA, NA, NA, 6.73e-5, 2.58e-5, 7.20e-5, 4.58e-5, 6.99e-5))
PREV[, `:=`(chr = as.integer(sub(":.*", "", SNP)), bp = as.integer(sub("^[0-9]+:([0-9]+):.*", "\\1", SNP)))]
p_bonf <- p_cut / nrow(PREV)
f_snps <- file.path(wd, "snps.txt"); writeLines(PREV$SNP, f_snps)

vm <- fread(f_in[["vmiss"]])[ID %in% PREV$SNP, .(SNP = ID, raw_fmiss = F_MISS)]
af <- fread(f_in[["afreq"]])[ID %in% PREV$SNP, .(SNP = ID, raw_alt_freq = ALT_FREQS)]
dup <- readLines(f_in[["dup"]])
in05 <- suppressWarnings(system2("grep", c("-cwF", "-f", shQuote(f_snps), shQuote(f_in[["g05"]])), stdout = TRUE))
X <- merge(merge(PREV, vm, by = "SNP", all.x = TRUE), af, by = "SNP", all.x = TRUE)
X[, qc_removed := fifelse(SNP %in% dup, "01-1 같은 위치 중복 (다중 대립)", fifelse(is.na(raw_fmiss), "원 자료에 없음",
                   fifelse(raw_fmiss > 0.05, sprintf("01-2 SNP 결측 > 5%% (%.0f%%)", 100 * raw_fmiss), "01-2 뒤 단계")))]
cat(sprintf("[%s] ① QC: 11개 중 01-2 통과 %s개 (0이어야 01-2에서 모두 빠진 것)\n", format(Sys.time(), "%T"), in05))

keep_all <- file.path(wd, "keep_all.txt"); fwrite(fread(f_in[["fam"]], header = FALSE)[, 1:2], keep_all, sep = " ", col.names = FALSE)
sub_bf <- file.path(wd, "prev11")
st <- system2(plink2, c("--dog", "--bfile", rel, "--extract", f_snps, "--keep", keep_all, "--make-bed", "--out", sub_bf, "--threads", "4"), stdout = FALSE, stderr = FALSE)
if (st != 0 || !file.exists(paste0(sub_bf, ".bed"))) stop("plink2 추출 실패")
bim <- fread(paste0(sub_bf, ".bim"), header = FALSE, col.names = c("chr", "SNP", "cm", "bp", "a1", "a2"))
cat(sprintf("[%s] ② 추출: SNP %d개 · 개 %d마리\n", format(Sys.time(), "%T"), nrow(bim), nrow(fread(paste0(sub_bf, ".fam"), header = FALSE))))
fam_all <- fread(f_in[["fam"]], header = FALSE)[, 1:2]; prev_ids <- fread(f_in[["prev_fam"]], header = FALSE)$V2
keep_new <- file.path(wd, "keep_new.txt"); keep_prev <- file.path(wd, "keep_prev.txt")
fwrite(fam_all[!V2 %in% prev_ids], keep_new, sep = " ", col.names = FALSE); fwrite(fam_all[V2 %in% prev_ids], keep_prev, sep = " ", col.names = FALSE)
n_new <- fam_all[!V2 %in% prev_ids, .N]; n_prev <- fam_all[V2 %in% prev_ids, .N]
cat(sprintf("[%s] ⑥ 이전 원고 표본에 없는 개 %d마리 · 있던 개 %d마리\n", format(Sys.time(), "%T"), n_new, n_prev))
runs <- data.table(run = c("main_rint", "main_raw", "purelike_rint", "new_rint", "prev_rint"),
                   pheno = c(f_in[["pheno"]], f_in[["pheno"]], f_in[["pl_pheno"]], f_in[["pheno"]], f_in[["pheno"]]), mpheno = c(1L, 2L, 1L, 1L, 1L),
                   qcovar = c(f_in[["qcovar"]], f_in[["qcovar"]], f_in[["pl_qcovar"]], f_in[["qcovar"]], f_in[["qcovar"]]),
                   covar = c(f_in[["covar"]], f_in[["covar"]], f_in[["pl_covar"]], f_in[["covar"]], f_in[["covar"]]),
                   keep = c("", "", f_in[["pl_keep"]], keep_new, keep_prev))
jobs <- CJ(run = runs$run, chr = sort(unique(bim$chr)))[runs, on = "run"]
jobs[, out := file.path(wd, sprintf("%s_chr%d", run, chr))]
st <- unlist(mclapply(seq_len(nrow(jobs)), function(i) { j <- jobs[i]
  a <- c("--bfile", sub_bf, "--autosome-num", "38", "--chr", j$chr, "--grm", file.path(grm_dir, "grm_all"), "--mlma-subtract-grm", file.path(grm_dir, paste0("grm_chr", j$chr)), "--mlma",
         "--pheno", j$pheno, "--mpheno", j$mpheno, "--qcovar", j$qcovar, "--covar", j$covar, "--thread-num", "2", "--out", j$out)
  if (nzchar(j$keep)) a <- c(a, "--keep", j$keep)
  system2(gcta, a, stdout = FALSE, stderr = FALSE) }, mc.cores = max(1, n_cores %/% 2), mc.preschedule = FALSE))
bad <- jobs[st != 0 | !file.exists(paste0(out, ".mlma"))]
if (nrow(bad)) stop("GCTA 실패: ", paste(basename(bad$out), collapse = ", "), " (로그: ", wd, ")")
D <- rbindlist(lapply(seq_len(nrow(jobs)), function(i) fread(paste0(jobs$out[i], ".mlma"))[, run := jobs$run[i]]))
D <- merge(D[, .(run, SNP, A1, A2, Freq, b, se, p)], PREV[, .(SNP, prev_A1, prev_beta)], by = "SNP")
D[, b_prev_allele := fifelse(A1 == prev_A1, b, fifelse(A2 == prev_A1, -b, NA_real_))]
D[, same_dir := sign(b_prev_allele) == sign(prev_beta)]
cat(sprintf("[%s] ② GCTA %d번 끝\n", format(Sys.time(), "%T"), nrow(jobs)))

system2(plink2, c("--dog", "--bfile", sub_bf, "--export", "A", "--out", file.path(wd, "dos")), stdout = FALSE, stderr = FALSE)
G <- fread(file.path(wd, "dos.raw"))
ph <- fread(f_in[["pheno"]], header = FALSE, col.names = c("FID", "IID", "rint", "raw"))
qc <- fread(f_in[["qcovar"]], header = FALSE); setnames(qc, c("FID", "IID", paste0("q", seq_len(ncol(qc) - 2))))
cv <- fread(f_in[["covar"]], header = FALSE); setnames(cv, c("FID", "IID", paste0("c", seq_len(ncol(cv) - 2))))
Z <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")), list(G[, !c("PAT", "MAT", "SEX", "PHENOTYPE")], ph, qc, cv))
cov_terms <- c(grep("^q", names(Z), value = TRUE), paste0("factor(", grep("^c", names(Z), value = TRUE), ")"))
gcols <- grep("^[0-9]+:", names(Z), value = TRUE)
CCR <- rbindlist(lapply(gcols, function(gc) {
  snp <- sub("_[ACGT]+$", "", gc); cnt <- sub(".*_", "", gc)
  d <- Z[, c(list(g = get(gc), rint = rint, age = q1), mget(sub("factor\\((.*)\\)", "\\1", cov_terms)))]
  n_obs <- sum(!is.na(d$g))
  fit <- tryCatch(summary(lm(reformulate(c("g", cov_terms), "rint"), data = d))$coefficients["g", ], error = function(e) rep(NA_real_, 4))
  d[, miss := as.integer(is.na(g))]
  mf <- tryCatch(summary(glm(miss ~ scale(rint) + scale(age), family = binomial, data = d))$coefficients["scale(rint)", ], error = function(e) rep(NA_real_, 4))
  data.table(SNP = snp, counted = cnt, n_obs = n_obs, cc_b = fit[1], cc_se = fit[2], cc_p = fit[4], miss_logOR_per_sd = mf[1], miss_p = mf[4])
}))
CCR <- merge(CCR, PREV[, .(SNP, prev_A1, prev_beta)], by = "SNP")
CCR[, cc_b_prev_allele := fifelse(counted == prev_A1, cc_b, -cc_b)][, cc_same_dir := sign(cc_b_prev_allele) == sign(prev_beta)]

read_win <- function(f) {
  aw <- paste(sprintf("($1 == %d && $3 >= %d && $3 <= %d)", PREV$chr, PREV$bp - win_kb * 1e3, PREV$bp + win_kb * 1e3), collapse = " || ")
  fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR == 1 || %s'", shQuote(f), aw))
}
RG <- rbindlist(lapply(c(main = "main", purelike = "pl_gwas"), function(k) {
  W <- read_win(f_in[[k]])[!is.na(p)]
  rbindlist(lapply(seq_len(nrow(PREV)), function(i) { w <- W[Chr == PREV$chr[i] & abs(bp - PREV$bp[i]) <= win_kb * 1e3]
    if (!nrow(w)) return(data.table(SNP = PREV$SNP[i], gwas = k, n_win = 0L))
    m <- w[which.min(p)]
    data.table(SNP = PREV$SNP[i], gwas = k, n_win = nrow(w), min_p = m$p, min_snp = m$SNP, min_dist_kb = (m$bp - PREV$bp[i]) / 1e3, local_bonf = p_cut / nrow(w)) }))
}), fill = TRUE)
RG[, gwas := fifelse(gwas == "pl_gwas", "purelike", gwas)]
RG[, region_signal := !is.na(min_p) & min_p < local_bonf]

gl <- fread(f_in[["gl"]], header = FALSE, col.names = c("GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"))[CHR %in% as.character(1:38)][, CHR := as.integer(CHR)]
orth <- unique(fread(f_in[["orth"]])[, .(GENE = dog, hs_symbol)], by = "GENE")
gr <- fread(f_in[["genes"]])[, .(GENE, NSNPS, ZSTAT, P)]
GN <- rbindlist(lapply(seq_len(nrow(PREV)), function(i) { g <- gl[CHR == PREV$chr[i]]
  ov <- g[STOP >= PREV$bp[i] - gene_kb * 1e3 & START <= PREV$bp[i] + gene_kb * 1e3]
  if (!nrow(ov)) { g[, dist := pmax(START - PREV$bp[i], PREV$bp[i] - STOP, 0)]; ov <- g[which.min(dist)]; how <- "가장 가까운" } else { ov[, dist := 0]; how <- "± 10 kb 겹침" }
  ov[, .(SNP = PREV$SNP[i], GENE, SYMBOL, map = how, dist_kb = dist / 1e3)] }))
GN <- merge(GN, orth, by = "GENE", all.x = TRUE)
GN <- merge(GN, gr, by = "GENE", all.x = TRUE)
GN[, label := fifelse(!is.na(hs_symbol) & hs_symbol != "", hs_symbol, SYMBOL)]
n_gene <- uniqueN(GN[!is.na(P), GENE])
GN[, gene_bonf := p_cut / n_gene]

dm <- D[run == "main_rint", .(SNP, v23_A1 = A1, v23_freq = Freq, v23_b = b, v23_se = se, v23_p = p, v23_b_prev_allele = b_prev_allele, v23_same_dir = same_dir)]
X <- Reduce(function(a, b) merge(a, b, by = "SNP", all.x = TRUE), list(X, dm,
       D[run == "main_raw", .(SNP, raw_p = p, raw_same_dir = same_dir)], D[run == "purelike_rint", .(SNP, purelike_p = p, purelike_same_dir = same_dir)],
       D[run == "new_rint", .(SNP, new_b_prev_allele = b_prev_allele, new_p = p, new_same_dir = same_dir)], D[run == "prev_rint", .(SNP, prevdogs_p = p, prevdogs_same_dir = same_dir)],
       CCR[, .(SNP, n_obs, cc_p, cc_same_dir, miss_logOR_per_sd, miss_p)],
       dcast(RG, SNP ~ gwas, value.var = c("min_p", "local_bonf", "region_signal")),
       GN[, .(genes = paste(sprintf("%s (p %.2g)", label, P), collapse = "; "), gene_min_p = suppressWarnings(min(P, na.rm = TRUE))), by = SNP]))
X[, verdict := fifelse(is.na(v23_p), "검정 불가", fifelse(v23_same_dir & v23_p < p_bonf, "재현", fifelse(v23_same_dir & v23_p < p_cut, "명목 같은 방향", "재현 안 됨")))]
X[, verdict_new := fifelse(is.na(new_p), "검정 불가", fifelse(new_same_dir & new_p < p_bonf, "독립 재현", fifelse(new_same_dir & new_p < p_cut, "명목 같은 방향", "재현 안 됨")))]
X[, SNP := factor(SNP, PREV$SNP)]; setorder(X, SNP)
fwrite(f3(X), paste0(out_prefix, "_snps.csv")); fwrite(f3(rbind(D, fill = TRUE)), paste0(out_prefix, "_direct.csv"))
fwrite(f3(RG), paste0(out_prefix, "_region.csv")); fwrite(f3(GN), paste0(out_prefix, "_genes.csv"))

sink(paste0(out_prefix, "_summary.txt"))
cat("13-13 이전 원고 (DAP 2022 순종 3,527마리, PLINK glm) locus의 v2023 재확인 (판정 규칙은 헤더; 같은 코호트 안 재분석, 독립 재현 아님)\n")
cat(sprintf("직접 검정 판정: 같은 방향 & p < %.4f (0.05 / %d) = 재현 · 같은 방향 & p < 0.05 = 명목 같은 방향\n", p_bonf, nrow(PREV)))
cat("\n■ ① QC · ② 직접 검정 (main RINT, GCTA MLMA-LOCO; 결측은 평균 대치) · 참고 열\n")
print(f3(X[, .(SNP, prev_table, prev_gene, prev_beta, raw_fmiss, raw_alt_freq, qc_removed, n_obs, v23_b_prev_allele, v23_p, verdict, raw_p, purelike_p, cc_p, cc_same_dir)]), row.names = FALSE)
cat(sprintf("\n■ ⑥ 독립 표본: 이전 원고 표본에 없는 v2023 개 %d마리 (new) · 있던 개 %d마리 (prevdogs, 참고)\n", n_new, n_prev))
print(f3(X[, .(SNP, prev_gene, v23_p, verdict, new_b_prev_allele, new_p, verdict_new, prevdogs_p, prevdogs_same_dir)]), row.names = FALSE)
cat("\n■ ③ 결측 여부 ~ RINT + 나이 (로지스틱; 표현형 1 SD당 log OR)\n")
print(f3(X[, .(SNP, raw_fmiss, miss_logOR_per_sd, miss_p)]), row.names = FALSE)
cat(sprintf("\n■ ④ 주변 ± %d kb 최소 p (QC 통과 SNP; 구간 Bonferroni)\n", win_kb)); print(f3(RG), row.names = FALSE)
cat(sprintf("\n■ ⑤ 유전자 수준 (v2023 main MAGMA; 유전자 %d개, Bonferroni %.4f)\n", n_gene, p_cut / n_gene)); print(f3(GN[, .(SNP, label, map, dist_kb, NSNPS, ZSTAT, P)]), row.names = FALSE)
cat("\n■ 판정 집계 (② 전체 · ⑥ 독립)\n"); print(X[, .N, by = .(prev_table, verdict, verdict_new)], row.names = FALSE)
sink()
unlink(wd, recursive = TRUE)
cat(sprintf("[%s] 끝 (%.1f분). 저장: %s_* (snps · direct · region · genes · summary)\n", format(Sys.time(), "%T"), as.numeric(difftime(Sys.time(), t0, units = "mins")), out_prefix))
print(f3(X[, .(SNP, prev_gene, raw_fmiss, n_obs, v23_p, verdict, new_p, verdict_new, prevdogs_p, min_p_main, gene_min_p)]), row.names = FALSE)

