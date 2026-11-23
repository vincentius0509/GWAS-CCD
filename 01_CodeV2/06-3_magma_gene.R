# 06-3_magma_gene.R - MAGMA gene-level analysis for each model, permutation null and gene-level heterogeneity. See README.md.
#        p < 1e-3 = "suggestive"
#        "Given the limited power (no genome-wide significant locus), gene-level associations at FDR q < 0.10 and p < 1 × 10⁻³ were reported
#         as exploratory/suggestive and were always compared with 20 phenotype-permuted GWAS; these require replication in an independent cohort."

library(data.table)
library(ggplot2)
library(parallel)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")    
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
out_prefix <- file.path(qc_dir, "06-3")
pic_prefix <- file.path(pic_dir, "06-3")
work_dir <- paste0(out_prefix, "_work")
ref_dir <- paste0(out_prefix, "_ref")
runs_tsv <- file.path(qc_dir, "05-1_runs.tsv")
gcta_dir <- file.path(qc_dir, "05-2_gwas")
regenie_dir <- file.path(qc_dir, "05-3_regenie")
het_prefix <- file.path(qc_dir, "06-2_het_")
bfile <- file.path(qc_dir, "03-3_maf0.01_hwe")
main_in <- list(pheno = file.path(qc_dir, "04-2_pheno.txt"), qcovar = file.path(qc_dir, "03-5_qcovar_PC5.txt"),
                covar = file.path(qc_dir, "03-5_covar.txt"))
v2022_dir <- file.path(base_dir, "GeneticData/analysis/03_GWAS_Pure")
gene_loc_file <- file.path(v2022_dir, "06-1_gene_loc.txt")                      # GeneID CHR START STOP STRAND SYMBOL (NCBI CanFam3.1)
gene_info_file <- file.path(v2022_dir, "06-1_gene_info.csv")
ortho_file <- file.path(v2022_dir, "06-1_ref", "gene_orthologs_human_dog.tsv")
dog_ft_file <- file.path(v2022_dir, "06-1_ref", "GCF_000002285.3_CanFam3.1_feature_table.txt.gz")
hs_info_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
magma <- "/data/00_tools/magma_v1.10/magma"
plink2 <- "/data/00_tools/plink2_linux_avx2_20231123/plink2"
for (f in c(runs_tsv, paste0(bfile, c(".bed", ".bim", ".fam")), unlist(main_in), gene_loc_file, gene_info_file, ortho_file, dog_ft_file, magma, plink2))
  if (!file.exists(f)) stop("파일 없음: ", f)
for (d in c(work_dir, ref_dir, pic_dir, file.path(work_dir, c("pval", "refbed")))) dir.create(d, showWarnings = FALSE, recursive = TRUE)

n_cores <- as.integer(Sys.getenv("N_CORES", "40"))
plink_threads <- as.integer(Sys.getenv("PLINK_THREADS", "8"))
n_null <- as.integer(Sys.getenv("N_NULL", "20"))
sens_all <- Sys.getenv("SENS_ALL", "1") == "1"
n_batch <- as.integer(Sys.getenv("N_BATCH", "40"))
chrs <- if (nzchar(Sys.getenv("CHRS"))) as.integer(strsplit(Sys.getenv("CHRS"), " ")[[1]]) else 1:38
chr_tag <- if (length(chrs) == 38) "" else paste0("_chr", paste(chrs, collapse = "-"))
res_dir <- paste0(out_prefix, "_magma", chr_tag)
batch_dir <- file.path(work_dir, paste0("batch", chr_tag))
for (d in c(res_dir, batch_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
run_ids_env <- if (nzchar(Sys.getenv("RUN_IDS"))) strsplit(Sys.getenv("RUN_IDS"), " ")[[1]] else NULL
win_main <- 10; win_sens <- 50                     # kb
trend_thr <- c("1e-5" = 1e-5, "1e-4" = 1e-4)
null_loci <- c("1e-5" = 13.5, "1e-4" = 105)
locus_kb <- 500
sexflag_min_freq <- 0.01

use06 <- c("main", "purelike", "mixlike", "int_age", "age8", "bin_age8", "young8", "sex_m", "sex_f", "wt_t1", "wt_t2", "wt_t3")
cmps <- list(sex = c("sex_m", "sex_f"), age = c("young8", "age8"), ancestry = c("purelike", "mixlike"), weight = c("wt_t1", "wt_t2", "wt_t3"))

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), axis.title = element_text(size = 13), axis.text = element_text(size = 11))
col_chr <- c("#2A78D6", "#898781")
p_to_lambda <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
clump_leads <- function(chr, bp, p, t) {
  h <- which(p < t); if (!length(h)) return(integer(0))
  h <- h[order(p[h])]; taken <- rep(FALSE, length(h)); out <- integer(0)
  for (i in seq_along(h)) if (!taken[i]) {
    out <- c(out, h[i]); taken <- taken | (chr[h] == chr[h[i]] & abs(bp[h] - bp[h[i]]) <= locus_kb * 1e3)
  }
  out
}
sh <- function(cmd, args, log) { st <- system2(cmd, args, stdout = log, stderr = log); if (st != 0) stop(basename(cmd), " 실패: ", log); invisible(st) }

runs <- fread(runs_tsv, na.strings = c("", "NA", "\"\""))
if (!all(use06 %in% runs$id)) stop("실행 표에 없는 06 모형: ", paste(setdiff(use06, runs$id), collapse = ", "))
bim <- fread(paste0(bfile, ".bim"), header = FALSE, select = c(1, 2, 4), col.names = c("Chr", "SNP", "bp"))
gloc <- fread(gene_loc_file, header = FALSE, col.names = c("GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"))
cat(sprintf("유전자 %s개, SNP %s개, 염색체 %s\n", format(nrow(gloc), big.mark = ","), format(nrow(bim), big.mark = ","),
            if (length(chrs) == 38) "1 ~ 38" else paste(chrs, collapse = " ")))

row_of <- function(id) { i <- which(runs$id == id)[1]; if (is.na(i)) stop("실행 표에 없음: ", id); as.list(runs[i]) }
read_snp_p <- function(id) {
  r <- row_of(id)
  if (r$engine == "existing") d <- fread(cmd = paste("zcat", shQuote(r$ref_result)), select = c("Chr", "SNP", "bp", "Freq", "p"))
  else if (r$engine == "gcta") d <- fread(cmd = paste("zcat", shQuote(file.path(gcta_dir, paste0(id, ".mlma.gz")))), select = c("Chr", "SNP", "bp", "Freq", "p"))
  else if (r$engine %in% c("regenie_bt", "regenie_int")) {
    x <- fread(cmd = paste("zcat", shQuote(file.path(regenie_dir, paste0(id, ".regenie.gz")))), select = c("CHROM", "GENPOS", "ID", "TEST", "LOG10P"))
    tsel <- grep(if (r$engine == "regenie_bt") "^ADD$" else "^ADD-INT_SNPx", unique(x$TEST), value = TRUE)[1]
    d <- x[TEST == tsel & !is.na(LOG10P), .(Chr = CHROM, SNP = ID, bp = GENPOS, Freq = NA_real_, p = 10^(-LOG10P))]
  } else stop(id, ": 엔진 ", r$engine)
  d[, p := suppressWarnings(as.numeric(p))]
  d
}
want <- function(run) is.null(run_ids_env) || run %in% run_ids_env || (grepl("^null", run) && "null" %in% run_ids_env) ||
  (grepl("^het_", run) && "het" %in% run_ids_env) || (grepl("_w50$", run) && sub("_w50$", "", run) %in% run_ids_env)

jobs <- list()
pval_file <- function(run) file.path(work_dir, "pval", paste0(sub("_w50$", "", run), ".txt"))
for (id in use06) {
  r <- row_of(id)
  keep <- if (r$n == nrow(fread(paste0(bfile, ".fam"), header = FALSE, select = 1))) NA_character_ else r$keep_file
  if (!is.na(keep) && !file.exists(keep)) stop(id, ": keep 파일 없음 ", keep)
  jobs[[id]] <- list(run = id, kind = "model", src = id, N = r$n, keep = keep, window = win_main, genes_only = FALSE)
}
for (id in if (sens_all) use06 else "main") {
  s50 <- jobs[[id]]; s50$run <- paste0(id, "_w50"); s50$kind <- "sens"; s50$window <- win_sens; jobs[[s50$run]] <- s50
}
for (cn in names(cmps)) jobs[[paste0("het_", cn)]] <- list(run = paste0("het_", cn), kind = "het", src = cn, N = sum(sapply(cmps[[cn]], function(i) row_of(i)$n)),
                                                       keep = NA_character_, window = win_main, genes_only = FALSE)
for (cn in names(cmps)) {
  h <- jobs[[paste0("het_", cn)]]; h$run <- paste0("het_", cn, "_w50"); h$kind <- "het_sens"; h$window <- win_sens; jobs[[h$run]] <- h
}
for (k in seq_len(n_null)) jobs[[paste0("null", k)]] <- list(run = paste0("null", k), kind = "null", src = k, N = row_of("main")$n, keep = NA_character_,
                                                          window = win_main, genes_only = FALSE)
jobs <- jobs[sapply(names(jobs), want)]
done <- function(run) all(file.exists(file.path(res_dir, paste0(run, c(".genes.out", ".genes.raw")))))
todo <- names(jobs)[!sapply(names(jobs), done)]
cat(sprintf("run %d개 (끝난 것 %d개는 건너뜀): %s\n", length(jobs), length(jobs) - length(todo), paste(todo, collapse = " ")))

for (run in todo) {
  j <- jobs[[run]]; f <- pval_file(run); if (file.exists(f)) next
  if (j$kind %in% c("model", "sens")) d <- read_snp_p(j$src)[is.finite(p) & p > 0, .(SNP, P = p)]
  else if (j$kind %in% c("het", "het_sens")) d <- fread(cmd = paste("zcat", shQuote(paste0(het_prefix, j$src, ".tsv.gz"))), select = c("SNP", "p_het"))[is.finite(p_het) & p_het > 0, .(SNP, P = p_het)]
  else next
  fwrite(d, f, sep = "\t")
}
null_todo <- todo[grepl("^null", todo)]; null_todo <- null_todo[!file.exists(vapply(null_todo, pval_file, ""))]
if (length(null_todo)) {
  wd <- file.path(work_dir, "null"); dir.create(wd, showWarnings = FALSE)
  y <- fread(main_in$pheno, header = FALSE, colClasses = list(character = 1:2))[, 1:3]; setnames(y, c("#FID", "IID", "y"))
  q <- fread(main_in$qcovar, header = FALSE, colClasses = list(character = 1:2)); setnames(q, c("#FID", "IID", paste0("q", seq_len(ncol(q) - 2))))
  cv <- fread(main_in$covar, header = FALSE, colClasses = "character"); setnames(cv, c("#FID", "IID", paste0("c", seq_len(ncol(cv) - 2))))
  fwrite(merge(q, cv, by = c("#FID", "IID")), file.path(wd, "covar.txt"), sep = "\t")
  set.seed(20261112); pp <- y[, 1:2]; for (k in seq_len(n_null)) pp[[paste0("null", k)]] <- sample(y$y)
  fwrite(pp, file.path(wd, "pheno.txt"), sep = "\t")
  cat(sprintf("null 순열 GWAS %d개 (plink2, 주 분석 표본)\n", n_null))
  sh(plink2, c("--dog", "--bfile", bfile, "--pheno", file.path(wd, "pheno.txt"), "--covar", file.path(wd, "covar.txt"), "--covar-variance-standardize",
               "--glm", "hide-covar", "cols=p", "--threads", n_cores, "--out", file.path(wd, "null")), file.path(wd, "null.stdout"))
  for (k in seq_len(n_null)) {
    f <- file.path(wd, sprintf("null.null%d.glm.linear", k))
    x <- fread(f); setnames(x, sub("^#", "", names(x)))                       # plink2 cols=p: #ID A1 P
    fwrite(x[is.finite(P) & P > 0, .(SNP = ID, P)], pval_file(paste0("null", k)), sep = "\t"); file.remove(f)
  }
}

keep_key <- function(keep) if (is.na(keep)) "all" else {
  k <- fread(keep, header = FALSE, colClasses = "character")[, 2][[1]]; paste0("k", substr(digest_md5(paste(sort(k), collapse = ",")), 1, 10))
}
digest_md5 <- function(x) { tf <- tempfile(); writeLines(x, tf); h <- tools::md5sum(tf); unlink(tf); unname(h) }
ref_of <- list()
for (run in todo) {
  j <- jobs[[run]]; key <- keep_key(j$keep)
  pre <- if (key == "all") bfile else file.path(work_dir, "refbed", key)
  if (key != "all" && !file.exists(paste0(pre, ".bed"))) {
    cat(sprintf("LD 기준 bfile 만들기: %s (%s)\n", run, basename(j$keep)))
    sh(plink2, c("--dog", "--bfile", bfile, "--keep", j$keep, "--make-bed", "--threads", n_cores, "--out", pre), paste0(pre, ".stdout"))
  }
  ref_of[[run]] <- pre
}

annot_full <- function(w) file.path(work_dir, sprintf("annot_w%d", w))
annot_of <- function(w) paste0(annot_full(w), chr_tag)
for (w in unique(sapply(jobs, `[[`, "window"))) {
  if (!file.exists(paste0(annot_full(w), ".genes.annot")))
    sh(magma, c("--annotate", "nonhuman", sprintf("window=%d", w), "--snp-loc", paste0(bfile, ".bim"), "--gene-loc", gene_loc_file, "--out", annot_full(w)),
       paste0(annot_full(w), ".stdout"))
  if (nzchar(chr_tag) && !file.exists(paste0(annot_of(w), ".genes.annot"))) {
    a <- readLines(paste0(annot_full(w), ".genes.annot"))
    keep_l <- grepl("^#", a) | sub(":.*", "", sapply(strsplit(a, "[\t ]+"), `[`, 2)) %in% as.character(chrs)
    writeLines(a[keep_l], paste0(annot_of(w), ".genes.annot"))
  }
}
bj <- CJ(run = todo, b = seq_len(n_batch))
cat(sprintf("MAGMA batch %d개 (run %d × batch %d), 동시 %d개\n", nrow(bj), length(todo), n_batch, n_cores))
t0 <- Sys.time()
bres <- mclapply(seq_len(nrow(bj)), function(i) {
  run <- bj$run[i]; b <- bj$b[i]; j <- jobs[[run]]
  out <- file.path(batch_dir, run)
  if (length(list.files(batch_dir, pattern = sprintf("^%s\\.batch%d_[0-9]+\\.genes\\.out$", run, b)))) return("있음")
  args <- c("--bfile", ref_of[[run]], "--pval", pval_file(run), sprintf("N=%d", j$N), "--gene-annot", paste0(annot_of(j$window), ".genes.annot"),
            "--batch", b, n_batch, if (j$genes_only) "--genes-only", "--out", out)
  st <- system2(magma, args, stdout = FALSE, stderr = FALSE)
  if (st == 0) return("완료")
  lg <- sprintf("%s.batch%d_%d.log", out, b, n_batch)
  if (file.exists(lg) && any(grepl("is empty", readLines(lg)))) "빔" else paste("실패", run, b)
}, mc.cores = n_cores, mc.preschedule = FALSE)
bad <- unlist(bres)[!unlist(bres) %in% c("완료", "있음", "빔")]
if (length(bad)) stop("MAGMA batch 실패: ", paste(head(bad, 10), collapse = "; "))
cat(sprintf("MAGMA batch 끝 (%.1f분)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
for (run in todo) {
  out <- file.path(batch_dir, run)
  bf <- list.files(batch_dir, pattern = paste0("^", run, "\\.batch[0-9]+_[0-9]+\\.genes\\.out$"), full.names = TRUE)
  if (length(bf) == 1) {
    for (ext in c(".genes.out", ".genes.raw")) { f <- sub("\\.genes\\.out$", ext, bf); if (file.exists(f)) file.copy(f, file.path(res_dir, paste0(run, ext)), overwrite = TRUE) }
  } else sh(magma, c("--merge", out, "--out", file.path(res_dir, run)), file.path(batch_dir, paste0(run, ".merge.stdout")))
  file.remove(list.files(batch_dir, pattern = paste0("^", run, "\\.batch[0-9]+_[0-9]+\\."), full.names = TRUE))
}
used_by <- split(names(ref_of), unlist(ref_of))
for (pre in names(used_by)) if (pre != bfile && all(sapply(used_by[[pre]], done))) file.remove(paste0(pre, c(".bed", ".bim", ".fam")))

sx <- lapply(c("sex_m", "sex_f"), function(id) fread(cmd = paste("zcat", shQuote(file.path(gcta_dir, paste0(id, ".mlma.gz")))), select = c("SNP", "Freq", "b", "se")))
sxm <- merge(sx[[1]], sx[[2]], by = "SNP", suffixes = c("_m", "_f"))
sxm[, `:=`(mono_m = !(is.finite(b_m) & is.finite(se_m) & se_m > 0), mono_f = !(is.finite(b_f) & is.finite(se_f) & se_f > 0))]
maf <- function(f) pmin(f, 1 - f)
sexflag <- sxm[(mono_f & !mono_m & maf(Freq_m) >= sexflag_min_freq) | (mono_m & !mono_f & maf(Freq_f) >= sexflag_min_freq),
               .(SNP, freq_m = Freq_m, freq_f = Freq_f, mono = fifelse(mono_f, "female", "male"))]
sexflag <- merge(bim, sexflag, by = "SNP")[order(Chr, bp)]
rm(sx, sxm)
fwrite(sexflag, paste0(out_prefix, "_sexlinked_snps.txt"), sep = "\t")
ann <- readLines(paste0(annot_of(win_main), ".genes.annot")); ann <- ann[!grepl("^#", ann)]
ann_l <- strsplit(ann, "[\t ]+")
gene_flag <- rbindlist(lapply(ann_l, function(x) { s <- x[-(1:2)]; k <- sum(s %in% sexflag$SNP); if (k) data.table(GENE = as.integer(x[1]), n_snp = length(s), n_sexflag = k) }))
if (!nrow(gene_flag)) gene_flag <- data.table(GENE = integer(0), n_snp = integer(0), n_sexflag = integer(0))
hs_file <- file.path(ref_dir, "Homo_sapiens.gene_info.gz")
if (!file.exists(hs_file)) { message("사람 gene_info 내려받기: ", hs_info_url); download.file(hs_info_url, hs_file, mode = "wb", quiet = TRUE) }
hs <- fread(cmd = paste("zcat", shQuote(hs_file)), select = c("GeneID", "Symbol", "chromosome"), colClasses = "character")
setnames(hs, c("hs_gene", "hs_symbol", "hs_chr")); hs[, hs_gene := as.integer(hs_gene)]
orth <- fread(ortho_file, select = c(2, 4, 5), col.names = c("hs_gene", "other_tax", "GENE"))[other_tax == 9615]
orth <- merge(orth[, .(GENE, hs_gene)], hs, by = "hs_gene")
orth <- orth[, .(hs_symbol = paste(unique(hs_symbol), collapse = ";"), hs_chr = paste(unique(hs_chr), collapse = ";"),
                 hs_on_xy = any(grepl("^(X|Y)$|X\\|Y", hs_chr))), by = GENE]
dogx <- fread(cmd = paste("zcat", shQuote(dog_ft_file)), sep = "\t", quote = "", select = c(1, 6, 15), col.names = c("feature", "chromosome", "symbol"))[feature == "gene" & chromosome %in% c("X", "Y") & nzchar(symbol)]
ginfo <- fread(gene_info_file, select = c("GeneID", "SYMBOL", "NAME"))[, .(GENE = GeneID, SYMBOL, NAME)]
bg <- merge(gloc[, .(GENE, CHR, START, STOP)], orth, by = "GENE", all.x = TRUE)
bg <- merge(bg, ginfo, by = "GENE", all.x = TRUE)
bg[, `:=`(hs_on_xy = !is.na(hs_on_xy) & hs_on_xy, dog_chrx_same_symbol = !is.na(SYMBOL) & SYMBOL %in% dogx$symbol)]
bg <- merge(bg, gene_flag, by = "GENE", all.x = TRUE); bg[is.na(n_sexflag), n_sexflag := 0L]
fl <- bg[n_sexflag > 0][order(CHR, START)]
sf_near <- sexflag[, { x <- gloc[CHR == Chr[1]]
  dd <- rbindlist(lapply(bp, function(p) { e <- pmax(0, x$START - p, p - x$STOP); data.table(sym = x$SYMBOL[which.min(e)], kb = min(e) / 1e3) }))
  dd <- dd[, .(kb = min(kb)), by = sym][order(kb)]
  .(n = .N, from_Mb = round(min(bp) / 1e6, 2), to_Mb = round(max(bp) / 1e6, 2), freq_other_median = round(median(fifelse(mono == "female", freq_m, freq_f)), 3),
    nearest_genes = paste(head(sprintf("%s(%.0f kb)", dd$sym, dd$kb), 3), collapse = " ")) }, by = Chr][order(-n)]

read_genes <- function(run) {
  f <- file.path(res_dir, paste0(run, ".genes.out")); if (!file.exists(f)) return(NULL)
  g <- fread(f, skip = "GENE")[, .(GENE, CHR, START, STOP, NSNPS, N, ZSTAT, P)]
  g[, `:=`(run = run, Q = p.adjust(P, "BH"))]
}
all_runs <- names(jobs)[sapply(names(jobs), done)]
gr <- rbindlist(lapply(all_runs, read_genes))
gr <- merge(gr, bg[, .(GENE, SYMBOL, n_sexflag, hs_on_xy)], by = "GENE", all.x = TRUE)
kinds <- sapply(jobs, `[[`, "kind")
gr[, kind := kinds[run]]
fwrite(gr[kind != "null"][order(run, P)], paste0(out_prefix, "_gene_results.csv.gz"))
sm <- gr[, .(n_genes = .N, lambda = round(p_to_lambda(P), 4), p_1e3 = sum(P < 1e-3), p_1e4 = sum(P < 1e-4), bonf = sum(P < 0.05 / .N),
             fdr05 = sum(Q < 0.05), fdr10 = sum(Q < 0.10), top5 = paste(head(fifelse(is.na(SYMBOL), as.character(GENE), SYMBOL)[order(P)], 5), collapse = ","),
             top_flag = sum(head(n_sexflag[order(P)], 50) > 0)), by = .(run, kind)]
sm[, run := factor(run, levels = all_runs)]; setorder(sm, run)
nl <- sm[kind == "null", .(lambda = mean(lambda), p_1e3 = mean(p_1e3), p_1e3_min = min(p_1e3), p_1e3_max = max(p_1e3), p_1e4 = mean(p_1e4),
                           bonf = mean(bonf), fdr05 = mean(fdr05), fdr10 = mean(fdr10), any_fdr10 = mean(fdr10 > 0), n = .N)]
rx <- gr[kind != "null" & (P < 1e-3 | Q < 0.10), .(run, kind, GENE, SYMBOL, CHR, START, STOP, NSNPS, ZSTAT, P = signif(P, 3), Q = signif(Q, 3),
                                                  q_lt_0.10 = Q < 0.10, n_sexflag)][order(run, P)]
if (nrow(sm[kind == "null"])) rx[, `:=`(null_mean_p1e3 = round(nl$p_1e3, 1), null_frac_runs_any_q10 = round(nl$any_fdr10, 2))]
fwrite(rx, paste0(out_prefix, "_gene_relaxed.csv"), bom = TRUE)
fwrite(sm, paste0(out_prefix, "_magma_summary.csv"), bom = TRUE)
ws <- rbindlist(lapply(intersect(c(use06, paste0("het_", names(cmps))), gr$run), function(id) if (paste0(id, "_w50") %in% gr$run) {
  a <- merge(gr[run == id, .(GENE, z10 = ZSTAT, p10 = P)], gr[run == paste0(id, "_w50"), .(GENE, z50 = ZSTAT, p50 = P)], by = "GENE")
  data.table(model = id, z_cor = round(cor(a$z10, a$z50), 3), top50_overlap = length(intersect(head(a[order(p10), GENE], 50), head(a[order(p50), GENE], 50))),
             top1_w10 = a[which.min(p10), GENE], top1_w50 = a[which.min(p50), GENE])
}))
if (nrow(ws)) ws[, `:=`(top1_w10 = bg$SYMBOL[match(top1_w10, bg$GENE)], top1_w50 = bg$SYMBOL[match(top1_w50, bg$GENE)])]

near_genes <- function(chr, bp) {
  g <- gloc[CHR == chr]; if (!nrow(g)) return(list(win = "", nearest = NA_character_, ngene = NA_integer_, dist = NA_real_))
  d <- pmax(0, g$START - bp, bp - g$STOP)
  w <- g[d <= win_main * 1e3]; i <- which.min(d)
  list(win = paste(w$SYMBOL, collapse = ";"), nearest = g$SYMBOL[i], ngene = g$GENE[i], dist = d[i] / 1e3)
}
tl <- rbindlist(lapply(intersect(use06, all_runs), function(id) {
  d <- read_snp_p(id)[is.finite(p)]
  gz <- gr[run == id, .(GENE, magma_p = P)]
  rbindlist(lapply(names(trend_thr), function(tn) {
    li <- clump_leads(d$Chr, d$bp, d$p, trend_thr[tn]); if (!length(li)) return(NULL)
    x <- d[li, .(Chr, SNP, bp, p)]
    ng <- rbindlist(lapply(seq_len(nrow(x)), function(i) as.data.table(near_genes(x$Chr[i], x$bp[i]))))
    x <- cbind(x, ng)
    x[, magma_p_nearest := gz$magma_p[match(ngene, gz$GENE)]]
    x[, sexflag_500kb := sapply(seq_len(.N), function(i) sum(sexflag$Chr == Chr[i] & abs(sexflag$bp - bp[i]) <= locus_kb * 1e3))]
    x[, `:=`(model = id, thr = tn)]
  }))
}), fill = TRUE)
if (nrow(tl)) {
  setcolorder(tl, c("model", "thr", "Chr", "SNP", "bp", "p", "win", "nearest", "dist", "magma_p_nearest", "sexflag_500kb"))
  setnames(tl, c("win", "nearest", "dist"), c("genes_10kb", "nearest_gene", "nearest_kb")); tl[, ngene := NULL]
  fwrite(tl, paste0(out_prefix, "_trend_loci.csv"), bom = TRUE)
}
fwrite(fl[, .(GENE, SYMBOL, NAME, CHR, START, STOP, n_snp, n_sexflag, hs_symbol, hs_chr, hs_on_xy, dog_chrx_same_symbol)],
       paste0(out_prefix, "_sexlinked_genes.csv"), bom = TRUE)

sink(paste0(out_prefix, "_magma_summary.txt"))
cat(sprintf("06-3 MAGMA 유전자 수준 분석 (window ±%d kb, snp-wise=mean, 염색체 %s)\n", win_main, if (length(chrs) == 38) "1 ~ 38" else paste(chrs, collapse = " ")))
cat("λ: 유전자 p의 median χ² / 0.455 | bonf: p < 0.05 / 유전자 수 | fdr05·fdr10: BH q < 0.05 · < 0.10 (0.10은 탐색적) | top_flag: 상위 50 유전자 중 성 연관 표시 유전자 수\n\n")
print(sm[kind != "null"], row.names = FALSE)
if (nrow(sm[kind == "null"])) {
  cat(sprintf("\nnull (주 분석 표본 순열 GWAS %d개 → MAGMA) 평균: λ %.3f, p < 1e-3 %.1f개 (범위 %d ~ %d), p < 1e-4 %.1f개, Bonferroni %.1f개, FDR < 0.05 %.1f개, FDR < 0.10 %.2f개\n",
              nl$n, nl$lambda, nl$p_1e3, nl$p_1e3_min, nl$p_1e3_max, nl$p_1e4, nl$bonf, nl$fdr05, nl$fdr10))
  cat(sprintf("  q < 0.10 유전자가 하나라도 나온 null run: %d / %d개 (%.0f%%) → 모형에서 q < 0.10이 하나쯤 나오는 것은 우연으로도 흔함\n",
              round(nl$any_fdr10 * nl$n), nl$n, 100 * nl$any_fdr10))
  cat("  → 모형의 유전자 수를 이 null과 비교 (유전자 수준에서도 우연 수준인지)\n")
}
if (nrow(rx)) {
  cat("\n완화 기준 (탐색·경향성용; 근거 아님): q < 0.10 유전자 (모든 run)\n")
  if (nrow(rx[q_lt_0.10 == TRUE])) print(rx[q_lt_0.10 == TRUE, .(run, SYMBOL, CHR, START, P, Q)], row.names = FALSE) else cat("  없음\n")
  cat("suggestive (p < 1e-3) 유전자 (run별 개수와 이름; null 평균 약 18개; 전체 표는 _gene_relaxed.csv)\n")
  print(rx[P < 1e-3, .(n = .N, genes = strtrim(paste(fifelse(is.na(SYMBOL), as.character(GENE), SYMBOL), collapse = ","), 150)), by = .(run, kind)], row.names = FALSE)
}
if (nrow(ws)) { cat("\nwindow 민감도 (±10 kb vs ±50 kb): 유전자 z 상관, 상위 50 유전자 겹침, 1위 유전자\n"); print(ws, row.names = FALSE) }
cat(sprintf("\n성 연관 의심 SNP (한 성별 단형성 & 다른 성별 MAF ≥ %.2f): %d개 — 염색체별 %s\n", sexflag_min_freq, nrow(sexflag),
            paste(sexflag[, .N, by = Chr][order(-N), sprintf("chr%d %d", Chr, N)], collapse = ", ")))
cat(sprintf("그 SNP이 ±%d kb에 든 유전자 %d개 (사람 ortholog가 X·Y에 있는 것 %d개, 개 chrX에 같은 이름 유전자 %d개)\n",
            win_main, nrow(fl), sum(fl$hs_on_xy), sum(fl$dog_chrx_same_symbol)))
cat(sprintf("  참고: 전체 유전자 중 사람 ortholog가 X·Y에 있는 것 %d / %d개 (개 상염색체 유전자라 원래 드묾; 검정은 하지 않음)\n", sum(bg$hs_on_xy), nrow(bg)))
if (nrow(fl)) print(fl[, .(SYMBOL, CHR, START, n_snp, n_sexflag, hs_symbol, hs_chr, dog_chrx_same_symbol)], row.names = FALSE)
cat("표시 SNP 묶음(염색체)별 위치와 가장 가까운 유전자 (괄호 = 유전자 몸체까지 거리):\n"); print(sf_near, row.names = FALSE)
if (nrow(tl)) {
  cat(sprintf("\n경향성 파악용 시사적 locus (근거 아님; SNP locus 실측 null 평균 1e-5 %.1f개 · 1e-4 %.0f개)\n", null_loci["1e-5"], null_loci["1e-4"]))
  print(tl[, .(loci = .N, with_gene_10kb = sum(nzchar(genes_10kb)), magma_p_lt_0.05 = sum(magma_p_nearest < 0.05, na.rm = TRUE), sexflag = sum(sexflag_500kb > 0)),
           by = .(model, thr)], row.names = FALSE)
  cat("\n주 분석 p < 1e-5 locus:\n"); print(tl[model == "main" & thr == "1e-5", .(Chr, SNP, bp, p = signif(p, 3), genes_10kb, nearest_gene, nearest_kb = round(nearest_kb, 1),
                                                                                   magma_p_nearest = signif(magma_p_nearest, 3))], row.names = FALSE)
}
sink()
cat(readLines(paste0(out_prefix, "_magma_summary.txt")), sep = "\n")

qq_of <- function(p) { p <- sort(p[is.finite(p) & p > 0]); n <- length(p); data.table(exp = -log10((seq_len(n) - 0.5) / n), obs = -log10(p), i = seq_len(n), n = n) }
qd <- gr[, qq_of(P), by = .(run, kind)]
nulls <- qd[kind == "null"]
shown <- qd[!kind %in% c("null", "sens", "het_sens")]
shown[, run := factor(run, levels = intersect(all_runs, unique(run)))]
band <- shown[, .(exp, lo = -log10(qbeta(0.975, i, n - i + 1)), hi = -log10(qbeta(0.025, i, n - i + 1))), by = run]
ng <- if (nrow(nulls)) rbindlist(lapply(levels(shown$run), function(r) { n_r <- shown[run == r, n[1]]; pr <- (seq_len(n_r) - 0.5) / n_r
  rbindlist(lapply(unique(nulls$run), function(k) { pk <- sort(gr[run == k, P]); data.table(run = r, k = k, exp = -log10(pr), obs = -log10(quantile(pk, pr, type = 1, names = FALSE))) })) })) else NULL
if (!is.null(ng)) ng[, run := factor(run, levels = levels(shown$run))]
gq <- ggplot() +
  geom_ribbon(data = band, aes(exp, ymin = lo, ymax = hi), fill = "#E1E0D9") +
  geom_abline(slope = 1, intercept = 0, color = "#898781", linewidth = 0.4) +
  { if (!is.null(ng)) geom_line(data = ng, aes(exp, obs, group = k), color = "#B4B2A9", linewidth = 0.4) } +
  geom_point(data = shown, aes(exp, obs), color = "#2A78D6", size = 0.5) +
  facet_wrap(~ run, ncol = 5) +
  labs(title = "Gene-level QQ (MAGMA)", subtitle = "grey band = 95% null interval; grey lines = permuted-phenotype GWAS (main sample) through the same MAGMA",
       x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p))) +
  theme_qc + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
ggsave(paste0(pic_prefix, "_gene_qq.png"), gq, width = 16, height = 3.2 * ceiling(length(levels(shown$run)) / 5) + 1.2, dpi = 300, limitsize = FALSE)

if ("main" %in% gr$run) {
  gm <- gr[run == "main"]
  off <- gloc[, .(len = max(STOP)), by = CHR][order(CHR)][, off := cumsum(as.numeric(shift(len, fill = 0)))]
  gm[, xpos := (START + STOP) / 2 + off$off[match(CHR, off$CHR)]]
  axis_df <- off[, .(CHR, mid = off + len / 2)]
  top <- head(gm[order(P)], 10)
  bonf <- -log10(0.05 / nrow(gm))
  gmh <- ggplot(gm, aes(xpos, -log10(P))) +
    geom_point(aes(color = factor(CHR %% 2)), size = 0.9, show.legend = FALSE) +
    geom_point(data = gm[n_sexflag > 0], shape = 4, size = 2.6, stroke = 0.9, color = "#EB6834") +
    geom_hline(yintercept = bonf, color = "#52514E", linetype = "dashed", linewidth = 0.4) +
    annotate("label", x = max(gm$xpos), y = bonf, label = sprintf("Bonferroni (0.05 / %s genes)", format(nrow(gm), big.mark = ",")),
             hjust = 1, size = 3.1, color = "#52514E", fill = "white", label.size = 0) +
    geom_text(data = top, aes(label = SYMBOL), vjust = -0.7, size = 3, color = "#1F1E1C", check_overlap = TRUE) +
    scale_color_manual(values = c("0" = col_chr[1], "1" = col_chr[2])) +
    scale_x_continuous(breaks = axis_df$mid, labels = ifelse(axis_df$CHR %% 2 == 1 | axis_df$CHR <= 10, axis_df$CHR, ""), expand = c(0.005, 0)) +
    scale_y_continuous(limits = c(0, max(bonf, -log10(min(gm$P)), na.rm = TRUE) + 0.8), expand = c(0, 0)) +
    labs(title = sprintf("Gene-level Manhattan: main (MAGMA, +/-%d kb)", win_main),
         subtitle = sprintf("lambda %.3f; labels = top 10 genes; orange x = gene containing a sex-linked-suspect SNP (kept in analysis)", p_to_lambda(gm$P)),
         x = "Chromosome", y = expression(-log[10](p[gene]))) +
    coord_cartesian(clip = "off") +
    theme_qc + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(10, 30, 10, 10),
                     plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#52514E"))
  ggsave(paste0(pic_prefix, "_main_gene_manhattan.png"), gmh, width = 16, height = 5.5, dpi = 300)
}
cat(sprintf("저장: %s_magma/, %s_gene_results.csv.gz, %s_magma_summary.*, %s_trend_loci.csv, %s_sexlinked_*, 그림 %s_gene_qq.png · _main_gene_manhattan.png\n",
            out_prefix, out_prefix, out_prefix, out_prefix, out_prefix, pic_prefix))

#    - main
#      · λ 0.915
#      · young8 0.63; sex_m, sex_f 0.46; purelike 0.45; mixlike 0.44
#      · wt_t1 ~ t3 0.27 ~ 0.31; age8 0.31; bin_age8 0.04; int_age −0.01
#      · sex_m–sex_f 0.008; wt_t1–t2 −0.017; wt_t1–t3 −0.006; wt_t2–t3 0.041
#      · young8–age8 0.023; purelike–mixlike −0.002
#      · main: CFI
#      · int_age: ZEB1
#    - het_age
#    - het_sex
#      CIDEC/RPUSD3/TTLL3 (chr20), GRXCR1 (chr13), PACSIN2 (chr10)
#    - FDR q < 0.10
#    het_sex_w50       0.829  8          1          0                SPARCL1 → SPARCL1        0.763
#    het_age_w50       0.796  11         5          0                ADORA2B → TTC19          0.767
#    het_ancestry_w50  0.873  6          1          0                TMPRSS7 → TAGLN3         0.771
#    het_weight_w50    0.721  3          0          0                PAPOLA → TULP4           0.760
