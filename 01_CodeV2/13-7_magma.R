# 13-7_magma.R - MAGMA gene and gene-set analysis for age-floor, binary and SNP x body weight models. See README.md.
suppressMessages({ library(data.table); library(parallel) })
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
out_prefix <- file.path(qc_dir, "13-7")
magma <- "/data/00_tools/magma_v1.10/magma"; plink2 <- "/data/00_tools/plink2_linux_avx2_20231123/plink2"
bfile <- file.path(d0, "03-3_maf0.01_hwe")
res63 <- file.path(d0, "06-3_magma")
annot <- c(`10` = file.path(d0, "06-3_work", "annot_w10.genes.annot"), `50` = file.path(d0, "06-3_work", "annot_w50.genes.annot"))
gene_loc_file <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/06-1_gene_loc.txt"
set_file <- c(primary = file.path(qc_dir, "06-4_sets_primary.txt"), aux = file.path(qc_dir, "06-4_sets_aux.txt"))
int_wt_f <- Sys.getenv("INT_WT_FILE", file.path(qc_dir, "13-6_regenie", "int_wt.regenie.gz"))
loci_breed_f <- file.path(qc_dir, "13-4_loci_breed.csv")
n_cores <- as.integer(Sys.getenv("N_CORES", "40")); n_batch <- as.integer(Sys.getenv("N_BATCH", "40"))
chrs <- if (nzchar(Sys.getenv("CHRS"))) as.integer(strsplit(Sys.getenv("CHRS"), " ")[[1]]) else 1:38
chr_tag <- if (length(chrs) == 38) "" else paste0("_chr", paste(chrs, collapse = "-"))
res_dir <- paste0(out_prefix, "_magma", chr_tag); gsa_dir <- paste0(out_prefix, "_gsa"); work_dir <- paste0(out_prefix, "_work")
for (d in c(res_dir, gsa_dir, work_dir, file.path(work_dir, c("pval", "refbed", "batch")))) dir.create(d, showWarnings = FALSE, recursive = TRUE)
for (f in c(magma, plink2, paste0(bfile, ".bed"), annot, gene_loc_file, set_file, int_wt_f, loci_breed_f)) if (!file.exists(f)) stop("파일 없음: ", f)
hap_kb <- 500
old12 <- c("main", "purelike", "mixlike", "int_age", "age8", "bin_age8", "young8", "sex_m", "sex_f", "wt_t1", "wt_t2", "wt_t3")
new8 <- c("agefloor6", "agefloor7", "agefloor9", "age10", "agefloor11", "agefloor12", "bin_all", "int_wt")
all20 <- c("young8", "main", "agefloor6", "agefloor7", "age8", "agefloor9", "age10", "agefloor11", "agefloor12", "wt_t1", "wt_t2", "wt_t3",
           "sex_m", "sex_f", "purelike", "mixlike", "int_age", "int_wt", "bin_all", "bin_age8")
n_all <- nrow(fread(paste0(bfile, ".fam"), header = FALSE, select = 1))
sh <- function(cmd, args, log) { st <- system2(cmd, args, stdout = log, stderr = log); if (st != 0) stop(basename(cmd), " 실패: ", log); invisible(st) }

read_p <- function(id) {
  if (id == "int_wt" || id == "bin_all") {
    f <- if (id == "int_wt") int_wt_f else file.path(d0, "05-3_regenie", "bin_all.regenie.gz")
    x <- fread(cmd = paste("zcat", shQuote(f)), select = c("ID", "TEST", "LOG10P", "N"))
    tsel <- if (id == "int_wt") "ADD-INT_SNPxlog2_wt" else "ADD"
    if (!tsel %in% x$TEST) stop(id, ": TEST ", tsel, " 없음")
    x <- x[TEST == tsel & !is.na(LOG10P)]; list(p = x[, .(SNP = ID, P = 10^(-LOG10P))], N = n_all, keep = NA_character_)
  } else {
    x <- fread(cmd = paste("zcat", shQuote(file.path(d0, "05-2_gwas", paste0(id, ".mlma.gz")))), select = c("SNP", "p"))
    kf <- file.path(d0, "05-1_inputs", paste0(id, "_keep.txt"))
    list(p = x[is.finite(p) & p > 0, .(SNP, P = p)], N = length(readLines(kf)), keep = kf)
  }
}
jobs <- rbindlist(lapply(new8, function(id) data.table(run = c(id, paste0(id, "_w50")), id = id, window = c(10L, 50L))))
done <- function(run) all(file.exists(file.path(res_dir, paste0(run, c(".genes.out", ".genes.raw")))))
todo <- jobs[!sapply(run, done)]
cat(sprintf("[%s] MAGMA 유전자 run %d개 (끝난 것 %d개는 건너뜀)%s\n", format(Sys.time(), "%T"), nrow(jobs), nrow(jobs) - nrow(todo),
            if (nzchar(chr_tag)) paste0(" — 시험: 염색체 ", paste(chrs, collapse = " ")) else ""))
info <- list()
for (id in unique(todo$id)) {
  r <- read_p(id); pf <- file.path(work_dir, "pval", paste0(id, ".txt")); fwrite(r$p, pf, sep = "\t")
  ref <- bfile
  if (!is.na(r$keep)) {
    ref <- file.path(work_dir, "refbed", paste0(id, chr_tag))
    if (!file.exists(paste0(ref, ".bed"))) sh(plink2, c("--dog", "--bfile", bfile, "--keep", r$keep, if (nzchar(chr_tag)) c("--chr", chrs), "--make-bed",
                                                        "--threads", min(n_cores, 16), "--out", ref), paste0(ref, ".stdout"))
  }
  info[[id]] <- list(pf = pf, N = r$N, ref = ref, nsnp = nrow(r$p))
  cat(sprintf("  %s: SNP %s개 · N %d · LD 기준 %s\n", id, format(nrow(r$p), big.mark = ","), r$N, if (is.na(r$keep)) "전체" else basename(r$keep)))
}
annot_of <- function(w) {
  if (!nzchar(chr_tag)) return(annot[[as.character(w)]])
  f <- file.path(work_dir, sprintf("annot_w%d%s.genes.annot", w, chr_tag))
  if (!file.exists(f)) { a <- readLines(annot[[as.character(w)]])
    keep_l <- grepl("^#", a) | sub(":.*", "", sapply(strsplit(a, "[\t ]+"), `[`, 2)) %in% as.character(chrs); writeLines(a[keep_l], f) }
  f
}

if (nrow(todo)) {
  bj <- CJ(run = todo$run, b = seq_len(n_batch)); bdir <- file.path(work_dir, "batch")
  t0 <- Sys.time()
  bres <- mclapply(seq_len(nrow(bj)), function(i) {
    run <- bj$run[i]; b <- bj$b[i]; j <- todo[run == bj$run[i]]; inf <- info[[j$id]]; out <- file.path(bdir, paste0(run, chr_tag))
    if (length(list.files(bdir, pattern = sprintf("^%s\\.batch%d_[0-9]+\\.genes\\.out$", paste0(run, chr_tag), b)))) return("있음")
    st <- system2(magma, c("--bfile", inf$ref, "--pval", inf$pf, sprintf("N=%d", inf$N), "--gene-annot", annot_of(j$window), "--batch", b, n_batch, "--out", out),
                  stdout = FALSE, stderr = FALSE)
    if (st == 0) return("완료")
    lg <- sprintf("%s.batch%d_%d.log", out, b, n_batch); if (file.exists(lg) && any(grepl("is empty", readLines(lg)))) "빔" else paste("실패", run, b)
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  bad <- unlist(bres)[!unlist(bres) %in% c("완료", "있음", "빔")]; if (length(bad)) stop("MAGMA batch 실패: ", paste(head(bad, 10), collapse = "; "))
  cat(sprintf("[%s] MAGMA batch 끝 (%.1f분)\n", format(Sys.time(), "%T"), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  for (run in todo$run) {
    out <- file.path(bdir, paste0(run, chr_tag)); bf <- list.files(bdir, pattern = paste0("^", paste0(run, chr_tag), "\\.batch[0-9]+_[0-9]+\\.genes\\.out$"), full.names = TRUE)
    if (length(bf) == 1) { for (ext in c(".genes.out", ".genes.raw")) file.copy(sub("\\.genes\\.out$", ext, bf), file.path(res_dir, paste0(run, ext)), overwrite = TRUE)
    } else sh(magma, c("--merge", out, "--out", file.path(res_dir, run)), file.path(bdir, paste0(run, ".merge.stdout")))
    file.remove(list.files(bdir, pattern = paste0("^", paste0(run, chr_tag), "\\.batch[0-9]+_[0-9]+\\."), full.names = TRUE))
  }
}
for (id in new8) { ref <- file.path(work_dir, "refbed", paste0(id, chr_tag)); if (file.exists(paste0(ref, ".bed")) && all(sapply(c(id, paste0(id, "_w50")), done))) file.remove(paste0(ref, c(".bed", ".bim", ".fam"))) }
raw_of <- function(run) { base <- sub("_w50$", "", run); file.path(if (base %in% old12) res63 else res_dir, paste0(run, ".genes.raw")) }
GS <- rbindlist(lapply(c(all20, paste0(all20, "_w50")), function(run) {
  f <- sub("\\.genes\\.raw$", ".genes.out", raw_of(run)); if (!file.exists(f)) return(NULL)
  g <- fread(f, skip = "GENE"); data.table(run = run, source = if (sub("_w50$", "", run) %in% old12) "06-3" else "13-7", n_gene = nrow(g),
                                           lambda = median(qchisq(g$P, 1, lower.tail = FALSE)) / qchisq(0.5, 1), n_p1e3 = sum(g$P < 1e-3), n_p1e4 = sum(g$P < 1e-4),
                                           bonf = sum(g$P < 0.05 / nrow(g)), top5 = paste(head(g[order(P), GENE], 5), collapse = ","))
}))
if (nzchar(chr_tag)) {
  out_prefix <- paste0(out_prefix, chr_tag); gsa_dir <- paste0(gsa_dir, chr_tag); dir.create(gsa_dir, showWarnings = FALSE)
  cat("*** 시험: 새 8개 run은 일부 염색체만 → 세트 분석 결과 해석 금지 ***\n")
}

lb <- unique(fread(loci_breed_f)[, .(SNP, chr, bp, breed_conc, top_breed, top_breed_far2_kb)], by = "SNP")
hap <- lb[breed_conc %in% TRUE & !is.na(top_breed_far2_kb) & top_breed_far2_kb >= hap_kb]
gl <- fread(gene_loc_file, header = FALSE, col.names = c("GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"))
ex <- rbindlist(lapply(seq_len(nrow(hap)), function(i) { h <- hap[i]; w <- h$top_breed_far2_kb * 1e3
  gl[CHR == h$chr & STOP >= h$bp - w & START <= h$bp + w, .(GENE, SYMBOL, SNP = h$SNP, top_breed = h$top_breed, span_kb = h$top_breed_far2_kb)] }))
fwrite(ex, paste0(out_prefix, "_nobreed_genes.csv"))
sp <- strsplit(readLines(set_file[["primary"]]), "[\t ]+")
nb <- vapply(sp, function(x) paste(c(x[1], setdiff(x[-1], as.character(ex$GENE))), collapse = " "), "")
nb_f <- paste0(out_prefix, "_sets_primary_nobreed.txt"); writeLines(nb, nb_f)
n_rm <- data.table(set = sapply(sp, `[`, 1), n_before = lengths(sp) - 1L, n_removed = vapply(sp, function(x) sum(x[-1] %in% as.character(ex$GENE)), 0L))
cat(sprintf("품종 하플로타입 후보 locus %d개 → 주변 유전자 %d개를 주 가설 세트에서 뺌 (세트에서 빠진 유전자 합 %d개)\n", nrow(hap), uniqueN(ex$GENE), sum(n_rm$n_removed)))

sfile <- c(primary = set_file[["primary"]], aux = set_file[["aux"]], nobreed = nb_f)
GJ <- rbind(data.table(run = all20, variant = "primary"), data.table(run = paste0(all20, "_w50"), variant = "primary"),
            data.table(run = all20, variant = "aux"), data.table(run = all20, variant = "nobreed"))
GJ[, `:=`(raw = vapply(run, raw_of, ""), out = file.path(gsa_dir, paste0(run, "__", variant)))]
if (length(miss <- GJ[!file.exists(raw), unique(raw)])) stop("genes.raw 없음: ", paste(basename(miss), collapse = ", "))
GJ[, done := file.exists(paste0(out, ".gsa.out")) & file.mtime(paste0(out, ".gsa.out")) > pmax(file.mtime(sfile[variant]), file.mtime(raw))]
cat(sprintf("[%s] MAGMA 세트 분석 %d개 (끝난 것 %d개 건너뜀)\n", format(Sys.time(), "%T"), nrow(GJ), sum(GJ$done)))
td <- GJ[done == FALSE]
st <- unlist(mclapply(seq_len(nrow(td)), function(i) system2(magma, c("--gene-results", shQuote(td$raw[i]), "--set-annot", shQuote(sfile[[td$variant[i]]]),
                                                                       "--out", shQuote(td$out[i])), stdout = FALSE, stderr = FALSE), mc.cores = n_cores, mc.preschedule = FALSE))
if (length(st) && any(st != 0)) stop("MAGMA 세트 분석 실패: ", paste(basename(td$out[st != 0]), collapse = ", "))
read_gsa <- function(f) { d <- fread(f, skip = "VARIABLE"); if ("FULL_NAME" %in% names(d)) d[, VARIABLE := FULL_NAME]; d[, .(set = VARIABLE, NGENES, BETA, BETA_STD, SE, P)] }
GA <- rbindlist(lapply(seq_len(nrow(GJ)), function(i) { x <- read_gsa(paste0(GJ$out[i], ".gsa.out")); x[, `:=`(run = GJ$run[i], variant = GJ$variant[i])]; x }))
GA[, `:=`(model = sub("_w50$", "", run), window = fifelse(grepl("_w50$", run), 50L, 10L))]
fwrite(GA, paste0(out_prefix, "_gsa_all.csv.gz")); fwrite(GS, paste0(out_prefix, "_gene_summary.csv"))

sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("13-7 MAGMA 유전자 · 세트 분석 (계산만; 판정은 13-8)\n")
cat("■ 유전자 수준 run 요약 (source 06-3 = 기존 결과 재사용, 13-7 = 이 파일)\n"); print(GS[, lapply(.SD, function(v) if (is.numeric(v) && !is.integer(v)) signif(v, 3) else v)], nrows = 50)
cat("\n■ E 민감도 세트: 세트별로 뺀 유전자 수\n"); print(n_rm)
cat("\n■ 주 가설 세트 p < 0.05 (±10 kb, run별 개수; 판정 아님)\n")
print(dcast(GA[variant == "primary" & window == 10L, .(n_p05 = sum(P < 0.05)), by = .(model)], . ~ model, value.var = "n_p05"))
sink()
cat(sprintf("저장: %s_* → 13-8\n", out_prefix))

