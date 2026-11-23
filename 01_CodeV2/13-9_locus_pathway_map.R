# 13-9_locus_pathway_map.R - Maps lead loci to genes and pathways for each model (input to 13-10). See README.md.
#   04_GWAS_2023/06-4_ref/palmer2021_brain_signature.tsv · gse304125_cortex_normcounts.tsv.gz · 13-9_ref/rna_tissue_consensus.tsv.gz (HPA, hpa_version.txt)

suppressMessages({ library(data.table); library(ggplot2); library(parallel); library(Matrix) })
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
pic_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/05_GWAS_2023"
test_mode <- Sys.getenv("TEST", "0") == "1"; tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("13-9", tag))
pre <- function(x) file.path(qc_dir, x)
n_cores <- as.integer(Sys.getenv("N_CORES", "8"))
plink1 <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
bfile <- file.path(d0, "03-3_maf0.01_hwe")
g22 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure"
f_in <- c(cmp = pre("13-1_compare.csv"), l33 = pre("13-3_loci.csv"), cs = pre("13-3_susie_cs.csv"), css = pre("13-3_susie_cs_snps.csv"),
          l35 = pre("13-5_loci_cov.csv"), sp = pre("06-4_sets_primary.txt"), sa = pre("06-4_sets_aux.txt"), si = pre("06-4_set_info.csv"),
          ai = pre("06-4_aux_info.csv"), om = pre("06-4_ortholog_map.csv"), gl = file.path(g22, "06-1_gene_loc.txt"))
f_ex <- c(pal = file.path(d0, "06-4_ref", "palmer2021_brain_signature.tsv"), ctx = file.path(d0, "06-4_ref", "gse304125_cortex_normcounts.tsv.gz"),
          hpa = file.path(d0, "13-9_ref", "rna_tissue_consensus.tsv.gz"))
for (f in c(f_in, f_ex, paste0(bfile, ".bed"), plink1)) if (!file.exists(f)) stop("파일 없음: ", f)
int_wt_f <- Sys.getenv("INT_WT_FILE", pre("13-6_regenie/int_wt.regenie.gz")); has_iw <- file.exists(int_wt_f)
res63 <- file.path(d0, "06-3_magma"); res137 <- Sys.getenv("MAGMA13_DIR", pre("13-7_magma"))
r2_min <- 0.5; win_kb <- 10; ld_kb <- 1000; near_max_kb <- 500; gap_kb <- 500; n_perm <- if (test_mode) 50L else 1000L
set_min <- 10; set_max <- 500
age_series <- c("young8", "main", "agefloor6", "agefloor7", "age8", "agefloor9", "age10", "agefloor11", "agefloor12")
wt_series <- c("wt_t1", "wt_t2", "wt_t3")
run_order <- c(age_series, wt_series, "sex_m", "sex_f", "purelike", "mixlike", "int_age", "int_wt", "bin_all", "bin_age8")
old12 <- c("main", "purelike", "mixlike", "int_age", "age8", "bin_age8", "young8", "sex_m", "sex_f", "wt_t1", "wt_t2", "wt_t3")
series_of <- function(m) fifelse(m %in% age_series, "나이", fifelse(m %in% wt_series, "체중", "기타"))
ld_grp <- function(m) fifelse(m %in% c("main", "int_age", "bin_all", "int_wt"), "pooled", fifelse(m == "bin_age8", "age8", m))
keep_of <- function(g) if (g == "pooled") file.path(d0, "04-1_grm", "grm_all.grm.id") else file.path(d0, "05-1_inputs", paste0(g, "_keep.txt"))
wd <- tempfile("r13_9_"); dir.create(wd)
f3 <- function(x) { x <- copy(x); for (j in names(x)) if (is.numeric(x[[j]]) && !is.integer(x[[j]])) set(x, j = j, value = signif(x[[j]], 3)); x }
col_of <- function(x, nm) intersect(nm, names(x))[1]
sc_names <- c(paste0("spec_", "class"), paste0("spec", "ass"))
t0 <- Sys.time(); set.seed(20261116)

read_sets <- function(f) { s <- strsplit(readLines(f), "[\t ]+"); setNames(lapply(s, function(x) as.integer(x[-1])), sapply(s, `[`, 1)) }
SP <- read_sets(f_in[["sp"]]); SA <- read_sets(f_in[["sa"]]); SA <- SA[lengths(SA) >= set_min & lengths(SA) <= set_max & !names(SA) %in% names(SP)]
si <- fread(f_in[["si"]]); ai <- fread(f_in[["ai"]])
set_name <- c(setNames(si$name, si$set), setNames(ai$name, ai$set)); fam_of <- setNames(si$family, si$set)
gl <- fread(f_in[["gl"]], header = FALSE, col.names = c("GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"))
gl <- gl[CHR %in% as.character(1:38)][, CHR := as.integer(CHR)]
om <- fread(f_in[["om"]])[, .(hs = paste(unique(hs_symbol), collapse = "/")), by = .(GENE = dog)]
gl <- merge(gl, om, by = "GENE", all.x = TRUE)
gl[, lab := fifelse(!is.na(hs) & toupper(hs) != toupper(SYMBOL), paste0(SYMBOL, "/", hs), SYMBOL)]
setkey(gl, GENE)
all_sets <- c(SP, SA); set_src <- c(setNames(rep("주 가설", length(SP)), names(SP)),
                                     setNames(fifelse(grepl("^KEGG", names(SA)), "KEGG", "GO"), names(SA)))
gidx <- setNames(seq_len(nrow(gl)), gl$GENE)
GS <- sparseMatrix(i = unlist(lapply(all_sets, function(s) gidx[as.character(s[as.character(s) %in% names(gidx)])])),
                   j = rep(seq_along(all_sets), vapply(all_sets, function(s) sum(as.character(s) %in% names(gidx)), 0L)),
                   x = 1, dims = c(nrow(gl), length(all_sets)), dimnames = list(gl$GENE, names(all_sets)))
GS@x[] <- 1
set_size <- colSums(GS)
cat(sprintf("[%s] 유전자 %d개 · 경로 %d개 (주 가설 %d · KEGG %d · GO %d; 보조는 유전자 %d ~ %d개)\n", format(Sys.time(), "%T"), nrow(gl), ncol(GS),
            length(SP), sum(set_src == "KEGG"), sum(set_src == "GO"), set_min, set_max))

bimp <- fread(paste0(bfile, ".bim"), header = FALSE, select = c(1, 2, 4), col.names = c("chr", "SNP", "bp"))
cmp <- fread(f_in[["cmp"]])
L <- cmp[thr == "1e-4" | thr == 1e-4, .(model, SNP, chr, bp, p, b, maf, region_id, level = fifelse(has_1e5 %in% TRUE, "1e-5", "1e-4"), near_13_1 = genes)]
if (has_iw) {
  iw <- fread(cmd = paste("zcat", shQuote(int_wt_f)), select = c("CHROM", "GENPOS", "ID", "A1FREQ", "TEST", "BETA", "LOG10P"))
  iw <- iw[TEST == "ADD-INT_SNPxlog2_wt" & LOG10P > 4][order(CHROM, GENPOS)]
  iw[, grp := cumsum(c(TRUE, diff(GENPOS) > gap_kb * 1e3 | diff(CHROM) != 0))]
  IL <- iw[, { i <- which.max(LOG10P); .(SNP = ID[i], chr = CHROM[i], bp = GENPOS[i], p = 10^(-LOG10P[i]), b = BETA[i], maf = min(A1FREQ[i], 1 - A1FREQ[i]),
                                         level = if (max(LOG10P) > 5) "1e-5" else "1e-4") }, by = grp][, grp := NULL]
  RG <- cmp[, .(lo = min(bp) - 250e3, hi = max(bp) + 250e3), by = .(region_id, chr)]
  IL[, region_id := vapply(seq_len(.N), function(i) { r <- RG[chr == IL$chr[i] & lo <= IL$bp[i] & hi >= IL$bp[i], region_id]; if (length(r)) r[1] else paste0("chr", IL$chr[i], "_iw", i) }, "")]
  near3 <- function(c0, b0) { g <- gl[CHR == c0][, d := pmax(0, START - b0, b0 - STOP)][order(d)][1:3]
    paste(sprintf("%s(%s)", g$SYMBOL, fifelse(g$d == 0, "in", paste0(round(g$d / 1e3), "kb"))), collapse = ", ") }
  IL[, near_13_1 := mapply(near3, chr, bp)]
  L <- rbind(L, IL[, .(model = "int_wt", SNP, chr, bp, p, b, maf, region_id, level, near_13_1)])
  cat(sprintf("[%s] int_wt locus: 1e-4 %d개 (1e-5 %d개)\n", format(Sys.time(), "%T"), nrow(IL), sum(IL$level == "1e-5")))
} else cat("int_wt 결과 없음 (", int_wt_f, ") → int_wt 건너뜀\n")
if (test_mode) L <- L[model %in% c("main", "wt_t1", if (has_iw) "int_wt" else "age8")]
L[, series := series_of(model)]
L[, n_models_region := uniqueN(model), by = region_id]
l33 <- fread(f_in[["l33"]]); l35 <- fread(f_in[["l35"]])
A3 <- l33[thr == "1e-4" | thr == 1e-4, .(model, SNP, n_cs, lead_pip, spec_grp = get(col_of(l33, sc_names)), shared)]
A5 <- l35[thr == "1e-4" | thr == 1e-4, .(model, SNP, breed_conc, top_breed, hap_cand, own_kind, wt_source)]
L <- merge(merge(L, unique(A3, by = c("model", "SNP")), by = c("model", "SNP"), all.x = TRUE), unique(A5, by = c("model", "SNP")), by = c("model", "SNP"), all.x = TRUE)
L[, lid := .I]
cat(sprintf("[%s] locus %d개 (모형 %d개; 1e-5 표시 %d개)\n", format(Sys.time(), "%T"), nrow(L), uniqueN(L$model), sum(L$level == "1e-5")))

cs <- fread(f_in[["cs"]]); css <- fread(f_in[["css"]])
csl <- merge(css, cs[, .(region_id, model, cs, lead)], by = c("region_id", "model", "cs"))
cs_has <- csl[, .(has_lead = any(SNP == lead[1])), by = .(region_id, model, cs, lead)][has_lead == TRUE]
cs_snp <- merge(csl, cs_has[, .(region_id, model, cs)], by = c("region_id", "model", "cs"))[, .(model, SNP = lead, snp_in = SNP)]
L[, grp := ld_grp(model)]
L[, range_src := fifelse(paste(model, SNP) %in% cs_snp[, paste(model, SNP)], "CS", "r2")]
need <- unique(L[range_src == "r2", .(grp, SNP)])
R2 <- rbindlist(mclapply(unique(need$grp), function(g) {
  writeLines(need[grp == g, unique(SNP)], file.path(wd, paste0("idx_", g, ".txt")))
  out <- file.path(wd, paste0("ld_", g))
  system2(plink1, c("--bfile", bfile, "--dog", "--threads", "1", "--keep", keep_of(g), "--r2", "--ld-snp-list", file.path(wd, paste0("idx_", g, ".txt")),
                    "--ld-window-kb", ld_kb, "--ld-window", "999999", "--ld-window-r2", r2_min, "--out", out), stdout = FALSE, stderr = FALSE)
  if (!file.exists(paste0(out, ".ld"))) stop("plink LD 실패: ", g)
  x <- fread(paste0(out, ".ld"))[, .(grp = g, SNP = SNP_A, b_in = BP_B)]
  rbind(x, data.table(grp = g, SNP = need[grp == g, unique(SNP)], b_in = bimp[match(need[grp == g, unique(SNP)], bimp$SNP), bp]))
}, mc.cores = min(n_cores, uniqueN(need$grp))))
rng_r2 <- merge(L[range_src == "r2", .(lid, grp, SNP)], R2, by = c("grp", "SNP"), allow.cartesian = TRUE)[, .(lid, b_in)]
rng_cs <- merge(L[range_src == "CS", .(lid, model, SNP)], cs_snp, by = c("model", "SNP"), allow.cartesian = TRUE)
rng_cs[, b_in := bimp[match(snp_in, bimp$SNP), bp]]
RNG <- rbind(rng_r2, rng_cs[, .(lid, b_in)])[, .(lo = min(b_in), hi = max(b_in), n_snp_range = .N), by = lid]
L <- merge(L, RNG, by = "lid", all.x = TRUE)
L[is.na(lo), `:=`(lo = bp, hi = bp, n_snp_range = 1L)]
L[, span_kb := (hi - lo) / 1e3]
cat(sprintf("[%s] SNP 범위: CS %d · r² %d, 넓이 중앙값 %.1f kb\n", format(Sys.time(), "%T"), sum(L$range_src == "CS"), sum(L$range_src == "r2"), median(L$span_kb)))

G2 <- gl[, .(chr = CHR, START, STOP, GENE)]; setkey(G2, chr, START, STOP)
map_genes <- function(D) {                                                      # D: id · chr · lo · hi → id · GENE · rule
  W <- D[, .(id, chr, s = lo - near_max_kb * 1e3, e = hi + near_max_kb * 1e3, lo, hi)]; setkey(W, chr, s, e)
  o <- foverlaps(G2, W, by.x = c("chr", "START", "STOP"), by.y = c("chr", "s", "e"), nomatch = NULL)
  o[, dist := pmax(0, START - hi, lo - STOP)]
  rg <- o[dist <= win_kb * 1e3, .(id, GENE, dist)]
  nr <- o[!id %in% rg$id][order(id, dist)][!duplicated(id), .(id, GENE, dist)]
  rbind(rg[, rule := "range"], nr[, rule := "nearest"])
}
LGo <- map_genes(L[, .(id = lid, chr, lo, hi)])
setnames(LGo, "id", "lid")
gp_of <- function(m) { f <- file.path(if (m %in% old12) res63 else res137, paste0(m, ".genes.out")); if (!file.exists(f)) return(NULL)
  fread(f, skip = "GENE", select = c("GENE", "P"))[, model := m] }
GP <- rbindlist(lapply(unique(L$model), gp_of))
if (!nrow(GP)) GP <- data.table(GENE = integer(), P = numeric(), model = character())
cat(sprintf("[%s] MAGMA 유전자 p: 모형 %d / %d개 있음%s\n", format(Sys.time(), "%T"), uniqueN(GP$model), uniqueN(L$model),
            if (uniqueN(GP$model) < uniqueN(L$model)) paste0(" (없음: ", paste(setdiff(unique(L$model), GP$model), collapse = " "), " → NA)") else ""))
LGo <- merge(LGo, L[, .(lid, model)], by = "lid")
LGo <- merge(LGo, GP, by = c("model", "GENE"), all.x = TRUE)
LGo[, lab := gl[.(LGo$GENE), lab]]
LGo[, sym := gl[.(LGo$GENE), SYMBOL]]

inc <- function(lg, ids) { j <- match(as.character(lg$GENE), rownames(GS)); ok <- !is.na(j)
  M <- sparseMatrix(i = match(lg$lid[ok], ids), j = j[ok], x = 1, dims = c(length(ids), nrow(GS))); M@x[] <- 1
  H <- M %*% GS; H@x[] <- as.numeric(H@x > 0); dimnames(H) <- list(ids, colnames(GS)); drop0(H) }
H_rg <- inc(LGo[rule == "range"], L$lid); H_nr <- inc(LGo, L$lid)
long_of <- function(rule_set, nm) {
  x <- LGo[rule %in% rule_set]; j <- match(as.character(x$GENE), rownames(GS))
  S <- summary(GS[j[!is.na(j)], , drop = FALSE]); xx <- x[!is.na(j)]
  data.table(lid = xx$lid[S$i], model = xx$model[S$i], GENE = xx$GENE[S$i], gene = xx$lab[S$i], gene_p = xx$P[S$i], gene_rule = xx$rule[S$i],
             set = colnames(GS)[S$j])[, `:=`(source = set_src[set], set_name = set_name[set], set_size = set_size[set])]
}
LGS <- long_of(c("range", "nearest"))
gstr <- LGo[order(lid, rule != "range", dist), .(genes_range = paste(lab[rule == "range"], collapse = ";"),
                 gene_p_range = paste(sprintf("%s(%s)", sym[rule == "range"], fifelse(is.na(P[rule == "range"]), "NA", formatC(P[rule == "range"], 1, format = "e"))), collapse = ";"),
                 min_gene_p = suppressWarnings(min(P[rule == "range"], na.rm = TRUE)),
                 nearest_gene = if (any(rule == "nearest")) lab[rule == "nearest"][1] else NA_character_,
                 nearest_kb = if (any(rule == "nearest")) round(dist[rule == "nearest"][1] / 1e3, 1) else NA_real_), by = lid]
gstr[is.infinite(min_gene_p), min_gene_p := NA]
sstr <- LGS[, .(primary_sets = paste(unique(set_name[source == "주 가설" & gene_rule == "range"]), collapse = " | "),
                primary_sets_nearest = paste(unique(set_name[source == "주 가설" & gene_rule == "nearest"]), collapse = " | "),
                kegg = { k <- unique(set_name[source == "KEGG"]); if (!length(k)) "" else paste0(if (all(gene_rule == "nearest")) "(가장 가까운 유전자) " else "", paste(k, collapse = " | ")) },
                n_go = uniqueN(set[source == "GO"])), by = lid]
L <- merge(L, gstr, by = "lid", all.x = TRUE); L <- merge(L, sstr, by = "lid", all.x = TRUE)
L[is.na(n_go), n_go := 0L]; for (j in c("genes_range", "primary_sets", "primary_sets_nearest", "kegg")) set(L, which(is.na(L[[j]])), j, "")
L[, n_genes_range := lengths(strsplit(genes_range, ";"))]
cat(sprintf("[%s] 유전자 연결: 범위 안 유전자 있음 %d / %d · 가장 가까운 유전자로만 %d · 없음 %d\n", format(Sys.time(), "%T"),
            sum(L$genes_range != ""), nrow(L), sum(L$genes_range == "" & !is.na(L$nearest_gene)), sum(L$genes_range == "" & is.na(L$nearest_gene))))

GSPEC <- data.table(GENE = as.integer(rownames(GS)), n_sets = as.integer(rowSums(GS)),
                    n_kegg = as.integer(rowSums(GS[, set_src[colnames(GS)] == "KEGG", drop = FALSE])),
                    n_go = as.integer(rowSums(GS[, set_src[colnames(GS)] == "GO", drop = FALSE])))
q_sp <- quantile(GSPEC[n_sets > 0, n_sets], c(0.25, 0.75), type = 1)
GSPEC[, path_spec := fifelse(n_sets == 0, "주석 없음", fifelse(n_sets <= q_sp[1], "특이", fifelse(n_sets >= q_sp[2], "범용", "중간")))]
GSPEC[, path_pct := fifelse(n_sets == 0, NA_real_, frank(fifelse(n_sets == 0, NA_integer_, n_sets), na.last = "keep", ties.method = "average") / sum(n_sets > 0))]
omr <- fread(f_in[["om"]])[, .(GENE = dog, hs_entrez = human, hs_sym = hs_symbol)]
pal <- fread(f_ex[["pal"]])[, .(hs_entrez = Entrez, palmer_dir = direction, palmer_cons = sprintf("%d/%d", NumSucc, NumTrials), palmer_p = PVal, palmer_q = QValue)]
pal <- merge(omr[, .(GENE, hs_entrez)], pal, by = "hs_entrez")[order(GENE, palmer_p)][!duplicated(GENE)][, hs_entrez := NULL]
brain_t <- c("amygdala", "basal ganglia", "cerebellum", "cerebral cortex", "choroid plexus", "hippocampal formation", "hypothalamus", "midbrain", "spinal cord")
hpa <- fread(cmd = paste("zcat", shQuote(f_ex[["hpa"]])), col.names = c("ens", "hs_sym", "tissue", "ntpm"))
hpa_br <- hpa[tissue %in% brain_t, .(brain_top = tissue[which.max(ntpm)]), by = hs_sym]
hpa[, t2 := fifelse(tissue %in% brain_t, "brain", tissue)]
hpa <- hpa[, .(ntpm = max(ntpm)), by = .(hs_sym, t2)]
hpa_assign <- function(t, v) {
  o <- order(-v); t <- t[o]; v <- v[o]; b <- v[t == "brain"]; ob <- v[t != "brain"]; n <- length(v)
  x <- log2(v + 1); tau <- if (max(x) > 0) sum(1 - x / max(x)) / (n - 1) else NA_real_
  cls <- if (!length(b) || b < 1) "뇌 미검출" else if (b >= 4 * max(ob)) "뇌 특이" else {
    grp <- FALSE; for (k in 2:5) if ("brain" %in% t[1:k] && mean(v[1:k]) >= 4 * v[k + 1]) { grp <- TRUE; break }
    if (grp) "뇌 포함 그룹 특이" else if (b >= 4 * mean(v)) "뇌 높음" else "뇌 발현" }
  list(hpa_grp = cls, hpa_tau = tau, hpa_top = t[1], hpa_brain_ntpm = if (length(b)) b else NA_real_)
}
HPA <- hpa[, hpa_assign(t2, ntpm), by = hs_sym]
HPA <- merge(HPA, hpa_br, by = "hs_sym", all.x = TRUE)
cls_rank <- c("뇌 특이" = 1, "뇌 포함 그룹 특이" = 2, "뇌 높음" = 3, "뇌 발현" = 4, "뇌 미검출" = 5)
HPA <- merge(omr[, .(GENE, hs_sym)], HPA, by = "hs_sym")[order(GENE, cls_rank[hpa_grp], -hpa_tau)][!duplicated(GENE)]
ctx <- fread(cmd = paste("zcat", shQuote(f_ex[["ctx"]])))
yc <- grep("^y", names(ctx), value = TRUE); oc <- grep("^o", names(ctx), value = TRUE)
ctx <- ctx[, .(SYMBOL = symbol, dog_ctx_mean = rowMeans(.SD[, c(yc, oc), with = FALSE]),
               dog_ctx_lfc_oy = log2((rowMeans(.SD[, oc, with = FALSE]) + 1) / (rowMeans(.SD[, yc, with = FALSE]) + 1)))][!duplicated(SYMBOL)]
GSPEC <- merge(GSPEC, gl[, .(GENE, SYMBOL)], by = "GENE")
GSPEC <- Reduce(function(a, b) merge(a, b, by = "GENE", all.x = TRUE), list(GSPEC, pal, HPA[, .(GENE, hpa_grp, hpa_tau, hpa_top, hpa_brain_ntpm, brain_top)]))
GSPEC <- merge(GSPEC, ctx, by = "SYMBOL", all.x = TRUE)
GSPEC[is.na(hpa_grp), hpa_grp := "자료 없음"]
GSPEC[, brain_elev := hpa_grp %in% c("뇌 특이", "뇌 포함 그룹 특이", "뇌 높음")]
GSPEC[, palmer_sig := !is.na(palmer_q) & palmer_q < 0.05]
GSPEC[, gtag := paste(c(if (path_spec == "특이") "s", if (brain_elev) "B", if (palmer_sig) paste0("A", fifelse(palmer_dir == "up", "+", "−"))), collapse = ","), by = GENE]
setkey(GSPEC, GENE)
cat(sprintf("[%s] 유전자 특이성: 경로 수 하위 25%% ≤ %d · 상위 25%% ≥ %d (주석 있는 유전자 %d개) · HPA 연결 %d (뇌 발현 높음 %d) · Palmer 연결 %d (q < 0.05 %d) · 개 대뇌피질 %d\n",
            format(Sys.time(), "%T"), q_sp[1], q_sp[2], GSPEC[n_sets > 0, .N], GSPEC[hpa_grp != "자료 없음", .N], sum(GSPEC$brain_elev),
            GSPEC[!is.na(palmer_dir), .N], sum(GSPEC$palmer_sig), GSPEC[!is.na(dog_ctx_mean), .N]))
gl_tag <- function(g) { x <- GSPEC[.(g)]; fifelse(x$gtag == "", x$SYMBOL, paste0(x$SYMBOL, "[", x$gtag, "]")) }
LGS[, `:=`(path_spec = GSPEC[.(LGS$GENE), path_spec], brain_elev = GSPEC[.(LGS$GENE), brain_elev], palmer_sig = GSPEC[.(LGS$GENE), palmer_sig],
           gtag = gl_tag(LGS$GENE))]
LGo[, `:=`(path_spec = GSPEC[.(LGo$GENE), path_spec], brain_elev = GSPEC[.(LGo$GENE), brain_elev], palmer_sig = GSPEC[.(LGo$GENE), palmer_sig])]
gsp <- LGo[rule == "range"][order(lid, dist)][, .(
  gene_spec_range = paste(vapply(GENE, function(g) { x <- GSPEC[.(g)]
    sprintf("%s[경로 %d·%s; %s%s; %s]", x$SYMBOL, x$n_sets, x$path_spec, x$hpa_grp, if (!is.na(x$hpa_top) && x$hpa_grp != "자료 없음") paste0("·top ", x$hpa_top) else "",
            if (is.na(x$palmer_dir)) "노화 -" else sprintf("노화 %s %s q%.1g", fifelse(x$palmer_dir == "up", "↑", "↓"), x$palmer_cons, x$palmer_q)) }, ""), collapse = "; "),
  n_genes_spec = sum(path_spec == "특이"), n_genes_brain = sum(brain_elev %in% TRUE), n_genes_palmer = sum(palmer_sig %in% TRUE)), by = lid]
L <- merge(L, gsp, by = "lid", all.x = TRUE)
for (j in c("n_genes_spec", "n_genes_brain", "n_genes_palmer")) set(L, which(is.na(L[[j]])), j, 0L)
L[is.na(gene_spec_range), gene_spec_range := ""]
gi_gs <- match(rownames(GS), GSPEC$GENE)
set_frac_spec <- as.numeric(crossprod(GS, as.numeric(GSPEC$path_spec[gi_gs] %in% "특이"))) / set_size
set_frac_brain <- as.numeric(crossprod(GS, as.numeric(GSPEC$brain_elev[gi_gs] %in% TRUE))) / set_size
names(set_frac_spec) <- names(set_frac_brain) <- colnames(GS)

bim_a <- bimp[chr %in% 1:38]
combos <- CJ(model = unique(L$model), lev = c("1e-4", "1e-5"), sorted = FALSE)
perm_one <- function(ids, H_obs, rule_set) {
  obs <- colSums(H_obs); keep <- which(obs > 0); if (!length(keep)) return(NULL)
  D <- L[match(ids, lid), .(dlo = lo - bp, dhi = hi - bp)]; n <- nrow(D)
  ge <- numeric(length(keep)); sm <- numeric(length(keep)); chunk <- 100L
  for (b0 in seq(1, n_perm, by = chunk)) {
    nb <- min(chunk, n_perm - b0 + 1); pick <- bim_a[sample.int(nrow(bim_a), n * nb, replace = TRUE)]
    P <- data.table(id = seq_len(n * nb), chr = pick$chr, lo = pick$bp + rep(D$dlo, nb), hi = pick$bp + rep(D$dhi, nb))
    mg <- map_genes(P)[rule %in% rule_set]; j <- match(as.character(mg$GENE), rownames(GS)); ok <- !is.na(j)
    M <- sparseMatrix(i = mg$id[ok], j = j[ok], x = 1, dims = c(n * nb, nrow(GS))); M@x[] <- 1
    Hn <- (M %*% GS[, keep, drop = FALSE]); Hn@x[] <- as.numeric(Hn@x > 0)
    PB <- sparseMatrix(i = rep(seq_len(nb), each = n), j = seq_len(n * nb), x = 1, dims = c(nb, n * nb))
    N <- as.matrix(PB %*% Hn)
    ge <- ge + colSums(sweep(N, 2, obs[keep], ">=")); sm <- sm + colSums(N)
  }
  data.table(set = colnames(GS)[keep], exp_ref = sm / n_perm, perm_p_ref = (1 + ge) / (1 + n_perm))
}
MAP <- rbindlist(mclapply(seq_len(nrow(combos)), function(k) {
  m <- combos$model[k]; lv <- combos$lev[k]
  ids <- L[model == m & (lv == "1e-4" | level == "1e-5"), lid]; if (!length(ids)) return(NULL)
  rbindlist(lapply(c("range", "nearest"), function(rl) {
    rule_set <- if (rl == "range") "range" else c("range", "nearest")
    H <- (if (rl == "range") H_rg else H_nr)[as.character(ids), , drop = FALSE]
    obs <- colSums(H); keep <- which(obs > 0); if (!length(keep)) return(NULL)
    pr <- perm_one(ids, H, rule_set)
    lg <- LGS[lid %in% ids & gene_rule %in% rule_set & set %in% colnames(GS)[keep]]
    lg <- merge(lg, L[, .(lid, SNP, lp = p)], by = "lid")
    desc <- lg[order(lp), .(genes = paste(unique(gene), collapse = ";"),
                            loci = paste(unique(sprintf("%s[%s]", SNP, paste(unique(gene), collapse = "/"))), collapse = " "),
                            loci_p = paste(unique(formatC(lp, 1, format = "e")), collapse = ";"), n_1e5 = sum(unique(lid) %in% L[level == "1e-5", lid])), by = set]
    spc <- lg[order(lp), .(n_genes_hit = uniqueN(GENE), n_spec_genes = uniqueN(GENE[path_spec %in% "특이"]), n_loci_spec = uniqueN(lid[path_spec %in% "특이"]),
                           n_brain_genes = uniqueN(GENE[brain_elev %in% TRUE]), n_loci_brain = uniqueN(lid[brain_elev %in% TRUE]),
                           n_palmer_genes = uniqueN(GENE[palmer_sig %in% TRUE]), n_loci_palmer = uniqueN(lid[palmer_sig %in% TRUE]),
                           genes_tag = paste(unique(gtag), collapse = ";")), by = set]
    x <- data.table(model = m, level = lv, rule = rl, n_loci_model = length(ids), set = colnames(GS)[keep], n_loci = as.integer(obs[keep]))
    x <- merge(merge(merge(x, desc, by = "set", all.x = TRUE), spc, by = "set", all.x = TRUE), pr, by = "set", all.x = TRUE)
    x[, `:=`(source = set_src[set], set_name = set_name[set], set_size = as.integer(set_size[set]), family = fam_of[set])]
    x
  }))
}, mc.cores = min(n_cores, nrow(combos)), mc.preschedule = FALSE))
MAP[, gathered := n_loci >= 2]
MAP[, `:=`(spec_gathered = n_loci_spec >= 2, brain_gathered = n_loci_brain >= 2, set_frac_spec_ref = set_frac_spec[set], set_frac_brain_ref = set_frac_brain[set])]
MAP[, model := factor(model, intersect(run_order, unique(model)))]
setorder(MAP, model, level, rule, -n_loci, set_size)
setcolorder(MAP, c("model", "level", "rule", "source", "set", "set_name", "family", "set_size", "n_loci_model", "n_loci", "n_1e5", "gathered", "genes", "loci", "loci_p", "exp_ref", "perm_p_ref",
                  "n_genes_hit", "n_spec_genes", "n_loci_spec", "spec_gathered", "n_brain_genes", "n_loci_brain", "brain_gathered", "n_palmer_genes", "n_loci_palmer",
                  "genes_tag", "set_frac_spec_ref", "set_frac_brain_ref"))
cat(sprintf("[%s] 지도: %d줄 (경로에 locus 2개 이상 모임 = 주 연결 · 1e-4 기준 %d줄)\n", format(Sys.time(), "%T"), nrow(MAP), MAP[rule == "range" & level == "1e-4" & gathered == TRUE, .N]))

L[, model := factor(model, intersect(run_order, unique(model)))]; setorder(L, model, p)
loc_cols <- c("model", "series", "level", "SNP", "chr", "bp", "p", "b", "maf", "region_id", "n_models_region", "spec_grp", "shared", "own_kind", "wt_source",
              "breed_conc", "top_breed", "hap_cand", "n_cs", "lead_pip", "range_src", "n_snp_range", "lo", "hi", "span_kb", "n_genes_range", "genes_range",
              "gene_p_range", "min_gene_p", "nearest_gene", "nearest_kb", "primary_sets", "primary_sets_nearest", "kegg", "n_go", "near_13_1",
              "n_genes_spec", "n_genes_brain", "n_genes_palmer", "gene_spec_range")
fwrite(f3(L[, ..loc_cols]), paste0(out_prefix, "_loci.csv"))
fwrite(f3(merge(LGS, L[, .(lid, SNP, level, p)], by = "lid")[, lid := NULL][order(model, p, gene, source, set)]), paste0(out_prefix, "_locus_gene_set.csv"))
fwrite(f3(MAP), paste0(out_prefix, "_pathway_map.csv"))
GU <- LGo[, .(models = paste(sort(unique(as.character(model))), collapse = " "), n_loci = uniqueN(lid), rule = paste(sort(unique(rule)), collapse = "/")), by = GENE]
GU <- merge(GU, GSPEC[, !"gtag"], by = "GENE")
GU[, `:=`(lab = gl$lab[match(GU$GENE, gl$GENE)], tag = gl_tag(GU$GENE))]
pcols <- colnames(GS)[set_src[colnames(GS)] == "주 가설"]
GU[, primary_sets := vapply(GU$GENE, function(g) { j <- match(as.character(g), rownames(GS)); if (is.na(j)) "" else
  paste(set_name[pcols[which(GS[j, pcols] > 0)]], collapse = " | ") }, "")]
setcolorder(GU, c("GENE", "SYMBOL", "lab", "tag", "models", "n_loci", "rule", "primary_sets"))
fwrite(f3(GU[order(-n_loci, SYMBOL)]), paste0(out_prefix, "_gene_spec.csv"))

dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
PM <- MAP[rule == "range" & level == "1e-4" & source == "주 가설"]
if (nrow(PM)) {
  PM[, lab := paste0(n_loci, "\n", vapply(strsplit(genes, ";"), function(g) { g <- sub("/.*", "", g); paste(c(head(g, 2), if (length(g) > 2) "…"), collapse = "\n") }, ""))]
  PM[, set_lab := factor(set_name, unique(set_name[order(family, set)]))]
  p1 <- ggplot(PM, aes(model, set_lab)) + geom_tile(aes(fill = n_loci), colour = "white") + geom_text(aes(label = lab), size = 1.9, lineheight = 0.8) +
    scale_fill_gradient(low = "#dde6f3", high = "#4a6fa5", name = "loci") + scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "1e-4 loci of each model mapped to the 14 pre-specified sets (genes within locus range +/- 10 kb; first 2 genes shown)") +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(pic_dir, paste0("13-9", tag, "_primary_map.png")), p1, width = 17, height = 9, dpi = 150)
}
if (nrow(PM)) {
  PM[, lab2 := sprintf("%d\ns %d · B %d\nA %d", n_loci, n_loci_spec, n_loci_brain, n_loci_palmer)]
  p3 <- ggplot(PM, aes(model, set_lab)) + geom_tile(aes(fill = n_loci_spec), colour = "white") + geom_text(aes(label = lab2), size = 2, lineheight = 0.8) +
    scale_fill_gradient(low = "#f3eadd", high = "#b0703a", name = "loci with\npathway-specific\ngene (s)") + scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Gene specificity of loci in the 14 pre-specified sets (top: all loci; s = pathway-specific gene, B = brain-elevated in HPA, A = brain-ageing gene, Palmer q < 0.05)") +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(pic_dir, paste0("13-9", tag, "_primary_spec_map.png")), p3, width = 17, height = 9, dpi = 150)
}
KM <- MAP[rule == "range" & level == "1e-4" & source == "KEGG"]
k_top <- KM[, .(mx = max(n_loci)), by = set][mx >= 3, set]
if (length(k_top)) {
  KM <- KM[set %in% k_top]; KM[, set_lab := factor(paste0(set_name, " (", set_size, ")"))]
  p2 <- ggplot(KM, aes(model, set_lab)) + geom_tile(aes(fill = n_loci), colour = "white") + geom_text(aes(label = n_loci), size = 2.4) +
    scale_fill_gradient(low = "#e8efe2", high = "#5b8a3c", name = "loci") + scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "KEGG pathways with >= 3 loci in at least one model (1e-4 loci, range +/- 10 kb)") +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(pic_dir, paste0("13-9", tag, "_kegg_map.png")), p2, width = 14, height = max(4, 0.25 * length(k_top) + 2), dpi = 150)
}

sink(paste0(out_prefix, "_summary.txt"))
cat("13-9 모형별 locus → 유전자 → 경로 지도 (서술용; exp_ref · perm_p_ref = 게놈 위치 순열 참고 열, 판정 · 거르기 · 순서에 안 씀)\n")
cat(sprintf("유전자 연결: 주 = locus SNP 범위 ± %d kb, 보조 = 없으면 가장 가까운 유전자 (≤ %d kb). 경로: 주 가설 14개 · KEGG · GO (유전자 %d ~ %d개). 순열 %d회\n",
            win_kb, near_max_kb, set_min, set_max, n_perm))
cat(sprintf("유전자 표시: s = 경로 특이 (들어 있는 경로 수 ≤ %d, 하위 25%%) · B = 뇌 발현 높음 (HPA 뇌 특이 · 뇌 포함 그룹 특이 · 뇌 높음) · A+ / A− = 사람 뇌 노화 증가 / 감소 (Palmer q < 0.05)\n", q_sp[1]))
cat("지도 괄호: [s locus · B locus · A locus] = 그런 유전자를 가진 locus 수\n")
cat("\n■ 모형별 locus 수\n")
print(L[, .(loci_1e4 = .N, loci_1e5 = sum(level == "1e-5"), gene_in_range = sum(genes_range != ""), nearest_only = sum(genes_range == "" & !is.na(nearest_gene)),
            primary_any = sum(primary_sets != ""), kegg_any = sum(kegg != "")), by = model])
for (m in levels(L$model)) {
  x <- L[model == m]
  cat(sprintf("\n================ %s (%s 계열) — 1e-4 locus %d개 (1e-5 %d개) ================\n", m, x$series[1], nrow(x), sum(x$level == "1e-5")))
  cat("● 1e-5 locus\n")
  y <- x[level == "1e-5"]
  if (nrow(y)) for (i in seq_len(nrow(y))) cat(sprintf("  %s  p %.1e  [%s]  유전자: %s%s\n      주 가설: %s | KEGG: %s\n", y$SNP[i], y$p[i], y$spec_grp[i],
        if (y$genes_range[i] != "") y$gene_p_range[i] else "(범위 안 없음)",
        if (y$genes_range[i] == "" && !is.na(y$nearest_gene[i])) sprintf(" → 가장 가까운 %s (%g kb)", y$nearest_gene[i], y$nearest_kb[i]) else "",
        fifelse(y$primary_sets[i] == "", "-", y$primary_sets[i]), fifelse(y$kegg[i] == "", "-", substr(y$kegg[i], 1, 150))))
  else cat("  없음\n")
  z <- MAP[model == m & level == "1e-4" & rule == "range" & source == "주 가설"]
  cat("● 주 가설 세트에 걸린 1e-4 locus (범위 ± 10 kb)\n")
  if (nrow(z)) for (i in seq_len(nrow(z))) cat(sprintf("  %-34s %d개%s [s %d · B %d · A %d]: %s  (참고 기대 %.2f · 순열 p %.2f)\n", z$set_name[i], z$n_loci[i], if (z$n_1e5[i] > 0) sprintf(" (1e-5 %d)", z$n_1e5[i]) else "",
                                                       z$n_loci_spec[i], z$n_loci_brain[i], z$n_loci_palmer[i], z$genes_tag[i], z$exp_ref[i], z$perm_p_ref[i]))
  else cat("  없음\n")
  zn <- MAP[model == m & level == "1e-4" & rule == "nearest" & source == "주 가설"]
  zn <- zn[n_loci > fcoalesce(z$n_loci[match(zn$set, z$set)], 0L)]
  if (nrow(zn)) cat(sprintf("  (가장 가까운 유전자까지 넣으면 더해지는 것: %s)\n", paste(sprintf("%s %d개 [%s]", zn$set_name, zn$n_loci, zn$genes), collapse = "; ")))
  for (src in c("KEGG", "GO")) {
    w <- MAP[model == m & level == "1e-4" & rule == "range" & source == src & gathered == TRUE]
    cat(sprintf("● %s 경로에 2개 이상 모인 locus (%d개 경로; 상위 %d개 = locus 수 많은 순 → 작은 경로 순)\n", src, nrow(w), min(nrow(w), 8)))
    if (nrow(w)) for (i in seq_len(min(nrow(w), 8))) cat(sprintf("  %-55s (%3d) %d개 [s %d · B %d · A %d]: %s  (참고 %.2f · %.3f)\n", substr(w$set_name[i], 1, 55), w$set_size[i], w$n_loci[i],
                                                                 w$n_loci_spec[i], w$n_loci_brain[i], w$n_loci_palmer[i], w$genes_tag[i], w$exp_ref[i], w$perm_p_ref[i]))
  }
}
cat("\n■ 여러 모형에 걸친 주 가설 세트 (1e-4, 범위 ± 10 kb; 모형별 locus 수)\n")
print(dcast(MAP[rule == "range" & level == "1e-4" & source == "주 가설"], set_name ~ model, value.var = "n_loci", fill = 0L, drop = c(TRUE, FALSE)))
cat("\n■ 1e-5 locus만의 주 가설 세트 (범위 ± 10 kb)\n")
print(MAP[rule == "range" & level == "1e-5" & source == "주 가설", .(model, set_name, n_loci, genes)])
cat("\n■ 유전자 특이성 분포: locus에 연결된 유전자 (범위 ± 10 kb) vs 게놈 전체 (참고)\n")
ru <- unique(LGo[rule == "range", GENE])
print(rbind(data.table(group = "locus 유전자", n = length(ru), GSPEC[.(ru), .(path_spec_frac = mean(path_spec == "특이"), path_broad_frac = mean(path_spec == "범용"),
                       no_annot_frac = mean(path_spec == "주석 없음"), brain_elev_frac = mean(brain_elev), palmer_frac = mean(palmer_sig))]),
            data.table(group = "게놈 전체", n = nrow(GSPEC), GSPEC[, .(path_spec_frac = mean(path_spec == "특이"), path_broad_frac = mean(path_spec == "범용"),
                       no_annot_frac = mean(path_spec == "주석 없음"), brain_elev_frac = mean(brain_elev), palmer_frac = mean(palmer_sig))])), digits = 3)
cat("\n■ 경로 특이 유전자로 locus 2개 이상 모인 경로 (spec_gathered; 1e-4 · 범위 ± 10 kb; 주 가설 전부 + KEGG · GO 모형마다 s locus 많은 순 → 작은 경로 순 5개)\n")
sg <- MAP[rule == "range" & level == "1e-4" & spec_gathered == TRUE][order(model, source != "주 가설", -n_loci_spec, set_size)]
sg[, rk := rowid(model, source)]; sg <- sg[source == "주 가설" | rk <= 5]
if (nrow(sg)) print(sg[, .(model, source, set_name = substr(set_name, 1, 50), set_size, n_loci, n_loci_spec, genes_tag = substr(genes_tag, 1, 90), frac_ref = round(set_frac_spec_ref, 2))], nrows = 500) else cat("  없음\n")
cat("\n■ 뇌 발현 높음 유전자로 locus 2개 이상 모인 경로 (brain_gathered; 같은 기준)\n")
bg <- MAP[rule == "range" & level == "1e-4" & brain_gathered == TRUE][order(model, source != "주 가설", -n_loci_brain, set_size)]
bg[, rk := rowid(model, source)]; bg <- bg[source == "주 가설" | rk <= 5]
if (nrow(bg)) print(bg[, .(model, source, set_name = substr(set_name, 1, 50), set_size, n_loci, n_loci_brain, genes_tag = substr(genes_tag, 1, 90), frac_ref = round(set_frac_brain_ref, 2))], nrows = 500) else cat("  없음\n")
cat("\n■ 주 가설 세트에 걸린 유전자의 특이성 (범위 ± 10 kb, 모든 모형)\n")
pg <- unique(LGS[source == "주 가설" & gene_rule == "range", .(GENE, set_name)])[, .(sets = paste(unique(set_name), collapse = " | ")), by = GENE]
pg <- merge(pg, GU[, .(GENE, SYMBOL, models, n_loci, n_sets, path_spec, hpa_grp, hpa_top, brain_top, palmer_dir, palmer_cons, palmer_q, dog_ctx_mean)], by = "GENE")
print(f3(pg[order(path_spec != "특이", !hpa_grp %in% c("뇌 특이", "뇌 포함 그룹 특이", "뇌 높음"), SYMBOL)][, GENE := NULL][, sets := substr(sets, 1, 60)]), nrows = 500)
sink()
unlink(wd, recursive = TRUE)
cat(sprintf("[%s] 끝 (%.1f분). 저장: %s_* (loci · locus_gene_set · pathway_map · gene_spec · summary) / 그림 13-9%s_*.png\n", format(Sys.time(), "%T"),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), out_prefix, tag))

