# 14-4_brain_expression.R - Brain-elevated expression (Human Protein Atlas) of the calcium-response set and of genes near loci. See README.md.
#   03_GWAS_Pure/06-1_gene_loc.txt

suppressMessages(library(data.table))
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
g22 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure"
f_in <- c(hpa = file.path(d0, "13-9_ref/rna_tissue_consensus.tsv.gz"), om = file.path(d0, "06-4_ortholog_map.csv"), sa = file.path(d0, "06-4_sets_aux.txt"),
          gs = file.path(d0, "13-9_gene_spec.csv"), loci = file.path(d0, "13-9_loci.csv"), gl = file.path(g22, "06-1_gene_loc.txt"))
for (f in f_in) if (!file.exists(f)) stop("파일 없음: ", f)
ca_set <- "GO_0071277"

brain_t <- c("amygdala", "basal ganglia", "cerebellum", "cerebral cortex", "choroid plexus", "hippocampal formation", "hypothalamus", "midbrain", "spinal cord")
hpa <- fread(cmd = paste("zcat", shQuote(f_in[["hpa"]])), col.names = c("ens", "hs_sym", "tissue", "ntpm"))
hpa[, t2 := fifelse(tissue %in% brain_t, "brain", tissue)]
hpa <- hpa[, .(ntpm = max(ntpm)), by = .(hs_sym, t2)]
hpa_assign <- function(t, v) {
  o <- order(-v); t <- t[o]; v <- v[o]; b <- v[t == "brain"]; ob <- v[t != "brain"]
  if (!length(b) || b < 1) "뇌 미검출" else if (b >= 4 * max(ob)) "뇌 특이" else {
    grp <- FALSE; for (k in 2:5) if ("brain" %in% t[1:k] && mean(v[1:k]) >= 4 * v[k + 1]) { grp <- TRUE; break }
    if (grp) "뇌 포함 그룹 특이" else if (b >= 4 * mean(v)) "뇌 높음" else "뇌 발현" }
}
HPA <- hpa[, .(hpa_grp = hpa_assign(t2, ntpm)), by = hs_sym]
cls_rank <- c("뇌 특이" = 1, "뇌 포함 그룹 특이" = 2, "뇌 높음" = 3, "뇌 발현" = 4, "뇌 미검출" = 5)
om <- fread(f_in[["om"]])[, .(GENE = dog, hs_sym = hs_symbol)]
HPA <- merge(om, HPA, by = "hs_sym")[order(GENE, cls_rank[hpa_grp])][!duplicated(GENE)]
gl <- fread(f_in[["gl"]], header = FALSE, col.names = c("GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"))[CHR %in% as.character(1:38)]
gl <- unique(gl[, .(GENE, SYMBOL, len_kb = (STOP - START + 1) / 1e3)], by = "GENE")
D <- merge(gl, HPA[, .(GENE, hpa_grp)], by = "GENE")
D[, brain := hpa_grp %in% c("뇌 특이", "뇌 포함 그룹 특이", "뇌 높음")]

GS <- fread(f_in[["gs"]])[grepl("range", rule)]
L <- fread(f_in[["loci"]])
sym_to_gene <- function(s) { s <- unique(unlist(strsplit(s[nzchar(s)], ";"))); gl[SYMBOL %in% s, GENE] }
sa <- strsplit(readLines(f_in[["sa"]]), "[\t ]+"); ca_genes <- as.integer(sa[[which(vapply(sa, `[`, "", 1) == ca_set)]][-1])
sets <- list(
  "Locus genes, all models (1e-4; main)" = GS$GENE,
  "Locus genes, all-dogs model only (1e-4)" = GS[grepl("(^| )main( |$)", models), GENE],
  "Locus genes, 1e-5 loci only" = sym_to_gene(L[level == "1e-05", genes_range]),
  "Cellular response to calcium ion (GO:0071277)" = ca_genes)

fit_or <- function(g) {
  d <- copy(D)[, inset := GENE %in% g]
  f0 <- glm(brain ~ inset, binomial, d); f1 <- glm(brain ~ inset + log(len_kb), binomial, d)
  c0 <- summary(f0)$coefficients["insetTRUE", ]; c1 <- summary(f1)$coefficients["insetTRUE", ]
  data.table(n_genes = d[inset == TRUE, .N], pct_brain = 100 * d[inset == TRUE, mean(brain)], pct_brain_rest = 100 * d[inset == FALSE, mean(brain)],
             len_median_kb = d[inset == TRUE, median(len_kb)], len_median_rest_kb = d[inset == FALSE, median(len_kb)],
             OR = exp(c0[1]), p = c0[4], OR_len = exp(c1[1]), OR_len_lo = exp(c1[1] - 1.96 * c1[2]), OR_len_hi = exp(c1[1] + 1.96 * c1[2]), p_len = c1[4])
}
R <- rbindlist(lapply(names(sets), function(k) cbind(comparison = k, fit_or(sets[[k]]))))
fwrite(R[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 3) else v)], file.path(d0, "14-4_brain_expression.csv"))

D[, len_decile := cut(len_kb, quantile(len_kb, 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
D[, locus := GENE %in% sets[[1]]]
B <- D[, .(len_kb_range = sprintf("%.1f-%.1f", min(len_kb), max(len_kb)), n_locus = sum(locus), pct_brain_locus = round(100 * mean(brain[locus]), 1),
           pct_brain_rest = round(100 * mean(brain[!locus]), 1)), by = len_decile][order(len_decile)]
fwrite(B, file.path(d0, "14-4_brain_by_length.csv"))

sink(file.path(d0, "14-4_brain_summary.txt"))
cat("14-4 locus 근처 유전자 · 칼슘 반응 세트의 HPA 뇌 발현 비율 (방법은 코드 머리말; p는 참고 — 유전자가 서로 독립이 아님)\n")
cat(sprintf("배경: HPA에 연결된 개 유전자 %d개 · 뇌 발현 높음 %d개 (%.1f%%)\n\n", nrow(D), sum(D$brain), 100 * mean(D$brain)))
print(R[, .(comparison, n_genes, pct_brain = round(pct_brain, 1), pct_brain_rest = round(pct_brain_rest, 1), len_kb = round(len_median_kb, 1),
            len_rest_kb = round(len_median_rest_kb, 1), OR = round(OR, 2), p = signif(p, 2), OR_len = round(OR_len, 2),
            CI = sprintf("%.2f-%.2f", OR_len_lo, OR_len_hi), p_len = signif(p_len, 2))], row.names = FALSE)
cat("\n■ 유전자 길이 10분위 안 비교 (locus 유전자, 주)\n"); print(B, row.names = FALSE)
sink()
cat(readLines(file.path(d0, "14-4_brain_summary.txt")), sep = "\n")
