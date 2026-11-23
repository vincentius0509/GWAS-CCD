# 13-10_pathway_test.R - Tests locus-nominated pathways after removing nominating loci, in complementary strata. See README.md.
suppressMessages({ library(data.table); library(ggplot2); library(parallel) })
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
pic_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/05_GWAS_2023"
test_mode <- Sys.getenv("TEST", "0") == "1"; tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("13-10", tag))
pre <- function(x) file.path(qc_dir, x)
n_cores <- as.integer(Sys.getenv("N_CORES", "8"))
magma <- "/data/00_tools/magma_v1.10/magma"
g22 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure"
f_in <- c(map = pre("13-9_pathway_map.csv"), loci = pre("13-9_loci.csv"), sa = pre("06-4_sets_aux.txt"), gl = file.path(g22, "06-1_gene_loc.txt"))
for (f in c(f_in, magma)) if (!file.exists(f)) stop("파일 없음: ", f)
res63 <- file.path(d0, "06-3_magma"); res137 <- Sys.getenv("MAGMA13_DIR", pre("13-7_magma"))
old12 <- c("main", "purelike", "mixlike", "int_age", "age8", "bin_age8", "young8", "sex_m", "sex_f", "wt_t1", "wt_t2", "wt_t3")
run_order <- c("young8", "main", "agefloor6", "agefloor7", "age8", "agefloor9", "age10", "agefloor11", "agefloor12", "wt_t1", "wt_t2", "wt_t3",
               "sex_m", "sex_f", "purelike", "mixlike", "int_age", "int_wt", "bin_all", "bin_age8")
set_min <- 10; set_max <- 300; loci_min <- 3; top_k <- 5; x_main <- 100; x_ref <- c(50, 10, 0); min_genes <- 10; p_cut <- 0.05
age_old <- c("age8", "agefloor9", "age10", "agefloor11", "agefloor12")
pair_of <- c(setNames(rep("young8", length(age_old)), age_old), young8 = "age8", sex_m = "sex_f", sex_f = "sex_m",
             purelike = "mixlike", mixlike = "purelike", wt_t1 = "wt_t3", wt_t3 = "wt_t1")
crit_of <- function(m) fifelse(m %in% age_old, "층 특이 (노령 > young8)", fifelse(m %in% names(pair_of), "공통 재현", "짝 없음"))
raw_of <- function(m, w50 = FALSE) file.path(if (m %in% old12) res63 else res137, paste0(m, if (w50) "_w50" else "", ".genes.raw"))
wd <- tempfile("r13_10_"); dir.create(wd)
f3 <- function(x) { x <- copy(x); for (j in names(x)) if (is.numeric(x[[j]]) && !is.integer(x[[j]])) set(x, j = j, value = signif(x[[j]], 3)); x }
t0 <- Sys.time()

M <- fread(f_in[["map"]], colClasses = list(character = "level"))
M <- M[rule == "range" & level == "1e-4" & source %in% c("KEGG", "GO") & set_size >= set_min & set_size <= set_max & n_loci >= loci_min]
if (test_mode) M <- M[model %in% c("main", "age8", "young8", "sex_f")]
M[, loci_snps := vapply(strsplit(loci, " "), function(x) paste(sort(unique(sub("\\[.*", "", x))), collapse = " "), "")]
C <- M[order(model, loci_snps, set_size, set)][!duplicated(paste(model, loci_snps))]
C <- merge(C, M[, .(n_same_loci = .N, same_loci_sets = paste(head(set_name[order(set_size)], 4), collapse = " | ")), by = .(model, loci_snps)], by = c("model", "loci_snps"))
setorder(C, model, -n_loci, set_size, set)
C[, rank := rowid(model)]; C[, tier := fifelse(rank <= top_k, "주 후보", "보조 후보")]
C[, cid := sprintf("c%04d", .I)]
C[, `:=`(pair = fifelse(model %in% names(pair_of), pair_of[model], NA_character_), crit = crit_of(model))]
cat(sprintf("[%s] 후보 %d개 (합치기 전 %d줄) · 주 후보 %d · 보조 %d · 모형 %d개\n", format(Sys.time(), "%T"), nrow(C), nrow(M), sum(C$tier == "주 후보"),
            sum(C$tier == "보조 후보"), uniqueN(C$model)))

read_sets <- function(f) { s <- strsplit(readLines(f), "[\t ]+"); setNames(lapply(s, function(x) as.integer(x[-1])), sapply(s, `[`, 1)) }
SA <- read_sets(f_in[["sa"]])
gl <- fread(f_in[["gl"]], header = FALSE, col.names = c("GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"))[CHR %in% as.character(1:38)][, CHR := as.integer(CHR)]
LO <- fread(f_in[["loci"]])[, .(model, SNP, chr, lo, hi)]
excl_genes <- function(m, snps, w) {
  r <- LO[model == m & SNP %in% snps]
  unique(unlist(lapply(seq_len(nrow(r)), function(i) gl[CHR == r$chr[i] & STOP >= r$lo[i] - w * 1e3 & START <= r$hi[i] + w * 1e3, GENE])))
}
SETS <- rbindlist(lapply(seq_len(nrow(C)), function(i) {
  g <- SA[[C$set[i]]]; snps <- strsplit(C$loci_snps[i], " ")[[1]]
  rbindlist(lapply(c(x_main, x_ref), function(w) { ex <- if (w == 0) integer() else excl_genes(C$model[i], snps, w)
    data.table(cid = C$cid[i], x = w, genes = list(setdiff(g, ex)), n_excl = sum(g %in% ex)) }))
}))
C <- merge(C, dcast(SETS, cid ~ paste0("n_excl_x", x), value.var = "n_excl"), by = "cid")
cat(sprintf("[%s] 제거 ± %d kb 후 남은 세트 유전자 중앙값 %g (원래 %g)\n", format(Sys.time(), "%T"), x_main,
            median(lengths(SETS[x == x_main, genes])), median(C$set_size)))

models <- intersect(run_order, unique(c(C$model, na.omit(C$pair))))
jobs <- rbindlist(lapply(models, function(m) {
  own <- SETS[cid %in% C[model == m, cid]][, nm := paste0(cid, "_x", x)]
  rep <- SETS[cid %in% C[pair == m, cid] & x == x_main][, nm := paste0(cid, "_r", x)]
  S <- rbind(own, rep); if (!nrow(S)) return(NULL)
  sf <- file.path(wd, paste0(m, ".sets")); writeLines(vapply(seq_len(nrow(S)), function(i) paste(c(S$nm[i], S$genes[[i]]), collapse = " "), ""), sf)
  data.table(model = m, w50 = c(FALSE, TRUE), sf = sf)
}))
jobs[, raw := mapply(raw_of, model, w50)]
miss <- jobs[!file.exists(raw)]
if (nrow(miss)) cat("MAGMA 유전자 결과 없음 → 건너뜀:", paste(basename(miss$raw), collapse = " "), "\n")
jobs <- jobs[file.exists(raw)]
jobs[, out := file.path(wd, paste0(model, fifelse(w50, "_w50", "")))]
st <- unlist(mclapply(seq_len(nrow(jobs)), function(i) system2(magma, c("--gene-results", shQuote(jobs$raw[i]), "--set-annot", shQuote(jobs$sf[i]), "--out", shQuote(jobs$out[i])),
                                                              stdout = FALSE, stderr = FALSE), mc.cores = n_cores, mc.preschedule = FALSE))
if (any(st != 0)) stop("MAGMA 실패: ", paste(basename(jobs$out[st != 0]), collapse = ", "))
read_gsa <- function(f) { d <- fread(f, skip = "VARIABLE"); if ("FULL_NAME" %in% names(d)) d[, VARIABLE := FULL_NAME]; d[, .(nm = VARIABLE, NGENES, BETA, SE, P)] }
GA <- rbindlist(lapply(seq_len(nrow(jobs)), function(i) read_gsa(paste0(jobs$out[i], ".gsa.out"))[, `:=`(run = jobs$model[i], w50 = jobs$w50[i])]))
GA[, `:=`(cid = sub("_.*", "", nm), ver = sub("^c[0-9]+_", "", nm))]
cat(sprintf("[%s] MAGMA 세트 분석 %d번 (세트 %d개 결과)\n", format(Sys.time(), "%T"), nrow(jobs), nrow(GA)))

own <- GA[w50 == FALSE & grepl("^x", ver)]
W <- dcast(own, cid ~ ver, value.var = c("NGENES", "BETA", "SE", "P"))
W50 <- GA[w50 == TRUE & ver == paste0("x", x_main), .(cid, P_w50_x100 = P)]
R <- GA[w50 == FALSE & ver == paste0("r", x_main), .(cid, pair_run = run, pair_NGENES = NGENES, pair_BETA = BETA, pair_SE = SE, pair_P = P)]
TT <- Reduce(function(a, b) merge(a, b, by = "cid", all.x = TRUE), list(C, W, W50, R))
mx <- paste0("x", x_main)
TT[, `:=`(n_left = get(paste0("NGENES_", mx)), beta = get(paste0("BETA_", mx)), se = get(paste0("SE_", mx)), p = get(paste0("P_", mx)))]
TT[, testable := !is.na(n_left) & n_left >= min_genes]
TT[, main_pass := testable & !is.na(p) & p < p_cut]
TT[, `:=`(diff_z = (beta - pair_BETA) / sqrt(se^2 + pair_SE^2))][, diff_p := 2 * pnorm(-abs(diff_z))]
TT[, pair_ok := fifelse(crit == "층 특이 (노령 > young8)", !is.na(diff_z) & diff_z > 0 & diff_p < p_cut,
                fifelse(crit == "공통 재현", !is.na(pair_P) & pair_P < p_cut, NA))]
TT[, verdict := fifelse(!testable, "검정 불가 (남은 유전자 < 10)", fifelse(!main_pass, "근거 없음",
                fifelse(pair_ok %in% TRUE, "견고", fifelse(crit == "짝 없음", "경로 수준 근거 (탐색적; 짝 없음)", "경로 수준 근거 (탐색적; 짝 미충족)"))))]
TT[, q_ref_tier := p.adjust(fifelse(testable, p, NA_real_), "BH"), by = tier]
TT[, q_ref_all := p.adjust(fifelse(testable, p, NA_real_), "BH")]
TT[, model := factor(model, intersect(run_order, unique(model)))]
setorder(TT, model, rank)

cand_cols <- c("cid", "model", "rank", "tier", "source", "set", "set_name", "set_size", "n_loci", "n_same_loci", "same_loci_sets", "loci_snps", "genes", "loci",
               "pair", "crit", paste0("n_excl_x", c(x_main, x_ref)))
fwrite(f3(TT[, ..cand_cols]), paste0(out_prefix, "_candidates.csv"))
test_cols <- c("cid", "model", "rank", "tier", "set", "set_name", "set_size", "n_loci", "genes", "n_left", "beta", "se", "p", "verdict", "testable", "main_pass",
               grep("^(NGENES|BETA|P)_x", names(TT), value = TRUE), "P_w50_x100", "pair", "crit", "pair_NGENES", "pair_BETA", "pair_P", "diff_z", "diff_p", "pair_ok", "q_ref_tier", "q_ref_all")
fwrite(f3(TT[, ..test_cols]), paste0(out_prefix, "_tests.csv"))

dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
TP <- TT[tier == "주 후보" & testable == TRUE]
if (nrow(TP)) {
  TP[, lab := fifelse(main_pass, substr(set_name, 1, 28), "")]
  TP[, vcol := factor(fifelse(verdict == "견고", "robust", fifelse(main_pass, "pathway-level (exploratory)", "not supported")),
                      c("robust", "pathway-level (exploratory)", "not supported"))]
  p1 <- ggplot(TP, aes(model, -log10(p))) + geom_hline(yintercept = -log10(p_cut), linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = vcol), size = 2.2, position = position_jitter(width = 0.15, height = 0, seed = 1)) +
    geom_text(aes(label = lab), size = 2, vjust = -0.8, check_overlap = TRUE) + scale_x_discrete(drop = FALSE) +
    scale_colour_manual(values = c("robust" = "#b03a2e", "pathway-level (exploratory)" = "#d68910", "not supported" = "grey60"), name = NULL, drop = FALSE) +
    labs(x = NULL, y = "-log10 p (leave-locus-out, genes within +/-100 kb of nominating loci removed)",
         title = "13-10: top-5 candidate pathways per model (13-9 loci >= 3), MAGMA competitive gene-set test without the nominating loci") +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom")
  ggsave(file.path(pic_dir, paste0("13-10", tag, "_primary.png")), p1, width = 13, height = 6, dpi = 150)
}

sink(paste0(out_prefix, "_summary.txt"))
cat("13-10 경로 수준 검정 (13-9 후보; leave-locus-out MAGMA + 겹치지 않는 짝)\n")
cat(sprintf("후보: KEGG · GO, 유전자 %d ~ %d개, 한 모형 locus ≥ %d, 같은 locus 집합만 합침. 주 후보 = 모형별 상위 %d개 (locus 수 → 작은 경로 순)\n", set_min, set_max, loci_min, top_k))
cat(sprintf("주 판정: 지명 locus ± %d kb 유전자 뺀 세트의 MAGMA 경쟁적 p < %.2f (남은 유전자 ≥ %d). 견고 = 주 판정 + 짝 기대 (나이 층: β 지명 > β young8, 차이 p < 0.05; 그 외: 짝 p < 0.05)\n",
            x_main, p_cut, min_genes))
cat("참고 열 (판정에 안 씀): ± 50 · ± 10 kb 제거, x0 (제거 없음 = 순환), MAGMA ± 50 kb 창, BH q. young8은 나이 짝의 공통 대조군으로 반복 사용 (독립 아님)\n")
cat("\n■ 후보 수 · 판정 (주 후보 / 보조 후보)\n")
print(dcast(TT[, .N, by = .(model, tier, verdict)], model + tier ~ verdict, value.var = "N", fill = 0L))
cat("\n■ 주 후보 결과 (모형별 상위 5개)\n")
print(f3(TT[tier == "주 후보", .(model, rank, set_name = substr(set_name, 1, 45), size = set_size, loci = n_loci, genes = substr(genes, 1, 45), left = n_left,
                                beta, p, p_x50 = get("P_x50"), p_x10 = get("P_x10"), p_x0 = get("P_x0"), p_w50 = P_w50_x100, pair, pair_P, diff_p, verdict)]), nrows = 500)
cat("\n■ 보조 후보 중 주 판정 p < 0.05 (순위 밖; \"보조 경향\"으로만 서술)\n")
b <- TT[tier == "보조 후보" & main_pass == TRUE]
if (nrow(b)) print(f3(b[, .(model, rank, set_name = substr(set_name, 1, 45), size = set_size, loci = n_loci, genes = substr(genes, 1, 50), left = n_left, p, pair, pair_P, diff_p, verdict)]), nrows = 500) else cat("  없음\n")
cat(sprintf("\n(참고) 주 후보 %d개 검정 중 p < 0.05 %d개 · 보조 후보 %d개 중 %d개. 주 후보 BH q < 0.05 %d개 · 전체 BH q < 0.05 %d개\n",
            TT[tier == "주 후보" & testable == TRUE, .N], TT[tier == "주 후보" & main_pass == TRUE, .N], TT[tier == "보조 후보" & testable == TRUE, .N],
            TT[tier == "보조 후보" & main_pass == TRUE, .N], TT[tier == "주 후보" & q_ref_tier < 0.05, .N], TT[q_ref_all < 0.05, .N]))
sink()
unlink(wd, recursive = TRUE)
cat(sprintf("[%s] 끝 (%.1f분). 저장: %s_* (candidates · tests · summary) / 그림 13-10%s_primary.png\n", format(Sys.time(), "%T"),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), out_prefix, tag))

