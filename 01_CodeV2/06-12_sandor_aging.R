# 06-12_sandor_aging.R - Ageing differential expression from 06-11 and enrichment of CSLB signals. See README.md.
#         04_GWAS_2023/06-11_ref/GCF_000002285.3_CanFam3.1_genomic.gtf.gz · 06-4_ref/gse304125_cortex_normcounts.tsv.gz

library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_mode <- Sys.getenv("TEST", "0") == "1"
tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("06-12", tag))
pic_prefix <- file.path(pic_dir, paste0("06-12", tag))
magma_out <- paste0(out_prefix, "_magma")
q11 <- if (test_mode) file.path(qc_dir, "06-11_test") else file.path(in_dir, "06-11")
star_dir <- paste0(q11, "_star")
gtf_f <- file.path(in_dir, "06-11_ref", "GCF_000002285.3_CanFam3.1_genomic.gtf.gz")
gse304 <- file.path(in_dir, "06-4_ref", "gse304125_cortex_normcounts.tsv.gz")
gse123 <- c(file.path(in_dir, "06-10_age_genes.csv"), file.path(qc_dir, "06-10_age_genes.csv"))
gse123 <- gse123[file.exists(gse123)][1]
magma <- "/data/00_tools/magma_v1.10/magma"
magma_dir <- file.path(in_dir, "06-3_magma")
cpm_min <- 1; fdr_cut <- 0.05
runs_main <- list(`1 (main)` = "main", `2 (age)` = c("int_age", "het_age"), `3 (explore)` = "age8")
runs_w50 <- c(main = "main_w50", int_age = "int_age_w50", het_age = "het_age_w50")
runs_null <- paste0("null", 1:20)
all_runs <- c(unlist(runs_main), runs_w50, runs_null)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(magma_out, showWarnings = FALSE)
for (f in c(magma, gtf_f, gse304, file.path(magma_dir, paste0(all_runs, ".genes.raw")))) if (!file.exists(f)) stop("파일 없음: ", f)

meta <- data.table(sample = paste0("CL_eto", 1:13),
                   group = rep(c("young", "old"), c(6, 7)),
                   age = c(4, 3, 3, 1, 1, 3, 14, 14, 17, 17, 14, 15, 14),
                   sex = c(rep("F", 11), "M", "M"),
                   breed = c("German shepherd dog", "Beagle", "Beagle", "Border collie", "Boxer", "Beagle", "Mixed (German shepherd dog)",
                             "Border collie", "Golden retriever", "Small Munsterlander", "Mixed breed", "Vizsla", "Mixed (border collie-mudi)"),
                   body_kg = c(20.2, 10.3, 11.2, 10.7, 20.1, 7.5, 29.7, 13.2, 22.0, 12.7, 25.4, 25.0, 21.5))
have <- meta$sample[file.exists(file.path(star_dir, meta$sample, "ReadsPerGene.out.tab"))]
if (!test_mode && length(have) != 13) stop("06-11 결과가 13개가 아님 (있는 것: ", paste(have, collapse = ", "), ")")
meta <- meta[sample %in% have]

rpg <- lapply(meta$sample, function(s) fread(file.path(star_dir, s, "ReadsPerGene.out.tab"), header = FALSE, skip = 4,
                                             col.names = c("gene_id", "unstr", "fwd", "rev")))
tot <- rbindlist(lapply(rpg, function(d) d[, .(unstr = sum(unstr), fwd = sum(fwd), rev = sum(rev))]))[, sample := meta$sample]
strand_col <- if (sum(tot$rev) > 3 * sum(tot$fwd)) "rev" else if (sum(tot$fwd) > 3 * sum(tot$rev)) "fwd" else "unstr"
C <- sapply(rpg, `[[`, strand_col); rownames(C) <- rpg[[1]]$gene_id; colnames(C) <- meta$sample
gl <- fread(cmd = sprintf("gzip -dc %s | awk -F'\\t' '$3==\"gene\"' | cut -f9", shQuote(gtf_f)), header = FALSE, sep = "\t", col.names = "attr")
gl[, `:=`(gene_id = sub('.*gene_id "([^"]+)".*', "\\1", attr), GENE = sub('.*db_xref "GeneID:([0-9]+)".*', "\\1", attr))]
gl <- unique(gl[grepl("GeneID:", attr), .(gene_id, GENE)])
magma_genes <- fread(file.path(magma_dir, "main.genes.out"), select = "GENE", colClasses = c(GENE = "character"))$GENE
map0 <- gl[gene_id %in% rownames(C)]
map <- map0[GENE %in% magma_genes]
n_map_magma <- nrow(map)
map <- map[!duplicated(gene_id) & !duplicated(GENE)]
cat(sprintf("유전자 대응: STAR 유전자 %s개 → GeneID 있는 것 %s개 → MAGMA 유전자 %s개 → gene_id·GeneID 중복 제거 뒤 %s개 (%d개 제거)\n",
            format(nrow(C), big.mark = ","), format(nrow(map0), big.mark = ","), format(n_map_magma, big.mark = ","),
            format(nrow(map), big.mark = ","), n_map_magma - nrow(map)))
C <- C[map$gene_id, , drop = FALSE]; rownames(C) <- map$GENE
lib <- colSums(C)
cat(sprintf("Sándor 2022: 표본 %d개, 가닥 = %s (전체 유전자 배정 read 정가닥 %s · 역가닥 %s), MAGMA 유전자 %s개, library 크기 %.1f ~ %.1f백만\n",
            ncol(C), strand_col, format(sum(tot$fwd), big.mark = ","), format(sum(tot$rev), big.mark = ","), format(nrow(C), big.mark = ","),
            min(lib) / 1e6, max(lib) / 1e6))

trigamma_inv <- function(x) { y <- 0.5 + 1 / x
  for (i in 1:50) { tri <- trigamma(y); dif <- tri * (1 - tri / x) / psigamma(y, 2); y <- y + dif; if (max(-dif / y) < 1e-8) break }; y }
de_fit <- function(samples, use_sex = TRUE, label) {
  m <- meta[match(samples, sample)]; Cs <- C[, samples, drop = FALSE]
  cpm <- t(t(Cs) / colSums(Cs)) * 1e6
  min_n <- min(table(m$group))
  keep <- rowSums(cpm >= cpm_min) >= min_n
  Cs <- Cs[keep, , drop = FALSE]
  lg <- log(Cs); gm <- rowMeans(lg); ok <- is.finite(gm)
  sf <- apply(lg[ok, , drop = FALSE] - gm[ok], 2, function(x) exp(median(x)))
  eff <- sf * exp(mean(log(colSums(Cs))))
  Y <- log2(t((t(Cs) + 0.5) / (eff + 1)) * 1e6)
  m[, group := factor(group, levels = c("young", "old"))]
  f <- if (use_sex && uniqueN(m$sex) > 1) ~ group + sex else ~ group
  D <- model.matrix(f, data = m)
  B <- t(qr.coef(qr(D), t(Y))); R <- Y - B %*% t(D)
  d <- nrow(m) - ncol(D); s2 <- rowSums(R^2) / d; v <- solve(crossprod(D))[2, 2]
  A <- rowMeans(Y)
  e <- log(s2) - digamma(d / 2) + log(d / 2)
  tr <- stats::lowess(A, e, f = 0.3); emean <- approx(tr$x, tr$y, xout = A, rule = 2, ties = mean)$y
  evar <- mean((e - emean)^2) * length(e) / (length(e) - 1) - trigamma(d / 2)
  if (evar > 0) { d0 <- 2 * trigamma_inv(evar); s02 <- exp(emean + digamma(d0 / 2) - log(d0 / 2)) } else { d0 <- Inf; s02 <- exp(emean) }
  s2p <- if (is.finite(d0)) (d0 * s02 + d * s2) / (d0 + d) else s02
  tt <- B[, 2] / sqrt(s2p * v)
  out <- data.table(GENE = rownames(Y), mean_expr = A, logFC = B[, 2], age_t = tt, p = 2 * pt(-abs(tt), d + min(d0, 1e6)))
  out[, fdr := p.adjust(p, "BH")]
  out[, `:=`(age_abs = abs(age_t), age_up = as.integer(fdr < fdr_cut & age_t > 0), age_down = as.integer(fdr < fdr_cut & age_t < 0))]
  attr(out, "info") <- data.table(analysis = label, n = nrow(m), n_young = sum(m$group == "young"), n_old = sum(m$group == "old"),
                                  model = paste(deparse(f), collapse = ""), genes = nrow(out), d0 = d0, up = sum(out$age_up), down = sum(out$age_down))
  attr(out, "Y") <- Y
  out[]
}
main_samples <- setdiff(meta$sample, if (test_mode) character(0) else "CL_eto1")
de <- de_fit(main_samples, TRUE, "main (CL_eto1 제외)")
sens <- list(all13 = de_fit(meta$sample, TRUE, "민감도: 13마리 전부"), nosex = de_fit(main_samples, FALSE, "민감도: 성별 공변량 없음"))
de_info <- rbindlist(c(list(attr(de, "info")), lapply(sens, attr, "info")))
print(de_info)
sym <- gl[, .(GENE, symbol = gene_id)][!duplicated(GENE)]
de <- merge(de, sym, by = "GENE", all.x = TRUE)
fwrite(de[order(p)], paste0(out_prefix, "_age_genes.csv"))
fwrite(rbindlist(lapply(names(sens), function(k) sens[[k]][, .(analysis = k, GENE, logFC, age_t, p, fdr)])), paste0(out_prefix, "_age_genes_sens.csv"))
fwrite(cbind(meta, used_in_main = meta$sample %in% main_samples, lib_size = lib[meta$sample]), paste0(out_prefix, "_samples.csv"))

g3 <- fread(cmd = sprintf("gzip -dc %s", shQuote(gse304)))
s3 <- names(g3)[-1]; cnt <- as.matrix(g3[, ..s3]); storage.mode(cnt) <- "double"
mg <- function(k) rowMeans(cnt[, substr(s3, 1, 1) == k, drop = FALSE])
e3 <- data.table(symbol = g3$symbol, m_y = mg("y"), m_o = mg("o"), m_c = mg("c"))
e3 <- e3[(m_y + m_o + m_c) / 3 >= 10][, symbol := toupper(symbol)][, if (.N == 1) .SD, by = symbol][, other := log2((m_o + 0.5) / (m_y + 0.5))]
gate_stats <- function(other, label) {
  x <- merge(de[!is.na(symbol)][, symbol := toupper(symbol)], other[, symbol := toupper(symbol)], by = "symbol")
  if (nrow(x) < 20) return(data.table(check = label, n_shared = nrow(x)))
  ct <- suppressWarnings(cor.test(x$age_t, x$other, method = "spearman", alternative = "greater"))
  sig <- x$fdr < fdr_cut; same <- sign(x$age_t) == sign(x$other)
  k <- sum(same[sig]); nn <- sum(sig); pb <- mean(same[!sig])
  data.table(check = label, n_shared = nrow(x), spearman = unname(ct$estimate), p_spearman = ct$p.value,
             n_sig = nn, same_dir = k, frac_same = k / nn, frac_background = pb,
             p_binom = if (nn > 0) binom.test(k, nn, pb, alternative = "greater")$p.value else NA_real_)
}
gate <- gate_stats(e3[, .(symbol, other)], "G1 GSE304125 old vs young (prefrontal, RNA-seq, 3 vs 3)")
g12 <- if (!is.na(gse123)) fread(gse123, colClasses = c(GENE = "character")) else NULL
if (!is.null(g12)) {
  g12s <- g12[symbol != "", .(symbol, other = age_t)][!duplicated(symbol)]
  gate <- rbind(gate, gate_stats(g12s, "G2 GSE12328 12 y vs 1 y (cerebral cortex, array, 6 vs 6; reference only)"), fill = TRUE)
}
G1 <- isTRUE(gate$p_spearman[1] < 0.05 && gate$p_binom[1] < 0.05)
fwrite(gate, paste0(out_prefix, "_gate.csv"))

cov_main <- paste0(out_prefix, "_covar.txt"); fwrite(de[, .(GENE, age_abs, age_up, age_down, mean_expr)], cov_main, sep = " ")
cov_cons <- NULL
if (!is.null(g12)) {
  cc <- merge(de[, .(GENE, age_t, fdr, mean_expr)], g12[, .(GENE, t12 = age_t, p12 = p)], by = "GENE")
  cc[, `:=`(cons_up = as.integer(fdr < fdr_cut & age_t > 0 & t12 > 0 & p12 < 0.05), cons_down = as.integer(fdr < fdr_cut & age_t < 0 & t12 < 0 & p12 < 0.05))]
  cov_cons <- paste0(out_prefix, "_covar_cons.txt"); fwrite(cc[, .(GENE, cons_up, cons_down, mean_expr)], cov_cons, sep = " ")
  cat(sprintf("합의 노화 유전자 (Sándor FDR < 0.05 & GSE12328 같은 방향 p < 0.05): 증가 %d · 감소 %d (두 자료 공통 배경 %d)\n",
              sum(cc$cons_up), sum(cc$cons_down), nrow(cc)))
}
run_magma <- function(run, covf, vars, tagm) {
  out <- file.path(magma_out, paste0(run, "_", tagm))
  if (!file.exists(paste0(out, ".gsa.out")))
    system2(magma, c("--gene-results", shQuote(file.path(magma_dir, paste0(run, ".genes.raw"))), "--gene-covar", shQuote(covf),
                     "--model", "condition-hide=mean_expr", "direction-covar=pos", paste0("analyse=list,", paste(vars, collapse = ",")),
                     "--out", shQuote(out)), stdout = paste0(out, ".stdout"), stderr = paste0(out, ".stdout"))
  f <- paste0(out, ".gsa.out"); if (!file.exists(f)) stop(run, ": MAGMA 결과 없음 (", out, ".stdout 확인)")
  g <- fread(cmd = sprintf("grep -v '^#' %s", shQuote(f)))
  g[, .(run = run, analysis = tagm, variable = VARIABLE, n_genes = NGENES, beta = BETA, beta_std = BETA_STD, se = SE, p = P)]
}
use_vars <- function(d, v) v[vapply(v, function(k) if (all(d[[k]] %in% 0:1)) sum(d[[k]]) >= 10 else TRUE, TRUE)]
v_main <- use_vars(de, c("age_abs", "age_up", "age_down"))
v_cons <- if (!is.null(cov_cons)) use_vars(cc, c("cons_up", "cons_down")) else character(0)
cat(sprintf("MAGMA 변수: 주 %s / 합의 %s (0/1 변수는 유전자 10개 이상일 때만)\n", paste(v_main, collapse = ", "),
            if (length(v_cons)) paste(v_cons, collapse = ", ") else "없음 (건너뜀)"))
gsa <- rbindlist(c(lapply(all_runs, run_magma, covf = cov_main, vars = v_main, tagm = "sandor"),
                   if (length(v_cons)) lapply(all_runs, run_magma, covf = cov_cons, vars = v_cons, tagm = "consensus")))
fwrite(gsa, paste0(out_prefix, "_gsa_all.csv"))
nul <- gsa[run %in% runs_null]
null_calib <- nul[, .(null_fp05 = mean(p < 0.05), null_min_p = min(p)), by = .(analysis, variable)]
stage_of <- function(st, an) {
  d <- gsa[analysis == an & run %in% runs_main[[st]], .(stage = st, run, variable, n_genes, beta_std, se, p)]
  if (!nrow(d)) return(NULL)
  if (st == "3 (explore)") d[, q := p.adjust(p, "BH"), by = run] else d[, q := p.adjust(p, "BH")]
  w <- gsa[analysis == an & run %in% runs_w50, .(run = names(runs_w50)[match(run, runs_w50)], variable, p_w50 = p)]
  d <- merge(d, w, by = c("run", "variable"), all.x = TRUE)
  d[, p_emp := vapply(seq_len(.N), function(i) (sum(nul[analysis == an & variable == d$variable[i], p] <= d$p[i]) + 1) / (length(runs_null) + 1), 0)]
  d <- merge(d, null_calib[analysis == an, .(variable, null_fp05)], by = "variable")
  d[, holds := p < 0.05 & (is.na(p_w50) & st == "3 (explore)" | !is.na(p_w50) & p_w50 < 0.05)]
  d[, pass := G1 & q < 0.05 & holds & null_fp05 <= 0.10]
  setorder(d, p)[]
}
stages <- lapply(setNames(names(runs_main), names(runs_main)), stage_of, an = "sandor")
for (st in names(stages)) fwrite(stages[[st]], paste0(out_prefix, "_stage", substr(st, 1, 1), ".csv"))
cons_tab <- if (length(v_cons)) rbindlist(lapply(names(runs_main), stage_of, an = "consensus")) else NULL
if (!is.null(cons_tab)) fwrite(cons_tab, paste0(out_prefix, "_sens_cons.csv"))
fwrite(null_calib, paste0(out_prefix, "_null_calib.csv"))

Yall <- attr(sens$all13, "Y"); top <- order(-apply(Yall, 1, var))[1:min(500, nrow(Yall))]
mds <- cmdscale(dist(t(Yall[top, ])), k = 2)
md <- cbind(meta[match(colnames(Yall), sample)], dim1 = mds[, 1], dim2 = mds[, 2])
g0 <- ggplot(md, aes(dim1, dim2, colour = group, label = sub("CL_eto", "", sample))) + geom_point(size = 2.5) + geom_text(nudge_y = 0.15, size = 3, colour = "grey20") +
  scale_colour_manual(values = c(young = "#2166ac", old = "#b2182b")) +
  labs(title = "Sándor 2022 prefrontal cortex: MDS (top 500 variable genes, log2 CPM)", subtitle = "number = CL_eto id (1 excluded in main analysis, as in the paper)") +
  theme_minimal(base_size = 10)
ggsave(paste0(pic_prefix, "_mds.png"), g0, width = 6.5, height = 5, dpi = 300)
xg <- merge(de[!is.na(symbol)][, symbol := toupper(symbol)], e3[, .(symbol, other)], by = "symbol")
g1 <- ggplot(xg, aes(age_t, other)) + geom_point(size = 0.5, alpha = 0.3, colour = "grey40") + geom_smooth(method = "lm", formula = y ~ x, colour = "#b2182b", se = FALSE) +
  labs(x = "Sándor 2022 old vs young (moderated t)", y = "GSE304125 old vs young (log2 FC)",
       title = "Gate G1: does the dog prefrontal ageing signature replicate?",
       subtitle = sprintf("Spearman %.3f (one-sided p %.2g); FDR < 0.05 genes same direction %.1f%% vs background %.1f%% (p %.2g)",
                          gate$spearman[1], gate$p_spearman[1], 100 * gate$frac_same[1], 100 * gate$frac_background[1], gate$p_binom[1])) +
  theme_minimal(base_size = 10)
ggsave(paste0(pic_prefix, "_gate.png"), g1, width = 7, height = 5, dpi = 300)
pd <- gsa[run %in% c(unlist(runs_main), runs_w50)][, `:=`(run = factor(run, levels = c(unlist(runs_main), runs_w50)), variable = paste(analysis, variable, sep = ": "))]
g2 <- ggplot(pd, aes(run, variable, fill = -log10(p))) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(p < 0.01, sprintf("%.1e", p), sprintf("%.2f", p))), size = 3, colour = "grey15") +
  scale_fill_gradient(low = "grey95", high = "#1f4e9c", limits = c(0, 3), oob = scales::squish) +
  labs(x = NULL, y = NULL, title = "Association with dog prefrontal ageing genes (MAGMA gene-property)",
       subtitle = sprintf("gate G1 %s; one-sided, conditioned on mean expression", if (G1) "PASSED" else "FAILED (do not interpret)")) +
  theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(pic_prefix, "_gsa.png"), g2, width = 7.5, height = 4, dpi = 300)

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("06-12 반려견 전전두엽 노화 유전자 (Sándor 2022) → 재현성 관문 (GSE304125) → CSLB 유전자 연관 (MAGMA gene-property)\n")
if (test_mode) cat("*** TEST=1: chr38 · 4마리 시험 결과. 해석하지 말 것 ***\n")
cat("탐색 분석 (06 결과를 본 뒤 추가; 설정·판정은 실행 전 헤더에 고정)\n\n")
cat(sprintf("■ 정량: 가닥 %s, MAGMA 유전자 %s개, library %.1f ~ %.1f백만 read\n", strand_col, format(nrow(C), big.mark = ","), min(lib) / 1e6, max(lib) / 1e6))
cat("\n■ 차등 발현 (노령 vs 어린; limma-trend 방식, BH FDR < 0.05)\n"); print(f3(de_info))
cat("  성별: 어린 개는 모두 암컷, 노령은 암 5 · 수 2. 주 모형 (~ 나이 + 성별)의 나이 효과 = 사실상 어린 암컷 vs 늙은 암컷 (나이 × 성별은 못 봄);\n",
    "  성별 없는 민감도 모형은 수컷 2마리가 나이 효과에 섞임. 품종: 어린 5마리 중 비글 3 → 품종 차이가 나이 효과에 섞일 수 있음 (n 12로 보정 불가)\n", sep = "")
cat(sprintf("  주 분석과 민감도의 나이 t 상관: 13마리 %.3f, 성별 없음 %.3f\n",
            cor(de$age_t, sens$all13$age_t[match(de$GENE, sens$all13$GENE)], use = "complete.obs"),
            cor(de$age_t, sens$nosex$age_t[match(de$GENE, sens$nosex$GENE)], use = "complete.obs")))
cat("  논문 보고: 3,436개 DEG (CL_eto1 제외, 논문의 기준·도구; 개수는 참고)\n")
cat("  상위 10개\n"); print(f3(de[order(p)][1:10, .(symbol, GENE, logFC, age_t, p, fdr, mean_expr)]))
cat(sprintf("\n■ 관문 (G1이 MAGMA 해석 조건): G1 %s\n", if (G1) "통과" else "실패 → MAGMA 결과는 참고 표 (해석하지 않음)")); print(f3(gate))
cat(sprintf("  대응된 유전자 수 n_shared: G1 %s개%s (수천 개보다 훨씬 적으면 기호 대응 문제 → 관문 검정력이 더 낮아짐)\n",
            format(gate$n_shared[1], big.mark = ","), if (nrow(gate) > 1) sprintf(", G2 %s개", format(gate$n_shared[2], big.mark = ",")) else ""))
cat("  G1 실패의 두 해석을 구분할 것: (가) Sándor 노화 signature가 다른 반려견 전전두엽에서 재현되지 않음,\n",
    "  (나) GSE304125가 집단당 3마리라 확인할 검정력이 없음. G2 (GSE12328 6 vs 6)도 실패하면 (가) 쪽 근거가 늘고,\n",
    "  G2만 통과하면 (나) 가능성이 남는다. 어느 경우든 G1 실패면 MAGMA 결과는 해석하지 않는다.\n", sep = "")
cat("\n■ null 보정 (표현형 순열 GWAS 20개)\n"); print(f3(null_calib))
for (st in names(stages)) if (!is.null(stages[[st]])) {
  cat(sprintf("\n■ %s단계 (%s): pass = G1 & q < 0.05 & holds & null_fp05 ≤ 0.10\n", substr(st, 1, 1), paste(runs_main[[st]], collapse = " + ")))
  print(f3(stages[[st]][, .(run, variable, n_genes, beta_std, se, p, q, p_w50, p_emp, null_fp05, holds, pass)]))
}
if (!is.null(cons_tab)) { cat("\n■ 민감도 (사용자 ③): Sándor ∩ GSE12328 합의 노화 유전자 (두 자료 공통 배경)\n")
  print(f3(cons_tab[, .(stage, run, variable, n_genes, beta_std, p, q, p_w50, p_emp, null_fp05, holds, pass)])) }
cat("\n■ 판정\n")
cat(sprintf("  관문 G1: %s\n", if (G1) "통과" else "실패"))
for (st in names(stages)) if (!is.null(stages[[st]])) cat(sprintf("  %s단계: pass %d개 (q < 0.05 %d, 명목 p < 0.05 %d / %d)\n", substr(st, 1, 1),
  sum(stages[[st]]$pass), sum(stages[[st]]$q < 0.05), sum(stages[[st]]$p < 0.05), nrow(stages[[st]])))
cat("  통과해도 \"탐색적 근거; 독립 표본 확인 필요\"로만 씀\n")
sink()
cat(sprintf("저장: %s_samples.csv · _age_genes.csv · _age_genes_sens.csv · _gate.csv · _gsa_all.csv · _stage1~3.csv · _sens_cons.csv · _null_calib.csv · _summary.txt / 그림 %s_mds.png · _gate.png · _gsa.png\n",
            out_prefix, pic_prefix))

