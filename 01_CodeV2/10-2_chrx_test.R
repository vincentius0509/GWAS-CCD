# 10-2_chrx_test.R - X-chromosome variance test (not testable: almost no heterozygous females after imputation). See README.md.
library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_mode <- Sys.getenv("TEST", "0") == "1"
tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("10-2", tag)); pic_prefix <- file.path(pic_dir, paste0("10-2", tag))
work_dir <- paste0(out_prefix, "_work")
p101 <- if (test_mode) file.path(qc_dir, "10-1_test") else file.path(in_dir, "10-1")
xbf <- paste0(p101, "_chrX"); grmx <- paste0(p101, "_grmX"); qual_f <- paste0(p101, "_quality.csv")
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
grm_a <- file.path(in_dir, "04-1_grm", "grm_all")
inp <- function(run) setNames(file.path(in_dir, "05-1_inputs", paste0(run, "_", c("pheno", "qcovar", "covar", "keep"), ".txt")), c("pheno", "qcovar", "covar", "keep"))
sets <- list(all = c(pheno = file.path(in_dir, "04-2_pheno.txt"), qcovar = file.path(in_dir, "03-5_qcovar_PC5.txt"), covar = file.path(in_dir, "03-5_covar.txt"), keep = NA),
             sex_m = inp("sex_m"), sex_f = inp("sex_f"), age8 = inp("age8"))
for (f in c(paste0(xbf, ".bed"), paste0(grmx, ".grm.bin"), gcta, paste0(grm_a, ".grm.bin"), na.omit(unlist(sets)))) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
xid <- fread(paste0(grmx, ".grm.id"), header = FALSE, colClasses = "character")
xbf39 <- file.path(work_dir, "chrX_as39")
if (!file.exists(paste0(xbf39, ".bim"))) {
  b <- fread(paste0(xbf, ".bim"), header = FALSE, colClasses = list(character = 1)); b[, V1 := "39"]
  fwrite(b, paste0(xbf39, ".bim"), sep = "\t", col.names = FALSE)
  for (e in c(".bed", ".fam")) file.symlink(normalizePath(paste0(xbf, e)), paste0(xbf39, e))
}
kx <- file.path(work_dir, "x_ids.txt"); fwrite(xid, kx, sep = " ", col.names = FALSE)
run_g <- function(args, out) system2(gcta, c(args, "--autosome-num", "38", "--out", out, "--threads", n_cores), stdout = paste0(out, ".stdout"), stderr = paste0(out, ".stdout"))
keep_args <- function(s) { k <- sets[[s]]["keep"]; if (is.na(k)) c("--keep", kx) else { kk <- file.path(work_dir, paste0("keep_", s, ".txt"))
  fwrite(merge(fread(k, header = FALSE, colClasses = "character"), xid, by = c("V1", "V2")), kk, sep = " ", col.names = FALSE); c("--keep", kk) } }
mpheno <- function(s) if (s == "all") c("--mpheno", "1") else character(0)

mg <- file.path(work_dir, "mgrm.txt"); writeLines(c(grm_a, grmx), mg)
greml <- function(name, s, lrt) {
  out <- file.path(work_dir, name)
  if (!file.exists(paste0(out, ".hsq")))
    run_g(c("--reml", "--mgrm", mg, keep_args(s), "--pheno", sets[[s]]["pheno"], mpheno(s), "--qcovar", sets[[s]]["qcovar"], "--covar", sets[[s]]["covar"],
            if (lrt) c("--reml-lrt", "2") else "--reml-no-lrt", "--reml-no-constrain", "--reml-maxit", "200"), out)
  if (!file.exists(paste0(out, ".hsq"))) { warning(name, ": GREML 실패 (", out, ".stdout)"); return(NULL) }
  h <- fread(paste0(out, ".hsq"), fill = TRUE); g <- function(t) suppressWarnings(as.numeric(h[Source == t, ]))
  data.table(model = name, n = g("n")[2], h2_auto = g("V(G1)/Vp")[2], se_auto = g("V(G1)/Vp")[3], h2_X = g("V(G2)/Vp")[2], se_X = g("V(G2)/Vp")[3],
             p_X = if (lrt) g("Pval")[2] else pnorm(g("V(G2)/Vp")[2] / g("V(G2)/Vp")[3], lower.tail = FALSE), test = if (lrt) "LRT" else "Wald")
}
gr <- rbindlist(list(greml("P_all_autoX", "all", TRUE), greml("S4_age8_autoX", "age8", FALSE)), fill = TRUE)
fwrite(gr, paste0(out_prefix, "_greml.csv"))

mlma <- function(s) {
  out <- file.path(work_dir, paste0("mlma_", s))
  if (!file.exists(paste0(out, ".mlma")))
    run_g(c("--mlma", "--bfile", xbf39, "--grm", grm_a, keep_args(s), "--pheno", sets[[s]]["pheno"], mpheno(s), "--qcovar", sets[[s]]["qcovar"], "--covar", sets[[s]]["covar"]), out)
  if (!file.exists(paste0(out, ".mlma"))) { warning(s, ": mlma 실패"); return(NULL) }
  fread(paste0(out, ".mlma"))[, run := s]
}
ml <- rbindlist(lapply(c("all", "sex_m", "sex_f", "age8"), mlma), fill = TRUE)
lam <- function(p) median(qchisq(p, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
msum <- ml[, .(n_snp = .N, lambda = lam(p), min_p = min(p, na.rm = TRUE), top_snp = SNP[which.min(p)], n_5e8 = sum(p < 5e-8, na.rm = TRUE),
               n_1e5 = sum(p < 1e-5, na.rm = TRUE), n_1e4 = sum(p < 1e-4, na.rm = TRUE)), by = run]
fwrite(msum, paste0(out_prefix, "_mlma_summary.csv"))
sm <- ml[run == "sex_m", .(SNP, b_m = b, se_m = se, p_m = p)]; sf <- ml[run == "sex_f", .(SNP, b_f = b, se_f = se, p_f = p)]
het <- merge(sm, sf, by = "SNP")[p_m < 1e-4 | p_f < 1e-4][, z_diff := (b_m - b_f) / sqrt(se_m^2 + se_f^2)][, p_diff := 2 * pnorm(-abs(z_diff))][order(p_diff)]
fwrite(het, paste0(out_prefix, "_sex_het.csv"))

pd <- ml[run == "all"]
if (nrow(pd)) {
  g1 <- ggplot(pd, aes(bp / 1e6, -log10(p))) + geom_point(size = 0.5, colour = "#1f4e9c") + geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "grey40") +
    labs(x = "chrX position (Mb; PAR excluded)", y = "-log10 p", title = "CSLB association on the X chromosome (all dogs, GCTA --mlma)") + theme_minimal(base_size = 10)
  ggsave(paste0(pic_prefix, "_manhattan.png"), g1, width = 8, height = 3.5, dpi = 300)
}
P <- if (nrow(gr)) gr[model == "P_all_autoX"] else data.table(h2_X = NA_real_, se_X = NA_real_, p_X = NA_real_)
if (!nrow(P)) P <- data.table(h2_X = NA_real_, se_X = NA_real_, p_X = NA_real_)
xq <- if (file.exists(qual_f)) fread(qual_f) else data.table(x_quality_ok = NA, female_het_O_over_E = NA)
H7 <- if (isFALSE(xq$x_quality_ok)) NA else isTRUE(P$p_X < 0.05)
f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("10-2 가설 H7: X 염색체와 CSLB (2026-11-13 추가, 판정 규칙은 실행 전 헤더에 고정)\n"); if (test_mode) cat("*** TEST=1 시험. 해석하지 말 것 ***\n")
cat("\n■ GREML (상염색체 GRM + X GRM)\n"); print(f3(gr))
cat("\n■ X 연관 검정 요약 (GCTA --mlma, 상염색체 GRM)\n"); print(f3(msum))
cat(sprintf("\n■ 성별 효과 차이 (수컷 · 암컷 중 하나에서 p < 1e-4인 SNP %d개; 가장 작은 p_diff %s)\n", nrow(het), if (nrow(het)) signif(min(het$p_diff), 3) else "NA"))
if (nrow(het)) print(f3(head(het, 10)))
cat(sprintf("\n■ X 유전형 품질 (10-1): 암컷 관측/기대 이형접합 비 %.3g, x_quality_ok = %s\n", xq$female_het_O_over_E, xq$x_quality_ok))
cat(sprintf("■ 판정 (주 검정 P_all_autoX): V(X)/Vp %.3g ± %.3g, LRT p %.3g → %s\n", P$h2_X, P$se_X, P$p_X,
            if (is.na(H7)) "H7 검정 불가 (자료 품질: 원본 X에 암컷 이형접합이 거의 없음). 위 결과는 탐색용 참고" else if (H7) "H7 지지 (Holm 보정 전)" else "H7 지지 안 됨"))
sink()
cat(sprintf("저장: %s_greml.csv · _mlma_summary.csv · _sex_het.csv · _summary.txt / 그림 %s_manhattan.png\n", out_prefix, pic_prefix))

