# 12-2_dominance_greml.R - GREML with additive and dominance components. See README.md.
#   12-2_greml.csv · 12-2_summary.txt · 12-2_work/

library(data.table)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_mode <- Sys.getenv("TEST", "0") == "1"
tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("12-2", tag)); work_dir <- paste0(out_prefix, "_work")
p121 <- if (test_mode) file.path(qc_dir, "12-1_test") else file.path(qc_dir, "12-1")
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
sets <- list(dense = c(A = paste0(p121, "_dense_A"), D = paste0(p121, "_dense_D.d")),
             pruned = c(A = if (test_mode) paste0(p121, "_pruned_A") else file.path(in_dir, "04-1_grm", "grm_all"), D = paste0(p121, "_pruned_D.d")))
froh_f <- file.path(in_dir, "07-1_froh.csv")
samp <- list(all = c(pheno = file.path(in_dir, "04-2_pheno.txt"), qcovar = file.path(in_dir, "03-5_qcovar_PC5.txt"), covar = file.path(in_dir, "03-5_covar.txt"), keep = NA),
             age8 = setNames(file.path(in_dir, "05-1_inputs", paste0("age8_", c("pheno", "qcovar", "covar", "keep"), ".txt")), c("pheno", "qcovar", "covar", "keep")))
for (f in c(gcta, froh_f, paste0(unlist(sets), ".grm.bin"), na.omit(unlist(samp)))) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
cor_gate <- 0.85; h2_06_6 <- c(h2 = 0.018, se = 0.014)
ids_grm <- fread(paste0(sets$dense["A"], ".grm.id"), header = FALSE, colClasses = "character")
fr <- fread(froh_f, select = c("FID", "IID", "F_ROH"), colClasses = c(FID = "character", IID = "character"))

prep <- function(s, with_f) {
  key <- paste0(s, if (with_f) "_F" else "_noF")
  kk <- file.path(work_dir, paste0("keep_", key, ".txt")); qq <- file.path(work_dir, paste0("qcovar_", key, ".txt"))
  if (!file.exists(qq)) {
    k <- if (is.na(samp[[s]]["keep"])) ids_grm else merge(fread(samp[[s]]["keep"], header = FALSE, colClasses = "character"), ids_grm, by = c("V1", "V2"))
    q <- fread(samp[[s]]["qcovar"], header = FALSE, colClasses = list(character = 1:2))
    if (with_f) { q <- merge(q, fr[, .(V1 = FID, V2 = IID, F_ROH)], by = c("V1", "V2")); k <- k[V2 %in% q$V2] }
    fwrite(k, kk, sep = " ", col.names = FALSE); fwrite(q, qq, sep = " ", col.names = FALSE)
  }
  c(keep = kk, qcovar = qq)
}
greml <- function(name, s, set, ad, with_f = TRUE) {
  out <- file.path(work_dir, name); inp <- prep(s, with_f)
  if (!file.exists(paste0(out, ".hsq"))) {
    g <- if (ad) { mg <- file.path(work_dir, paste0("mgrm_", set, ".txt")); writeLines(sets[[set]], mg); c("--mgrm", mg) } else c("--grm", sets[[set]]["A"])
    system2(gcta, c("--reml", g, "--keep", inp["keep"], "--pheno", samp[[s]]["pheno"], if (s == "all") c("--mpheno", "1"), "--qcovar", inp["qcovar"],
                    "--covar", samp[[s]]["covar"], "--reml-no-lrt", "--reml-no-constrain", "--reml-maxit", "500", "--out", out, "--threads", n_cores),
            stdout = paste0(out, ".stdout"), stderr = paste0(out, ".stdout"))
  }
  if (!file.exists(paste0(out, ".hsq"))) { warning(name, ": GREML 수렴 실패 (", out, ".stdout)"); return(data.table(model = name, converged = FALSE)) }
  h <- fread(paste0(out, ".hsq"), fill = TRUE); g <- function(t) suppressWarnings(as.numeric(h[Source == t, ]))
  if (ad) data.table(model = name, converged = TRUE, n = g("n")[2], logL = g("logL")[2], h2_A = g("V(G1)/Vp")[2], se_A = g("V(G1)/Vp")[3], d2 = g("V(G2)/Vp")[2], se_D = g("V(G2)/Vp")[3])
  else data.table(model = name, converged = TRUE, n = g("n")[2], logL = g("logL")[2], h2_A = g("V(G)/Vp")[2], se_A = g("V(G)/Vp")[3])
}
runs <- data.table(id = c("P_all_dense", "S1_age8_dense", "S2_all_pruned", "S3_all_dense_noF"),
                   s = c("all", "age8", "all", "all"), set = c("dense", "dense", "pruned", "dense"), with_f = c(TRUE, TRUE, TRUE, FALSE))
res <- rbindlist(lapply(seq_len(nrow(runs)), function(i) { r <- runs[i]
  a <- greml(paste0(r$id, "_A"), r$s, r$set, FALSE, r$with_f); ad <- greml(paste0(r$id, "_AD"), r$s, r$set, TRUE, r$with_f)
  o <- data.table(run = r$id, sample = r$s, grm_set = r$set, F_ROH_cov = r$with_f, ad_converged = isTRUE(ad$converged), n = ad$n,
                  h2_A_only = a$h2_A, se_A_only = a$se_A, h2_A_inAD = ad$h2_A, se_A_inAD = ad$se_A, d2 = ad$d2, se_D = ad$se_D)
  if (isTRUE(ad$converged) && isTRUE(a$converged)) {
    lrt <- max(0, 2 * (ad$logL - a$logL)); o[, `:=`(LRT = lrt, p_one = pnorm(-sign(d2) * sqrt(lrt)), p_two = pchisq(lrt, 1, lower.tail = FALSE), d2_upper95 = d2 + 1.96 * se_D)]
  }
  o }), fill = TRUE)

rd_grm <- function(p) { n <- nrow(fread(paste0(p, ".grm.id"), header = FALSE)); readBin(paste0(p, ".grm.bin"), "numeric", n = n * (n + 1) / 2, size = 4) }
offcor <- function(set, s) {
  ida <- fread(paste0(sets[[set]]["A"], ".grm.id"), header = FALSE, colClasses = "character")$V2
  idd <- fread(paste0(sets[[set]]["D"], ".grm.id"), header = FALSE, colClasses = "character")$V2
  keep <- fread(prep(s, TRUE)["keep"], header = FALSE, colClasses = "character")$V2
  ia <- which(ida %in% keep); id <- match(ida[ia], idd); if (anyNA(id)) return(NA_real_)
  pos <- function(i, j) { hi <- pmax(i, j); lo <- pmin(i, j); hi * (hi - 1) / 2 + lo }
  pr <- CJ(a = seq_along(ia), b = seq_along(ia))[a > b]
  if (nrow(pr) > 2e7) pr <- pr[sample(.N, 2e7)]
  xa <- rd_grm(sets[[set]]["A"]); xd <- rd_grm(sets[[set]]["D"])
  cor(xa[pos(ia[pr$a], ia[pr$b])], xd[pos(id[pr$a], id[pr$b])])
}
set.seed(1)
res[, cor_offdiag_A_D := mapply(offcor, grm_set, sample)]
res[, verdict := fcase(!ad_converged, "추정 실패 (AD 모형 수렴 안 됨)",
                       cor_offdiag_A_D > cor_gate, "구별 불가 (두 GRM 상관이 높아 분리 불안정)",
                       d2 > 0 & p_one < 0.05, "dominance 몫이 있다 (탐색적)",
                       d2 <= 0, "dominance 몫 없음 (음수 추정 = 표본 변동)",
                       default = "0과 구별되지 않음")]
fwrite(res, paste0(out_prefix, "_greml.csv"))

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]]) && j != "n") set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("12-2 가설 H9: 상가 + dominance GREML (추정 분석; 해석 규칙은 실행 전 헤더에 고정, 2026-11-14 개정)\n"); if (test_mode) cat("*** TEST=1 시험. 해석하지 말 것 ***\n")
print(f3(res[, .(run, n, h2_A_only, h2_A_inAD, d2, se_D, LRT, p_one, p_two, cor_offdiag_A_D, verdict)]))
P <- res[run == "P_all_dense"]
if (nrow(P)) {
  cat(sprintf("\n■ 판정 (주: 전체 · 전체 SNP GRM · F_ROH 공변량): δ² %.3g ± %.3g, 한쪽 p %.3g → %s", P$d2, P$se_D, P$p_one, P$verdict))
  if (P$verdict == "0과 구별되지 않음") cat(sprintf(" (95%% 상한 %.3g)", P$d2_upper95)); cat("\n")
  cat(sprintf("  상가 몫: 상가 단독 %.3g ± %.3g → dominance를 넣은 모형 %.3g ± %.3g (참고: 06-6 h2_main %.3g ± %.3g, pruning GRM · F_ROH 없음)\n",
              P$h2_A_only, P$se_A_only, P$h2_A_inAD, P$se_A_inAD, h2_06_6["h2"], h2_06_6["se"]))
}
sink()
cat(sprintf("저장: %s_greml.csv · _summary.txt\n", out_prefix))

