# 14-9_cfa38_breed_check.R - Same checks for the CFA38 locus. See README.md.
suppressMessages(library(data.table))
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
plink1 <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
snp <- "38:23706225:G:A"; chr <- "38"; minor <- "A"
leads <- data.table(locus = c("CFA19", "CFA38"), snp = c("19:41392616:G:T", snp), chr = c(19L, 38L), bp = c(41392616L, 23706225L))
inp <- file.path(d0, "05-1_inputs")
f_in <- c(ph = file.path(inp, "int_age_regenie_pheno.txt"), cv = file.path(inp, "int_age_regenie_covar.txt"), loco = file.path(d0, "13-6_regenie/step1/int_wt_1.loco"),
          reg = file.path(d0, "13-6_regenie/int_wt.regenie.gz"), anc = file.path(d0, "08-1_ancestry_wide.csv"), bed = file.path(d0, "03-3_maf0.01_hwe.bed"))
for (f in c(f_in, plink1)) if (!file.exists(f)) stop("파일 없음: ", f)
bfile <- sub("\\.bed$", "", f_in[["bed"]])
pure_cut <- 0.85; min_n <- 20; p_gw <- 5e-8; p_cut <- 0.05; p_lobo <- 1e-4; n_lobo <- 5
win_kb <- 250; r2_cut <- 0.5; p_nb <- 1e-5
wd <- tempfile("r14_9_"); dir.create(wd)

d <- merge(fread(f_in[["ph"]], colClasses = list(character = 1:2)), fread(f_in[["cv"]], colClasses = list(character = c("FID", "IID"))), by = c("FID", "IID"))
A <- fread(f_in[["anc"]]); A[, IID := as.character(IID)]; M <- as.matrix(A[, -(1:2)])
A[, `:=`(top = apply(M, 1, max), breed = colnames(M)[max.col(M, "first")])]
d <- merge(d, A[, .(IID, top, breed)], by = "IID")
writeLines(snp, file.path(wd, "snp.txt"))
st <- system2(plink1, c("--dog", "--bfile", bfile, "--extract", file.path(wd, "snp.txt"), "--recode", "A", "--out", file.path(wd, "g")), stdout = FALSE, stderr = FALSE)
if (st != 0) stop("plink 추출 실패")
G <- fread(file.path(wd, "g.raw"), colClasses = list(character = c("FID", "IID")))
gcol <- grep(paste0("^", chr, ":"), names(G), value = TRUE); counted <- sub(".*_", "", gcol)
d <- merge(d, G[, .(IID, g = get(gcol))], by = "IID")
setorder(d, IID)
if (nrow(d) != 7115 || anyNA(d[, .(CSLB_RINT, wt_cov, top)])) stop("개체 수 · 결측 확인: ", nrow(d))
n_miss <- sum(is.na(d$g)); d <- d[!is.na(g)]
if (counted != minor) d[, g := 2 - g]
d[, `:=`(lwt = log2(wt_cov), pure = top >= pure_cut)]
bn <- d[pure == TRUE, .N, by = breed][N >= min_n, breed]
d[, breed_fe := fifelse(pure & breed %in% bn, breed, "other")]
d[, wt3 := cut(wt_cov, quantile(wt_cov, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE, labels = c("Light", "Middle", "Heavy"))]
cov_q <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5)); cov_c <- c("Sex", "Intact", "Vision", "Cohort_BHS")
C0 <- model.matrix(as.formula(paste("~", paste(c(cov_q, cov_c), collapse = " + "))), data = d)
y_res <- residuals(lm.fit(C0, d$CSLB_RINT)); y_res <- (y_res - mean(y_res)) / sd(y_res)
L <- fread(f_in[["loco"]]); lr <- L[as.character(FID_IID) == chr]
if (nrow(lr) != 1) stop("LOCO에 CFA", chr, " 행 없음")
loco <- unlist(lr[, -1]); names(loco) <- sub("_.*$", "", names(L)[-1])
if (!all(d$IID %in% names(loco))) stop("LOCO에 없는 개체")
d[, y_adj := y_res - loco[IID]]
cat(sprintf("자료 %d마리 (유전형 결측 %d · 순종형 %d · 잡종형 %d) · %s 빈도 %.3f · 품종 고정효과 %d개 + 기타\n",
            nrow(d), n_miss, sum(d$pure), sum(!d$pure), minor, mean(d$g) / 2, length(bn)))

fit <- function(x, breed_fe = FALSE, label) {
  f <- paste("y_adj ~", paste(c(cov_q, cov_c, "lwt", if (breed_fe) "factor(breed_fe)", "g", "g:lwt"), collapse = " + "))
  X <- model.matrix(as.formula(f), data = x); X <- X[, colSums(abs(X)) > 0, drop = FALSE]
  qx <- qr(X); X <- X[, qx$pivot[seq_len(qx$rank)], drop = FALSE]
  b <- solve(crossprod(X), crossprod(X, x$y_adj)); e <- as.vector(x$y_adj - X %*% b); n <- nrow(X); k <- ncol(X)
  B <- solve(crossprod(X)); V <- B %*% crossprod(X * e) %*% B * n / (n - k)          # HC1
  ji <- match(c("g:lwt", "lwt:g"), colnames(X)); j <- c(match("g", colnames(X)), ji[!is.na(ji)][1])
  if (anyNA(j)) stop("계수 이름을 못 찾음: ", label); bb <- b[j]; VV <- V[j, j]
  w2 <- as.numeric(t(bb) %*% solve(VV, bb))
  data.table(model = label, n = n, n_minor = sum(x$g), b_g = bb[1], se_g = sqrt(VV[1, 1]), b_int = bb[2], se_int = sqrt(VV[2, 2]),
             p_int = 2 * pnorm(-abs(bb[2] / sqrt(VV[2, 2]))), p_2df = pchisq(w2, 2, lower.tail = FALSE))
}
R <- rbind(fit(d, FALSE, "M0 all dogs"), fit(d, TRUE, "S1 all dogs + breed fixed effects"), fit(d[pure == FALSE], FALSE, "S2 mixed-like only"),
           fit(d[pure == TRUE], TRUE, "S3 purebred-like only + breed fixed effects"))
BT <- d[pure == TRUE, .(n = .N, minor_count = sum(g), carriers = sum(g > 0), minor_freq = sum(g) / (2 * .N), wt_med = median(wt_cov)), by = breed][order(-minor_count)]
BT[, `:=`(share_minor = minor_count / sum(minor_count), share_dogs = n / sum(n))]
top_b <- head(BT$breed, n_lobo)
R <- rbind(R, rbindlist(lapply(top_b, function(b) fit(d[!(pure & breed == b)], FALSE, paste0("S4 without ", b)))))
fit_g <- function(x, label) {
  f <- paste("y_adj ~", paste(c(cov_q, cov_c, "lwt", "g"), collapse = " + "))
  m <- lm(as.formula(f), data = x); X <- model.matrix(m); e <- residuals(m); n <- nrow(X); k <- ncol(X)
  B <- solve(crossprod(X)); V <- B %*% crossprod(X * e) %*% B * n / (n - k); j <- match("g", colnames(X))
  data.table(model = label, n = n, n_minor = sum(x$g), b_g = coef(m)[["g"]], se_g = sqrt(V[j, j]), p_g = 2 * pnorm(-abs(coef(m)[["g"]] / sqrt(V[j, j]))))
}
S6 <- rbindlist(c(lapply(levels(d$wt3), function(w) fit_g(d[wt3 == w], paste0("All dogs, ", w))),
                  lapply(levels(d$wt3), function(w) fit_g(d[wt3 == w & pure == FALSE], paste0("Mixed-like, ", w)))))

awk_cond <- paste(sprintf("($1==%d && $2>=%d && $2<=%d)", leads$chr, leads$bp - win_kb * 1000L, leads$bp + win_kb * 1000L), collapse = " || ")
RG <- fread(cmd = sprintf("zcat %s | awk '$8==\"ADD-INT_2DF\" && (%s)'", shQuote(f_in[["reg"]]), awk_cond), header = FALSE,
            col.names = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA"))
if (!all(leads$snp %in% RG$ID)) stop("regenie에 lead SNP 없음")
writeLines(d$IID, file.path(wd, "iid1.txt")); fwrite(d[, .(FID, IID)], file.path(wd, "keep.txt"), sep = "\t", col.names = FALSE)
LD <- rbindlist(lapply(seq_len(nrow(leads)), function(i) {
  o <- file.path(wd, paste0("ld_", leads$locus[i]))
  st <- system2(plink1, c("--dog", "--bfile", bfile, "--chr", leads$chr[i], "--from-bp", leads$bp[i] - win_kb * 1000L, "--to-bp", leads$bp[i] + win_kb * 1000L,
                          "--r2", "--ld-snp", leads$snp[i], "--ld-window-kb", win_kb, "--ld-window", 99999, "--ld-window-r2", 0, "--out", o), stdout = FALSE, stderr = FALSE)
  if (st != 0) stop("plink LD 실패: ", leads$locus[i])
  x <- fread(paste0(o, ".ld"))[, .(ID = SNP_B, r2 = R2)]; x[, locus := leads$locus[i]]; x
}))
st <- system2(plink1, c("--dog", "--bfile", bfile, "--chr", "19,38", "--missing", "--out", file.path(wd, "miss")), stdout = FALSE, stderr = FALSE)
if (st != 0) stop("plink 결측 실패")
MS <- fread(file.path(wd, "miss.lmiss"))[, .(ID = SNP, f_miss = F_MISS)]
NB <- merge(LD, RG[, .(ID, bp = GENPOS, A1FREQ, N, p_2df = 10^-LOG10P)], by = "ID")[MS, on = "ID", nomatch = 0]
NB <- merge(NB, leads[, .(locus, lead_bp = bp)], by = "locus")[, dist_kb := (bp - lead_bp) / 1000][order(locus, -r2)]
SUP <- NB[, .(lead_p_2df = p_2df[dist_kb == 0], lead_f_miss = f_miss[dist_kb == 0], n_snp_window = .N, n_r2_05 = sum(r2 >= r2_cut & dist_kb != 0),
              n_r2_05_p1e5 = sum(r2 >= r2_cut & dist_kb != 0 & p_2df < p_nb), n_r2_05_p1e4 = sum(r2 >= r2_cut & dist_kb != 0 & p_2df < 1e-4),
              best_nb_p = suppressWarnings(min(p_2df[r2 >= r2_cut & dist_kb != 0]))), by = locus]
SUP[, verdict := fifelse(n_r2_05_p1e5 >= 1, "이웃 SNP 지지 있음", "이웃 SNP 지지 없음 (홀로 SNP)")]

f3 <- function(x) { x <- copy(x); for (j in names(x)) if (is.numeric(x[[j]]) && !is.integer(x[[j]])) set(x, j = j, value = signif(x[[j]], 3)); x }
fwrite(f3(rbind(R, S6, fill = TRUE)), file.path(d0, "14-9_cfa38_breed_check.csv"))
fwrite(f3(BT), file.path(d0, "14-9_cfa38_breeds.csv"))
fwrite(f3(NB[, .(locus, ID, bp, dist_kb, r2, A1FREQ, N, f_miss, p_2df)]), file.path(d0, "14-9_ld_support.csv"))

v1 <- fifelse(R[model == "S1 all dogs + breed fixed effects", p_2df] < p_gw, "품종 평균 차이로 설명되지 않음 (2자유도 p < 5e-8 유지)",
              sprintf("품종 보정 후 약해짐 (2자유도 p %.2g)", R[model == "S1 all dogs + breed fixed effects", p_2df]))
s2 <- R[model == "S2 mixed-like only"]; s0 <- R[model == "M0 all dogs"]
v2 <- fifelse(sign(s2$b_int) == sign(s0$b_int) & s2$p_int < p_cut, "잡종형에서도 재현 (같은 방향, p < 0.05)",
              sprintf("잡종형에서 재현 안 됨 (β %.3f, p %.2g)", s2$b_int, s2$p_int))
v3 <- fifelse(all(R[grepl("^S4", model), p_int] < p_lobo), "한 품종이 끌고 가지 않음 (어느 품종을 빼도 상호작용 p < 1e-4)",
              sprintf("특정 품종 의존 가능 (빼면 p ≥ 1e-4: %s)", paste(sub("S4 without ", "", R[grepl("^S4", model) & p_int >= p_lobo, model]), collapse = ", ")))
sink(file.path(d0, "14-9_cfa38_summary.txt"))
cat(sprintf("14-9 CFA38 %s 품종 확인 (판정 규칙은 코드 머리말, 14-6과 같음; %s = 소수 대립유전자, 빈도 %.3f · 유전형 결측 %d마리)\n", snp, minor, mean(d$g) / 2, n_miss))
cat(sprintf("판정: ① %s · ② %s · ③ %s\n\n", v1, v2, v3))
cat("■ 모형별 (b_int = 소수 대립유전자 하나당 · log2 체중 1 단위당 효과 변화; 음수 = 가벼운 개에서 효과가 큼)\n"); print(f3(R), row.names = FALSE)
cat("\n■ 체중 3분위 안 효과 (상호작용 항 없이)\n"); print(f3(S6), row.names = FALSE)
cat(sprintf("\n■ 순종형 %s 대립유전자 분포 (품종 %d개 중 가진 품종 %d개; 상위 10)\n", minor, nrow(BT), BT[minor_count > 0, .N])); print(f3(head(BT, 10)), row.names = FALSE)
cat(sprintf("\n■ 이웃 SNP 지지 (lead ±%d kb, r² ≥ %.1f인 SNP 중 2자유도 p < 1e-5가 1개 이상이면 지지)\n", win_kb, r2_cut)); print(f3(SUP), row.names = FALSE)
cat("\n  r² ≥ 0.5 이웃 상위 (locus별 p 순 5개)\n"); print(f3(NB[r2 >= r2_cut & dist_kb != 0][order(locus, p_2df)][, head(.SD, 5), by = locus][, .(locus, ID, dist_kb, r2, f_miss, p_2df)]), row.names = FALSE)
sink()
unlink(wd, recursive = TRUE)
cat(readLines(file.path(d0, "14-9_cfa38_summary.txt")), sep = "\n")

