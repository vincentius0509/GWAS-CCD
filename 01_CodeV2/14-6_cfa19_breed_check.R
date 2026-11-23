# 14-6_cfa19_breed_check.R - CFA19 SNP x body weight locus: breed adjustment, breed exclusion, neighbouring SNPs. See README.md.
suppressMessages(library(data.table))
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
plink1 <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
snp <- "19:41392616:G:T"; chr <- "19"
inp <- file.path(d0, "05-1_inputs")
f_in <- c(ph = file.path(inp, "int_age_regenie_pheno.txt"), cv = file.path(inp, "int_age_regenie_covar.txt"), loco = file.path(d0, "13-6_regenie/step1/int_wt_1.loco"),
          anc = file.path(d0, "08-1_ancestry_wide.csv"), bed = file.path(d0, "03-3_maf0.01_hwe.bed"))
for (f in c(f_in, plink1)) if (!file.exists(f)) stop("파일 없음: ", f)
pure_cut <- 0.85; min_n <- 20; p_gw <- 5e-8; p_cut <- 0.05; p_lobo <- 1e-4; n_lobo <- 5
wd <- tempfile("r14_6_"); dir.create(wd)

d <- merge(fread(f_in[["ph"]], colClasses = list(character = 1:2)), fread(f_in[["cv"]], colClasses = list(character = c("FID", "IID"))), by = c("FID", "IID"))
A <- fread(f_in[["anc"]]); A[, IID := as.character(IID)]; M <- as.matrix(A[, -(1:2)])
A[, `:=`(top = apply(M, 1, max), breed = colnames(M)[max.col(M, "first")])]
d <- merge(d, A[, .(IID, top, breed)], by = "IID")
writeLines(snp, file.path(wd, "snp.txt"))
st <- system2(plink1, c("--dog", "--bfile", sub("\\.bed$", "", f_in[["bed"]]), "--extract", file.path(wd, "snp.txt"), "--recode", "A", "--out", file.path(wd, "g")), stdout = FALSE, stderr = FALSE)
if (st != 0) stop("plink 추출 실패")
G <- fread(file.path(wd, "g.raw"), colClasses = list(character = c("FID", "IID")))
gcol <- grep("^19:", names(G), value = TRUE); counted <- sub(".*_", "", gcol)
d <- merge(d, G[, .(IID, g = get(gcol))], by = "IID")
setorder(d, IID)
if (nrow(d) != 7115 || anyNA(d[, .(CSLB_RINT, wt_cov, top)])) stop("개체 수 · 결측 확인: ", nrow(d))
d <- d[!is.na(g)]
if (counted != "T") d[, g := 2 - g]
d[, `:=`(lwt = log2(wt_cov), pure = top >= pure_cut)]
bn <- d[pure == TRUE, .N, by = breed][N >= min_n, breed]
d[, breed_fe := fifelse(pure & breed %in% bn, breed, "other")]
d[, wt3 := cut(wt_cov, quantile(wt_cov, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE, labels = c("Light", "Middle", "Heavy"))]
cov_q <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5)); cov_c <- c("Sex", "Intact", "Vision", "Cohort_BHS")
C0 <- model.matrix(as.formula(paste("~", paste(c(cov_q, cov_c), collapse = " + "))), data = d)
y_res <- residuals(lm.fit(C0, d$CSLB_RINT)); y_res <- (y_res - mean(y_res)) / sd(y_res)
L <- fread(f_in[["loco"]]); lr <- L[as.character(FID_IID) == chr]
if (nrow(lr) != 1) stop("LOCO에 CFA19 행 없음")
loco <- unlist(lr[, -1]); names(loco) <- sub("_.*$", "", names(L)[-1])
if (!all(d$IID %in% names(loco))) stop("LOCO에 없는 개체")
d[, y_adj := y_res - loco[IID]]
cat(sprintf("자료 %d마리 (순종형 %d · 잡종형 %d) · T 빈도 %.3f · 품종 고정효과 %d개 + 기타\n", nrow(d), sum(d$pure), sum(!d$pure), mean(d$g) / 2, length(bn)))

fit <- function(x, breed_fe = FALSE, label) {
  f <- paste("y_adj ~", paste(c(cov_q, cov_c, "lwt", if (breed_fe) "factor(breed_fe)", "g", "g:lwt"), collapse = " + "))
  X <- model.matrix(as.formula(f), data = x); X <- X[, colSums(abs(X)) > 0, drop = FALSE]
  qx <- qr(X); X <- X[, qx$pivot[seq_len(qx$rank)], drop = FALSE]
  b <- solve(crossprod(X), crossprod(X, x$y_adj)); e <- as.vector(x$y_adj - X %*% b); n <- nrow(X); k <- ncol(X)
  B <- solve(crossprod(X)); V <- B %*% crossprod(X * e) %*% B * n / (n - k)          # HC1
  ji <- match(c("g:lwt", "lwt:g"), colnames(X)); j <- c(match("g", colnames(X)), ji[!is.na(ji)][1])
  if (anyNA(j)) stop("계수 이름을 못 찾음: ", label); bb <- b[j]; VV <- V[j, j]
  w2 <- as.numeric(t(bb) %*% solve(VV, bb))
  data.table(model = label, n = n, n_T = sum(x$g), b_g = bb[1], se_g = sqrt(VV[1, 1]), b_int = bb[2], se_int = sqrt(VV[2, 2]),
             p_int = 2 * pnorm(-abs(bb[2] / sqrt(VV[2, 2]))), p_2df = pchisq(w2, 2, lower.tail = FALSE))
}
R <- rbind(fit(d, FALSE, "M0 all dogs"), fit(d, TRUE, "S1 all dogs + breed fixed effects"), fit(d[pure == FALSE], FALSE, "S2 mixed-like only"),
           fit(d[pure == TRUE], TRUE, "S3 purebred-like only + breed fixed effects"))

BT <- d[pure == TRUE, .(n = .N, T_count = sum(g), carriers = sum(g > 0), T_freq = sum(g) / (2 * .N), wt_med = median(wt_cov)), by = breed][order(-T_count)]
BT[, `:=`(share_T = T_count / sum(T_count), share_dogs = n / sum(n))]
top_b <- head(BT$breed, n_lobo)
R <- rbind(R, rbindlist(lapply(top_b, function(b) fit(d[!(pure & breed == b)], FALSE, paste0("S4 without ", b)))))
fit_g <- function(x, label) {
  f <- paste("y_adj ~", paste(c(cov_q, cov_c, "lwt", "g"), collapse = " + "))
  m <- lm(as.formula(f), data = x); X <- model.matrix(m); e <- residuals(m); n <- nrow(X); k <- ncol(X)
  B <- solve(crossprod(X)); V <- B %*% crossprod(X * e) %*% B * n / (n - k); j <- match("g", colnames(X))
  data.table(model = label, n = n, n_T = sum(x$g), b_g = coef(m)[["g"]], se_g = sqrt(V[j, j]), p_g = 2 * pnorm(-abs(coef(m)[["g"]] / sqrt(V[j, j]))))
}
S6 <- rbindlist(c(lapply(levels(d$wt3), function(w) fit_g(d[wt3 == w], paste0("All dogs, ", w))),
                  lapply(levels(d$wt3), function(w) fit_g(d[wt3 == w & pure == FALSE], paste0("Mixed-like, ", w)))))

f3 <- function(x) { x <- copy(x); for (j in names(x)) if (is.numeric(x[[j]]) && !is.integer(x[[j]])) set(x, j = j, value = signif(x[[j]], 3)); x }
fwrite(f3(rbind(R, S6, fill = TRUE)), file.path(d0, "14-6_cfa19_breed_check.csv"))
fwrite(f3(BT), file.path(d0, "14-6_cfa19_breeds.csv"))

v1 <- fifelse(R[model == "S1 all dogs + breed fixed effects", p_2df] < p_gw, "품종 평균 차이로 설명되지 않음 (2자유도 p < 5e-8 유지)",
              sprintf("품종 보정 후 약해짐 (2자유도 p %.2g)", R[model == "S1 all dogs + breed fixed effects", p_2df]))
s2 <- R[model == "S2 mixed-like only"]; s0 <- R[model == "M0 all dogs"]
v2 <- fifelse(sign(s2$b_int) == sign(s0$b_int) & s2$p_int < p_cut, "잡종형에서도 재현 (같은 방향, p < 0.05)", "잡종형에서 재현 안 됨")
v3 <- fifelse(all(R[grepl("^S4", model), p_int] < p_lobo), "한 품종이 끌고 가지 않음 (어느 품종을 빼도 상호작용 p < 1e-4)",
              sprintf("특정 품종 의존 가능 (빼면 p ≥ 1e-4: %s)", paste(sub("S4 without ", "", R[grepl("^S4", model) & p_int >= p_lobo, model]), collapse = ", ")))
sink(file.path(d0, "14-6_cfa19_summary.txt"))
cat(sprintf("14-6 CFA19 %s 품종 확인 (판정 규칙은 코드 머리말; T = 소수 대립유전자, 빈도 %.3f)\n", snp, mean(d$g) / 2))
cat(sprintf("판정: ① %s · ② %s · ③ %s\n\n", v1, v2, v3))
cat("■ 모형별 (b_int = T 하나당 · log2 체중 1 단위당 효과 변화; 음수 = 가벼운 개에서 효과가 큼)\n"); print(f3(R), row.names = FALSE)
cat("\n■ 체중 3분위 안 T 효과 (상호작용 항 없이)\n"); print(f3(S6), row.names = FALSE)
cat(sprintf("\n■ 순종형 T 대립유전자 분포 (품종 %d개 중 T를 가진 품종 %d개; 상위 10)\n", nrow(BT), BT[T_count > 0, .N])); print(f3(head(BT, 10)), row.names = FALSE)
sink()
unlink(wd, recursive = TRUE)
cat(readLines(file.path(d0, "14-6_cfa19_summary.txt")), sep = "\n")

