# 13-5_covariate_spec.R - Explains stratum-specific effects by age, body weight, sex and ancestry interactions. See README.md.
if (Sys.getenv("OPENBLAS_NUM_THREADS") != "1") {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"); quit(status = system2(file.path(R.home("bin"), "Rscript"), shQuote(f)))
}
suppressMessages({ library(data.table); library(ggplot2); library(parallel) })
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
pic_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/05_GWAS_2023"
test_mode <- Sys.getenv("TEST", "0") == "1"; tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("13-5", tag))
in13 <- file.path(qc_dir, "13-1"); in133 <- file.path(qc_dir, "13-3"); in134 <- file.path(qc_dir, "13-4")
n_cores <- as.integer(Sys.getenv("N_CORES", "8"))
plink1 <- "/data/00_tools/plink_linux_x86_64_20231018/plink"
bfile <- file.path(d0, "03-3_maf0.01_hwe")
inp <- file.path(d0, "05-1_inputs")
ph_f <- file.path(inp, "int_age_regenie_pheno.txt"); cv_f <- file.path(inp, "int_age_regenie_covar.txt")
loco_f <- file.path(d0, "05-3_regenie", "step1", "int_age_1.loco")
sens_f <- file.path(d0, "00_df_sensitivity.csv"); anc_f <- file.path(d0, "08-1_ancestry_wide.csv")
size_f <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/00_Breed_Count.csv"
for (f in c(ph_f, cv_f, loco_f, sens_f, anc_f, size_f, paste0(bfile, ".bed"), plink1, paste0(in13, "_compare.csv"), paste0(in133, "_loci.csv"), paste0(in134, "_loci_breed.csv")))
  if (!file.exists(f)) stop("파일 없음: ", f)
p_cut <- 0.05; pure_cut <- 0.85; min_n <- 20; min_mac_breed <- 5; min_breeds <- 8; hap_kb <- 500
mods <- c(age = "나이", wt = "체중", sex = "성별", anc = "유전 배경")
age_series <- c("young8", "main", "agefloor6", "agefloor7", "age8", "agefloor9", "age10", "agefloor11", "agefloor12")
wt_series <- c("wt_t1", "wt_t2", "wt_t3")
series_order <- c(age_series, wt_series, "sex_m", "sex_f", "purelike", "mixlike", "int_age", "bin_all", "bin_age8")
series_of <- function(m) fifelse(m %in% age_series, "나이", fifelse(m %in% wt_series, "체중", "기타"))
own_of <- function(m) fifelse(m %in% c(setdiff(age_series, "main"), "int_age", "bin_age8"), "age", fifelse(m %in% wt_series, "wt",
                       fifelse(m %in% c("sex_m", "sex_f"), "sex", fifelse(m %in% c("purelike", "mixlike"), "anc", "none"))))
keep_of <- function(m) if (m %in% c("main", "int_age", "bin_all")) NA_character_ else file.path(inp, paste0(if (m == "bin_age8") "age8" else m, "_keep.txt"))
wd <- tempfile("r13_5_"); dir.create(wd)

ph <- fread(ph_f, colClasses = list(character = 1:2)); cv <- fread(cv_f, colClasses = list(character = c("FID", "IID")))
d <- merge(ph, cv, by = c("FID", "IID"), sort = FALSE)
sens <- fread(sens_f, select = c("IID", "coi"), colClasses = list(character = "IID"))
d <- merge(d, sens, by = "IID", all.x = TRUE, sort = FALSE)
A <- fread(anc_f); A[, IID := as.character(IID)]; M <- as.matrix(A[, -(1:2)])
A[, `:=`(top = apply(M, 1, max), breed = colnames(M)[max.col(M, "first")])]
d <- merge(d, A[, .(IID, top, breed)], by = "IID", all.x = TRUE, sort = FALSE)
sz <- fread(size_f)[, .(breed = gsub("[^a-z]+", "_", tolower(Breed)), size_grp = get(paste0("size_", "class")))]
d <- merge(d, sz, by = "breed", all.x = TRUE, sort = FALSE)
setorder(d, IID)
if (anyNA(d[, .(CSLB_RINT, age_cslb, wt_cov, coi, top)])) stop("결측: 표현형 · 나이 · 체중 · COI · ancestry")
cat(sprintf("자료 %d마리 (순종형 %d · 잡종형 %d)\n", nrow(d), sum(d$top >= pure_cut), sum(d$top < pure_cut)))
if (nrow(d) != 7115 || anyDuplicated(d$IID)) stop("개체 수 · 중복 확인")
d[, `:=`(lwt = log2(wt_cov), lcoi = log(coi), fem = as.numeric(Sex == "Female"))]
d[, `:=`(m_age = age_cslb - mean(age_cslb), m_wt = lwt - mean(lwt), m_sex = fem - mean(fem), m_anc = lcoi - mean(lcoi))]
cov_q <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5)); cov_c <- c("Sex", "Intact", "Vision", "Cohort_BHS")
C0 <- model.matrix(as.formula(paste("~", paste(c(cov_q, cov_c), collapse = " + "))), data = d)
y_res <- residuals(lm.fit(C0, d$CSLB_RINT)); y_res <- (y_res - mean(y_res)) / sd(y_res)
L <- fread(loco_f); ids <- sub("_.*$", "", names(L)[-1])
LM <- as.matrix(L[, -1]); rownames(LM) <- L[[1]]; colnames(LM) <- ids
if (!all(d$IID %in% ids)) stop("LOCO에 없는 개체")
LM <- LM[, d$IID, drop = FALSE]
Cfull <- cbind(C0, lwt = d$lwt, lcoi = d$lcoi)
mixed <- d$top < pure_cut

# 13-1 locus · 13-3 · 13-4
cmp <- fread(paste0(in13, "_compare.csv"))[, .(model, SNP, thr, region_id, chr, bp, b, p)]
U <- unique(cmp[, .(SNP, chr, bp)])[order(chr, bp)]
if (test_mode) U <- U[round(seq(1, .N, length.out = 60))]
cmp <- cmp[SNP %in% U$SNP]

d[, size_pure := fifelse(top >= pure_cut, fifelse(is.na(size_grp), "unknown", size_grp), NA_character_)]
smd_vars <- c(age_cslb = "나이", lwt = "log2 체중", fem = "암컷", lcoi = "log COI", pure = "순종형")
d[, pure := as.numeric(top >= pure_cut)]
sd_all <- sapply(names(smd_vars), function(v) sd(d[[v]]))
CM <- rbindlist(lapply(series_order, function(m) {
  k <- keep_of(m); ii <- if (is.na(k)) d$IID else intersect(fread(k, header = FALSE, colClasses = "character")$V2, d$IID)
  x <- d[IID %in% ii]
  tb <- x[pure == 1, .N, by = breed][order(-N)]; sc <- x[pure == 1, .N, by = size_pure][, pct := 100 * N / sum(N)]
  row <- data.table(model = m, series = series_of(m), n = nrow(x), age_med = median(x$age_cslb), pct_age8 = 100 * mean(x$age_cslb >= 8),
                    pct_age11 = 100 * mean(x$age_cslb >= 11), wt_med = median(x$wt_cov), pct_female = 100 * mean(x$fem), pct_intact = 100 * mean(x$Intact == "intact"),
                    coi_med = median(x$coi), pct_pure = 100 * mean(x$pure), rint_mean = mean(x$CSLB_RINT),
                    top_breeds = paste(head(sprintf("%s %.0f%%", tb$breed, 100 * tb$N / sum(x$pure)), 5), collapse = ", "))
  for (s in c("Toy and Small", "Medium", "Standard", "Large", "Giant")) set(row, j = paste0("pure_", gsub(" ", "", s)), value = if (nrow(sc)) sum(sc[size_pure == s, pct]) else NA_real_)
  for (v in names(smd_vars)) set(row, j = paste0("smd_", v), value = (mean(x[[v]]) - mean(d[[v]])) / sd_all[[v]])
  row
}))

writeLines(U$SNP, file.path(wd, "snps.txt")); fwrite(d[, .(IID, IID)], file.path(wd, "keep.txt"), sep = "\t", col.names = FALSE)
st <- system2(plink1, c("--bfile", bfile, "--dog", "--keep-allele-order", "--keep", file.path(wd, "keep.txt"), "--extract", file.path(wd, "snps.txt"),
                        "--recode", "A", "--threads", min(n_cores, 8), "--out", file.path(wd, "g")), stdout = FALSE, stderr = FALSE)
if (st != 0) stop("plink 유전형 추출 실패: ", file.path(wd, "g.log"))
G <- fread(file.path(wd, "g.raw")); G[, IID := as.character(IID)]
bim <- fread(paste0(bfile, ".bim"), header = FALSE, select = c(2, 5), col.names = c("SNP", "A1"))[SNP %in% U$SNP]
gcols <- names(G)[-(1:6)]; gsnp <- sub("_[^_]+$", "", gcols); gcnt <- sub("^.*_", "", gcols)
if (!identical(gcnt, bim$A1[match(gsnp, bim$SNP)])) stop("유전형 개수 대립유전자가 bim A1과 다름 (--keep-allele-order 확인)")
G <- G[match(d$IID, IID)]; if (!identical(G$IID, d$IID)) stop("유전형 개체 순서")
GM <- as.matrix(G[, ..gcols]); colnames(GM) <- gsnp

hc1 <- function(X, y) {
  qx <- qr(X); b <- qr.coef(qx, y); if (anyNA(b)) return(NULL)
  e <- y - X %*% b; XtXi <- chol2inv(qr.R(qx))
  V <- XtXi %*% crossprod(X * as.vector(e)) %*% XtXi * nrow(X) / (nrow(X) - ncol(X))
  list(b = b, V = V, cn = colnames(X))
}
term_p <- function(f, nm) { k <- match(nm, f$cn); z <- f$b[k] / sqrt(f$V[k, k]); c(beta = f$b[k], se = sqrt(f$V[k, k]), p = 2 * pnorm(-abs(z))) }
fit_snp <- function(s) {
  chr <- U$chr[U$SNP == s]; g <- GM[, s]; ok <- !is.na(g)
  ya <- y_res - LM[as.character(chr), ]
  mm <- as.matrix(d[, .(m_age, m_wt, m_sex, m_anc)]); colnames(mm) <- names(mods)
  out <- list(SNP = s, n = sum(ok), maf = { f <- mean(g[ok]) / 2; min(f, 1 - f) })
  gi <- function(sel) { x <- g[sel] * mm[sel, , drop = FALSE]; colnames(x) <- paste0("gx_", names(mods)); x }
  for (k in names(mods)) {
    X <- cbind(Cfull[ok, ], g = g[ok], gi(ok)[, paste0("gx_", k), drop = FALSE]); f <- hc1(X, ya[ok])
    r <- if (is.null(f)) rep(NA_real_, 3) else term_p(f, paste0("gx_", k)); out[paste0("marg_", k, c("_b", "_se", "_p"))] <- as.list(r)
  }
  for (sub in c("all", "mixed")) {
    sel <- ok & (sub == "all" | mixed); X <- cbind(Cfull[sel, ], g = g[sel], gi(sel)); X <- X[, apply(X, 2, function(v) var(v) > 0) | colnames(X) == "(Intercept)"]
    f <- hc1(X, ya[sel]); pre <- if (sub == "all") "joint_" else "mixed_"
    for (k in names(mods)) { r <- if (is.null(f) || !paste0("gx_", k) %in% f$cn) rep(NA_real_, 3) else term_p(f, paste0("gx_", k))
                             out[paste0(pre, k, c("_b", "_se", "_p"))] <- as.list(r) }
    if (sub == "all" && !is.null(f)) {
      k4 <- match(paste0("gx_", names(mods)), f$cn); b4 <- f$b[k4]; W <- tryCatch(drop(t(b4) %*% solve(f$V[k4, k4], b4)), error = function(e) NA_real_)
      out$joint_4df_p <- pchisq(W, 4, lower.tail = FALSE); out$g_main_b <- f$b[match("g", f$cn)]
    }
  }
  pb <- d[, .I[top >= pure_cut]]; bb <- d$breed[pb]
  bt <- table(bb); brs <- names(bt)[bt >= min_n]
  out$breed <- rbindlist(lapply(brs, function(br) {
    i <- pb[bb == br]; i <- i[!is.na(g[i])]; gg <- g[i]; mac <- min(sum(gg), sum(2 - gg))
    if (mac < min_mac_breed || length(unique(gg)) < 2) return(NULL)
    db <- d[i]; vv <- c("age_cslb", "age_hinge", "pa_score", "owner_age", "lwt", "Sex", "Intact", "Vision", "Cohort_BHS")
    vv <- vv[sapply(vv, function(v) uniqueN(db[[v]]) > 1)]
    X <- model.matrix(as.formula(paste("~", paste(vv, collapse = " + "))), data = db)
    X <- X[, apply(X, 2, function(v) var(v) > 0) | colnames(X) == "(Intercept)", drop = FALSE]; X <- cbind(X, g = gg)
    qx <- qr(X); if (qx$rank < ncol(X)) X <- X[, qx$pivot[seq_len(qx$rank)], drop = FALSE]; if (!"g" %in% colnames(X)) return(NULL)
    qx <- qr(X); b <- qr.coef(qx, ya[i]); e <- ya[i] - X %*% b; s2 <- sum(e^2) / (nrow(X) - ncol(X)); V <- chol2inv(qr.R(qx)) * s2
    k <- match("g", colnames(X)); data.table(SNP = s, breed = br, n = length(i), mac = mac, beta = b[k], se = sqrt(V[k, k]))
  }))
  out
}
cat(sprintf("[%s] 대표 SNP %d개 × (단독 4 + 결합 2 + 품종 45)\n", format(Sys.time(), "%T"), nrow(U)))
FR <- mclapply(U$SNP, function(s) tryCatch(fit_snp(s), error = function(e) e), mc.cores = n_cores)
bad <- which(sapply(FR, inherits, "error")); if (length(bad)) stop("적합 실패: ", U$SNP[bad[1]], " — ", conditionMessage(FR[[bad[1]]]))
IA <- rbindlist(lapply(FR, function(x) as.data.table(x[setdiff(names(x), "breed")])), fill = TRUE)
BE <- rbindlist(lapply(FR, `[[`, "breed"))
bw <- d[top >= pure_cut, .(n_dogs = .N, wt_med = median(wt_cov)), by = breed][n_dogs >= min_n]
BE <- merge(BE, merge(bw, sz, by = "breed", all.x = TRUE), by = "breed")[order(SNP, breed)]
cat(sprintf("[%s] 적합 끝\n", format(Sys.time(), "%T")))

dl_meta <- function(b, se, x) {
  k <- length(b)
  na <- list(k_breeds = k, slope = NA_real_, slope_se = NA_real_, slope_p = NA_real_, Q_resid = NA_real_, Q_resid_p = NA_real_, I2_resid = NA_real_,
             tau2 = NA_real_, Q_total = NA_real_, Q_total_p = NA_real_)
  if (k < min_breeds || var(x) == 0) return(na)
  w <- 1 / se^2; X <- cbind(1, x)
  fe <- function(w) { XtW <- t(X * w); P <- solve(XtW %*% X); list(beta = P %*% XtW %*% b, P = P) }
  f0 <- fe(w); r <- b - X %*% f0$beta; QE <- sum(w * r^2)
  tr <- sum(w) - sum(diag(f0$P %*% t(X * w^2) %*% X))
  tau2 <- max(0, (QE - (k - 2)) / tr)
  f1 <- fe(1 / (se^2 + tau2)); sl <- f1$beta[2]; sl_se <- sqrt(f1$P[2, 2])
  w0 <- 1 / (se^2); mu <- sum(w0 * b) / sum(w0); Qt <- sum(w0 * (b - mu)^2)
  list(k_breeds = k, slope = sl, slope_se = sl_se, slope_p = 2 * pnorm(-abs(sl / sl_se)), Q_resid = QE, Q_resid_p = pchisq(QE, k - 2, lower.tail = FALSE),
       I2_resid = max(0, (QE - (k - 2)) / QE), tau2 = tau2, Q_total = Qt, Q_total_p = pchisq(Qt, k - 1, lower.tail = FALSE))
}
BM <- BE[, as.data.table(dl_meta(beta, se, log2(wt_med))), by = SNP]

LC <- merge(cmp, IA, by = "SNP", all.x = TRUE)
LC <- merge(LC, BM, by = "SNP", all.x = TRUE)
L33 <- fread(paste0(in133, "_loci.csv"))[, .(model, SNP, thr, specass, spec_axes, shared)]
L34 <- unique(fread(paste0(in134, "_loci_breed.csv"))[, .(SNP, breed_conc, top_breed, top_breed_far2_kb)], by = "SNP")
L34[, hap_cand := breed_conc %in% TRUE & !is.na(top_breed_far2_kb) & top_breed_far2_kb >= hap_kb]
LC <- merge(LC, L33, by = c("model", "SNP", "thr"), all.x = TRUE); LC <- merge(LC, L34, by = "SNP", all.x = TRUE)
LC[, `:=`(series = series_of(model), own = own_of(model), sel_bias = !model %in% c("main", "bin_all"))]
jp <- as.matrix(LC[, paste0("joint_", names(mods), "_p"), with = FALSE]); colnames(jp) <- names(mods)
LC[, joint_axes := apply(jp, 1, function(v) { h <- names(mods)[!is.na(v) & v < p_cut]; if (length(h)) paste(mods[h], collapse = "+") else "없음" })]
LC[, `:=`(own_marg_p = NA_real_, own_joint_p = NA_real_)]
for (ax in names(mods)) LC[own == ax, `:=`(own_marg_p = get(paste0("marg_", ax, "_p")), own_joint_p = get(paste0("joint_", ax, "_p")))]
other_hit <- function(own, i) { if (own == "none") return(""); o <- setdiff(names(mods), own); h <- o[!is.na(jp[i, o]) & jp[i, o] < p_cut]; paste(mods[h], collapse = "+") }
LC[, moved_to := mapply(other_hit, own, seq_len(.N))]
LC[, own_kind := fifelse(own == "none", "전체 표본 모형 (자기 축 없음)",
                 fifelse(is.na(own_marg_p), "계산 안 됨",
                 fifelse(own_marg_p >= p_cut, "자기 축 상호작용 없음",
                 fifelse(!is.na(own_joint_p) & own_joint_p < p_cut, "자기 축 유지",
                 fifelse(nzchar(moved_to), "다른 축으로 옮겨감", "구분 안 됨")))))]
LC[, wt_source := fifelse(is.na(marg_wt_p) | marg_wt_p >= p_cut, NA_character_,
                  fifelse(!is.na(slope_p) & slope_p < p_cut & !is.na(mixed_wt_p) & mixed_wt_p < p_cut, "둘 다",
                  fifelse(!is.na(slope_p) & slope_p < p_cut, "품종 체격 쪽",
                  fifelse(!is.na(mixed_wt_p) & mixed_wt_p < p_cut, "개체 체중 쪽", "구분 안 됨"))))]
LC[, model := factor(model, series_order)]; setorder(LC, model, -thr, chr, bp); LC[, model := as.character(model)]
summ <- function(x) x[, c(list(n_loci = .N, ref_chance = round(p_cut * .N, 1)),
                         setNames(lapply(c("자기 축 상호작용 없음", "자기 축 유지", "다른 축으로 옮겨감", "구분 안 됨"), function(v) sum(own_kind == v)),
                                  c("own_none", "own_kept", "own_moved", "own_mixed")),
                         list(moved_to_top = { t <- table(moved_to[own_kind == "다른 축으로 옮겨감"]); if (length(t)) paste(head(names(sort(t, decreasing = TRUE)), 2), collapse = ",") else "" },
                              wt_marg = sum(marg_wt_p < p_cut, na.rm = TRUE), wt_breed = sum(wt_source %in% "품종 체격 쪽"), wt_indiv = sum(wt_source %in% "개체 체중 쪽"),
                              wt_both = sum(wt_source %in% "둘 다"), joint4_p05 = sum(joint_4df_p < p_cut, na.rm = TRUE)))]
BYM <- rbind(LC[, summ(.SD), by = .(model, thr, shared)][, set := "전체"],
             LC[hap_cand %in% FALSE & breed_conc %in% FALSE, summ(.SD), by = .(model, thr, shared)][, set := "품종 집중 제외"])
BYM[, series := series_of(model)]; BYM[, model := factor(model, series_order)]; setorder(BYM, set, -thr, model, shared); BYM[, model := as.character(model)]

cmL <- melt(CM[, c("model", paste0("smd_", names(smd_vars))), with = FALSE], id.vars = "model")
cmL[, `:=`(model = factor(model, series_order), variable = factor(smd_vars[sub("^smd_", "", variable)], smd_vars))]
g1 <- ggplot(cmL, aes(variable, model, fill = value)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", value)), size = 2.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "SMD vs\nall dogs") + scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL, title = "Sample composition of each model (standardized difference from the full sample)") + theme_minimal(base_size = 9)
ggsave(file.path(pic_dir, paste0("13-5", tag, "_composition.png")), g1, width = 7, height = 6.5, dpi = 200)
ok_lev <- c("자기 축 유지", "다른 축으로 옮겨감", "구분 안 됨", "자기 축 상호작용 없음")
en_lev <- c("own axis kept", "moved to another axis", "not separable", "no own-axis interaction")
g2d <- LC[own != "none" & own_kind %in% ok_lev, .N, by = .(model, thr, shared = fifelse(shared %in% "공유 안 됨", "not shared", "shared / unclear"), own_kind)]
g2d[, `:=`(model = factor(model, series_order), thr = factor(format(thr, scientific = TRUE), c("1e-05", "1e-04")), own_kind = factor(en_lev[match(own_kind, ok_lev)], en_lev))]
g2 <- ggplot(g2d, aes(model, N, fill = own_kind)) + geom_col(position = "fill") + facet_grid(thr ~ shared) +
  scale_fill_manual(values = c("#B2182B", "#E17C05", "#BBBBBB", "#46ACC8"), name = NULL) + scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "share of loci", title = "Does the discovery-stratum interaction survive adjustment for the other three axes? (joint SNP x age/weight/sex/COI)",
       subtitle = "own axis = what the discovery model was stratified on; 'kept' is inflated by selection in that stratum and is not interpreted") +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom")
ggsave(file.path(pic_dir, paste0("13-5", tag, "_own_axis.png")), g2, width = 11, height = 7, dpi = 200)
jt <- melt(LC[, c("model", "thr", paste0("joint_", names(mods), "_p")), with = FALSE], id.vars = c("model", "thr"))[, .(pct = 100 * mean(value < p_cut, na.rm = TRUE)), by = .(model, thr, variable)]
jt[, `:=`(model = factor(model, series_order), axis = factor(c(age = "age", wt = "weight", sex = "sex", anc = "COI")[sub("joint_(.*)_p", "\\1", variable)], c("age", "weight", "sex", "COI")),
          thr = factor(format(thr, scientific = TRUE), c("1e-05", "1e-04")))]
g3 <- ggplot(jt, aes(model, pct, colour = axis, group = axis)) + geom_point() + geom_line() + facet_wrap(~thr, ncol = 1) +
  geom_hline(yintercept = 100 * p_cut, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(length(age_series), length(age_series) + length(wt_series)) + 0.5, linetype = "dotted", colour = "grey50") +
  labs(x = "age series | weight tertiles | others", y = "% loci with joint interaction p < 0.05", title = "Which axis carries each model's loci (joint model, full sample)",
       subtitle = "dashed = 5% (reference only)") + theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom")
ggsave(file.path(pic_dir, paste0("13-5", tag, "_joint_terms.png")), g3, width = 10, height = 7, dpi = 200)
ws <- unique(LC[!is.na(wt_source), .(SNP, wt_source, slope_p, mixed_wt_p, marg_wt_b, slope)], by = "SNP")
g4 <- ggplot(ws, aes(-log10(mixed_wt_p), -log10(slope_p), colour = wt_source)) + geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(p_cut), linetype = "dashed") + geom_vline(xintercept = -log10(p_cut), linetype = "dashed") +
  labs(x = "-log10 p, SNP x weight within mixed-breed dogs (individual weight)", y = "-log10 p, breed effect ~ breed median weight (breed size)",
       title = "Loci with a SNP x weight interaction: breed size or individual weight?", colour = NULL) + theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(pic_dir, paste0("13-5", tag, "_weight_source.png")), g4, width = 7, height = 6.5, dpi = 200)

f3 <- function(x) { x <- copy(x); for (j in names(x)) if (is.numeric(x[[j]]) && !is.integer(x[[j]])) set(x, j = j, value = signif(x[[j]], 3)); x }
fwrite(f3(CM), paste0(out_prefix, "_composition.csv")); fwrite(f3(IA), paste0(out_prefix, "_interaction.csv"))
fwrite(f3(BE), paste0(out_prefix, "_breed_effects.csv")); fwrite(f3(BM), paste0(out_prefix, "_breed_meta.csv"))
fwrite(f3(LC), paste0(out_prefix, "_loci_cov.csv")); fwrite(f3(BYM), paste0(out_prefix, "_by_model.csv"))
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("13-5 모형 특이성을 체중 · 품종 · 공변량으로 설명 (서술; 판정은 locus마다 명목 p < 0.05, 우연 기대 = 참고 열)\n"); if (test_mode) cat("*** TEST=1 시험. 해석하지 말 것 ***\n")
cat(sprintf("자료 %d마리 · 대표 SNP %d개 · locus %d개 · 품종 %d개 (순종형 ≥ %.2f, ≥ %d마리)\n", nrow(d), nrow(U), nrow(LC), nrow(bw), pure_cut, min_n))
cat("\n■ B 모형 구성 (smd_ = 전체 표본 대비 표준화 차이; pure_ = 그 모형 순종형의 체격 분류 %)\n")
print(f3(CM[, !c("top_breeds")]), nrows = 50); cat("상위 품종 (순종형 안 비율):\n"); print(CM[, .(model, top_breeds)], nrows = 50)
cat("\n■ C 결합 상호작용 — 대표 SNP 단위 (결합 모형에서 p < 0.05인 축 비율; 참고: 우연이면 축마다 5%)\n")
print(f3(IA[, .(n = .N, age = 100 * mean(joint_age_p < p_cut, na.rm = TRUE), weight = 100 * mean(joint_wt_p < p_cut, na.rm = TRUE),
                sex = 100 * mean(joint_sex_p < p_cut, na.rm = TRUE), coi = 100 * mean(joint_anc_p < p_cut, na.rm = TRUE),
                any_4df = 100 * mean(joint_4df_p < p_cut, na.rm = TRUE), weight_mixed = 100 * mean(mixed_wt_p < p_cut, na.rm = TRUE))]))
cat("단독 → 결합에서 각 축 상호작용 β가 남는 정도 (단독 p < 0.05인 SNP; β 비 중앙값 = 결합 β / 단독 β):\n")
print(f3(rbindlist(lapply(names(mods), function(k) IA[get(paste0("marg_", k, "_p")) < p_cut, .(axis = mods[[k]], n = .N,
      kept_p05 = sum(get(paste0("joint_", k, "_p")) < p_cut, na.rm = TRUE), beta_ratio_med = median(get(paste0("joint_", k, "_b")) / get(paste0("marg_", k, "_b")), na.rm = TRUE))]))))
cat("\n■ 자기 축 판정 (모형 × 기준 × 공유 여부; own_kept는 선택 편향으로 커져 해석 안 함 · own_moved = 다른 축으로 옮겨감 · moved_to_top = 가장 많이 옮겨간 축)\n")
for (s in c("전체", "품종 집중 제외")) for (sr in c("나이", "체중", "기타")) {
  cat(sprintf("[%s · %s 계열]\n", s, sr)); print(f3(BYM[set == s & series == sr, !c("set", "series")]), nrows = 200)
}
cat("\n■ 13-3 특이성 분류별 자기 축 판정 (1e-4, 층 모형)\n")
print(LC[thr < 5e-4 & thr > 5e-5 & own != "none", .N, by = .(specass, own_kind)][order(specass, -N)], nrows = 100)
cat("\n■ D 품종별 효과 메타회귀 (대표 SNP 단위; 품종 ≥ 8개) — slope = log2 품종 체중 1 증가당 효과 변화, Q_resid = 체격으로 설명 안 되는 품종 차이\n")
print(f3(BM[!is.na(slope_p), .(n_snp = .N, k_breeds_med = median(k_breeds), slope_p05 = sum(slope_p < p_cut), Qresid_p05 = sum(Q_resid_p < p_cut),
                               Qtotal_p05 = sum(Q_total_p < p_cut), I2_resid_med = median(I2_resid))]))
cat("체중 단독 상호작용이 있는 locus의 체중 해석 (대표 SNP 단위 · 모형 계열별):\n")
print(unique(LC[!is.na(wt_source), .(SNP, series, wt_source)])[, .N, by = .(series, wt_source)][order(series, -N)])
cat("\n■ E 품종 집중 · 하플로타입 후보와 결합 상호작용 (대표 SNP 단위)\n")
ue <- unique(LC[, .(SNP, breed_conc, hap_cand, joint_wt_p, joint_anc_p, slope_p, Q_resid_p)], by = "SNP")
print(f3(ue[, .(n = .N, weight_p05 = 100 * mean(joint_wt_p < p_cut, na.rm = TRUE), coi_p05 = 100 * mean(joint_anc_p < p_cut, na.rm = TRUE),
                breed_slope_p05 = 100 * mean(slope_p < p_cut, na.rm = TRUE), breed_Qresid_p05 = 100 * mean(Q_resid_p < p_cut, na.rm = TRUE)),
            by = .(group = fifelse(hap_cand %in% TRUE, "하플로타입 후보", fifelse(breed_conc %in% TRUE, "품종 집중 (나머지)", "품종 집중 아님")))]))
sink()
unlink(wd, recursive = TRUE)
cat(sprintf("저장: %s_* (composition · interaction · breed_effects · breed_meta · loci_cov · by_model · summary) / 그림 13-5%s_*.png\n", out_prefix, tag))

