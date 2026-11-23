# 11-3_premise_check.R - Checks whether the flagged conditions raise CSLB scores. See README.md.
#   · SurveyData/v2023/DAP_2023_{HLES,AFUS}_health_condition_v1.0.csv · DAP_2023_CSLB_v1.0.csv · DAP_2023_CODEBOOK_v1.1.csv
#   11-3_premise.csv (A) · _multi.csv (B) · _items.csv (C) · _screen_codes.csv · _screen_types.csv (D) · _screen_h2.csv (E) · _summary.txt

library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
sv <- file.path(base_dir, "SurveyData/v2023")
out_prefix <- file.path(qc_dir, "11-3"); pic_prefix <- file.path(pic_dir, "11-3"); work_dir <- paste0(out_prefix, "_work")
flags_f <- c(file.path(in_dir, "11-1_flags.csv"), file.path(qc_dir, "11-1_flags.csv")); flags_f <- flags_f[file.exists(flags_f)][1]
pheno_f <- file.path(in_dir, "04-2_pheno.txt"); qcov_f <- file.path(in_dir, "03-5_qcovar_PC5.txt"); cov_f <- file.path(in_dir, "03-5_covar.txt")
a8_keep <- file.path(in_dir, "05-1_inputs", "age8_keep.txt"); grm <- file.path(in_dir, "04-1_grm", "grm_all")
hles_f <- file.path(sv, "DAP_2023_HLES_health_condition_v1.0.csv"); afus_f <- file.path(sv, "DAP_2023_AFUS_health_condition_v1.0.csv")
cslb_f <- file.path(sv, "DAP_2023_CSLB_v1.0.csv"); cb_f <- file.path(base_dir, "DAP_2023_CODEBOOK_v1.1.csv")
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32")); no_greml <- Sys.getenv("NO_GREML", "0") == "1"; n_rand <- as.integer(Sys.getenv("RAND", "20"))
for (f in c(flags_f, pheno_f, qcov_f, cov_f, a8_keep, hles_f, afus_f, cslb_f, cb_f)) if (is.na(f) || !file.exists(f)) stop("파일 없음: ", f)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
pre_list <- c(101, 102, 115, 201, 205, 1215, 1216, 1210, 1298, 1201, 1202, 801, 1113, 1111, 1118, 1905, 1311, 1301, 1903, 511, 512)
min_dogs <- 20

qn <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5)); cn <- c("Sex", "Intact", "Vision", "Cohort_BHS")
ph <- fread(pheno_f, header = FALSE, colClasses = list(character = 1:2)); setnames(ph, c("FID", "IID", "rint0", "score"))
qc <- fread(qcov_f, header = FALSE, colClasses = list(character = 1:2)); setnames(qc, c("FID", "IID", qn))
cv <- fread(cov_f, header = FALSE, colClasses = "character"); setnames(cv, c("FID", "IID", cn))
fl <- fread(flags_f, colClasses = c(FID = "character", IID = "character"))
cats <- c("sensory_pre", "sens_severe_pre", "neuro_pre", "pain_pre", "systemic_pre")
if (!all(cats %in% names(fl))) stop("11-1 출력이 개정 전 형식: 11-1을 다시 실행 (", flags_f, ")")
D <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")), list(ph, qc, cv, fl[, c("FID", "IID", "excl_primary", cats), with = FALSE]))
D[, age8 := IID %in% fread(a8_keep, header = FALSE, colClasses = "character")$V2]

cs <- fread(cslb_f, integer64 = "double")[, dog_id := as.character(dog_id)][order(dog_id, -as.IDate(cslb_date))][, .SD[1], by = dog_id]
items <- c("cslb_pace", "cslb_stare", "cslb_stuck", "cslb_recognize", "cslb_walk_walls", "cslb_avoid", "cslb_find_food",
           "cslb_pace_6mo", "cslb_stare_6mo", "cslb_defecate_6mo", "cslb_food_6mo", "cslb_recognize_6mo", "cslb_active_6mo")
D <- merge(D, cs[, c("dog_id", "cslb_date", "cslb_score", items), with = FALSE], by.x = "IID", by.y = "dog_id")
stopifnot(all(D$score == D$cslb_score))

to_date <- function(y, m) {
  y <- suppressWarnings(as.integer(y)); m <- suppressWarnings(as.integer(m))
  y[!(y %in% 1990:2025)] <- NA; m[is.na(m) | m < 1 | m > 12] <- 6L
  as.IDate(ifelse(is.na(y), NA_character_, sprintf("%d-%02d-15", y, m)))
}
h <- fread(hles_f, select = c("dog_id", "hs_condition_type", "hs_condition", "hs_diagnosis_year", "hs_diagnosis_month"), integer64 = "double",
           colClasses = list(character = c("hs_diagnosis_year", "hs_diagnosis_month")))[
  , .(dog_id = as.character(dog_id), type = as.integer(hs_condition_type), code = as.integer(hs_condition), dx = to_date(hs_diagnosis_year, hs_diagnosis_month))]
a <- fread(afus_f, select = c("dog_id", "afus_hs_new_condition_type", "afus_hs_new_condition", "afus_hs_new_condition_diagnosis_year", "afus_hs_new_condition_diagnosis_month"),
           integer64 = "double", colClasses = list(character = c("afus_hs_new_condition_diagnosis_year", "afus_hs_new_condition_diagnosis_month")))[
  , .(dog_id = as.character(dog_id), type = as.integer(afus_hs_new_condition_type), code = as.integer(afus_hs_new_condition),
      dx = to_date(afus_hs_new_condition_diagnosis_year, afus_hs_new_condition_diagnosis_month))]
cond <- merge(unique(rbind(h, a)[!is.na(code) & dog_id %in% D$IID]), D[, .(dog_id = IID, cslb_date)], by = "dog_id")
cond <- cond[is.na(dx) | dx <= cslb_date]
cb <- fread(cb_f)
lab_of <- function(v) { r <- cb[Variable == v][1]; data.table(code = as.integer(strsplit(r$Values, " \\| ")[[1]]), label = strsplit(r$ValueLabels, " \\| ")[[1]]) }
labs <- lab_of("hs_condition"); tlabs <- lab_of("hs_condition_type")
D[, dementia_pre := IID %in% cond[code == 1205, dog_id]]
for (j in c(cats, "dementia_pre")) set(D, j = j, value = as.integer(D[[j]]))
nn <- cond[!(code %in% c(pre_list, 1205)), .(n = uniqueN(code)), by = dog_id]; D[, n_other := nn$n[match(IID, nn$dog_id)]]
D[is.na(n_other), n_other := 0L][, n_other_bin := factor(pmin(n_other, 5), levels = 0:5, labels = c("0", "1", "2", "3", "4", "5+"))]

cov_main <- c(qn, setdiff(cn, "Vision"))
fit <- function(s, terms, yv = "y", vision = FALSE, extra = NULL) {
  covs <- c(cov_main, if (vision) "Vision")
  covs <- covs[vapply(covs, function(j) uniqueN(s[[j]]) > 1, TRUE)]
  f <- as.formula(paste(yv, "~", paste(c(covs, terms, extra), collapse = " + ")))
  m <- lm(f, data = s); co <- summary(m)$coefficients; rn <- rownames(co)
  rbindlist(lapply(terms, function(t) { k <- rn[startsWith(rn, t)]
    if (!length(k)) return(data.table(term = t, beta = NA_real_, se = NA_real_, p = NA_real_))
    data.table(term = k, beta = co[k, 1], se = co[k, 2], p = co[k, 4]) }))
}
grp <- list(all = D, age8 = D[age8 == TRUE])
for (g in names(grp)) grp[[g]][, y := rint(score)]
n_of <- function(s, t) if (t %in% names(s)) sum(s[[t]]) else NA_integer_

A <- rbindlist(lapply(names(grp), function(g) { s <- grp[[g]]
  rbind(rbindlist(lapply(c(cats, "dementia_pre"), function(t) fit(s, t)[, `:=`(model = "따로", n_yes = n_of(s, t))])),
        rbindlist(lapply(c(cats, "dementia_pre"), function(t) fit(s, t, vision = TRUE)[, `:=`(model = "따로 + Vision 공변량", n_yes = n_of(s, t))])),
        fit(s, c("sensory_pre", "neuro_pre", "pain_pre", "systemic_pre"))[, `:=`(model = "네 범주 함께", n_yes = NA_integer_)])[, group := g] }))
fwrite(A, paste0(out_prefix, "_premise.csv"))

B_tab <- rbindlist(lapply(names(grp), function(g) { s <- copy(grp[[g]])
  s[, combo := paste0(fifelse(sensory_pre == 1, "S", "-"), fifelse(neuro_pre == 1, "N", "-"), fifelse(pain_pre == 1, "P", "-"), fifelse(systemic_pre == 1, "E", "-"))]
  s[, .(n = .N, cslb_mean = mean(score), y_mean = mean(y)), by = combo][order(-n)][, group := g] }))
B_fit <- rbindlist(lapply(names(grp), function(g) { s <- copy(grp[[g]])
  s[, n_cat := factor(pmin(sensory_pre + neuro_pre + pain_pre + systemic_pre, 3), levels = 0:3, labels = c("0", "1", "2", "3+"))]
  s[, n_cat_num := as.numeric(as.character(fifelse(n_cat == "3+", "3", as.character(n_cat))))]
  rbind(fit(s, "n_cat")[, model := "범주 수 (0 기준)"], fit(s, "n_cat_num")[, model := "범주 수 1개당 (추세)"],
        fit(s, "sensory_pre:neuro_pre", extra = c("sensory_pre", "neuro_pre", "pain_pre", "systemic_pre"))[, model := "감각 × 신경 상호작용"])[, group := g] }))
fwrite(rbind(B_fit, fill = TRUE), paste0(out_prefix, "_multi.csv")); fwrite(B_tab, paste0(out_prefix, "_multi_combos.csv"))

C <- rbindlist(lapply(names(grp), function(g) { s <- copy(grp[[g]])
  rbindlist(lapply(items, function(it) { s[, yi := as.numeric(scale(get(it)))]
    fit(s, c("sensory_pre", "neuro_pre", "pain_pre", "systemic_pre", "dementia_pre"), yv = "yi")[, item := sub("cslb_", "", it)] }))[, group := g] }))
fwrite(C, paste0(out_prefix, "_items.csv"))
gC <- ggplot(C[group == "age8"], aes(item, term, fill = beta)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(p < 0.05, sprintf("%.2f", beta), "")), size = 2.8) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "β (SD)") +
  labs(x = "CSLB item (z; higher = worse for every item)", y = NULL, title = "Dogs aged 8+: condition effects on each CSLB item (mutually adjusted; number shown if p < 0.05)") +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(paste0(pic_prefix, "_items.png"), gC, width = 10, height = 3.6, dpi = 300)

screen <- function(key, by_col, lab) {
  cnt <- cond[, .(n = uniqueN(dog_id)), by = c(by_col)][n >= min_dogs]
  rbindlist(lapply(cnt[[by_col]], function(k) {
    ids_k <- cond[get(by_col) == k, unique(dog_id)]
    rbindlist(lapply(names(grp), function(g) { s <- copy(grp[[g]])[, has := IID %in% ids_k]
      if (sum(s$has) < 5) return(NULL)
      fit(s, "has")[, `:=`(group = g, n_yes = sum(s$has))] }))[, (key) := k] }))
}
SC <- screen("code", "code")[, term := NULL]
SC <- merge(SC, labs, by = "code", all.x = TRUE)
SC[, fdr := p.adjust(p, "BH"), by = group][, `:=`(in_list = code %in% pre_list, is_dementia = code == 1205)]
setorder(SC, group, p); fwrite(SC, paste0(out_prefix, "_screen_codes.csv"))
D2 <- rbindlist(lapply(names(grp), function(g) { s <- grp[[g]]
  rbind(fit(s, "n_other_bin")[, model := "목록 밖 질환 수 (0 기준)"],
        fit(s, "log1p(n_other)")[, model := "log(1 + 목록 밖 질환 수)"],
        fit(s, c("sensory_pre", "neuro_pre", "pain_pre", "systemic_pre"), extra = "log1p(n_other)")[, model := "네 범주 함께 + 목록 밖 질환 수 보정"])[, group := g] }))
fwrite(D2, paste0(out_prefix, "_reporting.csv"))
ST <- screen("type", "type")[, term := NULL]; ST <- merge(ST, setnames(copy(tlabs), "code", "type"), by = "type", all.x = TRUE)
ST[, fdr := p.adjust(p, "BH"), by = group]; setorder(ST, group, p); fwrite(ST, paste0(out_prefix, "_screen_types.csv"))
cand <- SC[group == "all" & !in_list & !is_dementia & beta > 0 & fdr < 0.05]
gS <- ggplot(SC[group == "all"], aes(beta, -log10(p), colour = fifelse(is_dementia, "dementia (1205)", fifelse(in_list, "11-1 list", "not in list")))) +
  geom_point(size = 1.4) + geom_text(data = SC[group == "all" & (fdr < 0.05 | is_dementia)], aes(label = substr(label, 1, 28)), size = 2.4, vjust = -0.6, show.legend = FALSE) +
  labs(x = "β (CSLB RINT, SD)", y = "-log10 p", colour = NULL, title = sprintf("All dogs: every condition diagnosed before CSLB (>= %d dogs), one at a time", min_dogs)) +
  theme_minimal(base_size = 9)
ggsave(paste0(pic_prefix, "_screen.png"), gS, width = 9, height = 5, dpi = 300)

E <- NULL
if (nrow(cand) && !no_greml) {
  work <- work_dir; dir.create(work, showWarnings = FALSE, recursive = TRUE)
  cand_ids <- cond[code %in% cand$code, unique(dog_id)]
  wr <- function(d, f) { fwrite(d, f, sep = " ", col.names = FALSE, quote = FALSE); f }
  reml <- function(name, s) {
    out <- file.path(work, name)
    if (!file.exists(paste0(out, ".hsq"))) {
      s <- copy(s)[, y := rint(score)]
      wr(s[, .(FID, IID, y)], paste0(out, ".pheno")); wr(s[, .(FID, IID)], paste0(out, ".keep"))
      kq <- qn[vapply(qn, function(j) isTRUE(sd(s[[j]]) > 0), TRUE)]; kc <- cn[vapply(cn, function(j) uniqueN(s[[j]]) > 1, TRUE)]
      wr(s[, c("FID", "IID", kq), with = FALSE], paste0(out, ".qcovar")); wr(s[, c("FID", "IID", kc), with = FALSE], paste0(out, ".covar"))
      system2(gcta, c("--reml", "--grm", grm, "--keep", paste0(out, ".keep"), "--pheno", paste0(out, ".pheno"), "--qcovar", paste0(out, ".qcovar"),
                      "--covar", paste0(out, ".covar"), "--reml-maxit", "200", "--out", out, "--threads", n_cores), stdout = paste0(out, ".stdout"), stderr = paste0(out, ".stdout"))
    }
    if (!file.exists(paste0(out, ".hsq"))) return(NULL)
    hq <- fread(paste0(out, ".hsq"), fill = TRUE); gg <- function(t) suppressWarnings(as.numeric(hq[Source == t, ]))
    data.table(sample = name, n = gg("n")[2], h2 = gg("V(G)/Vp")[2], se = gg("V(G)/Vp")[3], p_LRT = gg("Pval")[2])
  }
  set.seed(20261114)
  runs <- list()
  for (g in c("all", "age8")) {
    base <- if (g == "all") D else D[age8 == TRUE]
    keep <- base[excl_primary == FALSE & !(IID %in% cand_ids)]
    runs[[paste0(g, "_screen")]] <- keep
    if (g == "age8") {
      base[, age_bin := pmin(floor(age_cslb), 15)]
      dn <- base[!(IID %in% keep$IID), .N, by = age_bin]
      for (r in seq_len(n_rand)) {
        out_ids <- unlist(lapply(seq_len(nrow(dn)), function(i) { pool <- base[age_bin == dn$age_bin[i], IID]; sample(pool, min(length(pool), dn$N[i])) }))
        runs[[sprintf("age8_screen_random_%02d", r)]] <- base[!(IID %in% out_ids)]
      }
    }
  }
  E <- rbindlist(Map(reml, names(runs), runs))
  fwrite(E, paste0(out_prefix, "_screen_h2.csv"))
}

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("11-3 정제의 전제 확인 (질환 → CSLB; OLS, 주 공변량, 진단은 쓰인 CSLB 이전)\n")
cat(sprintf("개: 전체 %d · 8세 이상 %d\n", nrow(grp$all), nrow(grp$age8)))
cat("\n■ A. 범주별 CSLB 상승 (β = CSLB RINT SD 단위)\n"); print(f3(A[, .(group, model, term, n_yes, beta, se, p)]))
a8 <- A[group == "age8" & model == "따로"]
for (t in c("sensory_pre", "neuro_pre", "pain_pre", "systemic_pre", "dementia_pre")) { r <- a8[term == t]
  cat(sprintf("  판정 (8세 이상) %-14s β %.3f (p %.3g) → %s\n", t, r$beta, r$p, if (isTRUE(r$beta > 0 && r$p < 0.05)) "CSLB를 올림 (전제 지지)" else "CSLB 상승 근거 없음")) }
cat("\n■ B. 다중 질환 (S 감각 · N 신경 · P 통증 · E 전신)\n"); print(f3(B_tab[n >= 5])); print(f3(B_fit[, .(group, model, term, beta, se, p)]))
cat("\n■ C. 문항별 (8세 이상, 네 범주 + 치매를 함께 넣은 모형, p < 0.05만)\n"); print(f3(C[group == "age8" & p < 0.05, .(term, item, beta, p)][order(term, p)]))
cat(sprintf("\n■ D. 전체 질환 선별 (CSLB 이전 진단 %d마리 이상인 코드 %d개, 전체 표본, FDR < 0.05)\n", min_dogs, SC[group == "all", .N]))
print(f3(SC[group == "all" & fdr < 0.05, .(code, label, n_yes, beta, p, fdr, in_list)]))
cat("\n■ D2. 보고 경향 (11-1 목록 · 치매 밖 질환 코드 수)\n"); print(f3(D2[, .(group, model, term, beta, se, p)]))
cat(sprintf("  목록 밖 질환 수 분포 (전체): %s\n", paste(names(table(D$n_other_bin)), table(D$n_other_bin), sep = ":", collapse = " ")))
cat("  몸 계통별 (FDR < 0.05):\n"); print(f3(ST[group == "all" & fdr < 0.05, .(type, label, n_yes, beta, p, fdr)]))
cat(sprintf("  11-1 목록 밖 추가 후보 (β > 0 · FDR < 0.05; 치매 제외): %s\n", if (nrow(cand)) paste(sprintf("%d %s", cand$code, cand$label), collapse = " · ") else "없음"))
if (!is.null(E)) {
  cat("\n■ E. 결과 기반 탐색 (사후): 주 정제 + 후보 코드 제외 표본의 h²\n"); print(f3(E[!grepl("random", sample)]))
  rr <- E[grepl("random", sample), h2]
  if (length(rr)) cat(sprintf("  8세 이상 나이 맞춘 무작위 제외 %d번: h² 중앙값 %.3f · 95 백분위 %.3f\n", length(rr), median(rr), quantile(rr, 0.95)))
}
sink()
cat(sprintf("저장: %s_premise.csv · _multi.csv · _multi_combos.csv · _items.csv · _screen_codes.csv · _screen_types.csv%s · _summary.txt / 그림 %s_items.png · _screen.png\n",
            out_prefix, if (!is.null(E)) " · _screen_h2.csv" else "", pic_prefix))

