# 06-5_magma_geneset.R - MAGMA competitive gene-set tests (prespecified and genome-wide). See README.md.
library(data.table)
library(ggplot2)
library(parallel)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
dry_run <- Sys.getenv("DRY_RUN") == "1"
out_prefix <- file.path(qc_dir, paste0("06-5", if (dry_run) "_dry" else ""))
pic_prefix <- file.path(pic_dir, paste0("06-5", if (dry_run) "_dry" else ""))
set_prefix <- file.path(qc_dir, "06-4")
res63_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023/06-3_magma")
gsa_dir <- paste0(out_prefix, "_gsa")
magma <- "/data/00_tools/magma_v1.10/magma"
n_cores <- as.integer(Sys.getenv("N_CORES", "16"))
dir.create(gsa_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)

models <- c("main", "purelike", "mixlike", "int_age", "age8", "bin_age8", "young8", "sex_m", "sex_f", "wt_t1", "wt_t2", "wt_t3")
hets <- c("het_sex", "het_age", "het_ancestry", "het_weight")
nulls <- paste0("null", 1:20)
stage1_run <- "main"; stage2_runs <- c("int_age", "het_age"); age_ref <- c("age8", "young8")
stage3_runs <- c(setdiff(models, c(stage1_run, stage2_runs)), setdiff(hets, stage2_runs))
q_cut <- 0.05; p_nom <- 0.05

set_files <- c(primary = "_sets_primary.txt", reported = "_sets_primary_reported.txt", p1e6 = "_sets_primary_p1e6.txt",
               age_up = "_sets_age_up.txt", aux = "_sets_aux.txt")
set_files <- setNames(paste0(set_prefix, set_files), names(set_files))
plan_file <- paste0(set_prefix, "_analysis_plan.txt")
setinfo_file <- paste0(set_prefix, "_set_info.csv")
for (f in c(set_files, plan_file, setinfo_file, magma)) if (!file.exists(f)) stop("파일 없음: ", f)

plan <- readLines(plan_file)
md5_plan <- regmatches(plan, regexec("^  ([0-9a-f]{32})  (\\S+)$", plan))
md5_plan <- rbindlist(lapply(md5_plan[lengths(md5_plan) == 3], function(x) data.table(md5 = x[2], file = x[3])))
md5_now <- data.table(file = basename(set_files), md5_now = unname(tools::md5sum(set_files)))
md5_chk <- merge(md5_now, md5_plan, by = "file", all.x = TRUE)
if (md5_chk[, any(is.na(md5) | md5 != md5_now)]) {
  print(md5_chk); stop("세트 파일이 06-4_analysis_plan.txt의 md5와 다름 → 06-4를 다시 돌려 계획을 새로 남긴 뒤 실행")
}
cat(sprintf("세트 파일 md5가 계획과 같음 (%s)\n", sub("^06-5 분석 계획 \\((.*?);.*", "\\1", plan[1])))

si <- fread(setinfo_file)
primary_sets <- si[kept == TRUE, set]
family_of <- setNames(si$family, si$set)
if (length(primary_sets) != 14) warning("주 가설 세트가 14개가 아님: ", length(primary_sets))

null_for <- function(run) {
  base <- sub("_w50$", "", run); all_runs <- c(models, hets)
  k <- (match(base, all_runs) - 1) %% 20 + 1
  paste0("null", k)
}
raw_of <- function(run) {
  r <- if (dry_run && !run %in% nulls) null_for(run) else run
  file.path(res63_dir, paste0(r, ".genes.raw"))
}
jobs <- rbind(
  data.table(run = c(models, paste0(models, "_w50"), hets, paste0(hets, "_w50"), nulls), variant = "primary"),
  CJ(run = c(models, hets), variant = c("reported", "p1e6", "age_up", "aux")),
  data.table(run = nulls, variant = "aux")
)
jobs[, kind := fifelse(run %in% nulls, "null", fifelse(grepl("_w50$", run), "w50", fifelse(run %in% hets, "het", "model")))]
jobs[, `:=`(raw = vapply(run, raw_of, ""), out = file.path(gsa_dir, paste0(run, "__", variant)))]
if (length(miss <- jobs[!file.exists(raw), unique(raw)])) stop("06-3 genes.raw 없음: ", paste(basename(miss), collapse = ", "))
jobs[, done := file.exists(paste0(out, ".gsa.out")) & file.mtime(paste0(out, ".gsa.out")) > file.mtime(set_files[variant])]
cat(sprintf("MAGMA 세트 검정 %d개 (끝난 것 %d개는 건너뜀), 동시 %d개%s\n", nrow(jobs), sum(jobs$done), n_cores,
            if (dry_run) " — DRY_RUN: 실제 run 대신 null genes.raw 사용" else ""))
t0 <- Sys.time()
todo <- jobs[done == FALSE]
st <- unlist(mclapply(seq_len(nrow(todo)), function(i) {
  system2(magma, c("--gene-results", shQuote(todo$raw[i]), "--set-annot", shQuote(set_files[[todo$variant[i]]]), "--out", shQuote(todo$out[i])),
          stdout = FALSE, stderr = FALSE)
}, mc.cores = n_cores, mc.preschedule = FALSE))
if (length(st) && any(st != 0)) stop("MAGMA 실패: ", paste(basename(todo$out[st != 0]), collapse = ", "), " (각 .log 확인)")
cat(sprintf("MAGMA 끝 (%.1f분)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

read_gsa <- function(f) {
  d <- fread(f, skip = "VARIABLE")
  if ("FULL_NAME" %in% names(d)) d[, VARIABLE := FULL_NAME]
  d[, .(set = VARIABLE, NGENES, BETA, BETA_STD, SE, P)]
}
gsa <- rbindlist(lapply(seq_len(nrow(jobs)), function(i) {
  d <- read_gsa(paste0(jobs$out[i], ".gsa.out")); d[, `:=`(run = jobs$run[i], variant = jobs$variant[i], kind = jobs$kind[i])]; d
}))
aux <- gsa[variant == "aux"]
gsa <- gsa[variant != "aux"]
gsa[, base_set := sub("_p1e6$", "", set)]
gsa[variant == "age_up" & set == "AGE_brain_up_Palmer2021", base_set := "AGE_brain_Palmer2021"]
if (gsa[variant != "primary", any(!base_set %in% primary_sets)]) stop("민감도 세트 중 주 가설 세트와 짝이 없는 것: ",
  paste(gsa[variant != "primary" & !base_set %in% primary_sets, unique(set)], collapse = ", "))
gsa[, family := family_of[base_set]]
pr <- gsa[variant == "primary"]
if (pr[, uniqueN(set), by = run][, any(V1 != length(primary_sets))]) stop("run마다 주 가설 세트 수가 다름")

nl <- pr[kind == "null"]
nl[, q_run := p.adjust(P, "BH"), by = run]
lam <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
null_calib <- nl[, .(n_null = .N, frac_p05 = mean(P < p_nom), min_p = min(P), lambda = round(lam(P), 3)), by = set]
null_overall <- list(
  lambda = lam(nl$P), frac_p05 = mean(nl$P < p_nom),
  stage1_fp = nl[, any(q_run < q_cut), by = run][, mean(V1)],
  stage2_fp = {
    pairs <- split(nulls, ceiling(seq_along(nulls) / 2))
    mean(vapply(pairs, function(pp) any(p.adjust(nl[run %in% pp, P], "BH") < q_cut), TRUE))
  })
emp_p <- function(s, p_obs) (1 + nl[set == s, sum(P <= p_obs)]) / (1 + nl[set == s, .N])

s1 <- pr[run == stage1_run, .(set, family, NGENES, BETA, BETA_STD, SE, P)]
s1[, q := p.adjust(P, "BH")]
s1[, `:=`(emp_p = mapply(emp_p, set, P), null_frac_p05 = null_calib$frac_p05[match(set, null_calib$set)])]
rep_tab <- pr[kind == "model" & run != stage1_run, .(n_models_p05 = sum(P < p_nom), models_p05 = paste(run[P < p_nom], collapse = ",")), by = set]
s1 <- merge(s1, rep_tab, by = "set", sort = FALSE)
s1[, `:=`(significant = q < q_cut, nominal = P < p_nom)]
setorder(s1, P)

s2 <- pr[run %in% stage2_runs, .(run, set, family, NGENES, BETA, BETA_STD, SE, P)]
s2[, q := p.adjust(P, "BH")]
s2[, `:=`(emp_p = mapply(emp_p, set, P), significant = q < q_cut, nominal = P < p_nom)]
setorder(s2, P)
ar <- dcast(pr[run %in% age_ref, .(run, set, BETA, SE, P)], set ~ run, value.var = c("BETA", "SE", "P"))
ar[, `:=`(z_diff = (BETA_age8 - BETA_young8) / sqrt(SE_age8^2 + SE_young8^2))]
ar[, p_diff := 2 * pnorm(-abs(z_diff))]
ar[, direction := fifelse(z_diff > 0, "age8 > young8", "age8 < young8")]
ar[, family := family_of[set]]
setorder(ar, p_diff)

s3 <- pr[run %in% stage3_runs, .(run, kind, set, family, NGENES, BETA, SE, P)]
s3[, q_run := p.adjust(P, "BH"), by = run]
s3[, `:=`(nominal = P < p_nom, null_frac_p05 = null_calib$frac_p05[match(set, null_calib$set)])]
s3_sum <- s3[, .(n_q05 = sum(q_run < q_cut), n_p05 = sum(P < p_nom), sets_p05 = paste(set[P < p_nom], collapse = ",")), by = .(run, kind)]

key_runs <- c(stage1_run, stage2_runs)
sens <- pr[run %in% key_runs, .(run, set, P_primary = P, BETA_primary = BETA)]
w50 <- pr[kind == "w50" & sub("_w50$", "", run) %in% key_runs, .(run = sub("_w50$", "", run), set, P_w50 = P, BETA_w50 = BETA)]
sens <- merge(sens, w50, by = c("run", "set"), all.x = TRUE)
for (v in c("reported", "p1e6", "age_up")) {
  d <- gsa[variant == v & run %in% key_runs, .(run, set = base_set, P)]
  if (v == "reported") d <- d[set %in% si[grepl("^GWAS_", set), set]]
  setnames(d, "P", paste0("P_", v))
  sens <- merge(sens, d, by = c("run", "set"), all.x = TRUE)
}
sens[, stage := fifelse(run == stage1_run, "1", "2")]
qs <- rbind(s1[, .(run = stage1_run, set, q_stage = q)], s2[, .(run, set, q_stage = q)])
sens <- merge(sens, qs, by = c("run", "set"))
sens[, n_sens := rowSums(!is.na(cbind(P_w50, P_reported, P_p1e6, P_age_up)))]
sens[, holds := {
  ps <- cbind(P_w50, P_reported, P_p1e6, P_age_up)
  apply(ps, 1, function(x) { x <- x[!is.na(x)]; if (!length(x)) NA else all(x < p_nom) })
}]
setorder(sens, stage, q_stage)

aux[, q_run := p.adjust(P, "BH"), by = run]
aux_null <- aux[kind == "null", .(n_q05 = sum(q_run < q_cut), lambda = lam(P), min_p = min(P)), by = run]
aux_null_overall <- list(frac_any_q05 = mean(aux_null$n_q05 > 0), mean_q05 = mean(aux_null$n_q05), lambda = median(aux_null$lambda),
                         min_p_median = median(aux_null$min_p))
aux_null_sets <- aux[kind == "null" & q_run < q_cut, .(run, set, NGENES, P, q_run)]
aux <- aux[kind != "null"]
aux_top <- aux[order(P), head(.SD, 20), by = run][, rank := seq_len(.N), by = run][, .(run, rank, set, NGENES, BETA, P, q_run)]
aux_sum <- aux[, .(n_sets = .N, n_q05 = sum(q_run < q_cut), lambda = round(lam(P), 3), top = paste(head(set[order(P)], 5), collapse = ",")), by = run]
aux_info_file <- paste0(set_prefix, "_aux_info.csv")
if (file.exists(aux_info_file)) {
  ai <- fread(aux_info_file, select = c("set", "name"))
  aux_top[, name := ai$name[match(set, ai$set)]]
}

fwrite(gsa, paste0(out_prefix, "_gsa_all.csv.gz"))
fwrite(s1, paste0(out_prefix, "_stage1.csv"))
fwrite(s2, paste0(out_prefix, "_stage2.csv"))
fwrite(ar, paste0(out_prefix, "_age_ref.csv"))
fwrite(s3, paste0(out_prefix, "_stage3.csv"))
fwrite(null_calib, paste0(out_prefix, "_null_calib.csv"))
fwrite(sens, paste0(out_prefix, "_sensitivity.csv"))
fwrite(aux, paste0(out_prefix, "_aux.csv.gz"))
fwrite(aux_top, paste0(out_prefix, "_aux_top.csv"))
fwrite(aux_null, paste0(out_prefix, "_aux_null.csv"))

run_order <- c(stage1_run, stage2_runs, age_ref, setdiff(stage3_runs, age_ref))
hm <- pr[run %in% run_order]
hm[, run := factor(run, levels = run_order)]
hm[, set := factor(set, levels = rev(primary_sets))]
hm[, mlp := pmin(-log10(P), 4)]
qmark <- rbind(s1[, .(run = stage1_run, set, q)], s2[, .(run, set, q)], s3[, .(run, set, q = q_run)])
hm <- merge(hm, qmark, by = c("run", "set"))
hm[, `:=`(run = factor(as.character(run), levels = run_order), set = factor(as.character(set), levels = rev(primary_sets)))]
hm[, lab := fifelse(q < q_cut, "*", fifelse(P < p_nom, ".", ""))]
hm[, stage := factor(fifelse(run == stage1_run, "Stage 1", fifelse(run %in% stage2_runs, "Stage 2 (age)", "Stage 3 (exploratory)")),
                     levels = c("Stage 1", "Stage 2 (age)", "Stage 3 (exploratory)"))]
g1 <- ggplot(hm, aes(run, set, fill = mlp)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lab), size = 5, vjust = 0.75) +
  scale_fill_gradient(low = "#f4f6fb", high = "#1f4e9c", limits = c(0, 4), name = "-log10 p\n(one-sided,\ncapped at 4)") +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = "MAGMA competitive gene-set test: 14 pre-specified sets",
       subtitle = "* = BH q < 0.05 within the stage family (stage 1: 14 tests; stage 2: 28; stage 3: 14 per run)   . = nominal p < 0.05 (not evidence)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), strip.text = element_text(face = "bold"))
ggsave(paste0(pic_prefix, "_set_heatmap.png"), g1, width = 12, height = 6.5, dpi = 300)

ab <- pr[run %in% age_ref, .(run, set, BETA, SE)]
ab[, set := factor(set, levels = rev(primary_sets))]
g2 <- ggplot(ab, aes(BETA, set, colour = run)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = BETA - 1.96 * SE, xmax = BETA + 1.96 * SE), position = position_dodge(width = 0.6), size = 0.35) +
  scale_colour_manual(values = c(age8 = "#b4462b", young8 = "#2b6cb4"), labels = c(age8 = "age >= 8 (n 2,249)", young8 = "age < 8 (n 4,866)"), name = NULL) +
  labs(x = "MAGMA gene-set beta (95% CI)", y = NULL, title = "Stage 2 reference: gene-set effect by age stratum",
       subtitle = "Descriptive only (not used for inference); strata are non-overlapping dogs") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(paste0(pic_prefix, "_age_beta.png"), g2, width = 8, height = 6.5, dpi = 300)

qq <- rbind(nl[, .(grp = "null (20 permuted GWAS)", P)], pr[run == stage1_run, .(grp = "main", P)], pr[run %in% stage2_runs, .(grp = "int_age + het_age", P)])
qq <- qq[order(P), .(obs = -log10(P), exp = -log10(ppoints(.N))), by = grp]
g3 <- ggplot(qq, aes(exp, obs, colour = grp)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") + geom_point(size = 1.6) +
  scale_colour_manual(values = c("null (20 permuted GWAS)" = "grey55", main = "#1f4e9c", "int_age + het_age" = "#b4462b"), name = NULL) +
  labs(x = "Expected -log10 p", y = "Observed -log10 p", title = "Gene-set p QQ (14 pre-specified sets)",
       subtitle = sprintf("null lambda = %.2f; null fraction p < 0.05 = %.3f", null_overall$lambda, null_overall$frac_p05)) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(paste0(pic_prefix, "_null_qq.png"), g3, width = 6, height = 5.5, dpi = 300)

fmt <- function(d) { d <- copy(d); for (j in intersect(c("BETA", "BETA_STD", "SE", "BETA_age8", "BETA_young8", "SE_age8", "SE_young8", "z_diff"), names(d))) set(d, j = j, value = signif(d[[j]], 3))
                     for (j in intersect(c("P", "q", "q_run", "emp_p", "p_diff", "P_primary", "P_w50", "P_reported", "P_p1e6", "P_age_up", "q_stage", "P_age8", "P_young8"), names(d))) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("06-5 MAGMA 유전자 세트 검정 (competitive, 한쪽 β > 0; 06-4_analysis_plan.txt 계획 그대로)\n")
if (dry_run) cat("*** DRY_RUN: 모든 실제 run 대신 null genes.raw를 넣은 시험 결과. 해석하지 말 것 ***\n")
cat(sprintf("세트: 주 가설 %d개 (%s)\n\n", length(primary_sets),
            paste(sprintf("%s %d", names(table(family_of[primary_sets])), table(family_of[primary_sets])), collapse = " · ")))

cat("■ null 보정 확인 (06-3 null 20개 × 주 가설 세트)\n")
cat(sprintf("  세트 p의 λ %.3f (기대 1), p < 0.05 비율 %.3f (기대 0.05)\n", null_overall$lambda, null_overall$frac_p05))
cat(sprintf("  1단계 판정 (14개 BH q < 0.05)을 null에 적용 → 하나라도 걸린 null run %.0f%% (기대 ≤ 5%%)\n", 100 * null_overall$stage1_fp))
cat(sprintf("  2단계 판정 (null 2개 짝 → 28개 BH) → 하나라도 걸린 짝 %.0f%% (짝 10개; null은 OLS 연관 GWAS라 int_age·het_age에는 근사)\n", 100 * null_overall$stage2_fp))
cat("  경험적 p (emp_p)는 null 20개라 가능한 값이 0.048, 0.095, ...뿐 → 참고용, 판정에 쓰지 않음\n")
print(fmt(null_calib[order(-frac_p05)]))

cat("\n■ 1단계 — 주 결론: main × 14개 (14개 안 BH)\n")
print(fmt(s1[, .(set, family, NGENES, BETA, SE, P, q, emp_p, null_frac_p05, n_models_p05, significant)]))
cat(sprintf("  → q < %.2f 세트: %s\n", q_cut, if (s1[, any(significant)]) paste(s1[significant == TRUE, set], collapse = ", ") else "없음"))
if (s1[, any(significant)]) { cat("  재현 확인 (다른 11개 모형에서 p < 0.05):\n"); print(s1[significant == TRUE, .(set, n_models_p05, models_p05)]) }
cat(sprintf("  명목 p < 0.05 (근거 아님): %s\n", if (s1[, any(nominal & !significant)]) paste(s1[nominal & !significant, sprintf("%s (p %.3g, null에서 이 세트 p < 0.05 비율 %.2f)", set, P, null_frac_p05)], collapse = "; ") else "없음"))

cat("\n■ 2단계 — 나이 가설: int_age · het_age × 14개 = 28개 (한 family BH)\n")
print(fmt(s2[, .(run, set, family, NGENES, BETA, SE, P, q, emp_p, significant)]))
cat(sprintf("  → q < %.2f: %s\n", q_cut, if (s2[, any(significant)]) paste(s2[significant == TRUE, paste0(run, ":", set)], collapse = ", ") else "없음"))
cat("  참고 표 (판정에 쓰지 않음): age8 vs young8 세트 β 차이\n")
print(fmt(ar[, .(set, family, BETA_age8, SE_age8, P_age8, BETA_young8, SE_young8, P_young8, z_diff, p_diff, direction)]))
cat("  (p_diff는 양쪽 검정; 판정에 쓰지 않음)\n")

cat("\n■ 3단계 — 탐색 (결론 아님): run마다 14개 안 BH\n")
print(s3_sum)

cat("\n■ 민감도 (1·2단계 run; holds = 있는 민감도에서 모두 p < 0.05)\n")
cat("  P_w50 = ±50 kb, P_reported = GWAS 세트 저자 유전자, P_p1e6 = GWAS 세트 p ≤ 1e-6, P_age_up = 뇌 노화 증가 유전자만\n")
print(fmt(sens[, .(stage, run, set, P_primary, q_stage, P_w50, P_reported, P_p1e6, P_age_up, n_sens, holds)]))
cat("  n_sens = 이 세트에 있는 민감도 수 (circadian·스플라이싱은 ±50 kb만). holds는 1·2단계에서 유의한 세트에만 의미가 있음\n")

cat("\n■ 보조 GO/KEGG (설명용; run마다 따로 BH)\n")
cat(sprintf("  null 20개에 같은 보조 세트: q < 0.05 세트가 하나라도 나온 null run %.0f%%, run당 평균 %.2f개, λ 중앙값 %.2f, 최소 p 중앙값 %.2g\n",
            100 * aux_null_overall$frac_any_q05, aux_null_overall$mean_q05, aux_null_overall$lambda, aux_null_overall$min_p_median))
if (nrow(aux_null_sets)) { cat("  null에서 q < 0.05가 된 보조 세트 (작은 세트·LD 묶음에서 잘 생김):\n"); print(aux_null_sets) }
cat("  → 실제 run의 보조 세트 q < 0.05는 이 null 비율과 비교해서만 읽는다.\n")
print(aux_sum)
cat("  main · int_age · het_age 상위 10:\n")
print(fmt(aux_top[run %in% key_runs & rank <= 10, .(run, rank, set, name = if ("name" %in% names(aux_top)) strtrim(name, 50) else "", NGENES, P, q_run)]))
sink()
cat(sprintf("저장: %s_gsa/ · _gsa_all.csv.gz · _stage1.csv · _stage2.csv · _age_ref.csv · _stage3.csv · _null_calib.csv · _sensitivity.csv · _aux.csv.gz · _aux_top.csv · _summary.txt\n그림: %s_set_heatmap.png · _age_beta.png · _null_qq.png\n",
            out_prefix, pic_prefix))

#      age8: PLEKHA1, ANK3, CYYR1, TMEM106B; bin_age8: HLA-DRB1 (DLA), WDR12, ANK3, APOE)
