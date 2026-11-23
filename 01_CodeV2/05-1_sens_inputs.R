# 05-1_sens_inputs.R - Defines the stratified and alternative covariate models and writes their input files. See README.md.
library(data.table)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")   
gwas_045 <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023/04-5_gwas")
out_prefix <- file.path(qc_dir, "05-1")
in_dir <- paste0(out_prefix, "_inputs")
dir.create(in_dir, showWarnings = FALSE, recursive = TRUE)

main_pc <- 5
pc_all_sets <- c(0, 5, 10, 15, 20, 25, 30)
age_knot <- 11
bin_cut <- 50
cont_main <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age")
cat_main <- c("Sex", "Intact", "Vision", "Cohort_BHS")

rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))

fam <- fread(file.path(qc_dir, "03-3_maf0.01_hwe.fam"), header = FALSE, select = 1:2, col.names = c("FID", "IID"))
fam[, `:=`(FID = as.character(FID), IID = as.character(IID))]
ph <- fread(file.path(qc_dir, "04-2_pheno.txt"), header = FALSE, col.names = c("FID", "IID", "rint_all", "raw"), colClasses = list(character = 1:2))
qc <- fread(file.path(qc_dir, sprintf("03-5_qcovar_PC%d.txt", max(main_pc, 5))), header = FALSE, colClasses = list(character = 1:2))
setnames(qc, c("FID", "IID", cont_main, paste0("PC", seq_len(ncol(qc) - 7))))
cv <- fread(file.path(qc_dir, "03-5_covar.txt"), header = FALSE, col.names = c("FID", "IID", cat_main), colClasses = list(character = 1:2))
pc_alt <- lapply(c(all = "all", qmat = "qmat"), function(m) {
  x <- fread(file.path(qc_dir, sprintf("03-5_qcovar_%s_PC10.txt", m)), header = FALSE, colClasses = list(character = 1:2))
  setnames(x, c("FID", "IID", cont_main, paste0(m, "_PC", 1:10)))[, c("IID", paste0(m, "_PC", 1:10)), with = FALSE]
})
sens <- fread(file.path(qc_dir, "00_df_sensitivity.csv"), colClasses = list(character = c("FID", "IID")))
keep_sets <- c("age8", "age10", "agefloor6", "agefloor7", "agefloor9", "agefloor11", "agefloor12",
               "purelike", "mixlike", "nomimic", "dob_exact", "afus", "ac", "no_bhs")
keeps <- lapply(setNames(keep_sets, keep_sets), function(s)
  fread(file.path(qc_dir, sprintf("00_keep_%s.txt", s)), header = FALSE, colClasses = "character")$V2)

d <- Reduce(function(a, b) merge(a, b, by = "IID", all.x = TRUE, sort = FALSE),
            list(fam, ph[, .(IID, rint_all, raw)], qc[, -"FID"], cv[, -"FID"], pc_alt$all, pc_alt$qmat, sens[, -"FID"]))
d <- d[match(fam$IID, IID)]
stopifnot(identical(d$IID, fam$IID), !anyNA(d[, c("raw", cont_main, cat_main), with = FALSE]))
stopifnot(max(abs(rint(d$raw) - d$rint_all)) < 1e-5)
for (m in c("all", "qmat")) {
  x <- fread(file.path(qc_dir, sprintf("03-5_qcovar_%s_PC10.txt", m)), header = FALSE, colClasses = list(character = 1:2))
  stopifnot(identical(x$V2, qc$IID), max(abs(as.matrix(x[, 3:7]) - as.matrix(qc[, 3:7]))) < 1e-9)
}
d[, Vision_narrow := fifelse(vision_cov_narrow == 1, "vision_problem", "normal")]
d[, bin50 := as.integer(raw >= bin_cut)]

keeps$young8 <- d[!IID %in% keeps$age8, IID]
stopifnot(all(d[IID %in% keeps$age8, age_cslb] >= 8), all(d[IID %in% keeps$young8, age_cslb] < 8))
stopifnot(setequal(unique(d$Sex), c("Male", "Female")))
keeps$sex_m <- d[Sex == "Male", IID]
keeps$sex_f <- d[Sex == "Female", IID]
wt_cut <- unname(quantile(d$wt_cov, c(1, 2) / 3))
d[, wt_tert := findInterval(wt_cov, wt_cut, left.open = TRUE) + 1L]
for (k in 1:3) keeps[[paste0("wt_t", k)]] <- d[wt_tert == k, IID]
strata <- list(age = c("age8", "young8"), sex = c("sex_m", "sex_f"), wt = paste0("wt_t", 1:3))
for (s in strata) {
  v <- unlist(lapply(s, function(k) intersect(keeps[[k]], d$IID)))
  stopifnot(!anyDuplicated(v), setequal(v, d$IID))
}
cat(sprintf("층: young8 %d / age8 %d | 수컷 %d / 암컷 %d | 체중 3분위 경계 %.1f·%.1f kg → %s마리\n",
            length(keeps$young8), sum(d$IID %in% keeps$age8), length(keeps$sex_m), length(keeps$sex_f), wt_cut[1], wt_cut[2],
            paste(d[, .N, keyby = wt_tert]$N, collapse = " / ")))
cat(sprintf("자료: %d마리, CSLB ≥ %d: %d마리\n", nrow(d), bin_cut, sum(d$bin50)))

pc_main <- paste0("PC", seq_len(main_pc))
spec <- function(id, group, desc, engine = "gcta", y = "rint", keep = NULL, add_q = NULL, add_c = NULL,
                 drop_q = NULL, drop_c = NULL, pcs = pc_main, swap_c = NULL, exclude = NULL)
  list(id = id, group = group, desc = desc, engine = engine, y = y, keep = keep, add_q = add_q, add_c = add_c,
       drop_q = drop_q, drop_c = drop_c, pcs = pcs, swap_c = swap_c, exclude = exclude)
runs <- list(
  spec("main_check", "check", "주 분석과 같은 설정 (05-2 CHECK=1에서 chr38 일부로 04-5 결과와 같은지 확인용, 전체 실행 안 함)", engine = "check"),
  spec("raw", "phenotype", "CSLB 원점수 (주: RINT)", y = "raw"),
  spec("bin_all", "phenotype", sprintf("CSLB ≥ %d case-control, 전체 (REGENIE Firth)", bin_cut), engine = "regenie_bt", y = "bin50"),
  spec("bin_age8", "phenotype", sprintf("CSLB ≥ %d case-control, 8세 이상 (REGENIE Firth)", bin_cut), engine = "regenie_bt", y = "bin50", keep = "age8"),
  spec("age8", "age", "8세 이상", keep = "age8"),
  spec("age10", "age", "10세 이상", keep = "age10"),
  spec("agefloor6", "age", "6세 이상", keep = "agefloor6"),
  spec("agefloor7", "age", "7세 이상", keep = "agefloor7"),
  spec("agefloor9", "age", "9세 이상", keep = "agefloor9"),
  spec("agefloor11", "age", "11세 이상", keep = "agefloor11"),
  spec("agefloor12", "age", "12세 이상", keep = "agefloor12"),
  spec("int_age", "age", "SNP × age_cslb 상호작용 (REGENIE --interaction)", engine = "regenie_int"),
  spec("purelike", "genetic", "순종형 (log COI 혼합분포)", keep = "purelike"),
  spec("mixlike", "genetic", "잡종형 (log COI 혼합분포)", keep = "mixlike"),
  spec("young8", "age", "8세 미만 (age8의 여집합, age8과 짝 비교)", keep = "young8"),
  spec("sex_m", "strata", "수컷만 (sex_f와 짝 비교)", keep = "sex_m"),
  spec("sex_f", "strata", "암컷만 (sex_m과 짝 비교)", keep = "sex_f"),
  spec("wt_t1", "strata", sprintf("체중 하위 1/3 (≤ %.1f kg)", wt_cut[1]), keep = "wt_t1"),
  spec("wt_t2", "strata", sprintf("체중 중간 1/3 (%.1f ~ %.1f kg)", wt_cut[1], wt_cut[2]), keep = "wt_t2"),
  spec("wt_t3", "strata", sprintf("체중 상위 1/3 (> %.1f kg)", wt_cut[2]), keep = "wt_t3"),
  spec("drop_age", "drop_cov", "age_cslb + age_hinge 빼기", drop_q = c("age_cslb", "age_hinge")),
  spec("drop_pa", "drop_cov", "pa_score 빼기", drop_q = "pa_score"),
  spec("drop_wt", "drop_cov", "wt_cov 빼기", drop_q = "wt_cov"),
  spec("drop_owner_age", "drop_cov", "owner_age 빼기", drop_q = "owner_age"),
  spec("drop_sex", "drop_cov", "Sex 빼기", drop_c = "Sex"),
  spec("drop_intact", "drop_cov", "Intact 빼기", drop_c = "Intact"),
  spec("drop_vision", "drop_cov", "Vision 빼기", drop_c = "Vision"),
  spec("drop_bhs", "drop_cov", "Cohort_BHS 빼기", drop_c = "Cohort_BHS"),
  spec("add_ancestry", "add_cov", "top_pct + COI 추가", add_q = c("top_pct", "coi")),
  spec("add_mimic_dx", "add_cov", "청력·발작·전정·뇌 질환 진단 이력 추가", add_q = c("hearing_cov", "seizure_cov", "vestibular_cov", "brain_dx")),
  spec("add_oa", "add_cov", "골관절염 진단 이력 추가", add_q = "oa_cov"),
  spec("add_ac", "add_cov", "AFUS 관찰 증상(시력·청력·이동 통증·이동 약화) 추가, AFUS 응답 개체", keep = "ac",
       add_q = c("ac_blind", "ac_deaf", "ac_mobility_pain", "ac_mobility_weak")),
  spec("add_health", "add_cov", "전반적 건강 평가 추가", add_q = "general_health"),
  spec("add_household", "add_cov", "가구 내 개·사람 수 + 거주 지역 추가", add_q = c("household_dogs", "household_people"), add_c = "home_area"),
  spec("add_ses", "add_cov", "보호자 소득 + 학력 추가 (소득 결측 제외)", add_q = c("owner_income", "owner_edu")),
  spec("vision_narrow", "add_cov", "Vision을 좁은 정의(실명·백내장)로 교체", swap_c = c(Vision = "Vision_narrow")),
  spec("nomimic", "subset", "CCD 유사 질환 진단 이력 제외", keep = "nomimic"),
  spec("dob_exact", "subset", "생년월일 정확한 개만", keep = "dob_exact"),
  spec("afus", "subset", "AFUS 1회 이상", keep = "afus"),
  spec("no_bhs", "subset", "Brain Health Study 제외", keep = "no_bhs"),
  spec("pc_all", "qc", sprintf("PC 방법 all PCA (PC %d개)", main_pc), pcs = paste0("all_PC", seq_len(main_pc))),
  spec("pc_qmat", "qc", sprintf("PC 방법 Q-matrix PCA (PC %d개)", main_pc), pcs = paste0("qmat_PC", seq_len(main_pc)))
)
ids <- vapply(runs, `[[`, "", "id")
stopifnot(!anyDuplicated(ids), grepl("^[A-Za-z0-9_]+$", ids))
if (main_pc > 10) stop("all·qmat PC 파일은 PC 10개까지만 있음")

check_design <- function(x, qn, cn) {
  dropped <- character()
  for (v in c(qn, cn)) if (uniqueN(x[[v]]) < 2) dropped <- c(dropped, v)
  qn <- setdiff(qn, dropped); cn <- setdiff(cn, dropped)
  repeat {
    f <- as.formula(paste("~", paste(c(qn, if (length(cn)) paste0("factor(", cn, ")")), collapse = " + ")))
    mm <- model.matrix(f, data = x)
    q <- qr(mm, tol = 1e-7)
    if (q$rank == ncol(mm)) break
    bad <- q$pivot[(q$rank + 1):ncol(mm)][1]
    term <- attr(mm, "assign")[bad]
    if (term == 0) stop("절편이 선형 종속으로 나옴 (확인 필요)")
    v <- c(qn, cn)[term]
    dropped <- c(dropped, v); qn <- setdiff(qn, v); cn <- setdiff(cn, v)
  }
  list(qn = qn, cn = cn, dropped = dropped)
}
wr <- function(x, f, header = FALSE) fwrite(x, f, sep = " ", col.names = header, quote = FALSE, na = "NA")

need <- unique(unlist(lapply(runs, function(r) c(r$pcs, r$add_q, r$add_c, unname(r$swap_c), if (r$y != "rint") r$y))))
if (length(setdiff(need, names(d)))) stop("자료에 없는 변수: ", paste(setdiff(need, names(d)), collapse = ", "))

tab <- list()
for (r in runs) {
  x <- copy(d)
  if (!is.null(r$keep)) x <- x[IID %in% keeps[[r$keep]]]
  if (!is.null(r$exclude)) x <- x[!IID %in% r$exclude]
  qn <- c(setdiff(cont_main, r$drop_q), r$pcs, r$add_q)
  cn <- setdiff(cat_main, r$drop_c)
  if (!is.null(r$swap_c)) cn[cn %in% names(r$swap_c)] <- r$swap_c[cn[cn %in% names(r$swap_c)]]
  cn <- c(cn, r$add_c)
  ycol <- if (r$y == "rint") "raw" else r$y
  n0 <- nrow(x)
  x <- x[complete.cases(x[, c(ycol, qn, cn), with = FALSE])]
  if (nrow(x) < 100) stop(r$id, ": 분석 개체가 너무 적음 (", nrow(x), ")")
  x[, y := if (r$y == "rint") rint(raw) else as.numeric(get(r$y))]
  ck <- check_design(x, qn, cn)
  base <- file.path(in_dir, r$id)
  wr(x[, .(FID, IID)], paste0(base, "_keep.txt"))
  if (r$engine %in% c("gcta", "check")) {
    wr(x[, .(FID, IID, y = round(y, 6))], paste0(base, "_pheno.txt"))
    wr(x[, c("FID", "IID", ck$qn), with = FALSE], paste0(base, "_qcovar.txt"))
    if (length(ck$cn)) wr(x[, c("FID", "IID", ck$cn), with = FALSE], paste0(base, "_covar.txt"))
  } else {
    pname <- if (r$engine == "regenie_bt") "CCD50" else "CSLB_RINT"
    p <- x[, .(FID, IID, y = if (r$engine == "regenie_bt") as.integer(y) else round(y, 6))]
    setnames(p, "y", pname)
    wr(p, paste0(base, "_regenie_pheno.txt"), header = TRUE)
    wr(x[, c("FID", "IID", ck$qn, ck$cn), with = FALSE], paste0(base, "_regenie_covar.txt"), header = TRUE)
  }
  tab[[r$id]] <- data.table(
    id = r$id, group = r$group, engine = r$engine, desc = r$desc, y = r$y,
    n = nrow(x), n_removed_missing = n0 - nrow(x),
    n_case = if (r$y == "bin50") sum(x$y) else NA_integer_,
    qcovar = paste(ck$qn, collapse = ","), covar = paste(ck$cn, collapse = ","), n_cov = length(ck$qn) + length(ck$cn),
    dropped = paste(ck$dropped, collapse = ","),
    pheno_file = if (r$engine %in% c("gcta", "check")) paste0(base, "_pheno.txt") else paste0(base, "_regenie_pheno.txt"),
    qcovar_file = if (r$engine %in% c("gcta", "check")) paste0(base, "_qcovar.txt") else paste0(base, "_regenie_covar.txt"),
    covar_file = if (r$engine %in% c("gcta", "check") && length(ck$cn)) paste0(base, "_covar.txt") else "",
    keep_file = paste0(base, "_keep.txt"), ref_result = "")
}
runs_tab <- rbindlist(tab)

ex <- rbindlist(lapply(pc_all_sets, function(k) data.table(
  id = if (k == main_pc) "main" else sprintf("pc%d", k), group = if (k == main_pc) "main" else "qc", engine = "existing",
  desc = if (k == main_pc) sprintf("주 분석 (PC-AiR PC %d개)", k) else sprintf("PC %d개 (04-5)", k), y = "rint", n = nrow(d),
  ref_result = file.path(gwas_045, sprintf("rint_maf0.01_PC%d.mlma.gz", k)))), fill = TRUE)
mf <- rbindlist(lapply(c("0.05", "0.1"), function(m) data.table(
  id = paste0("maf", m), group = "qc", engine = "maf_filter", desc = sprintf("MAF %s 세트 SNP로 주 분석 결과 거름", m), y = "rint",
  n = nrow(d), ref_result = file.path(qc_dir, sprintf("03-3_maf%s_hwe.bim", m)))), fill = TRUE)
runs_tab <- rbind(ex, runs_tab, mf, fill = TRUE)
for (cc in names(runs_tab)) if (is.character(runs_tab[[cc]])) set(runs_tab, which(is.na(runs_tab[[cc]])), cc, "")
missing_ref <- runs_tab[ref_result != "" & !file.exists(ref_result)]
if (nrow(missing_ref)) { print(missing_ref[, .(id, ref_result)]); warning("참고 결과 파일이 없음 (04-5 PC 개수 결과·MAF 세트 bim 확인)") }
fwrite(runs_tab, paste0(out_prefix, "_runs.tsv"), sep = "\t")

head(runs_tab)

sink(paste0(out_prefix, "_summary.txt"))
cat(sprintf("05-1 민감도 분석 모형: %d개 (GCTA %d, REGENIE %d, 기존 결과 %d, MAF 필터 %d) / 주 분석 PC %d개\n\n",
            nrow(runs_tab), sum(runs_tab$engine == "gcta"), sum(grepl("regenie", runs_tab$engine)),
            sum(runs_tab$engine == "existing"), sum(runs_tab$engine == "maf_filter"), main_pc))
print(runs_tab[, .(id, group, engine, n, n_case, n_removed_missing, n_cov, dropped, desc)], row.names = FALSE, nrows = 200)
sink()
cat(readLines(paste0(out_prefix, "_summary.txt")), sep = "\n")
