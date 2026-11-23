# 05-5_model_catalog.R - Catalogue of models (what differs from the primary model) and quality checks. See README.md.
library(data.table)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
runs_tsv <- file.path(qc_dir, "05-1_runs.tsv")
summ_csv <- file.path(qc_dir, "05-4_run_summary.csv")
gcta_dir <- file.path(qc_dir, "05-2_gwas")
regenie_dir <- file.path(qc_dir, "05-3_regenie")
gwas_045 <- file.path(qc_dir, "04-5_gwas")
out_prefix <- file.path(qc_dir, "05-5")
for (f in c(runs_tsv, summ_csv)) if (!file.exists(f)) stop("파일 없음: ", f)

lam_lo <- 0.95; lam_hi <- 1.05
tail_1e5 <- 2; tail_1e4 <- 1.5

runs <- fread(runs_tsv, na.strings = c("", "NA"))
summ <- fread(summ_csv, na.strings = c("", "NA"))
runs <- runs[id != "main_check"]
if (!setequal(runs$id, summ$id)) stop("05-1 실행 표와 05-4 요약의 모형이 다름: ", paste(symdiff <- union(setdiff(runs$id, summ$id), setdiff(summ$id, runs$id)), collapse = ", "))

main_row <- fread(runs_tsv, na.strings = c("", "NA"))[id == "main_check"]
split_v <- function(x) if (is.na(x) || x == "") character(0) else strsplit(x, ",", fixed = TRUE)[[1]]
main_q <- split_v(main_row$qcovar); main_c <- split_v(main_row$covar)
main_pc <- grep("^PC[0-9]+$", main_q, value = TRUE)

pc_k <- function(id) as.integer(sub("^pc", "", id))
for (i in which(runs$engine == "existing" & runs$id != "main")) {
  k <- pc_k(runs$id[i])
  runs$qcovar[i] <- paste(c(setdiff(main_q, main_pc), if (k > 0) paste0("PC", seq_len(k))), collapse = ",")
  runs$covar[i] <- paste(main_c, collapse = ",")
}
for (i in which(runs$id == "main" | runs$engine == "maf_filter")) {
  runs$qcovar[i] <- paste(main_q, collapse = ","); runs$covar[i] <- paste(main_c, collapse = ",")
}

pc_set <- function(v) grep("PC[0-9]+$", v, value = TRUE)
desc_pc <- function(v) {
  p <- pc_set(v); if (!length(p)) return("PC 없음")
  pre <- unique(sub("PC[0-9]+$", "", p)); lab <- setNames(c("PC-AiR", "all PCA", "Q-matrix PCA"), c("none", "all_", "qmat_"))[ifelse(pre == "", "none", pre)]
  sprintf("%s PC %d개", ifelse(is.na(lab), pre, lab), length(p))
}
cov_diff <- function(q, cv) {
  v <- c(split_v(q), split_v(cv)); m <- c(main_q, main_c)
  add <- setdiff(setdiff(v, pc_set(v)), m); rem <- setdiff(setdiff(m, pc_set(m)), v)
  pc_change <- !setequal(pc_set(v), pc_set(m))
  repl <- character(0)
  for (a in add) { r <- rem[startsWith(a, rem)]; if (length(r) == 1) { repl <- c(repl, paste0(r, " → ", a)); add <- setdiff(add, a); rem <- setdiff(rem, r) } }
  list(add = paste(add, collapse = ", "), rem = paste(rem, collapse = ", "),
       repl = paste(c(repl, if (pc_change) paste0(desc_pc(m), " → ", desc_pc(v))), collapse = ", "))
}
d <- runs[, cov_diff(qcovar, covar), by = id]
cat_tab <- merge(runs[, .(id, group, engine, y, n, n_case, n_removed_missing, dropped, desc)], d, by = "id", sort = FALSE)

N_main <- runs[id == "main", n]
y_lab <- c(rint = "CSLB RINT (주 분석과 같음)", raw = "CSLB 원점수", bin50 = "CSLB ≥ 50 이분형")
eng_lab <- c(existing = "GCTA MLMA-LOCO (04-5 결과)", gcta = "GCTA MLMA-LOCO", regenie_bt = "REGENIE Firth 로지스틱",
             regenie_int = "REGENIE SNP × age_cslb 상호작용", maf_filter = "주 분석 결과를 SNP 세트로 거름")
cat_tab[, `:=`(
  phenotype = y_lab[y],
  test = eng_lab[engine],
  n_excluded = N_main - n,
  subset = fcase(group %in% c("age", "genetic", "strata", "subset") & n < N_main, desc,
                 id == "bin_age8", "8세 이상", id == "add_ac", "AFUS 응답 개체 (AC 공변량 있는 개체)",
                 id == "add_ses", "소득 결측 제외", default = "")
)]
cat_tab[id == "maf0.05", snp_set := "MAF ≥ 0.05 SNP만"]
cat_tab[id == "maf0.1", snp_set := "MAF ≥ 0.1 SNP만"]
cat_tab[is.na(snp_set), snp_set := ""]
cat_tab[is.na(dropped) | dropped %in% c("", "\"\""), dropped := ""]
cat_tab[, auto_dropped := fifelse(dropped == "", "", paste0(dropped, " (하위집합에서 상수·공선성 → 05-1이 자동 제외)"))]
cat_tab[dropped != "", rem := vapply(seq_len(.N), function(i) paste(setdiff(strsplit(rem[i], ", ", fixed = TRUE)[[1]], strsplit(dropped[i], ",", fixed = TRUE)[[1]]), collapse = ", "), character(1))]

s <- summ[, .(id, n_snp, lambda, n_p5e8, n_p1e5, n_p1e4, z_cor_main, z_cor_expect, ref_1e5, ref_1e4, n_new_loci)]
cat_tab <- merge(cat_tab, s, by = "id", sort = FALSE)
cat_tab[, `:=`(tail_ratio_1e5 = n_p1e5 / (n_snp * 1e-5), tail_ratio_1e4 = n_p1e4 / (n_snp * 1e-4))]
cat_tab[, lambda_gate := !is.na(lambda) & lambda >= lam_lo & lambda <= lam_hi]
cat_tab[, tail_flag := tail_ratio_1e5 > tail_1e5 | tail_ratio_1e4 > tail_1e4]

count_constrained <- function(files) sum(vapply(files, function(f) any(grepl("constrained", readLines(f, warn = FALSE), ignore.case = TRUE)), logical(1)))
snp_line <- function(f) { x <- readLines(f, n = 3, warn = FALSE); m <- regmatches(x, regexpr("결과 SNP [0-9]+ / bim [0-9]+", x)); if (length(m)) m else NA_character_ }
chk <- rbindlist(lapply(cat_tab$id, function(rid) {
  eng <- cat_tab[id == rid, engine]
  if (eng == "existing") {
    k <- if (rid == "main") 5L else pc_k(rid); stem <- sprintf("rint_maf0.01_PC%d", k)
    lg <- Sys.glob(file.path(gwas_045, "log", paste0(stem, "_chr*.log"))); sm <- file.path(gwas_045, paste0(stem, "_summary.txt"))
    data.table(id = rid, n_chr_log = length(lg), n_constrained = count_constrained(lg), snp_check = snp_line(sm), check_from = "04-5")
  } else if (eng == "gcta") {
    lg <- Sys.glob(file.path(gcta_dir, "log", paste0(rid, "_chr*.log"))); sm <- file.path(gcta_dir, paste0(rid, "_summary.txt"))
    data.table(id = rid, n_chr_log = length(lg), n_constrained = count_constrained(lg), snp_check = snp_line(sm), check_from = "05-2")
  } else if (eng %in% c("regenie_bt", "regenie_int")) {
    lg <- Sys.glob(file.path(regenie_dir, "log", paste0(rid, "_chr*.log")))
    n_used <- vapply(lg, function(f) { x <- grep("number of individuals used in analysis", readLines(f, warn = FALSE), value = TRUE); if (length(x)) as.integer(sub(".*= *", "", x[length(x)])) else NA_integer_ }, integer(1))
    ok_n <- all(n_used == cat_tab[id == rid, n])
    data.table(id = rid, n_chr_log = length(lg), n_constrained = NA_integer_,
               snp_check = sprintf("염색체 %d개 분석 개체 = n: %s", length(lg), if (isTRUE(ok_n)) "예" else "아니오"), check_from = "05-3")
  } else data.table(id = rid, n_chr_log = NA_integer_, n_constrained = NA_integer_, snp_check = "주 분석 결과를 거름 (실행 없음)", check_from = "04-5")
}))
cat_tab <- merge(cat_tab, chk, by = "id", sort = FALSE)
cat_tab[, run_check := fcase(
  engine == "maf_filter", TRUE,
  engine %in% c("regenie_bt", "regenie_int"), n_chr_log == 38 & grepl(": 예$", snp_check),
  default = n_chr_log == 38 & !is.na(snp_check) & {
    m <- regmatches(snp_check, gregexpr("[0-9]+", snp_check)); vapply(m, function(z) length(z) == 2 && z[1] == z[2], logical(1)) })]
cat_tab[, gate := lambda_gate & run_check]

cat_tab[, caution := {
  x <- character(.N)
  x <- fifelse(tail_flag, paste0(x, sprintf("QQ 꼬리 부풀림 (p<1e-5 %.1f배, p<1e-4 %.1f배); ", tail_ratio_1e5, tail_ratio_1e4)), x)
  x <- fifelse(!is.na(n_constrained) & n_constrained > 0, paste0(x, sprintf("V(G) 경계값 고정 %d/38 염색체 (GRM 없는 일반 회귀); ", n_constrained)), x)
  x <- fifelse(!is.na(lambda) & lambda_gate & (lambda - lam_lo < 0.005 | lam_hi - lambda < 0.005), paste0(x, "λ가 관문 경계 근처; "), x)
  x <- fifelse(engine == "regenie_bt", paste0(x, sprintf("case %s개 (검정력 낮음); ", n_case)), x)
  sub("; $", "", x)
}]

cat_tab[, role := fcase(
  id == "main", "주 분석",
  group == "drop_cov", "진단용",
  group %in% c("qc", "add_cov", "subset") | id == "raw", "강건성 확인",
  group %in% c("phenotype", "age", "genetic", "strata"), "다른 질문",
  default = "기타")]
cat_tab[, role_note := fcase(
  role == "진단용", "일부러 보정을 뺀 모형 → 결과를 해석에 쓰지 않음 (drop_age는 나이 교란 포함)",
  role == "강건성 확인", "주 분석과 같은 질문 → 결과가 흔들리지 않는지 확인용",
  id %in% c("bin_all", "bin_age8"), "임상 기준(CSLB ≥ 50) 정의",
  id == "int_age", "나이에 따라 효과가 달라지는지 직접 검정",
  id == "young8", "8세 미만 (age8과 겹치지 않는 짝 → 노령 특이 효과 이질성 검정)",
  group == "age", "나이 든 개로 제한 (나이 의존 효과)",
  group == "strata", "성별·체중 층 (층끼리 겹치지 않음 → 효과 이질성 검정)",
  group == "genetic", "혈통 구성 (순종형·잡종형)",
  default = "")]

role_ord <- c("주 분석", "다른 질문", "강건성 확인", "진단용", "기타")
cat_tab[, role := factor(role, levels = role_ord)]
setorder(cat_tab, role, group)
out <- cat_tab[, .(id, role, group, desc, phenotype, test, n, n_excluded, n_case, subset, snp_set,
                   cov_added = add, cov_removed = rem, cov_replaced = repl, auto_dropped,
                   lambda = round(lambda, 4), lambda_gate, run_check, gate, caution,
                   tail_ratio_1e5 = round(tail_ratio_1e5, 2), tail_ratio_1e4 = round(tail_ratio_1e4, 2),
                   n_constrained, z_cor_main = round(z_cor_main, 3), z_cor_expect = round(z_cor_expect, 3),
                   n_p1e5, n_p1e4, ref_1e5, ref_1e4, n_new_loci, check_from, snp_check, role_note)]
fwrite(out, paste0(out_prefix, "_model_catalog.csv"), bom = TRUE)

change_txt <- function(r) {
  x <- c(if (r$phenotype != y_lab["rint"]) paste0("표현형: ", r$phenotype),
         if (!grepl("GCTA", r$test)) paste0("검정: ", r$test),
         if (r$subset != "") sprintf("대상: %s (n %s, 주 분석보다 %s마리 적음)", r$subset, format(r$n, big.mark = ","), format(r$n_excluded, big.mark = ",")),
         if (r$snp_set != "") paste0("SNP: ", r$snp_set),
         if (r$cov_added != "") paste0("추가: ", r$cov_added),
         if (r$cov_removed != "") paste0("제외: ", r$cov_removed),
         if (r$cov_replaced != "") paste0("교체: ", r$cov_replaced),
         if (r$auto_dropped != "") paste0("자동 제외: ", r$auto_dropped))
  if (!length(x)) "주 분석 설정 그대로" else paste(x, collapse = " / ")
}
sink(paste0(out_prefix, "_model_catalog.txt"))
cat("05-5 모형 목록 (주 분석: GCTA MLMA-LOCO, CSLB RINT, n 7,115, PC-AiR PC 5개)\n")
cat("주 분석 공변량: 연속형", paste(main_q, collapse = ", "), "/ 범주형", paste(main_c, collapse = ", "), "\n")
cat(sprintf("품질 관문: λ %.2f ~ %.2f (경계 포함) + 실행 검사 통과 → 통과 %d / %d 모형\n", lam_lo, lam_hi, sum(out$gate), nrow(out)))
cat(sprintf("주의 표시: QQ 꼬리 (p<1e-5 관측/기대 > %.1f 또는 p<1e-4 > %.1f), V(G) 경계값 고정, λ 경계 ±0.005, 이분형 case 수\n", tail_1e5, tail_1e4))
for (rl in role_ord) {
  sub_o <- out[role == rl]; if (!nrow(sub_o)) next
  cat(sprintf("\n== %s (%d개) ==\n", rl, nrow(sub_o)))
  for (i in seq_len(nrow(sub_o))) {
    r <- sub_o[i]
    cat(sprintf("\n[%s] %s  %s\n", if (r$gate) "통과" else "탈락", r$id, r$desc))
    cat("  변경:", change_txt(r), "\n")
    cat(sprintf("  λ %.3f%s | p<1e-5 %d개 (기대의 %.2f배) | p<1e-4 %d개 (%.2f배)%s\n", r$lambda, if (r$lambda_gate) "" else " ← λ 관문 밖",
                r$n_p1e5, r$tail_ratio_1e5, r$n_p1e4, r$tail_ratio_1e4,
                if (!is.na(r$z_cor_main)) sprintf(" | z 상관 %.3f (기대 %s)", r$z_cor_main, ifelse(is.na(r$z_cor_expect), "-", sprintf("%.3f", r$z_cor_expect))) else ""))
    cat(sprintf("  실행 검사 (%s): %s — %s\n", r$check_from, if (r$run_check) "통과" else "확인 필요", r$snp_check))
    if (r$caution != "") cat("  주의:", r$caution, "\n")
    if (r$role_note != "") cat("  역할:", r$role_note, "\n")
  }
}
cat("\n== 관문 탈락 모형 ==\n")
fail <- out[gate == FALSE]
if (!nrow(fail)) cat("없음\n") else print(fail[, .(id, role, desc, lambda, lambda_gate, run_check)], row.names = FALSE)
cat("\n== 관문 통과 + 다른 질문 (06 후보) ==\n")
print(out[gate == TRUE & role == "다른 질문", .(id, desc, n, lambda, caution)], row.names = FALSE)
sink()

cat(sprintf("저장: %s_model_catalog.csv / .txt  (관문 통과 %d / %d)\n", out_prefix, sum(out$gate), nrow(out)))

