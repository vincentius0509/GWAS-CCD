# 09-2_size_age.R - Tests body-size polygenic score x age in a GREML mixed model. See README.md.
library(data.table)
library(ggplot2)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_keep <- as.integer(Sys.getenv("TEST_KEEP", "0"))
tag <- if (test_keep > 0) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("09-2", tag)); pic_prefix <- file.path(pic_dir, paste0("09-2", tag))
work_dir <- paste0(out_prefix, "_work")
pgs_f <- c(file.path(in_dir, "09-1_size_pgs.csv"), file.path(qc_dir, "09-1_size_pgs.csv")); pgs_f <- pgs_f[file.exists(pgs_f)][1]
apc_f <- c(file.path(in_dir, "08-1_ancestry_pcs.txt"), file.path(qc_dir, "08-1_ancestry_pcs.txt")); apc_f <- apc_f[file.exists(apc_f)][1]
cslb_file <- file.path(base_dir, "SurveyData/v2023/DAP_2023_CSLB_v1.0.csv"); dog_file <- file.path(base_dir, "DogOverview/v2023/DAP_2023_DogOverview_v1.1.csv")
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
grm_full <- file.path(in_dir, "04-1_grm", "grm_all")
main_f <- c(pheno = file.path(in_dir, "04-2_pheno.txt"), qcovar = file.path(in_dir, "03-5_qcovar_PC5.txt"), covar = file.path(in_dir, "03-5_covar.txt"))
a8_pheno <- file.path(in_dir, "05-1_inputs", "age8_pheno.txt")
for (f in c(pgs_f, apc_f, cslb_file, dog_file, gcta, paste0(grm_full, ".grm.bin"), main_f, a8_pheno)) if (is.na(f) || !file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
qn <- c("age_cslb", "age_hinge", "pa_score", "wt_cov", "owner_age", paste0("PC", 1:5)); cn <- c("Sex", "Intact", "Vision", "Cohort_BHS")
q_nowt <- setdiff(qn, "wt_cov")
wr <- function(d, f) { fwrite(d, f, sep = " ", col.names = FALSE, quote = FALSE); f }
rd <- function(f, nm) { d <- fread(f, header = FALSE, colClasses = list(character = 1:2)); setnames(d, c("FID", "IID", nm)); d }

grm <- grm_full
if (test_keep > 0) {
  ids <- fread(paste0(grm_full, ".grm.id"), header = FALSE, colClasses = "character")
  kf <- wr(ids[seq(1, .N, by = test_keep)], file.path(work_dir, "test_keep.txt")); grm <- file.path(work_dir, "grm_test")
  if (!file.exists(paste0(grm, ".grm.bin"))) system2(gcta, c("--grm", grm_full, "--keep", kf, "--make-grm", "--out", grm), stdout = FALSE)
}
gid <- fread(paste0(grm, ".grm.id"), header = FALSE, colClasses = "character")$V2
M <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")), list(rd(main_f["pheno"], c("y", "score"))[, score := NULL], rd(main_f["qcovar"], qn),
            fread(main_f["covar"], header = FALSE, colClasses = "character", col.names = c("FID", "IID", cn))))
pg <- fread(pgs_f, colClasses = c(IID = "character", FID = "character"))
M <- merge(M, pg[, .(IID, PGS_auto, PGS_withX)], by = "IID")[IID %in% gid]
M[, `:=`(PGS_z = as.numeric(scale(PGS_auto)), PGSX_z = as.numeric(scale(PGS_withX)), wt_z = as.numeric(scale(wt_cov)), age_c = age_cslb - mean(age_cslb))]
M[, `:=`(PGS_x_age = PGS_z * age_c, PGSX_x_age = PGSX_z * age_c, wt_x_age = wt_z * age_c)]
apc <- rd(apc_f, paste0("aPC", 1:20)); M <- merge(M, apc[, -"FID"], by = "IID")
for (k in 1:5) M[, (paste0("aPC", k, "_x_age")) := get(paste0("aPC", k)) * age_c]
a8 <- rd(a8_pheno, "y8"); A <- merge(M, a8[, .(IID, y8)], by = "IID")[, y := y8]
cat(sprintf("전체 %d마리 · 8세 이상 %d마리. PGS_auto–체중 상관 %.3f\n", nrow(M), nrow(A), cor(M$PGS_auto, M$wt_cov)))

greml_fix <- function(name, d, qbase, extra) {
  out <- file.path(work_dir, name)
  kq <- qbase[vapply(qbase, function(j) isTRUE(sd(d[[j]]) > 0), TRUE)]; kc <- cn[vapply(cn, function(j) uniqueN(d[[j]]) > 1, TRUE)]
  qcols <- c(kq, extra)
  if (!file.exists(paste0(out, ".stdout")) || !any(grepl("Estimates of fixed effects", readLines(paste0(out, ".stdout"), warn = FALSE)))) {
    wr(d[, .(FID, IID, y)], paste0(out, ".pheno")); wr(d[, c("FID", "IID", qcols), with = FALSE], paste0(out, ".qcovar"))
    wr(d[, c("FID", "IID", kc), with = FALSE], paste0(out, ".covar")); wr(d[, .(FID, IID)], paste0(out, ".keep"))
    system2(gcta, c("--reml", "--grm", grm, "--keep", paste0(out, ".keep"), "--pheno", paste0(out, ".pheno"), "--qcovar", paste0(out, ".qcovar"),
                    "--covar", paste0(out, ".covar"), "--reml-est-fix", "--reml-no-lrt", "--reml-maxit", "200", "--out", out, "--threads", n_cores),
            stdout = paste0(out, ".stdout"), stderr = paste0(out, ".stdout"))
  }
  l <- readLines(paste0(out, ".stdout"), warn = FALSE); i <- grep("^Estimates of fixed effects", l)
  if (!length(i)) { warning(name, ": GREML 실패 (", out, ".stdout)"); return(NULL) }
  fx <- fread(text = l[(i + 2):length(l)], header = TRUE, fill = TRUE)[grepl("^X_", Source) | Source == "mean"]
  pos <- match(extra, qcols) + 1L
  data.table(model = name, term = extra, beta = as.numeric(fx[match(paste0("X_", pos), Source), Estimate]),
             se = as.numeric(fx[match(paste0("X_", pos), Source), SE]), n = nrow(d), grm = TRUE)
}
ols_fix <- function(name, d, qbase, extra) {
  kc <- cn[vapply(cn, function(j) uniqueN(d[[j]]) > 1, TRUE)]
  s <- summary(lm(as.formula(paste("y ~", paste(c(qbase, kc, extra), collapse = " + "))), data = d))$coefficients
  data.table(model = name, term = extra, beta = s[extra, 1], se = s[extra, 2], n = nrow(d), grm = FALSE)
}
res <- rbindlist(list(
  greml_fix("P_all_PGSxAge", M, q_nowt, c("PGS_z", "PGS_x_age")),
  greml_fix("S1_all_PGSxAge_plusWeight", M, qn, c("PGS_z", "PGS_x_age")),
  greml_fix("S2_all_WeightxAge", M, q_nowt, c("wt_z", "wt_x_age")),
  greml_fix("S3_under12_PGSxAge", M[age_cslb < 12], q_nowt, c("PGS_z", "PGS_x_age")),
  greml_fix("S4_all_PGSwithXxAge", M, q_nowt, c("PGSX_z", "PGSX_x_age")),
  greml_fix("S5_age8_PGS", A, q_nowt, "PGS_z"),
  ols_fix("S6_all_PGSxAge_noGRM", M, q_nowt, c("PGS_z", "PGS_x_age")),
  greml_fix("S7a_all_PGSxAge_plusAncPC", M, c(q_nowt, paste0("aPC", 1:20)), c("PGS_z", "PGS_x_age")),
  greml_fix("S7b_all_PGSxAge_plusAncPCxAge", M, c(q_nowt, paste0("aPC", 1:20), paste0("aPC", 1:5, "_x_age")), c("PGS_z", "PGS_x_age"))), fill = TRUE)

cs <- fread(cslb_file, select = c("dog_id", "cslb_date", "cslb_score"), colClasses = c(dog_id = "character"))[dog_id %in% M$IID]
dob <- fread(dog_file, select = c("dog_id", "Estimated_DOB"), colClasses = c(dog_id = "character"))
cs <- merge(cs, dob, by = "dog_id")[, age := as.numeric(as.IDate(cslb_date) - as.IDate(Estimated_DOB)) / 365.25][!is.na(cslb_score) & !is.na(age)]
setorder(cs, dog_id, age)
ch <- cs[, .(n = .N, age_first = age[1], age_last = age[.N], s_first = cslb_score[1], s_last = cslb_score[.N]), by = .(IID = dog_id)][n >= 2 & age_last - age_first >= 0.75]
ch[, rate := (s_last - s_first) / (age_last - age_first)]
ch <- merge(ch, M[, c("IID", "PGS_z", "Sex", "Intact", paste0("PC", 1:5)), with = FALSE], by = "IID")
long_fit <- function(name, d) {
  if (nrow(d) < 50) return(NULL)
  s <- summary(lm(as.formula(paste("rate ~ PGS_z + age_first + s_first + Sex + Intact +", paste(paste0("PC", 1:5), collapse = " + "))), data = d))$coefficients
  data.table(model = name, term = "PGS_z (연간 변화량, 점/년)", beta = s["PGS_z", 1], se = s["PGS_z", 2], n = nrow(d), grm = FALSE)
}
res <- rbind(res, long_fit("S8_long_all", ch), long_fit("S8_long_first8plus", ch[age_first >= 8]), fill = TRUE)
res[, z := beta / se][, p := 2 * pnorm(-abs(z))]
fwrite(res, paste0(out_prefix, "_results.csv"))
P <- res[model == "P_all_PGSxAge" & term == "PGS_x_age"]; H5 <- isTRUE(P$beta > 0 && P$p < 0.05)

M[, pgs_t := cut(PGS_z, quantile(PGS_z, 0:3 / 3), include.lowest = TRUE, labels = c("small (low PGS)", "medium", "large (high PGS)"))]
M[, age_bin := pmin(floor(age_cslb), 15)]
pd <- M[, .(y = mean(y), se = sd(y) / sqrt(.N), n = .N), by = .(pgs_t, age_bin)][n >= 15]
g1 <- ggplot(pd, aes(age_bin, y, colour = pgs_t)) + geom_line() + geom_pointrange(aes(ymin = y - 1.96 * se, ymax = y + 1.96 * se), size = 0.3) +
  scale_colour_manual(values = c(`small (low PGS)` = "#2166ac", medium = "grey50", `large (high PGS)` = "#b2182b")) +
  labs(x = "age at CSLB (years; 15 = 15+)", y = "CSLB (RINT), mean ± 95% CI", colour = "body-size PGS tertile",
       title = "Does CSLB rise faster with age in genetically larger dogs?", subtitle = "descriptive; the test is the GREML PGS × age term") + theme_minimal(base_size = 10)
ggsave(paste0(pic_prefix, "_pgs_age.png"), g1, width = 7, height = 4.5, dpi = 300)

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("09-2 가설 H5: 체격 유전 점수 × 나이와 CSLB (2026-11-13 추가, 판정 규칙은 실행 전 헤더에 고정)\n"); if (test_keep > 0) cat("*** 시험. 해석하지 말 것 ***\n")
cat("\n■ 계수 (PGS_z = 체격 PGS 1 SD, *_x_age = 1 SD × 나이 1년; CSLB_rint 단위, p = 양측 Wald)\n"); print(f3(res))
cat(sprintf("\n■ 판정 (주 검정 P_all_PGSxAge의 PGS_x_age): β %.3g ± %.3g, p %.3g → %s\n", P$beta, P$se, P$p,
            if (H5) "H5 지지 (Holm 보정 전)" else "H5 지지 안 됨"))
cat("  해석 주의: 생존 선택 (큰 개는 일찍 죽어 노령 표본에 덜 남음) 때문에 관찰된 개에서의 연관으로만 해석 (S3 참고)\n",
    "  PGS는 품종 평균 체중 기반 가중치라 사실상 '크기별 품종 계통'을 담고 (S7로 분리 시도), LCORL 등 대형견 locus 3개가 없어 대형견 쪽을 덜 담음\n",
    "  주 검정은 단면 (최근 1회 점수 × 측정 나이); 같은 개 안 변화는 S8 (CSLB 원점수 단위, 점/년; 처음 점수 · 처음 나이 보정)\n",
    "  S1 (체중 포함)과 주 검정의 PGS_x_age 비교 = 체중 매개 정도. 인과 (MR) 해석은 하지 않음 (IGF1 · GHR 등 성장인자 유전자의 다면 효과)\n", sep = "")
sink()
cat(sprintf("저장: %s_results.csv · _summary.txt / 그림 %s_pgs_age.png\n", out_prefix, pic_prefix))

