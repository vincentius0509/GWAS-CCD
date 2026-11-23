# 13-11_pathway_by_model.R - Tabulates the 13-10 results by model (no new computation). See README.md.
suppressMessages(library(data.table))
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
pre <- function(x) file.path(qc_dir, x)
f_in <- c(tests = pre("13-10_tests.csv"), cand = pre("13-10_candidates.csv"))
for (f in f_in) if (!file.exists(f)) stop("파일 없음: ", f)
out_dir <- pre("13-11_by_model"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
run_order <- c("young8", "main", "agefloor6", "agefloor7", "age8", "agefloor9", "age10", "agefloor11", "agefloor12", "wt_t1", "wt_t2", "wt_t3",
               "sex_m", "sex_f", "purelike", "mixlike", "int_age", "int_wt", "bin_all", "bin_age8")
desc_of <- c(young8 = "8세 미만", main = "전체", agefloor6 = "6세 이상", agefloor7 = "7세 이상", age8 = "8세 이상", agefloor9 = "9세 이상", age10 = "10세 이상",
             agefloor11 = "11세 이상", agefloor12 = "12세 이상", wt_t1 = "체중 하위 1/3", wt_t2 = "체중 중간 1/3", wt_t3 = "체중 상위 1/3", sex_m = "수컷", sex_f = "암컷",
             purelike = "순종형", mixlike = "잡종형", int_age = "SNP × 나이", int_wt = "SNP × 체중", bin_all = "이진 전체", bin_age8 = "이진 8세 이상")

TT <- fread(f_in[["tests"]]); CC <- fread(f_in[["cand"]])
X <- merge(TT, CC[, .(cid, source, same_loci_sets, loci_snps)], by = "cid", all.x = TRUE)
X <- X[, .(model, tier, rank, pathway = set_name, set_id = set, source, set_size, n_loci, nominating_genes = genes, same_loci_pathways = same_loci_sets,
           genes_left = n_left, beta, p_main = p, sig = fifelse(main_pass %in% TRUE, "*", ""), verdict, pair_model = pair, pair_rule = crit, pair_p = pair_P, diff_p,
           p_x50 = P_x50, p_x10 = P_x10, p_x0_circular = P_x0, p_magma_w50 = P_w50_x100, q_ref_tier, q_ref_all, nominating_loci = loci_snps)]
X[, tier_o := fifelse(tier == "주 후보", 1L, 2L)]
setorder(X, model, tier_o, rank); X[, tier_o := NULL]
num <- names(X)[vapply(X, function(v) is.numeric(v) && !is.integer(v), TRUE)]
for (j in num) set(X, j = j, value = signif(X[[j]], 3))

ms <- intersect(run_order, unique(X$model))
for (i in seq_along(ms)) fwrite(X[model == ms[i]], file.path(out_dir, sprintf("%02d_%s.csv", i, ms[i])), bom = TRUE)
OV <- X[, .(model_desc = desc_of[model[1]], n_tested = .N, n_primary = sum(tier == "주 후보"), primary_p05 = sum(tier == "주 후보" & sig == "*"),
            secondary_p05 = sum(tier == "보조 후보" & sig == "*"), robust = sum(verdict == "견고"),
            primary_hits = paste(pathway[tier == "주 후보" & sig == "*"], collapse = " | "),
            secondary_hits = paste(pathway[tier == "보조 후보" & sig == "*"], collapse = " | ")), by = model]
OV <- OV[match(ms, model)]
fwrite(OV, file.path(out_dir, "00_overview.csv"), bom = TRUE)
SG <- X[sig == "*"][order(match(model, ms), tier, p_main)]
fwrite(SG, file.path(out_dir, "00_significant.csv"), bom = TRUE)
cat(sprintf("저장: %s/ (모형 %d개 파일 + 00_overview.csv + 00_significant.csv) · 검정한 경로 %d개 · 주 판정 p < 0.05 %d개 (주 후보 %d · 보조 %d) · 견고 %d개\n",
            out_dir, length(ms), nrow(X), nrow(SG), SG[tier == "주 후보", .N], SG[tier == "보조 후보", .N], X[verdict == "견고", .N]))
