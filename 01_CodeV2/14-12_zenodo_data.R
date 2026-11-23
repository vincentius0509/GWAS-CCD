# 14-12_zenodo_data.R - Builds the Zenodo data folder (summary statistics, gene and gene-set results, dog and variant IDs). See README.md.
suppressMessages(library(data.table))
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
out <- file.path(d0, "14-12_zenodo_upload")
f <- function(x) { p <- file.path(d0, x); if (!file.exists(p)) stop("파일 없음: ", p); p }
chk <- function(ok, msg) { if (!isTRUE(ok)) stop("점검 실패: ", msg); cat("  ✔", msg, "\n") }
hangul <- "[가-힣]"

if (dir.exists(out)) { stopifnot(basename(out) == "14-12_zenodo_upload"); unlink(out, recursive = TRUE) }
sub_dirs <- c("01_GWAS_GCTA", "02_GWAS_REGENIE", "03_heterogeneity", "04_MAGMA_genes", "05_MAGMA_sets", "06_IDs")
for (s in sub_dirs) dir.create(file.path(out, s), recursive = TRUE)
put <- function(src, dst) {
  to <- file.path(out, dst); ok <- suppressWarnings(file.link(src, to)); if (!ok) ok <- file.copy(src, to)
  if (!ok) stop("못 넣음: ", src); invisible(to) }
wcsv <- function(x, dst) fwrite(x, file.path(out, dst))

run4 <- fread(f("05-4_run_summary.csv"))[, .(id, n)]
chk(!anyDuplicated(run4$id), "05-4 모형 id 하나씩")
M <- rbindlist(list(
  data.table(model = "main", group = "primary", description = "All dogs; primary model (RINT CSLB score, 5 PC-AiR PCs)"),
  data.table(model = c("young8", "age8"), group = "stratum: age", description = c("Dogs aged < 8 years", "Dogs aged >= 8 years")),
  data.table(model = c("sex_m", "sex_f"), group = "stratum: sex", description = c("Males", "Females")),
  data.table(model = c("purelike", "mixlike"), group = "stratum: genomic inbreeding",
             description = c("Purebred-like dogs (genomic inbreeding COI >= 0.1507)", "Mixed-like dogs (COI < 0.1507)")),
  data.table(model = c("wt_t1", "wt_t2", "wt_t3"), group = "stratum: body weight",
             description = c("Lightest body-weight tertile (<= 15.9 kg)", "Middle tertile (15.9-29.0 kg)", "Heaviest tertile (> 29.0 kg)")),
  data.table(model = c("agefloor6", "agefloor7", "agefloor9", "age10", "agefloor11", "agefloor12"), group = "age floor",
             description = sprintf("Dogs aged >= %d years", c(6, 7, 9, 10, 11, 12))),
  data.table(model = sprintf("pc%d", c(0, 10, 15, 20, 25, 30)), group = "covariate model",
             description = sprintf("All dogs; as primary but with %d PCs", c(0, 10, 15, 20, 25, 30))),
  data.table(model = c("pc_all", "pc_qmat", "add_ancestry", "add_mimic_dx", "add_oa", "add_health", "add_household", "vision_narrow"),
             group = "covariate model",
             description = c("All dogs; 5 PCs computed from all unrelated dogs instead of PC-AiR",
                             "All dogs; 5 PCs computed from breed-ancestry proportions instead of genotypes",
                             "All dogs; primary + largest breed-ancestry proportion and genomic inbreeding (COI)",
                             "All dogs; primary + diagnoses that can mimic cognitive signs (hearing, seizure, vestibular, or brain disease)",
                             "All dogs; primary + osteoarthritis diagnosis",
                             "All dogs; primary + overall health rating",
                             "All dogs; primary + household composition (numbers of dogs and people) and region",
                             "All dogs; vision-impairment covariate replaced by a narrower definition")),
  data.table(model = c("int_age", "int_wt"), group = "interaction (REGENIE)",
             description = c("All dogs; SNP x age interaction", "All dogs; SNP x log2 body weight interaction")),
  data.table(model = c("bin_all", "bin_age8"), group = "binary trait (REGENIE)",
             description = c("All dogs; CSLB score >= 50 vs. < 50 (Firth logistic regression)", "Dogs aged >= 8 years; CSLB score >= 50 vs. < 50"))))
M[, software := fifelse(model %in% c("int_age", "int_wt", "bin_all", "bin_age8"), "REGENIE v4.1.3", "GCTA 1.95.4 (MLMA-LOCO)")]
M[run4, on = .(model = id), n := i.n]
M[model == "int_wt", n := 7115L]
chk(M[, all(!is.na(n))], "모형 34개 표본 수 모두 있음")
M[, in_fig2a := model %in% c("main", "young8", "age8", "sex_m", "sex_f", "purelike", "mixlike", "wt_t1", "wt_t2", "wt_t3", "int_age", "bin_age8")]

cat("■ summary statistics\n")
src_gcta <- function(m) if (m == "main") f("04-5_gwas/rint_maf0.01_PC5.mlma.gz") else if (m %like% "^pc[0-9]+$")
  f(sprintf("04-5_gwas/rint_maf0.01_PC%s.mlma.gz", sub("pc", "", m))) else f(sprintf("05-2_gwas/%s.mlma.gz", m))
src_reg <- c(int_age = "05-3_regenie/int_age.regenie.gz", int_wt = "13-6_regenie/int_wt.regenie.gz",
             bin_all = "05-3_regenie/bin_all.regenie.gz", bin_age8 = "05-3_regenie/bin_age8.regenie.gz")
for (m in M[software %like% "GCTA", model]) put(src_gcta(m), file.path("01_GWAS_GCTA", paste0(m, ".mlma.gz")))
for (m in names(src_reg)) put(f(src_reg[[m]]), file.path("02_GWAS_REGENIE", paste0(m, ".regenie.gz")))
M[, file := fifelse(software %like% "GCTA", file.path("01_GWAS_GCTA", paste0(model, ".mlma.gz")), file.path("02_GWAS_REGENIE", paste0(model, ".regenie.gz")))]
hd <- fread(cmd = paste("zcat", shQuote(file.path(out, "01_GWAS_GCTA/main.mlma.gz")), "| head -3"))
chk(identical(names(hd), c("Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p")), "GCTA 열 = Chr SNP bp A1 A2 Freq b se p (개체 단위 열 없음)")
hr <- fread(cmd = paste("zcat", shQuote(file.path(out, "02_GWAS_REGENIE/int_age.regenie.gz")), "| head -3"))
chk(identical(names(hr)[1:12], c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P")), "REGENIE 열 표준 형식")

for (h in c("sex", "age", "ancestry", "weight")) put(f(sprintf("06-2_het_%s.tsv.gz", h)), file.path("03_heterogeneity", sprintf("het_%s.tsv.gz", h)))

cat("■ MAGMA genes\n")
sym <- unique(fread(cmd = paste("zcat", shQuote(f("06-3_gene_results.csv.gz"))))[, .(GENE, SYMBOL)])
chk(!anyDuplicated(sym$GENE), "GeneID → 기호 1:1")
mag_models <- c(M[in_fig2a == TRUE, model], "agefloor6", "agefloor7", "agefloor9", "age10", "agefloor11", "agefloor12", "int_wt", "bin_all")
g_src <- function(m) if (m %in% c("agefloor6", "agefloor7", "agefloor9", "age10", "agefloor11", "agefloor12", "int_wt", "bin_all"))
  f(sprintf("13-7_magma/%s.genes.out", m)) else f(sprintf("06-3_magma/%s.genes.out", m))
rd_genes <- function(p) fread(p, skip = "GENE")
G <- rbindlist(lapply(mag_models, function(m) rd_genes(g_src(m))[, model := m]))
GN <- rbindlist(lapply(1:20, function(i) rd_genes(f(sprintf("06-3_magma/null%d.genes.out", i)))[, model := paste0("null", i)]))
G <- rbind(G, GN)
G[sym, on = "GENE", SYMBOL := i.SYMBOL]
chk(G[, mean(!is.na(SYMBOL))] == 1, "모든 유전자에 기호 붙음")
G <- G[, .(model, GENE_ID = GENE, SYMBOL, CHR, START, STOP, NSNPS, N, ZSTAT, P)]
chk(G[model == "main", .N] == 18372, "main 유전자 18,372개")
wcsv(G, "04_MAGMA_genes/magma_gene_results.csv.gz")

cat("■ MAGMA gene sets\n")
SI <- fread(f("06-4_set_info.csv"))[kept == TRUE, .(set, family, source, name, n_genes = n_primary, dog_symbols)]
chk(nrow(SI) == 14, "사전 세트 14개")
AI <- fread(f("06-4_aux_info.csv"))
aux_ids <- sub("\\s.*", "", readLines(f("06-4_sets_aux.txt")))
AI <- AI[set %in% aux_ids, .(set, source, name, n_genes)]
chk(nrow(AI) == 5558 && !anyDuplicated(AI$set), "보조 세트 5,558개 이름표")
extra <- c("agefloor6", "agefloor7", "agefloor9", "age10", "agefloor11", "agefloor12", "int_wt", "bin_all")
fig_models <- M[in_fig2a == TRUE, model]
P5 <- fread(cmd = paste("zcat", shQuote(f("06-5_gsa_all.csv.gz"))))[variant == "primary" & kind %in% c("model", "null")]
chk(setequal(P5[kind == "model", unique(run)], fig_models) && P5[kind == "null", uniqueN(run)] == 20, "06-5 사전 세트: Fig. 2a 모형 12개 + null 20개")
P7 <- fread(cmd = paste("zcat", shQuote(f("13-7_gsa_all.csv.gz"))))
chk(all(P7$window %in% c(10, 50)), "13-7 window 값 = 10 · 50")
PS <- rbind(P5[, .(model = run, set, NGENES, BETA, BETA_STD, SE, P)],
            P7[variant == "primary" & window == 10 & run %in% extra, .(model = run, set, NGENES, BETA, BETA_STD, SE, P)])
PS <- merge(PS, SI[, .(set, set_name = name, family)], by = "set")
chk(PS[, .N, by = model][, all(N == 14)] && uniqueN(PS$model) == 40, "사전 세트: 모형 20 + null 20, 각 14줄")
chk(PS[model == "main" & set == "GO_0007623", .N] == 1, "main · 일주기 GO 한 줄")
wcsv(PS[, .(model, set, set_name, family, NGENES, BETA, BETA_STD, SE, P)][order(model, family, set)], "05_MAGMA_sets/prespecified_sets.csv")
A5 <- fread(cmd = paste("zcat", shQuote(f("06-5_aux.csv.gz"))))[variant == "aux" & kind == "model"]
chk(setequal(unique(A5$run), fig_models), "06-5 보조 세트: Fig. 2a 모형 12개")
AS <- rbind(A5[, .(model = run, set, NGENES, BETA, BETA_STD, SE, P)],
            P7[variant == "aux" & run %in% extra, .(model = run, set, NGENES, BETA, BETA_STD, SE, P)])
AS <- merge(AS, AI[, .(set, set_name = name)], by = "set", all.x = TRUE)
chk(AS[, .N, by = model][, all(N == 5558)] && uniqueN(AS$model) == 20 && !anyNA(AS$set_name), "보조 세트: 모형 20개 × 5,558")
wcsv(AS[, .(model, set, set_name, NGENES, BETA, BETA_STD, SE, P)][order(model, set)], "05_MAGMA_sets/supporting_sets.csv.gz")
T10 <- fread(f("13-10_tests.csv"))
tr_v <- c("견고" = "associated; supported in comparison stratum", "경로 수준 근거 (탐색적; 짝 없음)" = "associated; no comparison stratum",
          "경로 수준 근거 (탐색적; 짝 미충족)" = "associated; not supported in comparison stratum", "근거 없음" = "not associated (p >= 0.05)",
          "검정 불가 (남은 유전자 < 10)" = "not tested (< 10 genes left)")
tr_c <- c("공통 재현" = "supported if p < 0.05 in the complementary stratum",
          "층 특이 (노령 > young8)" = "supported if beta is larger than in dogs < 8 years with difference p < 0.05", "짝 없음" = "no comparison stratum")
chk(all(T10$verdict %in% names(tr_v)) && all(T10$crit %in% names(tr_c)), "13-10 판정 · 기준 값이 번역표에 모두 있음")
T10 <- T10[, .(model, set, set_name, set_size, n_loci, nominating_genes = genes, n_genes_left = n_left, BETA = beta, SE = se, P = p,
               comparison_stratum = pair, comparison_criterion = tr_c[crit], comparison_NGENES = pair_NGENES, comparison_BETA = pair_BETA,
               comparison_P = pair_P, diff_z, diff_p, result = tr_v[verdict])]
chk(T10[set == "GO_0071277" & model %in% c("sex_m", "sex_f"), all(result == tr_v[["견고"]])], "칼슘 반응 세트: 두 성별 모두 supported")
wcsv(T10[order(model, set)], "05_MAGMA_sets/locus_nominated_sets.csv")
lab <- setNames(M$description, M$model)
C1 <- fread(f("13-12_tests.csv"))[ver == "U100", .(model, sample = lab[model], NGENES, BETA, SE, P)]
chk(nrow(C1) == 10 && all(C1$NGENES == 76) && !anyNA(C1$sample), "칼슘 축소 세트 (76개) 모형 10개")
LG <- fread(f("13-12_logo.csv"))[, .(model, removed_GENE_ID = GENE, removed_SYMBOL = symbol, removed_gene_z = z, NGENES, BETA, P)]
TK <- fread(f("13-12_topk.csv"))[, .(model, k, removed_SYMBOL = removed, removed_gene_z = z, NGENES, BETA, P)]
U <- unique(LG[model == "sex_m", .(GENE_ID = removed_GENE_ID, SYMBOL = removed_SYMBOL)])
chk(nrow(U) == 76 && LG[, .N, by = model][, all(N == 76)], "LOGO: 모형 3개 × 유전자 76개")
wcsv(C1, "05_MAGMA_sets/calcium_set_reduced_tests.csv")
wcsv(LG, "05_MAGMA_sets/calcium_set_leave_one_gene_out.csv")
wcsv(TK, "05_MAGMA_sets/calcium_set_top_k_removal.csv")
wcsv(U[order(SYMBOL)], "05_MAGMA_sets/calcium_set_reduced_genes.csv")
S5 <- fread(f("14-5_age_pathway_scan.csv"))
S6c <- fread(f("14-11_TableS6c_full.csv"))
chk(nrow(S5) == 5558 && all.equal(S5[, names(S6c), with = FALSE], S6c, check.attributes = FALSE) == TRUE, "14-5 = Table S6c 전체판과 같은 값")
chk(S5[!is.na(diff_z) & abs(diff_z) > 0.01, all(sign(diff_z) == sign(BETA_age8 - BETA_young8))], "diff_z 방향 = β(≥ 8 y) − β(< 8 y) (|z| > 0.01; 저장된 β 반올림 때문)")
wcsv(S5[, .(set, set_name = name, n_genes, P_young8, P_agefloor6, P_agefloor7, P_age8, P_agefloor9, P_age10, P_agefloor11, P_agefloor12,
            P_int_age, BETA_young8, BETA_age8, diff_z, diff_p)], "05_MAGMA_sets/age_comparison_supporting_sets.csv")
put(f("06-4_sets_primary.txt"), "05_MAGMA_sets/set_definitions_prespecified.txt")
put(f("06-4_sets_aux.txt"), "05_MAGMA_sets/set_definitions_supporting.txt")
wcsv(SI, "05_MAGMA_sets/set_info_prespecified.csv")
wcsv(AI, "05_MAGMA_sets/set_info_supporting.csv")

cat("■ IDs\n")
ph <- fread(f("04-2_pheno.txt"), header = FALSE)
chk(nrow(ph) == 7115 && !anyDuplicated(ph$V2) && all(ph$V1 == ph$V2), "분석한 개 7,115마리 (FID = IID = dog_id)")
writeLines(c("dog_id", as.character(sort(ph$V2))), file.path(out, "06_IDs/dog_ids.txt"))
pv <- fread(f("03-3_maf0.01_hwe.pvar"), skip = "#CHROM")
setnames(pv, 1:5, c("CHROM", "POS", "ID", "REF", "ALT"))
chk(nrow(pv) == 4842736 && ncol(pv) == 5, "variant 4,842,736개, 열 5개")
fwrite(pv, file.path(out, "06_IDs/variants.tsv.gz"), sep = "\t")
mm <- fread(cmd = paste("zcat", shQuote(file.path(out, "01_GWAS_GCTA/main.mlma.gz")), "| head -100001"))
mm <- merge(mm[, .(ID = SNP, A1)], pv[, .(ID, ALT)], by = "ID", all.x = TRUE)
chk(!anyNA(mm$ALT) && all(mm$A1 == mm$ALT), "main GWAS SNP이 variant 목록에 있고 A1 = ALT (앞 10만 개; README에 적은 내용)")
rr <- fread(cmd = paste("zcat", shQuote(file.path(out, "02_GWAS_REGENIE/bin_all.regenie.gz")), "| head -100001"))
rr <- merge(rr[, .(ID, ALLELE1)], pv[, .(ID, ALT)], by = "ID")
chk(nrow(rr) == 1e5 && all(rr$ALLELE1 == rr$ALT), "REGENIE ALLELE1 = ALT (앞 10만 개)")

put(f("14-12_zenodo_README.md"), "README.md")
setcolorder(M, c("model", "group", "description", "n", "software", "in_fig2a", "file"))
wcsv(M, "model_key.csv")
cat("■ 전체 점검\n")
fl <- list.files(out, recursive = TRUE)
txt <- fl[!grepl("\\.(gz)$", fl)]
chk(!any(vapply(txt, function(x) any(grepl(hangul, readLines(file.path(out, x), warn = FALSE), perl = TRUE)), logical(1))), "평문 파일에 한글 없음")
gz_small <- fl[grepl("\\.csv\\.gz$", fl)]
chk(!any(vapply(gz_small, function(x) any(grepl(hangul, unlist(fread(cmd = paste("zcat", shQuote(file.path(out, x))), colClasses = "character")), perl = TRUE)), logical(1))), "csv.gz에 한글 없음")
chk(!any(vapply(txt, function(x) any(grepl("/data/|/home/", readLines(file.path(out, x), warn = FALSE))), logical(1))), "서버 경로 없음")
chk(!any(grepl("keep|pheno|covar|psam|fam|grm|log$", fl, ignore.case = TRUE)), "개체 단위 · 로그 파일 없음")
md5 <- tools::md5sum(file.path(out, fl))
writeLines(sprintf("%s  %s", md5, fl), file.path(out, "MD5SUMS.txt"))
sz <- file.size(file.path(out, fl))
cat(sprintf("\n파일 %d개, 합계 %.2f GB → %s\n", length(fl), sum(sz) / 1e9, out))
print(data.table(dir = dirname(fl), size_MB = sz / 1e6)[, .(files = .N, MB = round(sum(size_MB))), by = dir])
