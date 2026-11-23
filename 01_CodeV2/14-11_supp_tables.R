# 14-11_supp_tables.R - Assembles Supplementary Tables S1-S10 into a Word file (pandoc). See README.md.
#   S3  06-6_h2.csv · 06-8_h2.csv · 06-8_gxe.csv
#   S6  05-4_run_summary.csv · 13-10_tests.csv · 14-5_age_pathway_scan.csv
#   S7  13-12_topk.csv · 13-12_tests.csv (U100) · 14-4_brain_expression.csv
#   S9  14-3_power_table.csv · 14-3_power_design_grid.csv
#   S10 14-6_cfa19_breed_check.csv · 14-9_cfa38_breed_check.csv · 14-6_cfa19_breeds.csv · 14-9_cfa38_breeds.csv · 14-9_ld_support.csv

suppressMessages(library(data.table))
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
f <- function(x) { p <- file.path(d0, x); if (!file.exists(p)) stop("파일 없음: ", p); p }
out_md <- file.path(d0, "14-11_supp_tables.md"); out_docx <- file.path(d0, "14-11_supp_tables.docx")
if (Sys.which("pandoc") == "") stop("pandoc 없음")

sup <- c("-" = "⁻", "0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴", "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹")
to_sup <- function(s) paste0(sup[strsplit(s, "")[[1]]], collapse = "")
fp1 <- function(p) {
  if (is.na(p)) return("–")
  if (p < 0.001) { e <- floor(log10(p)); m <- p / 10^e; if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
    return(sprintf("%.1f×10%s", m, to_sup(as.character(e)))) }
  formatC(p, digits = 2, format = "fg", flag = "#") }
fp <- function(p) vapply(p, fp1, "")
fn <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
fb <- function(x, d = 3) ifelse(is.na(x), "–", formatC(x, digits = d, format = "f"))
fbs <- function(b, s, d = 3) ifelse(is.na(b), "–", sprintf("%s (%s)", fb(b, d), fb(s, d)))
it <- function(g) ifelse(is.na(g) | g == "", "", paste0("*", g, "*"))
breed_nm <- function(x) vapply(strsplit(gsub("_", " ", x), " "), function(w) paste0(toupper(substring(w, 1, 1)), substring(w, 2), collapse = " "), "")
genes_fmt <- function(s) vapply(strsplit(s, ", "), function(v) paste(vapply(v, function(g) {   # "CFI(in), GAR1(18kb)" → "*CFI* (within), *GAR1* (18 kb)"
  m <- regmatches(g, regexec("^(.+)\\((.+)\\)$", g))[[1]]
  if (length(m) != 3) return(it(g))
  paste0(it(m[2]), " (", if (m[3] == "in") "within" else sub("kb$", " kb", m[3]), ")") }, ""), collapse = ", "), "")
mname <- c(main = "All dogs", young8 = "< 8 y", age8 = "≥ 8 y", age10 = "≥ 10 y", agefloor6 = "≥ 6 y", agefloor7 = "≥ 7 y", agefloor9 = "≥ 9 y",
           agefloor11 = "≥ 11 y", agefloor12 = "≥ 12 y", sex_m = "Males", sex_f = "Females", purelike = "Purebred-like", mixlike = "Mixed-like",
           wt_t1 = "Light", wt_t2 = "Middle", wt_t3 = "Heavy", int_age = "SNP × age", bin_age8 = "≥ 8 y, binary", bin_all = "All dogs, binary", int_wt = "SNP × body weight")
snp_fmt <- function(id) sub("^(\\d+):(\\d+):(\\w+):(\\w+)$", "\\1:\\2 (\\3/\\4)", id)
md_tab <- function(D, align = NULL) {                  # data.table → Markdown pipe table
  D <- as.data.frame(D); D[] <- lapply(D, function(x) { x <- trimws(as.character(x)); x[is.na(x)] <- "–"; x <- gsub("\\|", "/", x)
    x <- gsub(">=", "≥", x, fixed = TRUE); x <- gsub("(\\d)-(\\d+ y)", "\\1–\\2", x)
    gsub("(^|[ (])-(\\d)", "\\1\u2212\\2", x) })
  if (is.null(align)) align <- c("l", rep("r", ncol(D) - 1))
  sep <- vapply(align, function(a) switch(a, l = ":---", r = "---:", c = ":---:"), "")
  c(paste0("| ", paste(names(D), collapse = " | "), " |"), paste0("|", paste(sep, collapse = "|"), "|"),
    apply(D, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "") }
md <- character(0)
add <- function(...) md <<- c(md, ...)
tab_no <- c(S1 = 1, S2 = 2, S3 = 4, S4 = 5, S5 = 3, S6 = 6, S7 = 8, S8 = 9, S9 = 10, S10 = 7)
chunk <- integer(0)
title <- function(id, txt) { new <- paste0("Table S", tab_no[[sub("^Table ", "", id)]]); chunk[new] <<- length(md) + 1L
  add(paste0("## ", new), "", paste0("**", new, "** ", txt), "") }
note <- function(txt) add(paste0("*", txt, "*"), "")
chk <- function(cond, msg) if (!isTRUE(cond)) stop("점검 실패: ", msg)

add("# Supplementary tables", "")

ids <- fread(f("04-2_pheno.txt"), header = FALSE)$V2
A <- fread(f("00_all_total.csv"), select = c("IID", "Breed_Status", "Breed_Size_Class_at_HLES", "top_pop", "top_pct", "coi", "genetic_type"))[IID %in% ids]
chk(nrow(A) == 7115 && A[genetic_type == "purelike", .N] == 3690, "S1 분석 개 7,115 · 순종형 3,690")
chk(A[genetic_type == "purelike", min(coi)] >= 0.1507 && A[genetic_type == "mixlike", max(coi)] < 0.1507, "순종형 = COI ≥ 0.1507 (각주 숫자)")
title("Table S1", "Breed composition of the 7,115 dogs analysed.")
add("**(a) Genetic classification by owner-reported breed status**", "")
ct <- dcast(A[, .N, by = .(genetic_type, Breed_Status)], genetic_type ~ Breed_Status, value.var = "N", fill = 0)
ct <- ct[match(c("purelike", "mixlike"), genetic_type)]
S1a <- data.table(`Genetic classification` = c("Purebred-like", "Mixed-like", "Total"),
                  `Owner-reported purebred` = fn(c(ct$Purebred, sum(ct$Purebred))),
                  `Owner-reported mixed breed` = fn(c(ct$`Mixed breed`, sum(ct$`Mixed breed`))),
                  Total = fn(c(ct$Purebred + ct$`Mixed breed`, nrow(A))))
add(md_tab(S1a))
add("**(b) Purebred-like dogs by breed (largest genetic ancestry component)**", "")
size_word <- function(s) sub("^(Toy and Small|Medium|Standard|Large|Giant).*$", "\\1", s)
P <- A[genetic_type == "purelike"]
B <- P[, .(n = .N, med = median(top_pct), size = { z <- table(size_word(Breed_Size_Class_at_HLES)); names(z)[which.max(z)] }), by = top_pop][order(-n, top_pop)]
big <- B[n >= 20]; small <- B[n < 20]
S1b <- rbind(big[, .(Breed = breed_nm(top_pop), Dogs = fn(n), `Median largest ancestry proportion` = fb(med, 2), `Most common owner-reported size class` = size)],
             data.table(Breed = sprintf("Other breeds (%d breeds, < 20 dogs each)", nrow(small)), Dogs = fn(sum(small$n)),
                        `Median largest ancestry proportion` = fb(P[top_pop %in% small$top_pop, median(top_pct)], 2), `Most common owner-reported size class` = "–"),
             data.table(Breed = "Total purebred-like", Dogs = fn(nrow(P)), `Median largest ancestry proportion` = fb(median(P$top_pct), 2), `Most common owner-reported size class` = "–"))
add(md_tab(S1b, c("l", "r", "r", "l")))
note("Purebred-like: genomic inbreeding coefficient (COI) ≥ 0.1507, the boundary between the low- and high-inbreeding modes of the COI distribution (two-component mixture model fitted to log COI). Breed = reference-panel breed with the largest genetic ancestry proportion. Size class from the owner-reported Health and Life Experience Survey (HLES).")

logn <- function(file, pat) { x <- readLines(f(file)); m <- regmatches(x, regexpr(pat, x)); as.numeric(sub("^(\\d+).*$", "\\1", m[1])) }
n_raw <- logn("01-2_geno0.05.log", "^\\d+ variants loaded"); n_geno <- logn("01-2_geno0.05.log", "^\\d+ variants remaining after main")
n_gm <- logn("01-4_groupmiss.log", "^\\d+ variants remaining after main"); n_mind_rm <- logn("01-4_mind0.05.log", "^\\d+ samples removed")
n_dog0 <- logn("01-2_geno0.05.log", "^\\d+ samples"); n_dog1 <- logn("01-4_mind0.05.log", "^\\d+ samples \\(\\d+ females, \\d+ males; \\d+ founders\\) remaining")
n_maf <- logn("02-1_maf0.01.log", "^\\d+ variants remaining after main"); n_hwe <- logn("03-3_maf0.01_hwe.log", "^\\d+ variants remaining after main")
chk(n_raw == 22397820 && n_hwe == 4842736 && n_dog1 == 7115, "S2 처음 · 마지막 수")
G1 <- fread(f("01-1_geno_threshold_summary.csv")); G3 <- fread(f("01-3_group_miss_criteria.csv")); G4 <- fread(f("01-4_mind_threshold_summary.csv"))
G5 <- fread(f("02-2_maf_set_summary.csv")); G6 <- fread(f("03-2_hwe_criteria.csv"))
n_grp <- nrow(fread(f("03-2_hwe_by_group.csv")))
alt <- function(lab, n) paste(sprintf("%s: %s", lab, fn(n)), collapse = "; ")
title("Table S2", "Genotype quality control.")
S2 <- data.table(
  Step = c("Imputed genotypes (DAP 2023 release)", "Variant missingness", "Ancestry-group differential missingness", "Sample missingness",
           "Minor allele frequency (variant missingness re-applied after sample removal)", "Hardy–Weinberg equilibrium"),
  `Criterion used` = c("–", "missing in < 5% of dogs", "removed if, within one ancestry group, missingness > 10% and exceeds that of other dogs (binomial p < 1×10⁻⁶)",
                       "missing in < 5% of SNPs", "MAF ≥ 0.01; missing in < 5% of dogs", sprintf("removed if p < 1×10⁻⁶ in ≥ 2 of %d ancestry groups", n_grp)),
  `SNPs removed` = c("–", fn(n_raw - n_geno), fn(n_geno - n_gm), "–", fn(n_gm - n_maf), fn(n_maf - n_hwe)),
  `SNPs remaining` = fn(c(n_raw, n_geno, n_gm, n_gm, n_maf, n_hwe)),
  `Dogs removed` = c("–", "–", "–", fn(n_mind_rm), "–", "–"),
  `Dogs remaining` = fn(c(n_dog0, n_dog0, n_dog0, n_dog1, n_dog1, n_dog1)),
  `Alternative thresholds evaluated (remaining)` = c("–",
    alt(paste0(G1[geno != 0.05, geno * 100], "%"), G1[geno != 0.05, n_kept]),
    alt(c("missingness > 10% only", "binomial test only"), G3[1:2, n_keep]),
    alt(paste0(G4[mind != 0.05, mind * 100], "%"), G4[mind != 0.05, n_keep]),
    alt(paste0("MAF ≥ ", G5[maf != 0.01, maf]), G5[maf != 0.01, n_keep]),
    alt(c("pooled sample", "≥ 1 group"), G6[1:2, n_keep])))
add(md_tab(S2, c("l", "l", "r", "r", "r", "r", "l")))
note("Alternative thresholds: SNPs remaining (dogs remaining for sample missingness) had that threshold been used at that step. Applying the HWE test to the pooled sample would remove 99% of SNPs because of population structure.")

H6 <- fread(f("06-6_h2.csv")); H8 <- fread(f("06-8_h2.csv"))[source == "06-8"]; GX <- fread(f("06-8_gxe.csv"))
lab6 <- c(h2_main = "All dogs | Most recent score", h2_mean = "All dogs | Mean of repeated scores", h2_young8 = "< 8 y | Most recent score",
          h2_age8 = "≥ 8 y | Most recent score", h2_sex_m = "Males | Most recent score", h2_sex_f = "Females | Most recent score")
lab8 <- c(h2_last8 = "≥ 8 y | Most recent score (+ number of scores at ≥ 8 y as covariate)", h2_mean8 = "≥ 8 y | Mean of scores at ≥ 8 y",
          h2_last8_m2 = "≥ 8 y, ≥ 2 scores at ≥ 8 y | Most recent score", h2_mean8_m2 = "≥ 8 y, ≥ 2 scores at ≥ 8 y | Mean of scores at ≥ 8 y")
HH <- rbind(H6[match(names(lab6), model), .(model, lab = lab6[model], n, h2, se, p = p_LRT)], H8[match(names(lab8), model), .(model, lab = lab8[model], n, h2, se, p = p_LRT)])
chk(!anyNA(HH$h2) && nrow(HH) == 10, "S3 유전율 줄 10개")
title("Table S3", "SNP heritability and genotype-by-age variance of CSLB scores.")
add("**(a) SNP heritability (GREML)**", "")
S3a <- HH[, .(Sample = sub(" \\|.*$", "", lab), Phenotype = sub("^.*\\| ", "", lab), n = fn(n), `h²SNP (SE)` = fbs(h2, se),
              `95% CI` = sprintf("%s–%s", fb(pmax(0, h2 - 1.96 * se)), fb(h2 + 1.96 * se)), p = fp(p))]
add(md_tab(S3a, c("l", "l", "r", "r", "r", "r")))
add("**(b) Genotype-by-age model (all dogs; age groups < 8 and ≥ 8 years)**", "")
S3b <- GX[, .(Phenotype = fifelse(model == "gxe_main", "Rank-based inverse-normal transformation across all dogs", "Rank-based inverse-normal transformation within age group"),
              n = fn(n), `Shared genetic variance, V(G)/Vp (SE)` = fbs(h2_G, se_G), `G × age variance, V(G×age)/Vp (SE)` = fbs(h2_GxE, se_GxE), `p (G × age)` = fp(p_GxE))]
add(md_tab(S3b, c("l", "r", "r", "r", "r")))
note("p: one-sided likelihood-ratio test. 95% CI = estimate ± 1.96 SE, lower limit truncated at zero.")

C4 <- fread(f("13-1_compare.csv"))[model == "main" & as.numeric(thr) == 1e-5][order(p)]
chk(nrow(C4) == 13 && all(C4$rob_stable), "S4 주 분석 lead 13개 · 모두 안정")
title("Table S4", "Suggestive lead loci and their replication between complementary strata.")
add("**(a) Lead SNPs with p < 1×10⁻⁵ in the full sample**", "")
C4[, `:=`(a1 = sub("^\\d+:\\d+:(\\w+):\\w+$", "\\1", SNP), a2 = sub("^.*:(\\w+)$", "\\1", SNP))][, oa := fifelse(ea == a1, a2, a1)]
chk(all(C4$ea == C4$a1 | C4$ea == C4$a2), "S4 효과 대립유전자가 SNP ID 안에 있음")
S4a <- C4[, .(`Lead SNP, CFA:position (effect/other allele)` = sprintf("%s:%s (%s/%s)", sub(":.*$", "", SNP), fn(bp), ea, oa),
              MAF = fb(maf, 3), `β` = fb(b, 3), p = fp(p), `Nearest genes (distance)` = genes_fmt(genes),
              `Minimum β ratio, 14 alternative models` = fb(rob_min_ratio, 2),
              `p het, sex` = fp(het_sex_p), `p het, age` = fp(het_age_p), `p het, purebred-like vs. mixed-like` = fp(het_ancestry_p), `p het, body weight` = fp(het_weight_p),
              `p, SNP × age` = fp(int_age_p))]
add(md_tab(S4a, c("l", "r", "r", "r", "l", "r", "r", "r", "r", "r", "r")))
note("β: effect of the effect allele on rank-based inverse-normal transformed CSLB score (GCTA MLMA-LOCO). Minimum β ratio: smallest ratio of β in 14 models with altered covariates (number of PCs, added or dropped covariates) to β in the primary model; all 13 kept the same direction. p het: test for a difference in β between non-overlapping strata (age < 8 vs. ≥ 8 years; males vs. females; purebred-like vs. mixed-like by genomic inbreeding, Table S1; body-weight tertiles). p, SNP × age: SNP × age interaction with age as a continuous variable (REGENIE).")
add("**(b) Lead SNPs nominated in one stratum (p < 1×10⁻⁴) that replicated in the complementary stratum (same direction, p < 0.05)**", "")
P4 <- fread(f("13-1_priority.csv"))[tier == 1 & as.numeric(thr) == 1e-4][order(chr, bp, p)]
rep_fmt <- function(s) vapply(strsplit(s, "; "), function(v) paste(vapply(v, function(x) {
  m <- regmatches(x, regexec("^(\\S+) b=(\\S+) p=(\\S+)$", x))[[1]]; sprintf("%s: β %s, p %s", mname[m[2]], m[3], fp(as.numeric(m[4]))) }, ""), collapse = "; "), "")
rep_ok <- vapply(seq_len(nrow(P4)), function(i) any(vapply(strsplit(P4$rep_partners[i], "; ")[[1]], function(x) {
  m <- regmatches(x, regexec("^(\\S+) b=(\\S+) p=(\\S+)$", x))[[1]]; sign(as.numeric(m[3])) == sign(P4$b[i]) && as.numeric(m[4]) < 0.05 }, TRUE)), TRUE)
chk(nrow(P4) == 24 && all(rep_ok), "S4b 24줄 모두 지명 층과 같은 방향 · p < 0.05인 짝이 있음")
S4b <- P4[, .(`Nominating stratum` = mname[model], `Lead SNP, CFA:position` = sprintf("%s:%s", chr, fn(bp)), MAF = fb(maf, 3), `β` = fb(b, 3), p = fp(p),
              `Complementary stratum` = rep_fmt(rep_partners), `Nearest genes (distance)` = genes_fmt(genes))]
add(md_tab(S4b, c("l", "l", "r", "r", "r", "l", "l")))
note("MAF, β and p: in the nominating stratum; β in the complementary stratum is for the same effect allele. Strata: age (< 8 y vs. ≥ 8 y; age-floor samples ≥ 9 to ≥ 12 y were compared with < 8 y), sex, genomic inbreeding (purebred-like vs. mixed-like; Table S1), and body-weight tertiles (light, middle, heavy; a nomination in one tertile was compared with each of the other two). The same locus can appear more than once when nominated in several strata; age-floor samples overlap each other. The lead loci of the full sample near *CFI*, *LRMDA*, *DTNBP1* and *PRDM16* replicated in two different stratifications each.")

F2 <- fread(f("14-2_fig2_data.csv"))[panel == "a"]; SI <- fread(f("06-4_set_info.csv"))
pre <- SI$set; chk(length(pre) == 14 && all(pre %in% F2$set), "S5 사전 세트 14개가 Fig. 2a 자료에 있음")
F2 <- F2[set %in% pre]; F2[, fam := SI$family[match(set, SI$set)]]
fam_nm <- c(circadian = "Circadian", sleep_gwas = "Sleep GWAS", AD_cognition_gwas = "AD & cognition GWAS", brain_aging_expr = "Brain ageing", splicing = "Splicing")
src_nm <- function(s, set) fifelse(grepl("^KEGG", s), "KEGG pathway (2022)", fifelse(s == "GO_BP", paste0("Gene Ontology biological process (GO:", sub("^GO_", "", set), ")"),
                         fifelse(s == "GWAS Catalog", "GWAS Catalog genes", fifelse(s == "curated", "Curated list", "Palmer et al. (2021), Tables S4 and S8"))))
ordS <- unique(F2[, .(set, name, fam)])[order(match(fam, names(fam_nm)), name)]
title("Table S5", "Prespecified gene sets and their associations.")
add("**(a) Set definitions and association in all dogs**", "")
mainS <- F2[run == "main"][match(ordS$set, set)]
S5a <- mainS[, .(Family = fam_nm[fam], `Gene set` = name, Source = src_nm(SI$source[match(set, SI$set)], set), `Genes tested` = NGENES, `β (SE)` = fbs(BETA, SE), p = fp(P))]
add(md_tab(S5a, c("l", "l", "l", "r", "r", "r")))
add("**(b) p in each model shown in Fig. 2a**", "")
runs_b <- c("main", "young8", "age8", "sex_m", "sex_f", "purelike", "mixlike", "wt_t1", "wt_t2", "wt_t3", "int_age", "bin_age8")
W5 <- dcast(F2, set ~ run, value.var = "P")[match(ordS$set, set)]
S5b <- data.table(`Gene set` = ordS$name); for (r in runs_b) S5b[[mname[r]]] <- fp(W5[[r]])
add(md_tab(S5b))
note("MAGMA competitive gene-set test (genes ± 10 kb): whether genes in the set are more strongly associated with CSLB scores than all other genes. Body-weight tertiles: light ≤ 15.9 kg, middle 15.9–29.0 kg, heavy > 29.0 kg. ≥ 8 y, binary: CSLB ≥ 50 vs. < 50 in dogs aged ≥ 8 years. Sources of the GWAS sets: chronotype, Jones et al. (2019); sleep duration, Dashti et al. (2019); insomnia, Jansen et al. (2019); Alzheimer's disease, Bellenguez et al. (2022); general cognitive ability, Davies et al. (2018).")

R6 <- fread(f("05-4_run_summary.csv"))
ord6 <- c("young8", "agefloor6", "agefloor7", "age8", "agefloor9", "age10", "agefloor11", "agefloor12", "bin_age8", "int_age")
R6 <- R6[match(ord6, id)]; chk(!anyNA(R6$n), "S6a 모형 10개")
title("Table S6", "Analyses restricted to older dogs.")
add("**(a) Genome-wide association by minimum age**", "")
S6a <- R6[, .(Sample = mname[id], n = fn(n), `λ` = fb(lambda, 3), `SNPs p < 5×10⁻⁸` = n_p5e8, `SNPs p < 1×10⁻⁵` = n_p1e5, `SNPs p < 1×10⁻⁴` = n_p1e4,
              `Correlation of z with full sample (expected)` = ifelse(is.na(z_cor_main), "–", sprintf("%s (%s)", fb(z_cor_main, 2), ifelse(is.na(z_cor_expect), "–", fb(z_cor_expect, 2)))),
              `Strongest locus not among full-sample loci` = ifelse(best_new == "", "–", sub("^chr(\\d+):([0-9.]+)Mb p=(.+)$", "CFA\\1, \\2 Mb", best_new)),
              `p` = ifelse(best_new == "", "–", fp(as.numeric(sub("^.*p=", "", best_new)))))]
add(md_tab(S6a, c("l", "r", "r", "r", "r", "r", "r", "l", "r")))
note("< 8 y: complement of ≥ 8 y. ≥ 8 y, binary: CSLB ≥ 50 vs. < 50 (REGENIE with Firth correction). SNP × age: REGENIE interaction test in all dogs. Expected correlation: √(n/7,115), the correlation expected if effects were the same as in the full sample.")
T6 <- fread(f("13-10_tests.csv")); old <- c("age8", "age10", "agefloor6", "agefloor7", "agefloor9", "agefloor11", "agefloor12", "bin_age8")
T6o <- T6[model %in% old]
add(sprintf("**(b) Gene sets nominated by loci (p < 1×10⁻⁴) in older-dog models: %s set–model tests (%s primary and %s secondary candidates); sets with p < 0.05 after removing genes within 100 kb of the nominating loci**",
            fn(nrow(T6o)), fn(T6o[tier == "주 후보", .N]), fn(T6o[tier == "보조 후보", .N])), "")
pass6 <- T6o[main_pass == TRUE][order(match(model, old), p)]
S6b <- pass6[, .(Model = mname[model], Candidate = fifelse(tier == "주 후보", "Primary", "Secondary"), `Gene set` = sprintf("%s (%s)", set_name, sub("_", ":", set)),
                 `Genes after removal` = n_left, `β` = fb(beta, 3), p = fp(p),
                 `< 8 y: p` = ifelse(is.na(pair_P), "–", fp(pair_P)), `p, difference from < 8 y` = ifelse(is.na(diff_p), "–", fp(diff_p)),
                 `Older-dog-specific` = fifelse(is.na(pair_ok), "not tested (no paired sample)", fifelse(pair_ok, "yes", "no")))]
add(md_tab(S6b, c("l", "l", "l", "r", "r", "r", "r", "r", "l")))
note("Older-dog-specific: β in the older-dog model larger than in dogs aged < 8 years with p < 0.05 for the difference. Age-floor samples ≥ 6 and ≥ 7 y overlap with < 8 y and were not compared.")
S6c_all <- fread(f("14-5_age_pathway_scan.csv"))
fwrite(S6c_all[, .(set, name, n_genes, P_young8, P_age8, BETA_young8, BETA_age8, diff_z, diff_p)][order(diff_p)], file.path(d0, "14-11_TableS6c_full.csv"))
add(sprintf("**(c) Genome-wide comparison of %s KEGG and Gene Ontology sets between dogs aged < 8 and ≥ 8 years: 20 sets with the smallest p for the difference**", fn(nrow(S6c_all))), "")
S6c <- S6c_all[order(diff_p)][1:20][, .(`Gene set` = name, ID = sub("_", ":", set), Genes = n_genes, `p, < 8 y` = fp(P_young8), `p, ≥ 8 y` = fp(P_age8),
                                       `β, < 8 y` = fb(BETA_young8, 3), `β, ≥ 8 y` = fb(BETA_age8, 3), `z, difference` = fb(diff_z, 2), `p, difference` = fp(diff_p),
                                       `Larger in` = fifelse(diff_z > 0, "≥ 8 y", "< 8 y"))]
add(md_tab(S6c, c("l", "l", "r", "r", "r", "r", "r", "r", "r", "l")))
note(sprintf("Full results for all %s sets are provided at [Zenodo data DOI]. z, difference = (β ≥ 8 y − β < 8 y) / √(SE² ≥ 8 y + SE² < 8 y).", fn(nrow(S6c_all))))

K7 <- fread(f("13-12_topk.csv")); U7 <- unique(fread(f("13-12_tests.csv"))[ver == "U100" & model %in% c("sex_m", "sex_f", "main"), .(model, NGENES, BETA, P)])
title("Table S7", "Genes contributing to the association of the cellular response to calcium ion set (GO:0071277) and brain-elevated expression.")
add("**(a) Cumulative removal of genes in order of decreasing gene-level z-score (Fig. S3b)**", "")
K7a <- rbind(U7[, .(model, k = 0L, removed = "– (all 76 genes)", z = NA_real_, NGENES, BETA, P)], K7[, .(model, k, removed = it(removed), z, NGENES, BETA, P)])
K7a <- K7a[order(match(model, c("sex_m", "sex_f", "main")), k)]
S7a <- K7a[, .(Sample = mname[model], `Genes removed` = k, `Gene removed at this step` = removed, `Gene z` = fb(z, 2), `Genes left` = NGENES, `β` = fb(BETA, 3), p = fp(P))]
add(md_tab(S7a, c("l", "r", "l", "r", "r", "r", "r")))
add("**(b) Brain-elevated expression (Human Protein Atlas)**", "")
E7 <- fread(f("14-4_brain_expression.csv"))
E7 <- E7[c(4, 1, 2, 3)]
lab7 <- c("Cellular response to calcium ion (GO:0071277)" = "Cellular response to calcium ion set (GO:0071277)", "Locus genes, all models (1e-4; main)" = "Genes near loci (p < 1×10⁻⁴), all models",
          "Locus genes, all-dogs model only (1e-4)" = "Genes near loci (p < 1×10⁻⁴), all-dogs model", "Locus genes, 1e-5 loci only" = "Genes near loci (p < 1×10⁻⁵)")
chk(all(E7$comparison %in% names(lab7)), "S7b 비교 이름")
S7b <- E7[, .(`Gene group` = lab7[comparison], Genes = fn(n_genes), `Brain-elevated (%)` = fb(pct_brain, 1), `Brain-elevated among other genes (%)` = fb(pct_brain_rest, 1),
              `Median length, kb (other genes)` = sprintf("%s (%s)", fb(len_median_kb, 1), fb(len_median_rest_kb, 1)), OR = fb(OR, 2), p = fp(p),
              `Length-adjusted OR (95% CI)` = sprintf("%s (%s–%s)", fb(OR_len, 2), fb(OR_len_lo, 2), fb(OR_len_hi, 2)), `p, length-adjusted` = fp(p_len))]
add(md_tab(S7b, c("l", "r", "r", "r", "r", "r", "r", "r", "r")))
note("Gene z: MAGMA gene-level association z-score with CSLB scores. β and p: MAGMA competitive test of the remaining genes. OR: odds ratio for brain-elevated expression compared with all other genes; length-adjusted OR from logistic regression including log gene length.")

r07 <- fread(f("07-2_results.csv"))[model == "P_age8_FROH"]; r08 <- fread(f("08-2_results.csv"))[model == "P_age8_anc"]
r09 <- fread(f("09-2_results.csv")); r10 <- fread(f("10-4_tests.csv")); r11 <- fread(f("11-2_h2.csv")); r12 <- fread(f("12-2_greml.csv"))[run == "P_all_dense"]
pg <- r09[model == "P_all_PGSxAge"]; y8 <- r10[model == "B2_age8_Y_male"]; mt8 <- r10[model == "B1_age8_mtCN"]
chk(nrow(r07) == 1 && nrow(r08) == 1 && nrow(pg) == 2 && nrow(y8) == 1 && nrow(mt8) == 1 && nrow(r12) == 1, "S8 주 검정 줄")
h11 <- function(s) r11[sample == s]
title("Table S8", "Tests of alternative explanations for CSLB score variation.")
S8 <- rbindlist(list(
  data.table(h = "Genomic inbreeding", t = "F_ROH (genomic inbreeding coefficient), dogs ≥ 8 y", n = fn(r07$n), e = fbs(r07$beta, r07$se), p = fp(r07$p), r = "β per unit F_ROH"),
  data.table(h = "Breed ancestry", t = "Variance explained by breed-ancestry relationship matrix, dogs ≥ 8 y", n = fn(r08$n), e = fbs(r08$prop, r08$se), p = fp(r08$p), r = "proportion of phenotypic variance; LRT"),
  data.table(h = "Body size (genetic)", t = "Body-size polygenic score × age, all dogs", n = fn(pg[term == "PGS_x_age", n]), e = fbs(pg[term == "PGS_x_age", beta], pg[term == "PGS_x_age", se], 4), p = fp(pg[term == "PGS_x_age", p]), r = "β per SD of score per year"),
  data.table(h = "", t = "Body-size polygenic score, main effect (same model; not prespecified)", n = fn(pg[term == "PGS_z", n]), e = fbs(pg[term == "PGS_z", beta], pg[term == "PGS_z", se], 3), p = fp(pg[term == "PGS_z", p]), r = "β per SD of score"),
  data.table(h = "X chromosome", t = "X-chromosome genetic variance", n = "–", e = "–", p = "–", r = "not testable: imputed X genotypes contained almost no heterozygous females"),
  data.table(h = "Y chromosome", t = "Y-chromosome read dosage, males ≥ 8 y", n = fn(y8$n), e = fbs(y8$beta, y8$se, 3), p = fp(y8$p), r = "β per SD"),
  data.table(h = "Mitochondrial DNA", t = "mtDNA copy number, dogs ≥ 8 y", n = fn(mt8$n), e = fbs(mt8$beta, mt8$se, 3), p = fp(mt8$p), r = "β per SD"),
  data.table(h = "Sensory or neurological disease", t = "SNP heritability after excluding dogs diagnosed with a sensory or neurological condition before the CSLB survey, all dogs (before exclusion)", n = sprintf("%s (%s)", fn(h11("all_primary")$n), fn(h11("all_unpurified")$n)),
             e = sprintf("%s; before exclusion %s", fbs(h11("all_primary")$h2, h11("all_primary")$se), fb(h11("all_unpurified")$h2)), p = fp(h11("all_primary")$p_LRT), r = "h²SNP (SE); one-sided LRT"),
  data.table(h = "", t = "Same, dogs ≥ 8 y (before exclusion)", n = sprintf("%s (%s)", fn(h11("age8_primary")$n), fn(h11("age8_unpurified")$n)),
             e = sprintf("%s; before exclusion %s", fbs(h11("age8_primary")$h2, h11("age8_primary")$se), fb(h11("age8_unpurified")$h2)), p = fp(h11("age8_primary")$p_LRT), r = "h²SNP (SE); one-sided LRT"),
  data.table(h = "Dominance", t = "Dominance genetic variance, all dogs", n = fn(r12$n), e = fbs(r12$d2, r12$se_D), p = fp(r12$p_one), r = "δ² (SE); one-sided LRT")), fill = TRUE)
setnames(S8, c("Hypothesis", "Test", "n", "Estimate (SE)", "p", "Scale / note"))
add(md_tab(S8, c("l", "l", "r", "r", "r", "l")))
note("Each hypothesis had one prespecified primary test, listed first. GREML models included the genomic relationship matrix and the covariates of the primary GWAS.")

P9 <- fread(f("14-3_power_table.csv")); D9 <- fread(f("14-3_power_design_grid.csv"))
title("Table S9", "Sample sizes required to detect SNP heritability or a single-SNP effect with 80% power under alternative survey designs.")
add("**(a) Current design and dogs aged ≥ 12 years**", "")
S9a <- P9[, .(Target = fifelse(grepl("heritability", target), "SNP heritability > 0 (one-sided p < 0.05)", "Single SNP (p < 5×10⁻⁸)"), Design = design,
              Reliability = fb(reliability, 2), `True effect` = sprintf("%s%%", formatC(100 * effect, format = "fg", digits = 2)), `Observed effect` = sprintf("%s%%", fb(100 * effect_obs, 2)),
              `Dogs required` = fn(N_req), Range = ifelse(is.na(N_req_lo), "–", sprintf("%s–%s", fn(N_req_lo), fn(N_req_hi))))]
add(md_tab(S9a, c("l", "l", "r", "r", "r", "r", "r")))
add("**(b) Designs by age at enrolment and number of annual surveys**", "")
S9b <- D9[, .(`Age at enrolment` = start_age, `Annual surveys` = k, `Single-survey reliability` = fb(rel1, 3), `Reliability of mean` = fb(reliability, 2),
              `Dogs, h²SNP 5%` = fn(N_h2_5), `Dogs, h²SNP 10%` = fn(N_h2_10), `Dogs, SNP 1%` = fn(N_snp_1), `Dogs, SNP 2%` = fn(N_snp_2), `Surveys, h²SNP 5%` = fn(surveys_total_h2_5))]
add(md_tab(S9b, c("l", "r", "r", "r", "r", "r", "r", "r", "r")))
note("Observed effect = true effect × reliability. Reliability of the mean of k surveys from the Spearman–Brown formula. SNP heritability: GREML standard error scaled with sample size from the five models fitted in this study; range across these models. Single SNP: variance explained, two-sided α = 5×10⁻⁸.")

# Table S10  CFA19 · CFA38
c19 <- fread(f("14-6_cfa19_breed_check.csv")); c38 <- fread(f("14-9_cfa38_breed_check.csv")); setnames(c19, "n_T", "n_minor")
mlab <- function(m) { m <- sub("^M0 all dogs$", "All dogs", m); m <- sub("^S1 all dogs \\+ breed fixed effects$", "All dogs + breed fixed effects", m)
  m <- sub("^S2 mixed-like only$", "Dogs with largest breed ancestry < 0.85 only", m); m <- sub("^S3 purebred-like only \\+ breed fixed effects$", "Dogs with largest breed ancestry ≥ 0.85 only + breed fixed effects", m)
  m <- ifelse(grepl("^S4 without ", m), paste0("Excluding ", breed_nm(sub("^S4 without ", "", m))), m)
  m <- sub("^Mixed-like, ", "Largest breed ancestry < 0.85, ", m); sub("^(All dogs|Largest breed ancestry < 0.85), (Light|Middle|Heavy)$", "\\1, \\2 tertile", m) }
mk10 <- function(D, locus) D[, .(Locus = locus, Model = mlab(model), n = fn(n), `Minor alleles` = fn(n_minor),
                                 `β SNP (SE)` = fbs(b_g, se_g), `β SNP × log₂ weight (SE)` = fbs(b_int, se_int), `p, interaction` = fp(p_int), `p, joint 2-df` = fp(p_2df), `p, SNP (tertile)` = fp(p_g))]
title("Table S10", "Body-weight-dependent loci on CFA19 and CFA38.")
add("**(a) Interaction models**", "")
add(md_tab(rbind(mk10(c19, "CFA19 19:41392616 (near *NXPH2*)"), mk10(c38, "CFA38 38:23706225 (near *DISP1*)")), c("l", "l", rep("r", 7))))
note("Ordinary least squares with HC1 standard errors on the phenotype adjusted for the leave-one-chromosome-out polygenic prediction. Breed in this table: reference-panel breed with the largest genetic ancestry proportion; breed fixed effects for breeds with ≥ 20 dogs whose largest ancestry proportion is ≥ 0.85. Excluded breeds: the five breeds with the most minor alleles. Tertile rows: SNP effect within body-weight tertiles. The joint 2-df p in the main text (7.7×10⁻⁹ and 3.1×10⁻⁸) are from REGENIE.")
b19 <- fread(f("14-6_cfa19_breeds.csv")); b38 <- fread(f("14-9_cfa38_breeds.csv")); setnames(b19, c("T_count", "T_freq", "share_T"), c("minor_count", "minor_freq", "share_minor"))
mkb <- function(D, locus) D[order(-minor_count)][1:10][, .(Locus = locus, Breed = breed_nm(breed), Dogs = fn(n), `Minor alleles` = fn(minor_count), Carriers = fn(carriers),
                                                         `Minor-allele frequency` = fb(minor_freq, 3), `Median weight (kg)` = fb(wt_med, 1), `Share of minor alleles in these dogs` = fb(share_minor, 3))]
add(sprintf("**(b) Ten breeds with the most minor alleles (dogs with largest breed ancestry ≥ 0.85; the minor allele was present in %d of %d breeds at CFA19 and %d of %d at CFA38)**",
            b19[minor_count > 0, .N], nrow(b19), b38[minor_count > 0, .N], nrow(b38)), "")
add(md_tab(rbind(mkb(b19, "CFA19"), mkb(b38, "CFA38")), c("l", "l", rep("r", 6))))
L10 <- fread(f("14-9_ld_support.csv"))[p_2df < 1e-5][order(locus, dist_kb)]
add("**(c) SNPs within ± 250 kb in linkage disequilibrium with the lead SNP (r² ≥ 0.5) with joint-test p < 1×10⁻⁵**", "")
S10c <- L10[, .(Locus = locus, SNP = snp_fmt(ID), `Distance from lead (kb)` = fb(dist_kb, 1), `r² with lead` = fb(r2, 2), `Allele frequency` = fb(A1FREQ, 3),
                n = fn(N), `Genotype missingness` = sprintf("%s%%", fb(100 * f_miss, 1)), `p, joint 2-df (REGENIE)` = fp(p_2df))]
add(md_tab(S10c, c("l", "l", "r", "r", "r", "r", "r", "r")))
note("Rows with distance 0 are the lead SNPs.")

bnd <- c(chunk, length(md) + 1L)
parts <- setNames(lapply(seq_along(chunk), function(i) md[bnd[i]:(bnd[i + 1] - 1L)]), names(chunk))
chk(setequal(names(parts), paste0("Table S", 1:10)), "표 번호 S1 ~ S10이 하나씩")
md <- c(md[seq_len(chunk[1] - 1L)], unlist(parts[paste0("Table S", 1:10)], use.names = FALSE))
writeLines(md, out_md)
st <- system2("pandoc", c(shQuote(out_md), "-f", "markdown", "-t", "docx", "-o", shQuote(out_docx)), stdout = TRUE, stderr = TRUE)
if (!file.exists(out_docx)) stop("pandoc 실패: ", paste(st, collapse = "\n"))
cat(sprintf("저장: %s\n      %s\n      %s\n", out_docx, out_md, file.path(d0, "14-11_TableS6c_full.csv")))
cat(sprintf("S1 품종 ≥ 20마리 %d개 (%d마리) · 나머지 %d품종 %d마리\n", nrow(big), sum(big$n), nrow(small), sum(small$n)))
cat(sprintf("S2 단계: %s → %s → %s → %s → %s SNP · 개 %d → %d\n", fn(n_raw), fn(n_geno), fn(n_gm), fn(n_maf), fn(n_hwe), n_dog0, n_dog1))
cat(sprintf("S4b 재현 줄 %d · S6b 통과 %d · S10c 이웃 SNP 줄 %d\n", nrow(P4), nrow(pass6), nrow(L10)))
