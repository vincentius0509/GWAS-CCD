# 13-12_calcium_check.R - Calcium-response set: removal of nominating loci, leave-one-gene-out, cumulative removal, age comparison. See README.md.
suppressMessages({ library(data.table); library(ggplot2); library(parallel) })
options(width = 220)
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
qc_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
pic_dir <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/00_Picture/05_GWAS_2023"
out_prefix <- file.path(qc_dir, "13-12")
pre <- function(x) file.path(qc_dir, x)
n_cores <- as.integer(Sys.getenv("N_CORES", "8"))
magma <- "/data/00_tools/magma_v1.10/magma"
g22 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure"
f_in <- c(cand = pre("13-10_candidates.csv"), loci = pre("13-9_loci.csv"), sa = pre("06-4_sets_aux.txt"), orth = pre("06-4_ortholog_map.csv"),
          gl = file.path(g22, "06-1_gene_loc.txt"))
for (f in c(f_in, magma)) if (!file.exists(f)) stop("파일 없음: ", f)
res63 <- file.path(d0, "06-3_magma"); res137 <- Sys.getenv("MAGMA13_DIR", file.path(d0, "13-7_magma"))
tgt <- "GO_0071277"
m_sex <- c("sex_m", "sex_f"); m_logo <- c("sex_m", "sex_f", "main"); m_age <- c("young8", "age8"); m_ladder <- c("agefloor9", "age10", "agefloor11", "agefloor12")
models <- c(m_sex, "main", m_age, "int_age", m_ladder)
desc_of <- c(sex_m = "Males", sex_f = "Females", main = "All dogs", young8 = "< 8 y", age8 = ">= 8 y", int_age = "SNP x age",
             agefloor9 = ">= 9 y", age10 = ">= 10 y", agefloor11 = ">= 11 y", agefloor12 = ">= 12 y")
x_main <- 100; x_ref <- c(50, 10); p_cut <- 0.05; k_max <- 10
res_of <- function(m, w50 = FALSE, ext = ".genes.raw") file.path(if (m %in% m_ladder) res137 else res63, paste0(m, if (w50) "_w50" else "", ext))
wd <- tempfile("r13_12_"); dir.create(wd)
f3 <- function(x) { x <- copy(x); for (j in names(x)) if (is.numeric(x[[j]]) && !is.integer(x[[j]])) set(x, j = j, value = signif(x[[j]], 3)); x }
t0 <- Sys.time()

read_sets <- function(f) { s <- strsplit(readLines(f), "[\t ]+"); setNames(lapply(s, function(x) as.integer(x[-1])), sapply(s, `[`, 1)) }
S0 <- read_sets(f_in[["sa"]])[[tgt]]
if (is.null(S0)) stop("세트 없음: ", tgt)
gl <- fread(f_in[["gl"]], header = FALSE, col.names = c("GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"))[CHR %in% as.character(1:38)][, CHR := as.integer(CHR)]
orth <- unique(fread(f_in[["orth"]])[, .(GENE = dog, hs_symbol)], by = "GENE")
sym_of <- function(g) { s <- orth$hs_symbol[match(g, orth$GENE)]; fifelse(is.na(s) | s == "", gl$SYMBOL[match(g, gl$GENE)], s) }
CC <- fread(f_in[["cand"]])[set == tgt & model %in% m_sex]
if (nrow(CC) != 2) stop("13-10 후보에 두 성별의 ", tgt, "가 없음")
NOM <- CC[, .(SNP = strsplit(loci_snps, " ")[[1]]), by = model]
LO <- unique(fread(f_in[["loci"]])[, .(model, SNP, chr, lo, hi)], by = c("model", "SNP"))
NL <- LO[NOM, on = .(model, SNP)]
if (anyNA(NL$lo)) stop("13-9 locus 범위가 없는 지명 SNP: ", paste(NL[is.na(lo), SNP], collapse = " "))
near_genes <- function(R, w) unique(unlist(lapply(seq_len(nrow(R)), function(i) gl[CHR == R$chr[i] & STOP >= R$lo[i] - w * 1e3 & START <= R$hi[i] + w * 1e3, GENE])))
V <- list(x0 = S0)
for (w in c(x_main, x_ref)) V[[paste0("U", w)]] <- setdiff(S0, near_genes(NL, w))
for (m in m_sex) V[[paste0("own_", m)]] <- setdiff(S0, near_genes(NL[model == m], x_main))
U <- V[[paste0("U", x_main)]]
ex_U <- setdiff(S0, U)
cat(sprintf("[%s] %s: 유전자 %d개 · 지명 locus %d개 (sex_m %d · sex_f %d) · 합집합 ± %d kb 제거 %d개 (%s) → U%d %d개\n", format(Sys.time(), "%T"), tgt, length(S0),
            nrow(NL), NL[model == "sex_m", .N], NL[model == "sex_f", .N], x_main, length(ex_U), paste(sym_of(ex_U), collapse = ", "), x_main, length(U)))

GZ <- rbindlist(lapply(m_logo, function(m) fread(res_of(m, ext = ".genes.out"))[GENE %in% U, .(model = m, GENE, z = ZSTAT)]))

jobs <- rbindlist(lapply(models, function(m) {
  S <- c(V, setNames(list(setdiff(U, near_genes(LO[model == m], x_main))), paste0("U", x_main, "own")))
  if (m %in% m_logo) {
    S <- c(S, setNames(lapply(U, function(g) setdiff(U, g)), paste0("logo_", U)))
    zz <- GZ[model == m][order(-z)]$GENE
    S <- c(S, setNames(lapply(seq_len(min(k_max, length(zz))), function(k) setdiff(U, zz[seq_len(k)])), paste0("topk_", seq_len(min(k_max, length(zz))))))
  }
  sf <- file.path(wd, paste0(m, ".sets")); writeLines(vapply(names(S), function(n) paste(c(n, S[[n]]), collapse = " "), ""), sf)
  data.table(model = m, w50 = c(FALSE, TRUE), sf = sf)
}))
jobs[, raw := mapply(res_of, model, w50)]
miss <- jobs[!file.exists(raw)]
if (nrow(miss)) cat("MAGMA 유전자 결과 없음 → 건너뜀:", paste(basename(miss$raw), collapse = " "), "\n")
jobs <- jobs[file.exists(raw)]
jobs[, out := file.path(wd, paste0(model, fifelse(w50, "_w50", "")))]
st <- unlist(mclapply(seq_len(nrow(jobs)), function(i) system2(magma, c("--gene-results", shQuote(jobs$raw[i]), "--set-annot", shQuote(jobs$sf[i]), "--out", shQuote(jobs$out[i])),
                                                              stdout = FALSE, stderr = FALSE), mc.cores = n_cores, mc.preschedule = FALSE))
if (any(st != 0)) stop("MAGMA 실패: ", paste(basename(jobs$out[st != 0]), collapse = ", "))
read_gsa <- function(f) { d <- fread(f, skip = "VARIABLE"); if ("FULL_NAME" %in% names(d)) d[, VARIABLE := FULL_NAME]; d[, .(ver = VARIABLE, NGENES, BETA, SE, P)] }
GA <- rbindlist(lapply(seq_len(nrow(jobs)), function(i) read_gsa(paste0(jobs$out[i], ".gsa.out"))[, `:=`(model = jobs$model[i], w50 = jobs$w50[i])]))
cat(sprintf("[%s] MAGMA 세트 분석 %d번 (세트 판 %d개 결과)\n", format(Sys.time(), "%T"), nrow(jobs), nrow(GA)))

core <- c(paste0("U", x_main), paste0("U", x_ref), "x0", paste0("own_", m_sex), paste0("U", x_main, "own"))
TT <- GA[w50 == FALSE & ver %in% core]
TT <- merge(TT, GA[w50 == TRUE & ver %in% core, .(model, ver, P_w50 = P)], by = c("model", "ver"), all.x = TRUE)
TT[, `:=`(model = factor(model, models), ver = factor(ver, core))]; setorder(TT, ver, model)
TT[, desc := desc_of[as.character(model)]]
get1 <- function(m, v = paste0("U", x_main), col = "P") TT[model == m & ver == v][[col]][1]
dz <- function(a, b, v = paste0("U", x_main)) { z <- (get1(a, v, "BETA") - get1(b, v, "BETA")) / sqrt(get1(a, v, "SE")^2 + get1(b, v, "SE")^2); c(z = z, p = 2 * pnorm(-abs(z))) }

LG <- GA[w50 == FALSE & grepl("^logo_", ver), .(model, GENE = as.integer(sub("logo_", "", ver)), NGENES, BETA, P)]
LG[, symbol := sym_of(GENE)]; LG <- merge(LG, GZ, by = c("model", "GENE"), all.x = TRUE)
LG[, model := factor(model, m_logo)]; setorder(LG, model, -P)
TK <- GA[w50 == FALSE & grepl("^topk_", ver), .(model, k = as.integer(sub("topk_", "", ver)), NGENES, BETA, P)]
TK <- merge(TK, GZ[order(-z), .(k = seq_len(.N), removed = sym_of(GENE), z), by = model], by = c("model", "k"), all.x = TRUE)
TK[, model := factor(model, m_logo)]; setorder(TK, model, k)

pU <- setNames(sapply(models, get1), models)
v1 <- with(list(a = pU[["sex_m"]] < p_cut, b = pU[["sex_f"]] < p_cut), fifelse(a & b, "두 성별 지명 locus를 모두 빼도 양쪽 성별에서 유지",
                                                                          fifelse(a | b, "한쪽 성별만 유지", "지명 locus 의존")))
v2 <- fifelse(pU[["main"]] < p_cut, "전체 표본에서도 일관", "전체 표본에서 근거 없음")
lo_v <- LG[, .(max_p = max(P), worst = symbol[which.max(P)], n_lost = sum(P >= p_cut), lost = paste(symbol[P >= p_cut], collapse = " · ")), by = model]
lo_v[, verdict := fifelse(max_p < p_cut, "한 유전자에 의존하지 않음", paste0("그 유전자 의존: ", lost))]
k_lost <- TK[, .(k_lost = if (any(P >= p_cut)) min(k[P >= p_cut]) else NA_integer_), by = model]
d_age <- dz("age8", "young8"); d_sex <- dz("sex_m", "sex_f")
v4 <- fifelse(d_age[["z"]] > 0 & d_age[["p"]] < p_cut, "8세 이상에서 더 큼 (탐색)", "나이 차이 근거 없음")
VD <- rbind(
  data.table(test = "① 합집합 제거", models = "sex_m · sex_f", stat = sprintf("U%d p: 수컷 %.3g · 암컷 %.3g (성별 차이 p %.2g, 참고)", x_main, pU[["sex_m"]], pU[["sex_f"]], d_sex[["p"]]), verdict = v1),
  data.table(test = "② 전체 표본", models = "main", stat = sprintf("U%d p %.3g (main 자신의 locus도 뺀 판 %.3g, 참고)", x_main, pU[["main"]], get1("main", paste0("U", x_main, "own"))), verdict = v2),
  lo_v[, .(test = "③ 유전자 하나씩 빼기", models = as.character(model),
           stat = sprintf("유전자 %d개 각각 뺀 p의 최대 %.3g (%s 뺄 때)", length(U), max_p, worst), verdict)],
  data.table(test = "④ 나이 (탐색)", models = "age8 vs young8", stat = sprintf("β 8세+ %.3f · 8세 미만 %.3f, 차이 z %.2f p %.2g (참고: int_age p %.3g)",
             get1("age8", col = "BETA"), get1("young8", col = "BETA"), d_age[["z"]], d_age[["p"]], pU[["int_age"]]), verdict = v4))
VD <- merge(VD, k_lost[, .(models = as.character(model), k_lost_ref = k_lost)], by = "models", all.x = TRUE, sort = FALSE)

fwrite(f3(TT), paste0(out_prefix, "_tests.csv")); fwrite(f3(LG), paste0(out_prefix, "_logo.csv"))
fwrite(f3(TK), paste0(out_prefix, "_topk.csv")); fwrite(VD, paste0(out_prefix, "_verdict.csv"))

dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)
grp_of <- function(m) factor(fifelse(m %in% m_sex, "Sex (non-overlapping)", fifelse(m == "main", "All dogs", fifelse(m %in% m_age, "Age (non-overlapping)", "Age floor (overlapping)"))),
                             c("Sex (non-overlapping)", "All dogs", "Age (non-overlapping)", "Age floor (overlapping)"))
P1 <- TT[ver == paste0("U", x_main) & model != "int_age"]
P1[, `:=`(grp = grp_of(as.character(model)), lab = factor(desc, rev(desc_of[setdiff(models, "int_age")])), ptxt = sprintf("P = %.2g", P))]
p1 <- ggplot(P1, aes(BETA, lab, colour = grp)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = BETA - 1.96 * SE, xmax = BETA + 1.96 * SE), height = 0, linewidth = 0.6) + geom_point(size = 2.4) +
  geom_text(aes(label = ptxt), colour = "grey25", size = 2.8, hjust = 0, nudge_y = 0.32) +
  scale_colour_manual(values = c("Sex (non-overlapping)" = "#2f6f9f", "All dogs" = "#222222", "Age (non-overlapping)" = "#c26a1b", "Age floor (overlapping)" = "#b9a17a"), name = NULL) +
  labs(x = sprintf("MAGMA competitive beta (+/- 95%% CI); nominating loci of both sexes +/- %d kb removed", x_main), y = NULL,
       title = sprintf("Cellular response to calcium ion (GO:0071277, %d of %d genes)", length(U), length(S0))) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggsave(file.path(pic_dir, "13-12_strata.png"), p1, width = 7.5, height = 5, dpi = 300)

P2 <- copy(LG)[, desc := factor(desc_of[as.character(model)], desc_of[m_logo])]
P2[, ord := frank(-P, ties.method = "first"), by = model]
ref2 <- TT[ver == paste0("U", x_main) & model %in% m_logo, .(desc = factor(desc_of[as.character(model)], desc_of[m_logo]), P)]
note2 <- P2[ord <= 3, .(txt = paste0("Largest P when removed:\n", paste(sprintf("%s (P = %.2g)", symbol, P), collapse = "\n"))), by = desc]
p2 <- ggplot(P2, aes(ord, -log10(P))) + geom_hline(yintercept = -log10(p_cut), linetype = "dashed", colour = "grey45") +
  geom_hline(data = ref2, aes(yintercept = -log10(P)), colour = "#2f6f9f", linewidth = 0.5) +
  geom_point(size = 1.3, colour = "grey30") +
  geom_text(data = note2, aes(x = Inf, y = 0, label = txt), hjust = 1.05, vjust = -0.2, size = 2.6, colour = "grey15", lineheight = 0.95, inherit.aes = FALSE) +
  facet_wrap(~ desc, nrow = 1) + expand_limits(y = 0) +
  labs(x = "Genes, ordered by P after removing that gene", y = "-log10 P (set without one gene)",
       title = "Leave-one-gene-out: blue = full set, dashed = P 0.05") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
ggsave(file.path(pic_dir, "13-12_logo.png"), p2, width = 9, height = 3.8, dpi = 300)

sink(paste0(out_prefix, "_summary.txt"))
cat("13-12 칼슘 경로 확인 분석 (GO:0071277 cellular response to calcium ion; 판정 규칙은 실행 전 고정, 헤더 참고)\n")
cat(sprintf("세트 유전자 %d개 · 두 성별 지명 locus %d개 (%s) · 합집합 ± %d kb 제거 %d개 (%s) → U%d %d개\n", length(S0), nrow(NL), paste(NL$SNP, collapse = " "),
            x_main, length(ex_U), paste(sym_of(ex_U), collapse = ", "), x_main, length(U)))
cat("\n■ 판정\n"); print(VD, row.names = FALSE)
cat("\n■ 모형 × 세트 판 (U100 = 주; 나머지 참고)\n")
print(f3(dcast(TT, model + desc ~ ver, value.var = "P")), row.names = FALSE)
cat("\n■ U100 β · SE · ± 50 kb 창 p\n"); print(f3(TT[ver == paste0("U", x_main), .(model, desc, NGENES, BETA, SE, P, P_w50)]), row.names = FALSE)
cat("\n■ ③ 빼면 p가 가장 커지는 유전자 (모형별 상위 5)\n"); print(f3(LG[, head(.SD, 5), by = model][, .(model, symbol, z, P)]), row.names = FALSE)
cat("\n■ ③ 참고: 유전자 z 상위 k개를 차례로 뺄 때\n"); print(f3(dcast(TK, k ~ model, value.var = "P")), row.names = FALSE)
cat("  상위 k개 이름:\n"); print(TK[, .(top = paste(removed, collapse = " > ")), by = model], row.names = FALSE)
sink()
unlink(wd, recursive = TRUE)
cat(sprintf("[%s] 끝 (%.1f분). 저장: %s_* (tests · logo · topk · verdict · summary) / 그림 13-12_strata.png · 13-12_logo.png\n", format(Sys.time(), "%T"),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), out_prefix))
print(VD, row.names = FALSE)

