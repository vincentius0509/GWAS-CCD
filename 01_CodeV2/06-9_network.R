# 06-9_network.R - Network propagation of gene-level signals on a protein interaction network. See README.md.
# - null

library(data.table)
library(Matrix)
library(ggplot2)
options(width = 220)
set.seed(20261112)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
out_prefix <- file.path(qc_dir, "06-9")
pic_prefix <- file.path(pic_dir, "06-9")
magma_dir <- file.path(in_dir, "06-3_magma")
ref_dir <- file.path(in_dir, "06-9_ref")
links_f <- file.path(ref_dir, "9606.protein.links.v12.0.txt.gz")
info_f <- file.path(ref_dir, "9606.protein.info.v12.0.txt.gz")
orth_f <- file.path(in_dir, "06-4_ortholog_map.csv")
sets_f <- file.path(in_dir, "06-4_sets_primary.txt")
n_perm <- as.integer(Sys.getenv("N_PERM", "10000"))
score_min <- 700; close_bp <- 1e6; alpha_main <- 0.5; alpha_sens <- c(0.3, 0.7); set_min <- 10; n_deg_bin <- 10
runs_main <- list(`1 (main)` = "main", `2 (age)` = c("int_age", "het_age"), `3 (explore)` = "age8")
runs_w50 <- c(main = "main_w50", int_age = "int_age_w50", het_age = "het_age_w50")
runs_null <- paste0("null", 1:20)
all_runs <- c(unlist(runs_main), runs_w50, runs_null)
for (f in c(links_f, info_f, orth_f, sets_f, file.path(magma_dir, paste0(all_runs, ".genes.out"))))
  if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE); dir.create(pic_dir, showWarnings = FALSE, recursive = TRUE)

read_gz <- function(f, ...) fread(cmd = paste("gzip -dc", shQuote(f)), ...)
info <- read_gz(info_f, select = 1:2, col.names = c("id", "sym"))
ln <- read_gz(links_f)[combined_score >= score_min]
om <- fread(orth_f, colClasses = c(dog = "character"))
if (anyDuplicated(om$hs_symbol) || anyDuplicated(om$dog)) stop("ortholog 표가 1:1이 아님")
genes <- fread(file.path(magma_dir, "main.genes.out"), colClasses = c(GENE = "character"))
ln[, `:=`(d1 = om$dog[match(info$sym[match(protein1, info$id)], om$hs_symbol)], d2 = om$dog[match(info$sym[match(protein2, info$id)], om$hs_symbol)])]
e <- unique(ln[!is.na(d1) & !is.na(d2) & d1 != d2 & d1 %in% genes$GENE & d2 %in% genes$GENE, .(a = pmin(d1, d2), b = pmax(d1, d2))])
n_edge_mapped <- nrow(e)
pos <- genes[, .(GENE, CHR, mid = (START + STOP) / 2)]
e <- merge(merge(e, pos, by.x = "a", by.y = "GENE"), pos, by.x = "b", by.y = "GENE", suffixes = c(".a", ".b"))
nose <- e[CHR.a == CHR.b & abs(mid.a - mid.b) < close_bp, .N]
e <- e[!(CHR.a == CHR.b & abs(mid.a - mid.b) < close_bp)]
nodes <- sort(unique(c(e$a, e$b))); n <- length(nodes)
A <- sparseMatrix(i = match(e$a, nodes), j = match(e$b, nodes), x = 1, dims = c(n, n), symmetric = TRUE)
A <- as(A, "generalMatrix")
deg <- rowSums(A)
Dm <- Diagonal(x = 1 / sqrt(deg)); W <- Dm %*% A %*% Dm
Wrow <- Diagonal(x = 1 / deg) %*% A
deg_bin <- cut(rank(deg, ties.method = "first"), n_deg_bin, labels = FALSE)
bin_idx <- split(seq_len(n), deg_bin)
net_info <- data.table(string_edges_ge700 = nrow(ln), dog_edges = n_edge_mapped, close_edges_removed = nose, edges = nrow(e),
                       nodes = n, magma_genes = nrow(genes), pct_genes_in_net = 100 * n / nrow(genes),
                       deg_median = median(deg), deg_max = max(deg))
cat(sprintf("네트워크: STRING ≥ %d 연결 %s개 → 개 유전자 연결 %s개 (1 Mb 안 %s개 제외) → 유전자 %s개 (MAGMA 유전자의 %.1f%%), 차수 중앙값 %g · 최대 %g\n",
            score_min, format(nrow(ln), big.mark = ","), format(n_edge_mapped, big.mark = ","), format(nose, big.mark = ","),
            format(n, big.mark = ","), net_info$pct_genes_in_net, net_info$deg_median, net_info$deg_max))

propagate <- function(F0, alpha, tol = 1e-10, maxit = 1000) {
  F0 <- as.matrix(F0); Fk <- F0
  for (it in seq_len(maxit)) { Fn <- alpha * as.matrix(W %*% Fk) + (1 - alpha) * F0; if (max(abs(Fn - Fk)) < tol) return(Fn); Fk <- Fn }
  stop("전파가 수렴하지 않음")
}

sl <- strsplit(readLines(sets_f), "[ \t]+"); sets <- setNames(lapply(sl, `[`, -1), vapply(sl, `[`, "", 1))
set_tab <- data.table(set = names(sets), n_genes = lengths(sets), n_in_net = vapply(sets, function(s) sum(s %in% nodes), 0L))
use_sets <- set_tab[n_in_net >= set_min, set]
S <- sapply(use_sets, function(s) { v <- numeric(n); i <- match(intersect(sets[[s]], nodes), nodes); v[i] <- 1 / length(i); v })
U <- lapply(c(alpha_main, alpha_sens), function(a) propagate(S, a))
names(U) <- as.character(c(alpha_main, alpha_sens))
cat(sprintf("세트 %d개 중 네트워크 유전자 ≥ %d개인 세트 %d개 (빠진 세트: %s)\n", nrow(set_tab), set_min, length(use_sets),
            if (length(setdiff(set_tab$set, use_sets))) paste(setdiff(set_tab$set, use_sets), collapse = ", ") else "없음"))

read_z <- function(run) {
  g <- fread(file.path(magma_dir, paste0(run, ".genes.out")), colClasses = c(GENE = "character"))
  z <- g$ZSTAT[match(nodes, g$GENE)]
  if (anyNA(z)) { warning(run, ": 네트워크 유전자 ", sum(is.na(z)), "개의 z가 없어 0으로 둠"); z[is.na(z)] <- mean(z, na.rm = TRUE) }
  (z - mean(z)) / sd(z)
}
perm_index <- replicate(n_perm, { p <- seq_len(n); for (b in bin_idx) p[b] <- b[sample.int(length(b))]; p })
rnet <- function(Z) { Z <- as.matrix(Z); NZ <- as.matrix(Wrow %*% Z); colSums(scale(Z) * scale(NZ)) / (nrow(Z) - 1) }
analyse <- function(run) {
  z <- read_z(run)
  t_obs <- lapply(U, function(u) drop(crossprod(u, z)))
  Zp <- matrix(z[perm_index], n, n_perm)
  t_perm <- lapply(U, function(u) crossprod(u, Zp))
  p_set <- mapply(function(o, pm) (rowSums(pm >= o) + 1) / (n_perm + 1), t_obs, t_perm, SIMPLIFY = FALSE)
  sets_dt <- rbindlist(lapply(names(U), function(a) data.table(run = run, alpha = as.numeric(a), set = use_sets, t = t_obs[[a]],
                                                              t_null_mean = rowMeans(t_perm[[a]]), t_null_sd = apply(t_perm[[a]], 1, sd), p = p_set[[a]])))
  sets_dt[, z_perm := (t - t_null_mean) / t_null_sd]
  r_obs <- rnet(z)
  r_perm <- unlist(lapply(split(seq_len(n_perm), ceiling(seq_len(n_perm) / 1000)), function(k) rnet(Zp[, k, drop = FALSE])))
  glob <- data.table(run = run, r_net = r_obs, r_null_mean = mean(r_perm), r_null_sd = sd(r_perm), p = (sum(r_perm >= r_obs) + 1) / (n_perm + 1))
  top <- data.table(run = run, GENE = nodes, z = z, F = drop(propagate(z, alpha_main)), degree = deg)
  list(sets = sets_dt, glob = glob, top = top)
}
t0 <- Sys.time()
res <- lapply(setNames(all_runs, all_runs), function(r) { x <- analyse(r); cat(".", sep = ""); x }); cat("\n")
cat(sprintf("순열 %s회 × 모형 %d개 끝 (%.1f분)\n", format(n_perm, big.mark = ","), length(all_runs), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sets_all <- rbindlist(lapply(res, `[[`, "sets")); glob_all <- rbindlist(lapply(res, `[[`, "glob"))
sets_all <- merge(sets_all, set_tab, by = "set")

null_calib <- sets_all[run %in% runs_null & alpha == alpha_main, .(null_fp05 = mean(p < 0.05), null_min_p = min(p)), by = set]
null_glob <- glob_all[run %in% runs_null]
cat(sprintf("null 20개: 세트 p < 0.05 비율 %.3f (기대 0.05; 세트별 최대 %.2f), 전역 r_net %.4f ~ %.4f, 전역 p < 0.05 %d / 20\n",
            sets_all[run %in% runs_null & alpha == alpha_main, mean(p < 0.05)], max(null_calib$null_fp05),
            min(null_glob$r_net), max(null_glob$r_net), sum(null_glob$p < 0.05)))

stage_of <- function(st) {
  d <- sets_all[run %in% runs_main[[st]] & alpha == alpha_main, .(stage = st, run, set, n_genes, n_in_net, t, z_perm, p)]
  if (st == "3 (explore)") d[, q := p.adjust(p, "BH"), by = run] else d[, q := p.adjust(p, "BH")]
  sens <- rbindlist(list(
    sets_all[alpha %in% alpha_sens & run %in% runs_main[[st]], .(run, set, key = paste0("p_alpha", alpha), p)],
    sets_all[alpha == alpha_main & run %in% runs_w50, .(run = names(runs_w50)[match(run, runs_w50)], set, key = "p_w50", p)]))
  if (nrow(sens)) d <- merge(d, dcast(sens, run + set ~ key, value.var = "p"), by = c("run", "set"), all.x = TRUE)
  pc <- intersect(c("p_w50", "p_alpha0.3", "p_alpha0.7"), names(d))
  d[, holds := p < 0.05 & Reduce(`&`, lapply(pc, function(k) !is.na(get(k)) & get(k) < 0.05))]
  d <- merge(d, null_calib, by = "set", all.x = TRUE)
  d[, null_ok := null_fp05 <= 0.10]
  d[, pass := q < 0.05 & holds & null_ok]
  setorder(d, p)[]
}
stages <- lapply(setNames(names(runs_main), names(runs_main)), stage_of)

cal <- sets_all[run %in% runs_null & alpha == alpha_main, .(null_mean_z = mean(z_perm), null_sd_z = sd(z_perm)), by = set]
for (st in names(stages)) {
  d <- merge(stages[[st]], cal, by = "set", all.x = TRUE)
  d[, `:=`(z_cal = (z_perm - null_mean_z) / null_sd_z)][, p_cal := pnorm(z_cal, lower.tail = FALSE)]
  if (st == "3 (explore)") d[, q_cal := p.adjust(p_cal, "BH"), by = run] else d[, q_cal := p.adjust(p_cal, "BH")]
  stages[[st]] <- setorder(d, p)[]
}
decomp <- function(run, set) {
  z <- read_z(run); k <- match(set, use_sets); u <- U[[as.character(alpha_main)]][, k]; inset <- nodes %in% sets[[set]]
  Zp <- matrix(z[perm_index], n, n_perm)
  pp <- function(w) (sum(drop(crossprod(w, Zp)) >= sum(w * z)) + 1) / (n_perm + 1)
  cz <- u * z; ord <- order(-cz)[1:5]
  data.table(run = run, set = set, t_total = sum(u * z), t_self = sum((u * inset) * z), t_neighbour = sum((u * !inset) * z),
             weight_self = sum(u[inset]) / sum(u), p_self = pp(u * inset), p_neighbour = pp(u * !inset),
             top5_genes = paste(om$hs_symbol[match(nodes[ord], om$dog)], collapse = ","), top5_share = sum(cz[ord]) / sum(cz))
}
dec_rows <- rbindlist(lapply(stages, function(d) d[p < 0.05, .(run, set)]))
decomp_tab <- if (nrow(dec_rows)) rbindlist(Map(decomp, dec_rows$run, dec_rows$set)) else data.table()
nullg <- glob_all[run %in% runs_null]
glob_all[, `:=`(z_cal = (r_net - mean(nullg$r_net)) / sd(nullg$r_net))][, p_cal := pnorm(z_cal, lower.tail = FALSE)]
robust <- rbindlist(lapply(setdiff(glob_all$run, runs_null), function(r) {
  z <- read_z(r); zr <- qnorm((rank(z) - 0.5) / n); zc <- pmin(pmax(z, -2.5), 2.5)
  data.table(run = r, r_net_rank = rnet(zr), r_netip2.5 = rnet(zc))
}))
glob_all <- merge(glob_all, robust, by = "run", all.x = TRUE)

top_main <- res[["main"]]$top[order(-F)][1:30]
in_set <- vapply(top_main$GENE, function(g) paste(use_sets[vapply(use_sets, function(s) g %in% sets[[s]], TRUE)], collapse = ";"), "")
top_main[, `:=`(rank = .I, symbol = om$hs_symbol[match(GENE, om$dog)], in_primary_set = in_set, magma_p = 2 * pnorm(-abs(z)))]
top_main[, magma_p := NULL]
top_main <- merge(top_main, genes[, .(GENE, CHR, START, P_magma = P)], by = "GENE")[order(rank)]

fwrite(net_info, paste0(out_prefix, "_network_info.csv"))
fwrite(sets_all, paste0(out_prefix, "_sets.csv"))
for (st in names(stages)) fwrite(stages[[st]], paste0(out_prefix, "_stage", substr(st, 1, 1), ".csv"))
fwrite(glob_all, paste0(out_prefix, "_global.csv"))
fwrite(null_calib, paste0(out_prefix, "_null_calib.csv"))
fwrite(top_main, paste0(out_prefix, "_top_genes.csv"))
fwrite(decomp_tab, paste0(out_prefix, "_decomp.csv"))

hm <- sets_all[alpha == alpha_main & run %in% c(unlist(runs_main), runs_w50)]
hm[, run := factor(run, levels = c(unlist(runs_main), runs_w50))]
g1 <- ggplot(hm, aes(run, set, fill = -log10(p))) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(p < 0.01, sprintf("%.1e", p), sprintf("%.2f", p))), size = 2.6, colour = "grey15") +
  scale_fill_gradient(low = "grey95", high = "#1f4e9c", limits = c(0, 4), oob = scales::squish) +
  labs(x = NULL, y = NULL, title = "Network-propagated set scores (RWR, alpha 0.5, STRING >= 700)",
       subtitle = "cell = degree-preserving permutation p (one-sided); _w50 = +/-50 kb gene z") +
  theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(pic_prefix, "_sets.png"), g1, width = 8.5, height = 5.5, dpi = 300)
gd <- copy(glob_all)[, type := fifelse(run %in% runs_null, "null GWAS (20)", "model")]
g2 <- ggplot(gd, aes(r_net, reorder(run, r_net), colour = type)) + geom_point(size = 2) +
  geom_vline(xintercept = max(null_glob$r_net), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c(model = "#1f4e9c", `null GWAS (20)` = "grey60")) +
  labs(x = "network autocorrelation r_net = cor(z, neighbour mean z)", y = NULL, colour = NULL,
       title = "Do association signals cluster in the PPI network?", subtitle = "dashed = max of 20 null GWAS") +
  theme_minimal(base_size = 10)
ggsave(paste0(pic_prefix, "_global.png"), g2, width = 7, height = 7, dpi = 300)

f3 <- function(d) { d <- copy(d); for (j in names(d)) if (is.numeric(d[[j]]) && !is.integer(d[[j]])) set(d, j = j, value = signif(d[[j]], 3)); d }
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("06-9 네트워크 전파 (STRING v12 사람 ≥ 700 → 개 1:1 ortholog, RWR α 0.5, 차수 보존 순열", format(n_perm, big.mark = ","), "회)\n")
cat("탐색 분석 (06 결과를 본 뒤 추가; 설정·판정은 실행 전 헤더에 고정)\n\n■ 네트워크\n"); print(f3(net_info))
cat("\n■ 세트 (네트워크 유전자 수)\n"); print(set_tab)
cat("\n■ null 보정 (표현형 순열 GWAS 20개, α 0.5)\n"); print(f3(null_calib[order(-null_fp05)]))
cat(sprintf("  전체 세트 p < 0.05 비율 %.3f (기대 0.05). null_ok = null 20개 중 p < 0.05가 2개 이하 (10%% 이하)\n",
            sets_all[run %in% runs_null & alpha == alpha_main, mean(p < 0.05)]))
for (st in names(stages)) {
  cat(sprintf("\n■ %s단계 (%s): BH q, holds = ±50 kb와 α 0.3·0.7에서도 p < 0.05, pass = q < 0.05 & holds & null_ok\n", substr(st, 1, 1), paste(runs_main[[st]], collapse = " + ")))
  print(f3(stages[[st]][, intersect(c("run", "set", "n_in_net", "z_perm", "p", "q", "p_w50", "p_alpha0.3", "p_alpha0.7", "holds", "null_fp05", "null_ok", "pass", "z_cal", "p_cal", "q_cal"), names(stages[[st]])), with = FALSE]))
}
cat("  (실행 뒤 추가한 진단) z_cal·p_cal·q_cal = null GWAS 20개의 z_perm 평균·SD로 다시 표준화한 값과 그 단계 BH. pass 판정에는 쓰지 않음\n")
cat("\n■ (실행 뒤 추가한 진단) 명목 p < 0.05 세트의 자기·이웃 분해 (t_total = t_self + t_neighbour; 각각 차수 보존 순열 p)\n")
if (nrow(decomp_tab)) print(f3(decomp_tab)) else cat("  해당 없음\n")
cat("  p_neighbour가 크면 네트워크 이웃이 더한 신호가 없음 → 06-5 세트 검정과 같은 신호를 다른 null로 본 것\n")
cat("\n■ 전역: 연관 신호가 네트워크 이웃끼리 뭉치는가 (r_net; 차수 보존 p, null GWAS 20개 범위)\n")
print(f3(glob_all[!run %in% runs_null]))
cat(sprintf("  null GWAS 20개 r_net: 평균 %.4f, 범위 %.4f ~ %.4f, 차수 보존 p < 0.05가 %d / 20개\n",
            mean(null_glob$r_net), min(null_glob$r_net), max(null_glob$r_net), sum(null_glob$p < 0.05)))
cat("  판정: 모형의 r_net이 null 20개 최대값을 넘을 때만 \"뭉침\"으로 봄 (차수 보존 p는 LD 상관을 깨서 너무 작을 수 있음)\n")
cat("  (실행 뒤 추가한 진단) z_cal·p_cal = null 20개 평균·SD 기준. r_net_rank·r_netip2.5 = 순위 변환·±2.5 자름 (극단 유전자 영향 확인)\n")
cat("\n■ 주 분석 전파 점수 상위 30개 (기술용; 검정 아님. P_magma = 06-3 유전자 p)\n")
print(f3(top_main[, .(rank, symbol, GENE, CHR, z, F, degree, P_magma, in_primary_set)]))
cat("\n■ 판정\n")
for (st in names(stages)) cat(sprintf("  %s단계: pass %d개 (q < 0.05 %d개, 명목 p < 0.05 %d개 / %d)\n", substr(st, 1, 1),
                                      sum(stages[[st]]$pass, na.rm = TRUE), sum(stages[[st]]$q < 0.05), sum(stages[[st]]$p < 0.05), nrow(stages[[st]])))
clus <- glob_all[!run %in% runs_null & r_net > max(null_glob$r_net), run]
cat(sprintf("  전역 뭉침 (null 최대값 초과): %s\n", if (length(clus)) paste(clus, collapse = ", ") else "없음"))
cat("  통과해도 \"탐색적 근거; 독립 표본 확인 필요\"로만 씀\n")
sink()
cat(sprintf("저장: %s_network_info.csv · _sets.csv · _stage1~3.csv · _global.csv · _null_calib.csv · _top_genes.csv · _summary.txt / 그림 %s_sets.png · _global.png\n", out_prefix, pic_prefix))

