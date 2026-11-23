# 04-4_pc_choice.R - Summarises lambda of the simulated phenotypes to guide the number of PCs. See README.md.
library(data.table)
library(ggplot2)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")  
pic_dir <- file.path(base_dir, "GeneticData/analysis/00_Picture/05_GWAS_2023")
out_prefix <- file.path(qc_dir, "04-4")
pic_prefix <- file.path(pic_dir, "04-4")
lambda_cut <- 1.05

theme_qc <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

pc_sets <- c(0, seq(5, 30, by = 5)); n_sim <- 20
n_ind <- nrow(fread(file.path(qc_dir, "03-3_maf0.01_hwe.fam"), header = FALSE, select = 1))
n_test <- length(readLines(file.path(qc_dir, "04-3_sim", "test_snps.txt")))

lam <- fread(file.path(qc_dir, "04-3_sim", "lambda.tsv"), colClasses = "character", fill = TRUE)
stopifnot(identical(names(lam)[1:5], c("n_pc", "sim", "lambda", "vg_vp", "n_snp")))
bad <- lam[, which(!grepl("^[0-9]+$", n_pc) | !grepl("^[0-9]+$", sim) | !grepl("^[0-9.]+$", lambda) | !grepl("^[0-9]+$", n_snp))]
if (ncol(lam) > 5) bad <- union(bad, which(rowSums(lam[, -(1:5)] != "" & !is.na(lam[, -(1:5)])) > 0))
if (length(bad)) { print(lam[bad]); stop("lambda.tsv 형식 오류 ", length(bad), "줄 (위 출력) → 04-3 다시 실행 필요") }
lam <- lam[, .(n_pc = as.integer(n_pc), sim = as.integer(sim), lambda = as.numeric(lambda), n_snp = as.integer(n_snp))]
dup <- lam[duplicated(lam[, .(n_pc, sim)])]
if (nrow(dup)) { print(dup); stop("같은 n_pc·sim 조합이 여러 번 기록됨") }
miss <- CJ(n_pc = pc_sets, sim = seq_len(n_sim))[!lam, on = .(n_pc, sim)]
if (nrow(miss)) { print(miss); stop("빠진 조합 ", nrow(miss), "개 → RESUME=1 bash 04-3_Sim_MLMA.sh 로 채운 뒤 다시 실행") }
extra <- lam[!n_pc %in% pc_sets | !sim %in% seq_len(n_sim)]
if (nrow(extra)) { print(extra); stop("예상하지 않은 조합이 있음 (pc_sets·n_sim 설정 확인)") }
if (any(lam$n_snp < 0.99 * n_test)) { print(lam[n_snp < 0.99 * n_test]); warning("λ 계산에 쓴 SNP가 검정 SNP(", n_test, ")보다 1% 넘게 적은 실행이 있음") }
if (any(lam$lambda <= 0 | !is.finite(lam$lambda))) stop("λ가 0 이하이거나 숫자가 아님")

log_check <- function(k, s) {
  f <- file.path(qc_dir, "04-3_sim", "log", sprintf("PC%d_SIM%d.log", k, s))
  out <- list(n_common = NA_integer_, converged = FALSE, constrained = NA, vg_vp = NA_real_)
  if (!file.exists(f)) return(out)
  x <- readLines(f)
  nc <- grep("individuals are in common", x, value = TRUE)
  if (length(nc)) out$n_common <- as.integer(sub("^([0-9]+) .*", "\\1", nc[length(nc)]))
  out$converged <- any(grepl("Log-likelihood ratio converged", x))
  out$constrained <- any(grepl("constrained", x, ignore.case = TRUE))
  i <- grep("^Iter\\.", x)
  if (out$converged && length(i)) {
    r <- strsplit(trimws(x[(i[1] + 1):length(x)]), "\\s+")
    r <- Filter(function(v) length(v) >= 4 && grepl("^[0-9]+$", v[1]), r)
    if (length(r)) { v <- suppressWarnings(as.numeric(r[[length(r)]][3:4])); if (all(is.finite(v)) && sum(v) > 0) out$vg_vp <- v[1] / sum(v) }
  }
  out
}
lc <- rbindlist(mapply(function(k, s) as.data.table(log_check(k, s)), lam$n_pc, lam$sim, SIMPLIFY = FALSE))
lam <- cbind(lam, lc)
chk <- lam[is.na(n_common) | n_common != n_ind | !converged | is.na(vg_vp) | constrained %in% TRUE]
if (nrow(chk)) {
  print(chk)
  if (any(is.na(chk$n_common) | chk$n_common != n_ind)) stop("로그가 없거나 분석 개체 수가 ", n_ind, "가 아닌 실행이 있음 (위 출력)")
  warning("수렴 문구 없음·V(G)/Vp 계산 불가·경계값 고정 실행이 있음 (위 출력). λ는 유효하지만 vg_vp 해석에 주의")
}
cat(sprintf("로그 검사: %d개 실행, 개체 수 %d 일치, 수렴 %d, 경계값 고정 %d\n",
            nrow(lam), n_ind, sum(lam$converged), sum(lam$constrained %in% TRUE)))
sum_tab <- lam[, .(n_sim = .N, lambda_mean = round(mean(lambda), 3), lambda_sd = round(sd(lambda), 3),
                   lambda_max = round(max(lambda), 3), vg_vp_mean = round(mean(vg_vp, na.rm = TRUE), 3),
                   n_vg = sum(!is.na(vg_vp))), keyby = n_pc]
print(sum_tab)

#     n_pc n_sim lambda_mean lambda_sd lambda_max vg_vp_mean  n_vg
#    <int> <int>       <num>     <num>      <num>      <num> <int>
# 1:     0    20       1.001     0.017      1.038      0.404    20
# 2:     5    20       0.997     0.016      1.036      0.395    20
# 3:    10    20       0.995     0.016      1.032      0.388    20
# 4:    15    20       0.993     0.016      1.034      0.378    20
# 5:    20    20       0.991     0.016      1.023      0.362    20
# 6:    25    20       0.988     0.018      1.028      0.352    20
# 7:    30    20       0.985     0.017      1.016      0.333    20

lam_pc <- fread(file.path(qc_dir, "03-5_pca_lambda.csv"))[label %in% c("none", "pcair"), .(n_pc, lambda_pc_only = lambda_mean)]
if (!all(pc_sets %in% lam_pc$n_pc)) warning("03-5_pca_lambda.csv에 없는 PC 개수: ", paste(setdiff(pc_sets, lam_pc$n_pc), collapse = ", "))
sum_tab <- merge(sum_tab, lam_pc, by = "n_pc", all.x = TRUE)
print(sum_tab[, .(n_pc, lambda_pc_only, lambda_grm_pc = lambda_mean)])

ok <- sum_tab[lambda_mean <= lambda_cut]
n_pc_main <- if (nrow(ok)) min(ok$n_pc) else sum_tab[which.min(lambda_mean), n_pc]
cat(sprintf("주 분석 PC 개수: %d (규칙: 평균 λ ≤ %.2f인 최소 개수%s)\n", n_pc_main, lambda_cut,
            if (nrow(ok)) "" else " — 해당 없음, 평균 λ 최소 개수 사용"))
if (sum_tab[n_pc == n_pc_main, lambda_max] > 1.1) cat(sprintf("참고: PC %d개에서 λ 최대 %.3f (> 1.1)\n", n_pc_main, sum_tab[n_pc == n_pc_main, lambda_max]))
if (n_pc_main == 0) cat("PC 0개: 03-5에 PC0 공변량 파일이 없으므로 04-5가 PC 없이 연속형 공변량 5개만 담은 04-5_qcovar_PC0.txt를 만들어 사용\n")

p <- ggplot(melt(sum_tab, id.vars = "n_pc", measure.vars = c("lambda_pc_only", "lambda_mean"))[!is.na(value)],
            aes(x = n_pc, y = value, color = variable)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_hline(yintercept = c(1, lambda_cut), linetype = c("dashed", "dotted")) +
  scale_color_manual(values = c(lambda_pc_only = "grey50", lambda_mean = "firebrick"),
                     labels = c(lambda_pc_only = "PC only (linear regression, 03-4)", lambda_mean = "GRM + PC (GCTA MLMA)")) +
  labs(title = "Simulated ancestry-structured phenotype: lambda", x = "Number of PCs (PC-AiR)", y = "Mean lambda (20 sims)", color = NULL) +
  theme_qc + theme(legend.position = "bottom")
ggsave(paste0(pic_prefix, "_sim_lambda.png"), p, width = 8, height = 5.5, dpi = 300)

fwrite(sum_tab, paste0(out_prefix, "_sim_lambda_summary.csv"))
writeLines(as.character(n_pc_main), paste0(out_prefix, "_pc_choice.txt"))

