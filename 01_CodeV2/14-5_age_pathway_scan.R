# 14-5_age_pathway_scan.R - Genome-wide comparison of gene-set effects between dogs < 8 and >= 8 years. See README.md.
suppressMessages(library(data.table))
d0 <- "/data/02_ByProject/DogAgingProject/GeneticData/analysis/04_GWAS_2023"
floors <- c("agefloor6", "agefloor7", "agefloor9", "age10", "agefloor11", "agefloor12")
f_in <- c(aux = file.path(d0, "06-5_aux.csv.gz"), info = file.path(d0, "06-4_aux_info.csv"),
          setNames(file.path(d0, "13-7_gsa", paste0(floors, "__aux.gsa.out")), floors))
for (f in f_in) if (!file.exists(f)) stop("파일 없음: ", f)
old <- c("age8", "agefloor9", "age10", "agefloor11", "agefloor12")
p_cut <- 0.05

rd <- function(f, m) { d <- fread(f, skip = "VARIABLE"); if ("FULL_NAME" %in% names(d)) d[, VARIABLE := FULL_NAME]; d[, .(set = VARIABLE, NGENES, BETA, SE, P, model = m)] }
A <- rbind(fread(cmd = paste("zcat", shQuote(f_in[["aux"]])))[variant == "aux" & kind == "model" & run %in% c("young8", "age8", "int_age"), .(set, NGENES, BETA, SE, P, model = run)],
           rbindlist(lapply(floors, function(m) rd(f_in[[m]], m))))
n_set <- A[model == "age8", uniqueN(set)]
W <- dcast(A, set ~ model, value.var = c("BETA", "SE", "P"))
W[, `:=`(diff_z = (BETA_age8 - BETA_young8) / sqrt(SE_age8^2 + SE_young8^2))][, diff_p := 2 * pnorm(-abs(diff_z))]
W[, q_ref := p.adjust(diff_p, "BH")]
W[, n_old_p05 := rowSums(.SD < p_cut, na.rm = TRUE), .SDcols = paste0("P_", old)]
W <- merge(W, fread(f_in[["info"]])[, .(set, name, n_genes)], by = "set", all.x = TRUE)
W[, side := fifelse(diff_p < p_cut, fifelse(diff_z > 0, "8세 이상 쪽이 큼", "8세 미만 쪽이 큼"), "차이 없음")]
out_cols <- c("set", "name", "n_genes", paste0("P_", c("young8", "agefloor6", "agefloor7", old, "int_age")), "BETA_young8", "BETA_age8", "diff_z", "diff_p", "q_ref", "side", "n_old_p05")
f3 <- function(x) { x <- copy(x); for (j in names(x)) if (is.numeric(x[[j]]) && !is.integer(x[[j]])) set(x, j = j, value = signif(x[[j]], 3)); x }
fwrite(f3(W[order(diff_p), ..out_cols]), file.path(d0, "14-5_age_pathway_scan.csv"))

fam <- c(GABA = "GABA|gamma-aminobutyric|GABAergic", Calcium = "calcium", "Metal ion" = "metal ion")
FM <- rbindlist(lapply(names(fam), function(k) W[grepl(fam[[k]], name, ignore.case = TRUE)][, family := k]))
sh <- function(d) d[, .(name = substr(name, 1, 55), n_genes, P_young8 = signif(P_young8, 2), P_age8 = signif(P_age8, 2), diff_z = round(diff_z, 2), diff_p = signif(diff_p, 2))]
sink(file.path(d0, "14-5_age_pathway_summary.txt"))
cat("14-5 게놈 전체 세트의 나이 비교 (탐색, 서술; 사전 기준 판정은 13-10 — 노령 > 8세 미만 경로 0개)\n")
cat(sprintf("세트 %d개. 8세 이상 vs 8세 미만 (겹치지 않음) 효과 차이 p < 0.05: 8세 이상 쪽이 큼 %d · 8세 미만 쪽이 큼 %d (참고: BH q < 0.05 %d개)\n",
            n_set, W[side == "8세 이상 쪽이 큼", .N], W[side == "8세 미만 쪽이 큼", .N], W[q_ref < 0.05, .N]))
cat("\n■ 8세 이상 쪽이 큰 세트 상위 10\n"); print(sh(W[diff_z > 0][order(diff_p)][1:10]), row.names = FALSE)
cat("\n■ 8세 미만 쪽이 큰 세트 상위 10\n"); print(sh(W[diff_z < 0][order(diff_p)][1:10]), row.names = FALSE)
cat(sprintf("\n■ 나이 하한 층 %d개 (8 ~ 12세 이상) 중 p < 0.05인 층 수 분포 (층끼리 겹쳐 독립 재현 아님)\n", length(old))); print(table(W$n_old_p05))
cat("\n  5개 층 모두 p < 0.05:\n"); print(W[n_old_p05 == length(old), .(name = substr(name, 1, 55), n_genes, P_young8 = signif(P_young8, 2), P_age8 = signif(P_age8, 2),
                                                                            P_ge12 = signif(P_agefloor12, 2), diff_p = signif(diff_p, 2))][order(P_age8)], row.names = FALSE)
cat("\n■ 관련 세트끼리 같은 방향인가 (이름 핵심어로 모음)\n")
for (k in names(fam)) { cat(sprintf("  [%s]\n", k)); print(sh(FM[family == k][order(diff_p)]), row.names = FALSE) }
sink()
cat(readLines(file.path(d0, "14-5_age_pathway_summary.txt")), sep = "\n")
