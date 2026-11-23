# 10-1_chrx_qc.R - X-chromosome genotype QC and X-specific GRM. See README.md.
library(data.table)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
in_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")
test_mode <- Sys.getenv("TEST", "0") == "1"
tag <- if (test_mode) "_test" else ""
out_prefix <- file.path(qc_dir, paste0("10-1", tag)); work_dir <- paste0(out_prefix, "_work")
plink2 <- "/data/00_tools/plink2_linux_avx2_20231123/plink2"
gcta <- "/data/00_tools/gcta-1.95.4/gcta-1.95.4-linux-x86_64/gcta64"
n_cores <- as.integer(Sys.getenv("N_CORES", "32"))
raw_b <- file.path(base_dir, "GeneticData/v2023/DogAgingProject_GeneticData_CuratedRelease_2023_ed")
id_f <- file.path(in_dir, "04-1_grm", "grm_all.grm.id")
for (f in c(plink2, gcta, paste0(raw_b, c(".bed", ".bim", ".fam")), id_f)) if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
par_end <- 6650000; geno_max <- 0.05; maf_min <- 0.01; hwe_p <- 1e-6
p2 <- function(args, log) { system2(plink2, c("--dog", "--threads", n_cores, args), stdout = log, stderr = log); invisible(0) }

ids <- fread(id_f, header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
if (test_mode) ids <- ids[seq(1, .N, by = 8)]
kf <- file.path(work_dir, "keep.txt"); fwrite(ids, kf, sep = " ", col.names = FALSE)
rng <- if (test_mode) c("--from-bp", "50000000", "--to-bp", "60000000") else c("--from-bp", as.character(par_end + 1), "--to-bp", "200000000")

x0 <- file.path(work_dir, "x0")
if (!file.exists(paste0(x0, ".bed"))) p2(c("--bfile", raw_b, "--chr", "X", rng, "--keep", kf, "--make-bed", "--out", x0), paste0(x0, ".stdout"))
fam <- fread(paste0(x0, ".fam"), header = FALSE); n_m <- sum(fam$V5 == 1); n_f <- sum(fam$V5 == 2)
bim0 <- fread(paste0(x0, ".bim"), header = FALSE)
cat(sprintf("X (PAR 제외%s): SNP %s개, 개 %d마리 (수컷 %d · 암컷 %d)\n", if (test_mode) ", 시험 50 ~ 60 Mb" else "", format(nrow(bim0), big.mark = ","), nrow(fam), n_m, n_f))

for (s in c("males", "females")) {
  o <- file.path(work_dir, paste0("miss_", s))
  if (!file.exists(paste0(o, ".vmiss"))) p2(c("--bfile", x0, paste0("--keep-", s), "--missing", "variant-only", "--out", o), paste0(o, ".stdout"))
}
mm <- fread(file.path(work_dir, "miss_males.vmiss")); mf <- fread(file.path(work_dir, "miss_females.vmiss"))
ok_call <- intersect(mm[F_MISS <= geno_max, ID], mf[F_MISS <= geno_max, ID])
cf <- file.path(work_dir, "call_ok.txt"); writeLines(ok_call, cf)
x1 <- file.path(work_dir, "x1")
if (!file.exists(paste0(x1, ".bed"))) p2(c("--bfile", x0, "--extract", cf, "--maf", maf_min, "--make-bed", "--out", x1), paste0(x1, ".stdout"))
grp <- fread(file.path(in_dir, "03-1_group_hwe", "group_pheno.txt"), colClasses = c(FID = "character", IID = "character"))
fem <- fam[V5 == 2, as.character(V2)]
grp <- grp[GRP != "NONE" & IID %in% fem]
glist <- grp[, .N, by = GRP][N >= 20, GRP]
fails <- rbindlist(lapply(glist, function(g) {
  o <- file.path(work_dir, paste0("hwe_", g)); gk <- file.path(work_dir, paste0("keep_", g, ".txt"))
  if (!file.exists(paste0(o, ".hardy.x"))) { fwrite(grp[GRP == g, .(FID, IID)], gk, sep = " ", col.names = FALSE)
    p2(c("--bfile", x1, "--keep", gk, "--hardy", "--out", o), paste0(o, ".stdout")) }
  h <- fread(paste0(o, ".hardy.x")); data.table(ID = h$ID, fail = h$P < hwe_p)
}))
n_fail <- fails[, .(n = sum(fail, na.rm = TRUE)), by = ID]
oe <- rbindlist(lapply(glist, function(g) { h <- fread(file.path(work_dir, paste0("hwe_", g, ".hardy.x")))
  data.table(group = g, n_female = grp[GRP == g, .N], O_het = mean(h[["O(FEMALE_HET_A1)"]], na.rm = TRUE), E_het = mean(h[["E(FEMALE_HET_A1)"]], na.rm = TRUE),
             frac_p1e6 = mean(h$P < hwe_p, na.rm = TRUE)) }))
oe_ratio <- oe[, weighted.mean(O_het, n_female) / weighted.mean(E_het, n_female)]
x_quality_ok <- isTRUE(oe_ratio >= 0.5)
x1_ids <- fread(paste0(x1, ".bim"), header = FALSE)$V2
hw_ok <- if (x_quality_ok) setdiff(x1_ids, n_fail[n >= 2, ID]) else x1_ids
fwrite(oe, paste0(out_prefix, "_female_het.csv"))
fwrite(data.table(x_quality_ok = x_quality_ok, female_het_O_over_E = oe_ratio, hwe_filter_applied = x_quality_ok), paste0(out_prefix, "_quality.csv"))
hf <- file.path(work_dir, "hwe_ok.txt"); writeLines(hw_ok, hf)
xq <- paste0(out_prefix, "_chrX")
if (!file.exists(paste0(xq, ".bed"))) p2(c("--bfile", x1, "--extract", hf, "--make-bed", "--out", xq), paste0(xq, ".stdout"))
n1 <- nrow(fread(paste0(x1, ".bim"), header = FALSE)); nq <- nrow(fread(paste0(xq, ".bim"), header = FALSE))

pr <- paste0(out_prefix, "_chrX_prune")
if (!file.exists(paste0(pr, ".prune.in"))) p2(c("--bfile", xq, "--indep-pairwise", "500kb", "1", "0.2", "--out", pr), paste0(pr, ".stdout"))
gx <- paste0(out_prefix, "_grmX")
if (!file.exists(paste0(gx, ".grm.bin")))
  system2(gcta, c("--bfile", xq, "--extract", paste0(pr, ".prune.in"), "--autosome-num", "38", "--make-grm-xchr", "--out", gx, "--thread-num", n_cores),
          stdout = paste0(gx, ".stdout"), stderr = paste0(gx, ".stdout"))
if (!file.exists(paste0(gx, ".grm.bin"))) stop("X GRM 실패: ", gx, ".stdout")
n_pr <- length(readLines(paste0(pr, ".prune.in")))
gl <- readLines(paste0(gx, ".stdout")); n_gx <- sub(" valid SNPs.*", "", grep("valid SNPs are included", gl, value = TRUE))

counts <- data.table(step = c("X SNP (PAR 제외)", "수컷 · 암컷 call rate ≥ 0.95", "+ MAF ≥ 0.01", if (x_quality_ok) "+ 암컷 HWE W2 (QC 최종)" else "HWE 미적용 (X 품질 미달; QC 최종)", "LD pruning (GRM용)", "X GRM에 쓰인 SNP"),
                     n = c(nrow(bim0), length(ok_call), n1, nq, n_pr, as.integer(n_gx)))
fwrite(counts, paste0(out_prefix, "_qc_counts.csv"))
sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("10-1 X 염색체 QC (가설 H7 준비)\n"); if (test_mode) cat("*** TEST=1 시험 (일부 개 · 50 ~ 60 Mb). 해석하지 말 것 ***\n")
cat(sprintf("개 %d마리 (수컷 %d · 암컷 %d), PAR 제외 기준 %s bp. 암컷 HWE 집단 %d개 (암컷 20마리 이상)\n\n", nrow(fam), n_m, n_f,
            format(par_end, big.mark = ","), length(glist)))
print(counts)
cat(sprintf("\n■ X 유전형 품질: 집단 안 암컷 관측/기대 이형접합 비 %.3f → %s\n", oe_ratio,
            if (x_quality_ok) "정상 범위 (HWE W2 적용)" else "미달 (< 0.5): 원본 X에 암컷 이형접합이 거의 없음 → HWE 미적용, 10-2에서 H7은 '검정 불가 (자료 품질)'"))
print(oe[order(-n_female)][1:min(.N, 8)])
cat(sprintf("\n수컷 call rate 중앙값 %.3f · 암컷 %.3f (모든 X SNP; 수컷 이형접합은 결측으로 처리됨)\n", median(1 - mm$F_MISS), median(1 - mf$F_MISS)))
sink()
cat(sprintf("저장: %s_chrX.{bed,bim,fam} · _chrX_prune.prune.in · _grmX.grm.* · _qc_counts.csv · _female_het.csv · _quality.csv · _summary.txt\n", out_prefix))

