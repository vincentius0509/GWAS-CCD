#!/bin/bash
# 03-4_PCA.sh - Genotype principal components: PC-AiR (primary) and two alternatives (all unrelated dogs; ancestry-proportion PCA). See README.md.
set -euo pipefail

PLINK2=/data/00_tools/plink2_linux_avx2_20231123/plink2
GENO_DIR=/data/02_ByProject/DogAgingProject/GeneticData
OUT_DIR=${GENO_DIR}/analysis/04_GWAS_2023

IN_PREFIX=${OUT_DIR}/03-3_maf0.05_hwe
ANC_FILE=${GENO_DIR}/v2023/DogAgingProject_GeneticData_CuratedRelease_2023.ancestry.tsv
WORK_DIR=${OUT_DIR}/03-4_pca
PRUNE_WIN=500kb
PRUNE_R2=0.2
KING_CUT=0.177
PCAIR_KIN=0.02209709
PCAIR_DIV=-0.02209709
N_PC=30
PC_SETS=(5 10 15 20 25 30)
N_SIM=20
STRUCT_VAR=0.3
SIM_SEED=1
THREADS=32
MEM_MB=200000

for f in ${IN_PREFIX}.pgen ${IN_PREFIX}.pvar ${IN_PREFIX}.psam ${ANC_FILE}; do
  if [ ! -f "$f" ]; then
    echo "파일 없음: $f"
    exit 1
  fi
done
mkdir -p ${WORK_DIR}/pcair ${WORK_DIR}/all ${WORK_DIR}/qmat

$PLINK2 --pfile $IN_PREFIX --dog \
        --indep-pairwise ${PRUNE_WIN} ${PRUNE_R2} \
        --threads $THREADS --memory $MEM_MB \
        --out ${WORK_DIR}/prune > /dev/null
echo "pruning 후 SNP 수: $(wc -l < ${WORK_DIR}/prune.prune.in)"

$PLINK2 --pfile $IN_PREFIX --dog \
        --extract ${WORK_DIR}/prune.prune.in \
        --make-king square bin \
        --threads $THREADS --memory $MEM_MB \
        --out ${WORK_DIR}/king > /dev/null

$PLINK2 --pfile $IN_PREFIX --dog \
        --extract ${WORK_DIR}/prune.prune.in \
        --king-cutoff ${KING_CUT} \
        --threads $THREADS --memory $MEM_MB \
        --out ${WORK_DIR}/king > /dev/null
awk 'NR > 1 {print $1 "\t" $2}' ${WORK_DIR}/king.king.cutoff.in.id > ${WORK_DIR}/all/ref_select.txt
echo "all 기준 집단: $(wc -l < ${WORK_DIR}/all/ref_select.txt) / $(( $(wc -l < ${IN_PREFIX}.psam) - 1 ))"

Rscript - "${WORK_DIR}" "${PCAIR_KIN}" "${PCAIR_DIV}" <<'EOF'
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]; kin_thr <- as.numeric(args[2]); div_thr <- as.numeric(args[3])
ids <- fread(file.path(wd, "king.king.id"))   # 열: #FID IID
n <- nrow(ids)
# plink2는 행 우선(row-major)으로 저장 → byrow = TRUE로 명시 (대칭 행렬이라 결과는 같지만 가정을 드러내 둠)
K <- matrix(readBin(file.path(wd, "king.king.bin"), what = "double", n = n * n, size = 8), n, n, byrow = TRUE)
stopifnot(length(K) == n * n, !anyNA(K), isSymmetric(K, tol = 1e-6))
diag(K) <- 0
rel <- K > kin_thr                                  # 근연 쌍
div <- (K < div_thr) & !rel                         # divergent 쌍 (근연이면 제외)
n_rel <- rowSums(rel); n_div <- rowSums(div)
kin_sum <- rowSums(K * rel)                         # 근연 쌍 kinship 합
set <- ifelse(n_rel == 0, "unrelated", NA_character_)
nr <- n_rel                                         # 남은(미배정) 근연 개체 수
active <- is.na(set)
# nr[x] = x의 "아직 미배정인" 근연 개체 수. 어떤 개가 unrelated로 갈 때는 nr = 0, 즉 미배정 근연 개체가 없을 때뿐이므로
# 그 개 때문에 다른 개의 nr을 줄일 필요가 없음 (related로 가는 개만 이웃의 nr을 줄이면 불변식 유지)
while (any(active)) {
  cand <- which(active)
  # ① 남은 근연 수 최대 → ② divergent 쌍 최소 → ③ kinship 합 최소 → ④ 첫 번째
  i <- cand[order(-nr[cand], n_div[cand], kin_sum[cand])[1]]
  set[i] <- "related"; active[i] <- FALSE
  nb <- which(rel[i, ] & active)                   # 아직 미배정인 그 개의 근연 개체
  nr[nb] <- nr[nb] - 1
  done <- nb[nr[nb] == 0]
  set[done] <- "unrelated"; active[done] <- FALSE
}
stopifnot(!anyNA(set))
u <- which(set == "unrelated")
stopifnot(!any(rel[u, u]))                          # 기준 집단 안에는 근연 쌍이 없어야 함
# 참고: 탐욕(greedy) 분할이라 최대 집합은 아님 — related인데 기준 집단 안에 근연 개체가 하나도 없는 개(기준 집단에 더 넣을 수 있는 개) 수
cat("related 중 기준 집단에 근연 개체가 없는 개(추가 가능, GENESIS 규칙 그대로 두고 기록만):", sum(set == "related" & rowSums(rel[, u, drop = FALSE]) == 0), "\n")
part <- data.table(FID = ids[[1]], IID = ids[[2]], set, n_rel, n_div, kin_sum = round(kin_sum, 4))
fwrite(part, file.path(wd, "pcair_partition.tsv"), sep = "\t")
fwrite(part[set == "unrelated", .(FID, IID)], file.path(wd, "pcair", "ref_select.txt"), sep = "\t", col.names = FALSE)
cat(sprintf("PC-AiR 분할: 근연 쌍 %d, divergent 쌍 %d / 근연 개체가 없는 개 %d → 기준 집단(unrelated) %d, 투영(related) %d\n",
            sum(rel) / 2, sum(div) / 2, sum(n_rel == 0), length(u), sum(set == "related")))
EOF

for LAB in pcair all; do
  D=${WORK_DIR}/${LAB}
  $PLINK2 --pfile $IN_PREFIX --dog \
          --keep ${D}/ref_select.txt \
          --extract ${WORK_DIR}/prune.prune.in \
          --freq \
          --pca allele-wts ${N_PC} \
          --threads $THREADS --memory $MEM_MB \
          --out ${D}/ref_pca > /dev/null
  mv ${D}/ref_pca.afreq ${D}/ref_freq.afreq
  ALLELE_FILE=${D}/ref_pca.eigenvec.allele
  A1_COL=$(head -1 ${ALLELE_FILE} | tr '\t' '\n' | grep -n -x "A1" | cut -d: -f1)
  PC1_COL=$(head -1 ${ALLELE_FILE} | tr '\t' '\n' | grep -n -x "PC1" | cut -d: -f1)
  if [ -z "${A1_COL}" ] || [ -z "${PC1_COL}" ]; then
    echo "열 탐색 실패 (${ALLELE_FILE}): A1='${A1_COL}', PC1='${PC1_COL}'"; exit 1
  fi
  PCN_COL=$(( PC1_COL + N_PC - 1 ))
  $PLINK2 --pfile $IN_PREFIX --dog \
          --extract ${WORK_DIR}/prune.prune.in \
          --read-freq ${D}/ref_freq.afreq \
          --score ${ALLELE_FILE} 2 ${A1_COL} header-read no-mean-imputation variance-standardize \
          --score-col-nums ${PC1_COL}-${PCN_COL} \
          --threads $THREADS --memory $MEM_MB \
          --out ${D}/proj > /dev/null
  awk -v n=${N_PC} 'BEGIN {FS = OFS = "\t"}
       NR == 1 {for (i = 1; i <= NF; i++) if ($i == "PC1_AVG") c = i
                if (!c || $(c + n - 1) != "PC" n "_AVG") {print "열 탐색 실패: PC1_AVG ~ PC" n "_AVG" > "/dev/stderr"; exit 1}
                printf "#FID\tIID"; for (j = 1; j <= n; j++) printf "\tPC%d", j; printf "\n"; next}
       {printf "%s\t%s", $1, $2; for (j = 0; j < n; j++) printf "\t%s", $(c + j); printf "\n"}' \
      ${D}/proj.sscore > ${D}/covar_pc.txt
  echo "${LAB}: 기준 집단 $(wc -l < ${D}/ref_select.txt)마리 PCA → 투영 $(( $(wc -l < ${D}/covar_pc.txt) - 1 ))마리"
done

Rscript - "${IN_PREFIX}.psam" "${ANC_FILE}" "${WORK_DIR}" "${N_PC}" "${N_SIM}" "${STRUCT_VAR}" "${SIM_SEED}" <<'EOF'
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
psam <- fread(args[1], select = 1:2, col.names = c("FID", "IID"))
anc <- fread(args[2], select = c("dog", "pop", "pct"))[, `:=`(dog = as.character(dog), pct = as.numeric(pct))]   # pct가 문자로 읽힐 수 있음
A <- dcast(anc[dog %in% as.character(psam$IID)], dog ~ pop, value.var = "pct", fill = 0)
A <- A[match(as.character(psam$IID), dog)]
stopifnot(!anyNA(A$dog))
M <- as.matrix(A[, -1])
n_pc <- as.integer(args[4])
q <- prcomp(M, center = TRUE, scale. = FALSE)
qpc <- data.table(`#FID` = psam$FID, IID = psam$IID, round(q$x[, 1:n_pc], 6))
setnames(qpc, 3:(n_pc + 2), paste0("PC", 1:n_pc))
fwrite(qpc, file.path(args[3], "qmat", "covar_pc.txt"), sep = "\t")
fwrite(data.table(PC = 1:ncol(M), var = q$sdev^2), file.path(args[3], "qmat", "eigenval.txt"), sep = "\t")
cat(sprintf("Q-matrix PCA: 참조 집단 %d개, PC1~10 설명 비율 %.1f%%, PC1~30 %.1f%%\n", ncol(M),
            100 * sum(q$sdev[1:10]^2) / sum(q$sdev^2), 100 * sum(q$sdev[1:n_pc]^2) / sum(q$sdev^2)))
n_sim <- as.integer(args[5]); sv <- as.numeric(args[6]); set.seed(as.integer(args[7]))
pheno <- data.table(`#FID` = psam$FID, IID = psam$IID)
for (s in seq_len(n_sim)) {
  g <- as.vector(M %*% rnorm(ncol(M)))
  g <- (g - mean(g)) / sd(g) * sqrt(sv)
  pheno[, (paste0("SIM", s)) := round(g + rnorm(nrow(M), 0, sqrt(1 - sv)), 5)]
}
fwrite(pheno, file.path(args[3], "sim_pheno.txt"), sep = "\t")
cat("가상 표현형:", n_sim, "개, 참조 집단", ncol(M), "개, 구조 효과 분산 비율", sv, "\n")
EOF

LAMBDA_FILE=${WORK_DIR}/sim_lambda.tsv
GLM_DIR=${WORK_DIR}/glm_tmp
mkdir -p ${GLM_DIR}
echo -e "label\tn_pc\tsim\tlambda" > ${LAMBDA_FILE}
run_glm() {
  local label=$1 k=$2
  local prefix=${GLM_DIR}/${label}_PC${k}
  if [ "${k}" -eq 0 ]; then
    $PLINK2 --pfile $IN_PREFIX --dog --extract ${WORK_DIR}/prune.prune.in \
            --pheno ${WORK_DIR}/sim_pheno.txt \
            --glm allow-no-covars hide-covar \
            --threads $THREADS --memory $MEM_MB --out ${prefix} > /dev/null
  else
    $PLINK2 --pfile $IN_PREFIX --dog --extract ${WORK_DIR}/prune.prune.in \
            --pheno ${WORK_DIR}/sim_pheno.txt \
            --covar ${WORK_DIR}/${label}/covar_pc.txt --covar-name PC1-PC${k} --covar-variance-standardize \
            --glm hide-covar \
            --threads $THREADS --memory $MEM_MB --out ${prefix} > /dev/null
  fi
  for s in $(seq 1 ${N_SIM}); do
    f=${prefix}.SIM${s}.glm.linear
    lam=$(awk 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "T_STAT") c = i; if (!c) {print "열 탐색 실패: T_STAT" > "/dev/stderr"; exit 1}; next} $c != "NA" {print $c * $c}' $f \
          | sort -g | awk '{a[NR] = $1} END {m = (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2; printf "%.4f", m / 0.4549364}')
    echo -e "${label}\t${k}\t${s}\t${lam}" >> ${LAMBDA_FILE}
    rm -f $f
  done
}
run_glm none 0
for LAB in pcair all qmat; do
  for K in "${PC_SETS[@]}"; do
    run_glm ${LAB} ${K}
  done
done
rm -rf ${GLM_DIR}
echo "λ 계산 완료: $(( $(wc -l < ${LAMBDA_FILE}) - 1 ))줄"   # (1 + 3 x 6) x 20 = 380
