#!/bin/bash

PLINK2_PATH="/data/00_tools/plink2_linux_avx2_20231123/plink2"

INPUT_BFILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/v2022/DogAgingProject_GeneticData_CuratedRelease_2022"
OUTPUT_BFILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/12_Geno_Beta"

# GUCY1A2
# 5:25529170:A:T v

# CACNA1C
# 27:44801518:G:A v
# 27:44847773:C:T v
# 27:44860521:T:C v
# 27:44790338:C:T v
# 27:45065570:A:G v
# 27:44235256:G:C v
# 27:44941531:C:T v
# 27:44520608:T:G v
# 27:44929991:C:T v

# PLCB4
# 24:12787773:G:A v
# 24:12787779:T:C v
# 24:12787767:A:G => GWAS

# PLCB4
# 24:12787779:T:C v
# 24:12787773:G:A v
# 24:12787779:T:A => GWAS
# 24:12785743:C:T v
# 24:12785679:G:A v
# 24:12782819:C:T v
# 24:12782528:G:A v
# 24:12780920:G:A v
# 24:12780924:T:C v
# 24:12752625:A:T v
# 24:12724858:A:T v
# 24:12724860:A:T v
# 24:12818735:A:G v
# 24:12818736:T:A v
# 24:12753225:G:A v
# 24:12754140:A:G v

# lncRNA
# 23:32276785:G:A v
# 23:32366169:C:T v

# ---------- 1) 특정 SNP의 유전형 추출 (A-format) ----------
$PLINK2_PATH \
    --bfile $INPUT_BFILE \
    --snp 24:12787779:T:C \
    --recode A \
    --dog \
    --out $OUTPUT_BFILE/snp_24_12787779_genotype

# phenotype + covariate 파일 설정
PHENO_FILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_pheno001.txt"
COVAR_FILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/04_cov_001_1e6_10.txt"

# ---------- 2) 주변 영역(Chr5:25.52–25.54 Mb)에서 β 계산 ----------
$PLINK2_PATH \
    --bfile $INPUT_BFILE \
    --chr 24 \
    --from-bp 12737779 \
    --to-bp 12837779 \
    --glm cols=+a1freq hide-covar \
    --pheno $PHENO_FILE \
    --pheno-name phe1 \
    --covar $COVAR_FILE \
    --covar-variance-standardize \
    --dog \
    --out $OUTPUT_BFILE/snp_24_12787779_beta
