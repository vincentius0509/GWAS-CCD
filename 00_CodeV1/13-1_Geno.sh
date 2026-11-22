#!/bin/bash

PLINK2_PATH="/data/00_tools/plink2_linux_avx2_20231123/plink2"
INPUT_BFILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6"
OUTPUT_BFILE="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/12_Geno_Beta"

# 5:25529184:C:T v
# 23:32288824:G:A v
# 24:12787767:A:G v
# 24:12787779:T:C v
# 27:44543118:A:T

$PLINK2_PATH \
    --bfile $INPUT_BFILE \
    --snp 27:44543118:A:T \
    --recode A \
    --dog \
    --out $OUTPUT_BFILE/snp_27_44543118_A_T_genotype
