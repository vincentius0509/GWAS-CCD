#!/bin/bash

PLINK_PATH=/data/00_tools/plink_linux_x86_64_20231018/plink

HWE_THRESHOLDS=1e-6


declare -A INPUT_BFILES
# INPUT_BFILES[0.001]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_MAF_0.001_Filtered"
# INPUT_BFILES[0.005]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_MAF_0.005_Filtered"
INPUT_BFILES[0.01]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_MAF_0.01_Filtered"
INPUT_BFILES[0.05]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_MAF_0.05_Filtered"
INPUT_BFILES[0.1]="/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/02_MAF_0.1_Filtered"

for MAF in "${!INPUT_BFILES[@]}"; do
  INPUT_BFILE=${INPUT_BFILES[$MAF]}
  for HWE in "${HWE_THRESHOLDS[@]}"; do
    HWE_OUT=./analysis/03_GWAS_Pure/03_MAF_${MAF}_HWE_${HWE}

    $PLINK_PATH \
      --bfile "$INPUT_BFILE" \
      --dog \
      --hwe "$HWE" \
      --make-bed \
      --out "$HWE_OUT" \
      --threads 1
  done
  echo "HWE filtering done for MAF=${MAF}"
done
