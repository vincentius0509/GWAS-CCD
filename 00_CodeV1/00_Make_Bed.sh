#!/bin/bash

# 1st
PLINK_PATH=/data/00_tools/plink_linux_x86_64_20231018/plink
INPUT_BFILE=/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/DogAgingProject_GeneticData_CuratedRelease_2022
OUTPUT_BFILE=/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/00_purebred_data

$PLINK_PATH --bfile $INPUT_BFILE \
            --keep /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/00_purebred_id.txt \
            --make-bed \
            --dog \
            --out $OUTPUT_BFILE
