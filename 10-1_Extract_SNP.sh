PLINK2_PATH="/data/00_tools/plink2_linux_avx2_20231123/plink2"

# $PLINK2_PATH \
#   --bfile /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/03_MAF_0.01_HWE_1e-6 \
#   --extract /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/10_snplist.txt \
#   --export A \
#   --dog \
#   --out /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/10_snp_selected.txt

$PLINK2_PATH \
  --bfile /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/00_purebred_data \
  --extract /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/10_snplist.txt \
  --export A \
  --dog \
  --out /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Pure/10_snp_selected_total