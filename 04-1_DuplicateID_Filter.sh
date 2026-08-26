#!/bin/bash

# INPUT FILE
# /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Group/04_pheno01.txt
# /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Group/04_pheno005.txt
# /home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Group/04_pheno001.txt

# OUTPUT FILE: INPUT FILE 과 같은 경로, 파일 명은 원래 파일명에 _filtered 추가

# 입력 파일 목록
PHENO_FILES=(
"/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Group/04_pheno001.txt"
"/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Group/04_pheno005.txt"
"/home/yhSong/data/02_ByProject/DogAgingProject/GeneticData/analysis/03_GWAS_Group/04_pheno01.txt"
)

for FILE in "${PHENO_FILES[@]}"; do
  echo "▶ 처리 중: $FILE"

  OUTPUT_FILE="${FILE%.txt}_filtered.txt"

  awk '
  BEGIN {
    FS = OFS = "\t"
  }

  NR == 1 {
    print $0
    next
  }

  {
    key = $2  # IID 기준
    phe1_values[key] = $3  # phe1은 첫 번째 값을 사용 (혹은 평균 계산 원할 경우 배열로 쌓기)
    if ($4 == 2) phe2_max[key] = 2
    else if (!(key in phe2_max)) phe2_max[key] = 1
  }

  END {
    for (id in phe1_values) {
      print id, id, phe1_values[id], phe2_max[id]
    }
  }
  ' "$FILE" | sort -k2,2n > "$OUTPUT_FILE"

  echo "✔ 저장 완료: $OUTPUT_FILE"
done
