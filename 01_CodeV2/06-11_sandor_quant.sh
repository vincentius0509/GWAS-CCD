#!/usr/bin/env bash
# 06-11_sandor_quant.sh - Downloads and quantifies a public dog prefrontal-cortex ageing RNA-seq data set (SRA PRJNA815057). See README.md.

set -euo pipefail
BASE=/data/02_ByProject/DogAgingProject/GeneticData/analysis
OUT=$BASE/04_GWAS_2023
REFROOT=$BASE/04_GWAS_2023/06-11_ref
STAR=/data/00_tools/cellranger-7.2.0/lib/bin/STAR
N_CORES=${N_CORES:-32}
TEST=${TEST:-0}
REMOVE_FASTQ=${REMOVE_FASTQ:-1}
SRATK=/data/00_tools/sratoolkit.3.4.1-ubuntu64/bin
TAG=""; [ "$TEST" = "1" ] && TAG="_test"
PFX=$OUT/06-11$TAG
FQ=${PFX}_fastq
SRA=${PFX}_sra
ST=${PFX}_star
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.3_CanFam3.1
FNA=GCF_000002285.3_CanFam3.1_genomic.fna.gz
GTF=GCF_000002285.3_CanFam3.1_genomic.gtf.gz
mkdir -p "$OUT" "$REFROOT" "$FQ" "$ST" "$SRA"
[ -x "$STAR" ] || { echo "STAR 없음: $STAR"; exit 1; }
for t in fasterq-dump fastq-dump vdb-validate; do [ -x "$SRATK/$t" ] || { echo "sra-tools 없음: $SRATK/$t"; exit 1; }; done
log() { echo "[$(date '+%F %T')] $*"; }

cd "$REFROOT"
[ -s md5checksums.txt ] || curl -sSL --fail --retry 3 -o md5checksums.txt "$NCBI/md5checksums.txt"
for f in $FNA $GTF; do
  want=$(grep "/$f\$" md5checksums.txt | cut -d' ' -f1)
  if [ ! -s "$f" ] || [ "$(md5sum "$f" | cut -d' ' -f1)" != "$want" ]; then
    log "참조 내려받기: $f"; curl -sSL --fail --retry 3 -o "$f" "$NCBI/$f"
    [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want" ] || { echo "md5 불일치: $f"; exit 1; }
  fi
done
log "참조 md5 확인: $FNA · $GTF"

if [ "$TEST" = "1" ]; then
  IDX=$REFROOT/star_index_chr38_test
  if [ ! -s "$IDX/SA" ]; then
    mkdir -p "$IDX"
    gzip -dc $FNA | awk '/^>/{p=($1==">NC_006620.3")} p' > "$IDX/genome.fa"
    gzip -dc $GTF | awk -F'\t' '$1=="NC_006620.3"' > "$IDX/genes.gtf"
    log "STAR 색인 (시험, chr38)"
    "$STAR" --runMode genomeGenerate --runThreadN "$N_CORES" --genomeDir "$IDX" --genomeFastaFiles "$IDX/genome.fa" \
      --sjdbGTFfile "$IDX/genes.gtf" --sjdbOverhang 150 --genomeSAindexNbases 11 --outFileNamePrefix "$IDX/" > "$IDX/build.log" 2>&1
  fi
else
  IDX=$REFROOT/star_index
  if [ ! -s "$IDX/SA" ]; then
    mkdir -p "$IDX"
    [ -s "$REFROOT/genome.fa" ] || gzip -dc $FNA > "$REFROOT/genome.fa"
    [ -s "$REFROOT/genes.gtf" ] || gzip -dc $GTF > "$REFROOT/genes.gtf"
    log "STAR 색인 (전체 유전체; 약 1시간, 메모리 약 35 GB)"
    "$STAR" --runMode genomeGenerate --runThreadN "$N_CORES" --genomeDir "$IDX" --genomeFastaFiles "$REFROOT/genome.fa" \
      --sjdbGTFfile "$REFROOT/genes.gtf" --sjdbOverhang 150 --outFileNamePrefix "$IDX/" > "$IDX/build.log" 2>&1
    rm -f "${REFROOT:?}/genome.fa"
  fi
fi
[ -s "$IDX/SA" ] || { echo "STAR 색인 실패: $IDX/build.log 확인"; exit 1; }
log "STAR 색인: $IDX"

RUNS=${PFX}_runs.tsv
[ -s "$RUNS" ] || curl -sSL --fail --retry 3 -o "$RUNS" \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA815057&result=read_run&fields=run_accession,sample_alias,read_count,fastq_ftp,fastq_md5"
[ "$(tail -n +2 "$RUNS" | wc -l)" = "13" ] || { echo "ENA 표의 run 수가 13이 아님"; exit 1; }
if [ "$TEST" = "1" ]; then SAMPLES="CL_eto2 CL_eto7"; else SAMPLES=$(tail -n +2 "$RUNS" | cut -f2 | sort -V | tr '\n' ' '); fi
run_of() { awk -F'\t' -v s="$1" '$2==s{print $1}' "$RUNS"; }
nread_of() { awk -F'\t' -v s="$1" '$2==s{print $3}' "$RUNS"; }

get_sra() {
  local s=$1 r u f want have i
  r=$(awk -F'\t' -v s="$s" '$2==s{print $1}' "$RUNS"); [ -n "$r" ] || { echo "$s: run 없음"; return 1; }
  u=https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r; f=$SRA/$r.sra
  [ -s "$f.ok" ] && { echo "[$(date '+%F %T')] SRA 있음 (검증됨): $s ($r)"; return 0; }
  want=$(curl -sI --fail --retry 5 "$u" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}' | tail -1)
  [ -n "$want" ] || { echo "$s: SRA 크기를 알 수 없음 ($u)"; return 1; }
  for i in $(seq 1 30); do
    have=$(stat -c %s "$f" 2>/dev/null || echo 0)
    [ "$have" = "$want" ] && break
    [ "$have" -gt "$want" ] && rm -f "${f:?}"
    curl -sS --fail --retry 5 --retry-all-errors --connect-timeout 60 --speed-limit 100000 --speed-time 120 -C - -o "$f" "$u" || sleep 20
  done
  [ "$(stat -c %s "$f" 2>/dev/null || echo 0)" = "$want" ] || { echo "$s: SRA 내려받기 실패 (30번 시도)"; return 1; }
  if "$SRATK/vdb-validate" "$f" > "$f.validate.log" 2>&1 && ! grep -qi "err" "$f.validate.log"; then
    echo "$want" > "$f.ok"; echo "[$(date '+%F %T')] SRA 준비 (크기 $want, vdb-validate 통과): $s ($r)"
  else
    echo "$s: vdb-validate 실패 ($f.validate.log) → 파일을 지움. 다시 실행하면 새로 받음"; rm -f "${f:?}"; return 1
  fi
}
export -f get_sra; export RUNS SRA SRATK
if [ "$TEST" = "1" ]; then
  for s in $SAMPLES; do
    r=$(run_of "$s")
    [ -s "$FQ/${s}_R1.fastq" ] || { "$SRATK/fastq-dump" -X 100000 --split-files -O "$FQ" "$r" > /dev/null 2>&1; mv "$FQ/${r}_1.fastq" "$FQ/${s}_R1.fastq"; mv "$FQ/${r}_2.fastq" "$FQ/${s}_R2.fastq"; }
    nl=$(wc -l < "$FQ/${s}_R1.fastq"); [ "$nl" = "400000" ] || { echo "시험 FASTQ 이상 ($s: $nl줄)"; exit 1; }
  done
else
  log "SRA 내려받기 (동시 4개): $SAMPLES"
  if ! echo "$SAMPLES" | tr ' ' '\n' | grep -v '^$' | xargs -P 4 -I{} bash -c 'get_sra {}'; then
    echo "일부 표본의 SRA를 받지 못함 (위 메시지 확인). 다시 실행하면 받다 만 파일부터 이어받음"; exit 1
  fi
fi

to_fastq() {
  local s=$1 r n1
  [ "$TEST" = "1" ] && return 0
  r=$(run_of "$s")
  if [ -s "$FQ/${s}_R1.fastq" ] && [ -s "$FQ/${s}_R2.fastq" ]; then return 0; fi
  log "FASTQ 풀기: $s ($r)"
  rm -rf "${FQ:?}/tmp_$s"; mkdir -p "$FQ/tmp_$s"
  "$SRATK/fasterq-dump" --split-files --threads 8 --temp "$FQ/tmp_$s" -O "$FQ" "$SRA/$r.sra" > "$FQ/tmp_$s.log" 2>&1
  rm -rf "${FQ:?}/tmp_$s"
  mv "$FQ/${r}_1.fastq" "$FQ/${s}_R1.fastq"; mv "$FQ/${r}_2.fastq" "$FQ/${s}_R2.fastq"
  n1=$(( $(wc -l < "$FQ/${s}_R1.fastq") / 4 ))
  [ "$n1" = "$(nread_of "$s")" ] || { echo "$s: FASTQ read 수 $n1 ≠ ENA 표 $(nread_of "$s")"; exit 1; }
}

for s in $SAMPLES; do
  o=$ST/$s
  if [ -s "$o/ReadsPerGene.out.tab" ]; then log "정렬 건너뜀 (있음): $s"; continue; fi
  to_fastq "$s"
  mkdir -p "$o"; log "STAR 정렬: $s"
  "$STAR" --runThreadN "$N_CORES" --genomeDir "$IDX" --readFilesIn "$FQ/${s}_R1.fastq" "$FQ/${s}_R2.fastq" \
    --quantMode GeneCounts --outSAMtype None --outFileNamePrefix "$o/" > "$o/star.stdout" 2>&1
  [ -s "$o/ReadsPerGene.out.tab" ] || { echo "STAR 실패: $s ($o/Log.out 확인)"; exit 1; }
  if [ "$REMOVE_FASTQ" = "1" ] && [ "$TEST" != "1" ]; then rm -f "${FQ:?}/${s}"_R[12].fastq; fi
done

SUM=${PFX}_summary.tsv
printf "sample\tinput_reads\tunique_pct\tmulti_pct\tgene_unstranded\tgene_forward\tgene_reverse\n" > "$SUM"
for s in $SAMPLES; do
  o=$ST/$s
  inp=$(grep "Number of input reads" "$o/Log.final.out" | cut -f2)
  uni=$(grep "Uniquely mapped reads %" "$o/Log.final.out" | cut -f2)
  mul=$(grep "% of reads mapped to multiple loci" "$o/Log.final.out" | cut -f2)
  read -r c2 c3 c4 < <(awk 'NR>4{a+=$2;b+=$3;c+=$4} END{print a,b,c}' "$o/ReadsPerGene.out.tab")
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$inp" "$uni" "$mul" "$c2" "$c3" "$c4" >> "$SUM"
done
column -t "$SUM"
log "끝: $SUM (가닥: gene_reverse가 gene_forward보다 훨씬 크면 역가닥 라이브러리 → 06-12가 자동 판단)"

