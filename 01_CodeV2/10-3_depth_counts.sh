#!/usr/bin/env bash
# 10-3_depth_counts.sh - Aligns DAP raw reads (SRA PRJNA800779) and counts reads in 10 kb bins. See README.md.

set -euo pipefail
BASE=/data/02_ByProject/DogAgingProject/GeneticData/analysis
OUT=$BASE/04_GWAS_2023
IN=$BASE/04_GWAS_2023
REF=$IN/10-3_ref
REL=/data/02_ByProject/DogAgingProject/GeneticData/v2023/DogAgingProject_GeneticData_CuratedRelease_2023.tsv
GRMID=$IN/04-1_grm/grm_all.grm.id
FNA=$IN/06-11_ref/GCF_000002285.3_CanFam3.1_genomic.fna.gz
MM2=/data/00_tools/minimap2-2.28_x64-linux/minimap2
SAMT=/data/00_tools/samtools-1.19.2/samtools
SRATK=/data/00_tools/sratoolkit.3.4.1-ubuntu64/bin
YACC=NC_051844.1
N_JOBS=${N_JOBS:-24}; THREADS=${THREADS:-5}
TEST=${TEST:-0}; PILOT=${PILOT:-0}; MAX_SPOTS=${MAX_SPOTS:-2000000}; TOPUP=${TOPUP:-0}
TAG=""; [ "$TEST" = "1" ] && { TAG="_test"; MAX_SPOTS=500000; }
[ "$PILOT" != "0" ] && [ "$TEST" != "1" ] && TAG="_pilot"
SFX=""; [ "$TOPUP" = "1" ] && { SFX=".top"; TAG="${TAG}_topup"; }
PFX=$OUT/10-3$TAG
CNT=$OUT/10-3_counts
[ "$TEST" = "1" ] && CNT=$OUT/10-3_test_counts
TMPD=${PFX}_tmp
mkdir -p "$OUT" "$REF" "$CNT" "$TMPD"
for t in "$MM2" "$SAMT" "$SRATK/fastq-dump"; do [ -x "$t" ] || { echo "도구 없음: $t"; exit 1; }; done
for f in "$REL" "$GRMID" "$FNA"; do [ -s "$f" ] || { echo "파일 없음: $f"; exit 1; }; done
log() { echo "[$(date '+%F %T')] $*"; }

FA=$REF/canfam3_plusY.fa
if [ ! -s "$FA.ok" ]; then
  log "chrY 받기 ($YACC, NCBI efetch)"
  curl -sSL --fail --retry 5 -o "$REF/$YACC.fa" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$YACC&rettype=fasta&retmode=text"
  ylen=$(grep -v '^>' "$REF/$YACC.fa" | tr -d '\n' | wc -c)
  [ "$ylen" = "3937623" ] || { echo "chrY 길이 이상: $ylen (기대 3,937,623)"; exit 1; }
  log "합친 FASTA 만들기 (줄 길이 80으로 통일)"
  { gzip -dc "$FNA"; cat "$REF/$YACC.fa"; } | awk '/^>/{if(s!=""){for(i=1;i<=length(s);i+=80)print substr(s,i,80); s=""} print $1; next} {s=s $0; if(length(s)>=8000){n=int(length(s)/80)*80; for(i=1;i<=n;i+=80)print substr(s,i,80); s=substr(s,n+1)}} END{for(i=1;i<=length(s);i+=80)print substr(s,i,80)}' > "$FA"
  "$SAMT" faidx "$FA"
  echo ok > "$FA.ok"
fi
MAP=$REF/seq_map.tsv
[ -s "$MAP" ] || { gzip -dc "$FNA" | grep '^>' | sed 's/^>//' | awk '{acc=$1; if($0 ~ /mitochondrion/) print acc"\tMT"; else if(match($0,/chromosome [0-9X]+,/)){c=substr($0,RSTART+11,RLENGTH-12); if($0 !~ /unlocalized|unplaced/) print acc"\t"(c=="X"?"X":c)}}' > "$MAP"; echo -e "$YACC\tY" >> "$MAP"; }
[ "$(wc -l < "$MAP")" = "41" ] || { echo "염색체 이름표 줄 수가 41이 아님 ($MAP)"; exit 1; }
BINS=$REF/bins10k.tsv
if [ ! -s "$BINS" ]; then
  log "10 kb 구간 GC 표 계산"
  awk -v mapf="$MAP" 'BEGIN{while((getline l < mapf)>0){split(l,a,"\t"); m[a[1]]=a[2]}; OFS="\t"; print "seq","chr","bin","len","gc","n"}
    function flush(){ if(keep && len>0){ acgt=len-nn; print sq, m[sq], b, len, (acgt>0? gc/acgt : "NA"), nn/len } gc=0; nn=0; len=0 }
    /^>/{flush(); sq=substr($1,2); keep=(sq in m); b=0; next}
    keep{ s=$0; L=length(s); g=gsub(/[GCgc]/,"",s); s=$0; n=gsub(/[^ACGTacgt]/,"",s); gc+=g; nn+=n; len+=L; if(len>=10000){flush(); b++} }
    END{flush()}' "$FA" > "$BINS"
fi
IDX=$REF/canfam3_plusY.sr.mmi
if [ ! -s "$IDX" ]; then log "minimap2 색인 (-x sr; 수십 분)"; "$MM2" -x sr -t "$THREADS" -d "$IDX" "$FA" > "$REF/mmi.log" 2>&1; fi

SAMP=$OUT/10-3_samples.tsv
awk -F'\t' 'NR==FNR{k[$2]=1; next} FNR>1 && ($1 in k) && $6 ~ /^[SED]RR[0-9]+$/{print $1"\t"$6"\t"$7}' "$GRMID" "$REL" | sort -k1,1 > "$SAMP"
log "대상: 주 분석 개 중 SRR 있는 $(wc -l < "$SAMP")마리 (SRR 칸이 NA · 빈 칸인 개 제외; 2026-11-14 수정: 빈 칸 1마리 (76379) 때문에 xargs가 그 뒤 개 번호 · SRR을 한 칸씩 밀었음)"
RINFO=$OUT/10-3_runinfo.csv
[ -s "$RINFO" ] || curl -sSL --fail --retry 5 -o "$RINFO" "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA800779" || log "runinfo를 받지 못함 (10-4에서 기술 공변량 일부 없이 진행)"
LIST=${PFX}_todo.tsv
if [ "$TOPUP" = "1" ]; then
  awk -F'\t' '$3=="M"' "$SAMP" | { if [ "$TEST" = "1" ]; then awk 'c++ < 2'; else cat; fi; } | while IFS=$'\t' read -r d r x; do if [ -s "$CNT/$d.stats" ]; then printf '%s\t%s\t%s\n' "$d" "$r" "$x"; fi; done > "$LIST"
  log "수컷 추가 읽기 (조각 $((MAX_SPOTS + 1)) ~ $((2 * MAX_SPOTS))번째): 첫 실행이 끝난 수컷 $(wc -l < "$LIST")마리"
elif [ "$TEST" = "1" ]; then
  { awk -F'\t' '$3=="M" && c++ < 2' "$SAMP"; awk -F'\t' '$3=="F" && c++ < 2' "$SAMP"; } > "$LIST"
elif [ "$PILOT" != "0" ]; then
  h=$((PILOT / 2))
  { awk -F'\t' '$3=="M"' "$SAMP" | awk 'BEGIN{srand(1)}{print rand()"\t"$0}' | sort -k1,1g | awk -v h=$h 'NR <= h' | cut -f2-
    awk -F'\t' '$3=="F"' "$SAMP" | awk 'BEGIN{srand(1)}{print rand()"\t"$0}' | sort -k1,1g | awk -v h=$h 'NR <= h' | cut -f2-; } > "$LIST"
else
  cp "$SAMP" "$LIST"
fi
log "이번 실행 대상 $(wc -l < "$LIST")마리 (이미 끝난 개는 건너뜀: $(ls "$CNT" | grep -cE "^[0-9]+${SFX//./\\.}\.stats\$" || true)마리 완료)"

one_dog() {
  local dog=$1 r=$2 o t want got i ps spots rng
  [[ "$dog" =~ ^[0-9]+$ && "$r" =~ ^[SED]RR[0-9]+$ ]] || { echo "잘못된 입력 (개 번호 · SRR 어긋남): $dog $r"; return 1; }
  o=$CNT/$dog$SFX; [ -s "$o.stats" ] && return 0
  t=$TMPD/$dog; mkdir -p "$t"
  spots=$(awk -F, -v r="$r" '$1==r{print $4; exit}' "$RINFO" 2>/dev/null); [[ "$spots" =~ ^[0-9]+$ ]] || spots=0
  if [ "$TOPUP" = "1" ]; then
    if [ "$spots" -gt 0 ] && [ "$spots" -le "$MAX_SPOTS" ]; then
      printf 'total\t0\nmapped\t0\nq30\t0\nrun\t%s\nwant\t0\nnote\tno_extra_spots\n' "$r" > "$o.stats"; : | gzip -c > "$o.counts.gz"; rm -rf "${t:?}"; return 0
    fi
    want=$MAX_SPOTS; [ "$spots" -gt 0 ] && [ $((spots - MAX_SPOTS)) -lt "$MAX_SPOTS" ] && want=$((spots - MAX_SPOTS))
    rng=(-N $((MAX_SPOTS + 1)) -X $((2 * MAX_SPOTS)))
  else
    want=$spots
    if [ "$MAX_SPOTS" -gt 0 ] && { [ "$want" = 0 ] || [ "$want" -gt "$MAX_SPOTS" ]; }; then want=$MAX_SPOTS; fi
    rng=(); [ "$MAX_SPOTS" -gt 0 ] && rng=(-X "$MAX_SPOTS")
  fi
  for i in 1 2 3; do
    "$SRATK/fastq-dump" "${rng[@]}" --split-spot -Z "$r" 2> "$t/dump.log" \
      | "$MM2" -ax sr -t "$THREADS" "$IDX" /dev/stdin 2> "$t/mm2.log" \
      | awk -v so="$o.stats.tmp" 'BEGIN{OFS="\t"} /^@/{next}
          { fl=$2; if (and(fl,2304) || !and(fl,64)) next; tot++; if (and(fl,4)) next; map++;
            if ($5 >= 30 && !and(fl,1536)) { q++; c[$3 OFS int(($4-1)/10000)]++ } }
          END{ for (k in c) print k, c[k]; print "total",tot+0"\nmapped",map+0"\nq30",q+0 > so }' \
      | sort -k1,1 -k2,2n | gzip -c > "$o.counts.gz.tmp"
    ps="${PIPESTATUS[*]}"
    got=$(awk '$1=="total"{print $2}' "$o.stats.tmp" 2>/dev/null)
    if [ "$ps" = "0 0 0 0 0" ] && [ -n "$got" ] && { [ "$want" = 0 ] && [ "$got" -gt 0 ] || [ "$got" = "$want" ]; }; then
      mv "$o.counts.gz.tmp" "$o.counts.gz"; echo -e "run\t$r\nwant\t$want" >> "$o.stats.tmp"; mv "$o.stats.tmp" "$o.stats"; rm -rf "${t:?}"
      echo "[$(date '+%F %T')] 완료 $dog ($r): $(head -3 "$o.stats" | tr '\n' ' ')"; return 0
    fi
    echo "[$(date '+%F %T')] $dog ($r): 시도 $i 실패 (종료 코드 $ps, 조각 $got / 기대 $want) → 30초 뒤 다시"; rm -f "$o.counts.gz.tmp" "$o.stats.tmp"; sleep 30
  done
  echo "$dog: 3번 모두 실패 ($t/dump.log · mm2.log)"; return 1
}
[ -s "$LIST" ] || { log "이번 실행 대상이 없음 (TOPUP이면 첫 실행이 끝난 수컷이 있어야 함)"; exit 0; }
export -f one_dog; export CNT TMPD SRATK MM2 IDX THREADS MAX_SPOTS RINFO TOPUP SFX
set +e
awk -F'\t' '$1 ~ /^[0-9]+$/ && $2 ~ /^[SED]RR[0-9]+$/{print $1, $2}' "$LIST" | xargs -P "$N_JOBS" -n 2 bash -c 'one_dog "$0" "$1"'
st=$?
set -e
done_n=$(cut -f1 "$LIST" | while read -r d; do [ -s "$CNT/$d$SFX.stats" ] && echo "$d"; done | wc -l)
log "끝: 대상 $(wc -l < "$LIST")마리 중 $done_n마리 완료 (실패가 있으면 다시 실행하면 그 개만 다시 함)"
[ "$st" = "0" ] || log "일부 개 실패 (위 메시지)"
