#!/bin/bash
#SBATCH --job-name=k2_fungi
#SBATCH --cpus-per-task=48
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/k2_fungi_%j.out
# Kraken2 + Bracken of the 32 non-host read sets against the extended database (base libraries + 594 fungal genomes
# incl. Neocallimastigomycota and the three diet species), with the published parameters (confidence 0.2,
# minimum-hit-groups 2, Bracken read length 150). One job so the 136 GB hash is loaded once into the page cache.
set -uo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/taxonomy; export LC_ALL=C
U=$PROJECT/03.unified_analysis; DB=$REFDIR/kraken2_pluspf_fungiext_2026-11
O=$PROJECT/05.reanalysis_2026-11/fungal_profiling; mkdir -p $O/kraken $O/bracken
# guard: refuse to run against a database that lacks the added genomes (build 2 produced exactly that)
[ "$(wc -l < $DB/seqid2taxid.map)" -gt 500000 ] || { echo 'DB lacks the added genomes (seqid2taxid < 500k) — abort'; exit 1; }
echo "=== start $(date) $(hostname)"
for S in $(cat $U/samples_all32.txt); do
  [ -s $O/kraken/$S.report ] && continue
  kraken2 --db $DB --threads 48 --confidence 0.2 --minimum-hit-groups 2 --paired --gzip-compressed \
    --report $O/kraken/$S.report --output - $U/02_host_removal/$S/${S}_nonhost_R1.fq.gz $U/02_host_removal/$S/${S}_nonhost_R2.fq.gz > /dev/null 2> $O/kraken/$S.log
  echo "$S $(grep -m1 classified $O/kraken/$S.log | tr -s ' ') $(date +%H:%M)"
done
for S in $(cat $U/samples_all32.txt); do
  for L in S G; do bracken -d $DB -i $O/kraken/$S.report -o $O/bracken/$S.$L.bracken -w $O/bracken/$S.$L.report -r 150 -l $L -t 10 > /dev/null 2>&1; done
done
echo "=== fungal taxa in the reports (sum over samples):"
cat $O/kraken/*.report | awk -F'\t' '$6 ~ /Pleurotus eryngii|Pleurotus ostreatus|Hypsizygus marmoreus|Neocallimastigomycota|Piromyces|Anaeromyces|Neocallimastix|Fungi$/ {gsub(/^ +/,"",$6); c[$6]+=$2} END{for(k in c) print k"\t"c[k]}' | sort -k2,2nr
echo "=== DONE $(date)"
