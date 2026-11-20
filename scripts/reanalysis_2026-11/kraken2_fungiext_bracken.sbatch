#!/bin/bash
#SBATCH --job-name=k2_bracken
#SBATCH --cpus-per-task=60
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/kraken2_fungiext_bracken_%j.out
# Needs a node with ~500 GB RAM: the Kraken2 index is held in memory.
# Finish build 3 (a job failed only in Bracken): kraken2_fungiext_build2.sbatch removed the local build directory before
# bracken-build and had never copied seqid2taxid.map to NFS, so kmer2read_distr aborted. The kraken2 classification of the
# library (database.kraken, 102 GB) completed and is reused. Here: regenerate seqid2taxid.map from library/*/prelim_map.txt
# (all 570,994 entries are TAXID-type), assert the count, then run the remaining Bracken steps.
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/taxonomy; export LC_ALL=C
NFS=$REFDIR/kraken2_pluspf_fungiext_2026-11; cd $NFS
[ -s hash.k2d ] && [ -s database.kraken ] || { echo "hash.k2d / database.kraken missing"; exit 1; }
cat library/*/prelim_map.txt | awk -F'\t' '$1=="TAXID"{print $2"\t"$3}' > seqid2taxid.map
n=$(wc -l < seqid2taxid.map); echo "seqid2taxid.map lines: $n"; [ "$n" -eq 570994 ] || { echo "expected 570994"; exit 1; }
grep -c 'kraken:taxid|5323' seqid2taxid.map
echo "=== bracken start $(date) $(hostname)"
bracken-build -d $NFS -t 60 -k 35 -l 150
ls -la $NFS/database150mers.kmer_distrib $NFS/database150mers.kraken 2>/dev/null || true
echo "=== bracken done $(date)"
kraken2-inspect --db $NFS --threads 16 | grep -iE "Pleurotus ostreatus|Pleurotus eryngii|Hypsizygus|Neocallimastigomycota|Piromyces" | head
echo "=== ALL DONE $(date)"
