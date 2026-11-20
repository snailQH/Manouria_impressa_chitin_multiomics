#!/bin/bash
#SBATCH --job-name=k2_build
#SBATCH --cpus-per-task=60
#SBATCH --mem=400G
#SBATCH --time=72:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/kraken2_fungiext_build_%j.out
# Step 2: build the extended Kraken2 DB and the Bracken 150-mer distribution (480 GB node).
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/taxonomy
DB=$REFDIR/kraken2_pluspf_fungiext_2026-11
echo "=== build start $(date) $(hostname)"; kraken2 --version | head -1
kraken2-build --build --db $DB --threads 60
echo "=== kraken2 build done $(date)"; ls -la $DB/*.k2d
bracken-build -d $DB -t 60 -k 35 -l 150
echo "=== bracken done $(date)"
kraken2-inspect --db $DB --skip-counts | head -5
kraken2-inspect --db $DB --threads 16 | grep -iE "Pleurotus ostreatus|Pleurotus eryngii|Hypsizygus|Neocallimastigomycota|Piromyces" | head
echo "=== ALL DONE $(date)"
