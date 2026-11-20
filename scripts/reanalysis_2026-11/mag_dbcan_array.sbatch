#!/bin/bash
#SBATCH --job-name=mag_dbcan
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --time=6:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/mag_dbcan_%A_%a.out
# dbCAN3 (prok mode) per dereplicated MAG, as an array (the sequential loop in mag_postbin.sbatch ran at ~4.7 min/MAG).
# Skips MAGs that already have overview.tsv.
set -uo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/dbcan
P=$PROJECT/05.reanalysis_2026-11/mags/postbin
g=$(ls $P/drep/dereplicated_genomes/*.fa | sed -n "${SLURM_ARRAY_TASK_ID}p"); b=$(basename $g .fa)
[ -s $P/dbcan/$b/overview.tsv ] && { echo "$b done"; exit 0; }
rm -rf $P/dbcan/$b
run_dbcan CAZyme_annotation --input_raw_data $g --mode prok --db_dir $REFDIR/dbcan --output_dir $P/dbcan/$b --threads 4 --methods diamond,hmm,dbCANsub > /dev/null 2>&1 || echo "dbcan failed $b"
echo "$b $(date)"
