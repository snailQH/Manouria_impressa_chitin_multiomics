#!/bin/bash
#SBATCH --job-name=raw2mzml
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/raw2mzml_%A_%a.out
# Thermo .raw -> indexed mzML (ThermoRawFileParser 2.0, bioconda; DIA-NN 1.8.1 on Linux cannot read .raw)
set -uo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/proteomics
R=$PROJECT/05.reanalysis_2026-11; F=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $R/raw_list.txt); B=$(basename $F .raw)
[ -s $R/mzml/$B.mzML ] && { echo "$B exists"; exit 0; }
echo "=== $B $(date) $(hostname)"; ThermoRawFileParser -i $F -o $R/mzml -f 2 2>&1 | tail -3   # -f 2 = indexed mzML; -z = gzip? no: keep plain for DIA-NN
ls -la $R/mzml/$B.mzML* ; echo "=== done $(date)"
