#!/bin/bash
#SBATCH --job-name=dbcan_u32
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/dbcan_u32_%A_%a.out
# dbCAN3 (run_dbcan v5; HMMER dbCAN-HMM + DIAMOND CAZy + dbCAN-sub) on the unified 32-sample
# non-redundant protein catalogue, split into 8 chunks (array 1-8). STP.hmm is absent
# (403 at source) - only needed for CGC/STP annotation, which we do not run.
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/dbcan
R=$PROJECT/05.reanalysis_2026-11/dbcan
C=$PROJECT/03.unified_analysis/15_unified32/catalog/nr_proteins.faa
i=${SLURM_ARRAY_TASK_ID}; mkdir -p $R/chunks $R/out
if [ ! -s $R/chunks/split.done ]; then
  # one task does the split; others wait for the flag
  if [ "$i" = "1" ]; then
    seqkit split2 -p 8 -O $R/chunks -j 4 $C >/dev/null 2>&1 \
      || python3 - "$C" "$R/chunks" <<'PY'
import sys
src,out=sys.argv[1],sys.argv[2]; n=8
fh=[open(f"{out}/nr_proteins.part_{k:03d}.faa","w") for k in range(1,n+1)]
k=-1
for line in open(src):
    if line.startswith(">"): k=(k+1)%n
    fh[k].write(line)
for f in fh: f.close
PY
    touch $R/chunks/split.done
  else
    while [ ! -s $R/chunks/split.done ] && [ ! -e $R/chunks/split.done ]; do sleep 60; done
  fi
fi
F=$(ls $R/chunks/nr_proteins.part_*.faa | sed -n "${i}p")
echo "=== chunk $i: $F  $(grep -c '>' $F) proteins  $(date) $(hostname)"
run_dbcan CAZyme_annotation --input_raw_data $F --mode protein --db_dir $REFDIR/dbcan \
  --output_dir $R/out/part_$i --threads 16 --methods diamond,hmm,dbCANsub
echo "=== done $(date)"; ls -la $R/out/part_$i
