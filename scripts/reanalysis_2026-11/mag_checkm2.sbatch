#!/bin/bash
#SBATCH --job-name=mag_checkm2
#SBATCH --cpus-per-task=36
#SBATCH --mem=70G
#SBATCH --time=48:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/mag_checkm2_%j.out
# Needs a CPU with AVX-512; CheckM2's TensorFlow model aborts without it.
# CheckM2 on all bins. Runs on the analysis server (Xeon Gold, AVX-512): CheckM2's TensorFlow model dies with
# 'Illegal instruction' on the Sandy Bridge compute nodes (the high-memory compute nodes). Writes postbin/bin_quality.tsv.
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/binning
R=$PROJECT/05.reanalysis_2026-11/mags; P=$R/postbin; mkdir -p $P/all_bins
n=$(ls $R/*/bins.done | wc -l); [ "$n" = 32 ] || { echo "only $n/32 samples binned"; exit 1; }
for f in $R/*/bins/*.fa; do ln -sf $f $P/all_bins/; done
echo "=== $(ls $P/all_bins | wc -l) bins  $(date) $(hostname)"
export CHECKM2DB=$REFDIR/checkm2/CheckM2_database/uniref100.KO.1.dmnd
[ -s $P/checkm2/quality_report.tsv ] || checkm2 predict --threads 36 --input $P/all_bins -x fa --output-directory $P/checkm2 --database_path $CHECKM2DB
python3 - <<'PY'
import pandas as pd
P='$PROJECT/05.reanalysis_2026-11/mags/postbin'
q=pd.read_csv(f'{P}/checkm2/quality_report.tsv',sep='\t')
q['quality']=['high' if c>=90 and x<=5 else 'medium' if c>=50 and x<=10 else 'low' for c,x in zip(q.Completeness,q.Contamination)]
q.to_csv(f'{P}/bin_quality.tsv',sep='\t',index=False); print(q.quality.value_counts.to_dict)
keep=q[q.quality!='low']
open(f'{P}/mq_bins.txt','w').write('\n'.join(f'{P}/all_bins/{n}.fa' for n in keep.Name)+'\n')
pd.DataFrame({'genome':keep.Name+'.fa','completeness':keep.Completeness,'contamination':keep.Contamination}).to_csv(f'{P}/genomeInfo.csv',index=False)
PY
echo "=== checkm2 done $(date)"
