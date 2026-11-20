#!/bin/bash
#SBATCH --job-name=mag_depth
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/mag_depth_%j.out
# per-sample depth of every (pre-dereplication) bin from the binning depth files (length-weighted mean contig depth)
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/binning
python3 - <<'PY'
import pandas as pd, glob, os
R='$PROJECT/05.reanalysis_2026-11/mags'; rows=[]
for d in sorted(glob.glob(f'{R}/*/depth.txt')):
    s=d.split('/')[-2]; dep=pd.read_csv(d,sep='\t').set_index('contigName')
    for b in glob.glob(f'{R}/{s}/bins/*.fa'):
        ids=[l[1:].split[0] for l in open(b) if l.startswith('>')]; sub=dep.loc[dep.index.intersection(ids)]
        if len(sub)==0: continue
        w=sub.contigLen; rows.append(dict(sample=s,bin=os.path.basename(b)[:-3],n_contigs=len(sub),size_bp=int(w.sum),mean_depth=float((sub.totalAvgDepth*w).sum/w.sum)))
pd.DataFrame(rows).to_csv(f'{R}/postbin/bin_depth_per_sample.tsv',sep='\t',index=False); print(len(rows),'bins with depth')
PY
