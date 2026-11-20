#!/bin/bash
#SBATCH --job-name=mag_post
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/mag_post_%j.out
# After per-sample binning (mag_bin_sample.sbatch, 32/32):
#  1. (mag_checkm2.sbatch, on the analysis server because of AVX) CheckM2 -> keep MIMAG medium+ (completeness >= 50 %, contamination <= 10 %)
#  2. dRep dereplication at 95 % ANI (species-level MAG set)
#  3. GTDB-Tk R226 classify_wf on the dereplicated set
#  4. dbCAN3 (prok mode) on every dereplicated MAG -> GH18/GH19/GH20 per MAG
#  5. per-sample MAG abundance: reads already mapped to each sample's assembly -> bin coverage (depth.txt from binning)
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh
R=$PROJECT/05.reanalysis_2026-11/mags; P=$R/postbin; mkdir -p $P/all_bins; T=32
[ -s $P/bin_quality.tsv ] || { echo "run mag_checkm2.sbatch (the analysis server) first"; exit 1; }
conda activate $CONDA_ROOT/envs/binning
echo "=== $(wc -l < $P/mq_bins.txt) medium+ bins  $(date) $(hostname)"
[ -d $P/drep/dereplicated_genomes ] || dRep dereplicate $P/drep -g $P/mq_bins.txt --genomeInfo $P/genomeInfo.csv -sa 0.95 -nc 0.30 -p $T -comp 50 -con 10
echo "=== dRep done: $(ls $P/drep/dereplicated_genomes | wc -l) species-level MAGs  $(date)"
export GTDBTK_DATA_PATH=$REFDIR/gtdbtk/release226
[ -s $P/gtdbtk/gtdbtk.bac120.summary.tsv ] || gtdbtk classify_wf --genome_dir $P/drep/dereplicated_genomes -x fa --out_dir $P/gtdbtk --cpus $T --pplacer_cpus 8 --skip_ani_screen
echo "=== gtdbtk done $(date)"
conda activate $CONDA_ROOT/envs/dbcan
mkdir -p $P/dbcan
for g in $P/drep/dereplicated_genomes/*.fa; do b=$(basename $g .fa); [ -s $P/dbcan/$b/overview.tsv ] || run_dbcan CAZyme_annotation --input_raw_data $g --mode prok --db_dir $REFDIR/dbcan --output_dir $P/dbcan/$b --threads $T --methods diamond,hmm,dbCANsub > /dev/null 2>&1 || echo "dbcan failed $b"; done
echo "=== dbcan done $(date)"
# per-sample abundance of every (pre-dereplication) bin from the binning depth files: mean contig depth weighted by length
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
echo "=== ALL DONE $(date)"
