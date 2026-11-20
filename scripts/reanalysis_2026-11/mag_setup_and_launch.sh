#!/bin/bash
# Build the `binning` conda env, launch the 32-sample binning array,
# then fetch the CheckM2 and GTDB-Tk (r226) reference databases into $REFDIR.
set -u
R=$PROJECT/05.reanalysis_2026-11
source $CONDA_ROOT/etc/profile.d/conda.sh
if [ ! -x $CONDA_ROOT/envs/binning/bin/metabat2 ]; then
  conda create -y -q -p $CONDA_ROOT/envs/binning -c conda-forge -c bioconda \
    metabat2 checkm2 gtdbtk=2.4.0 drep > $R/logs/conda_binning.log 2>&1
fi
for t in metabat2 jgi_summarize_bam_contig_depths checkm2 gtdbtk dRep; do
  [ -x $CONDA_ROOT/envs/binning/bin/$t ] && echo "ok $t" || echo "MISSING $t"
done
sbatch --array=1-32%6 $R/mags/mag_bin_sample.sbatch
conda activate $CONDA_ROOT/envs/binning
mkdir -p $REFDIR/checkm2 $REFDIR/gtdbtk
if ! ls $REFDIR/checkm2/CheckM2_database/*.dmnd >/dev/null 2>&1; then
  checkm2 database --download --path $REFDIR/checkm2 > $R/logs/checkm2_db.log 2>&1
fi
echo "checkm2 db: $(ls $REFDIR/checkm2/CheckM2_database/ 2>/dev/null)"
cd $REFDIR/gtdbtk
if [ ! -d release226 ]; then
  URL=https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
  curl -sfL -o gtdbtk_r226_data.tar.gz "$URL" > $R/logs/gtdbtk_db.log 2>&1 \
    && tar xzf gtdbtk_r226_data.tar.gz >> $R/logs/gtdbtk_db.log 2>&1 && rm -f gtdbtk_r226_data.tar.gz
fi
echo "gtdbtk db: $(ls $REFDIR/gtdbtk/)"
echo "SETUP DONE $(date)"
