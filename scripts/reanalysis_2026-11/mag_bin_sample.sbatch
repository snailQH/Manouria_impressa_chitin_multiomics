#!/bin/bash
#SBATCH --job-name=mag_bin
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/mag_bin_%A_%a.out
# Per-sample MAG binning (array 1-32 over samples_all32.txt):
#   Bowtie2 map the sample's non-host reads to its own MEGAHIT assembly ->
#   jgi_summarize_bam_contig_depths -> MetaBAT2 2.15 (min contig 1,500 bp, seed 1).
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh
U=$PROJECT/03.unified_analysis
R=$PROJECT/05.reanalysis_2026-11/mags
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $U/samples_all32.txt); O=$R/$S; mkdir -p $O; T=20
[ -s $O/bins.done ] && { echo "$S already done"; exit 0; }
echo "=== $S $(date) $(hostname)"
A=$U/03_assembly/$S/final.contigs.fa
R1=$U/02_host_removal/$S/${S}_nonhost_R1.fq.gz; R2=$U/02_host_removal/$S/${S}_nonhost_R2.fq.gz
conda activate $CONDA_ROOT/envs/chitin_meta
if [ ! -s $O/$S.bam ]; then
  bowtie2-build --threads $T $A $O/asm > $O/build.log 2>&1
  bowtie2 -x $O/asm -1 $R1 -2 $R2 -p $T --no-unal 2> $O/bowtie2.log \
    | samtools view -@ 4 -bS - | samtools sort -@ 4 -m 2G -o $O/$S.bam -
  samtools index $O/$S.bam
fi
grep overall $O/bowtie2.log
# MetaBAT2 2.15 build: the 2.18 bioconda build (env binning) needs AVX2 and dies with 'Illegal instruction'
# on the Sandy Bridge compute nodes (E5-4650). CheckM2/dRep/GTDB-Tk stay in env binning.
conda activate $CONDA_ROOT/envs/metabat215
jgi_summarize_bam_contig_depths --outputDepth $O/depth.txt $O/$S.bam
metabat2 -i $A -a $O/depth.txt -o $O/bins/$S.bin -m 1500 -t $T --seed 1 > $O/metabat2.log 2>&1
n=$(ls $O/bins/*.fa 2>/dev/null | wc -l); echo "$S bins: $n"
rm -f $O/asm.*.bt2*
touch $O/bins.done; echo "=== done $(date)"
