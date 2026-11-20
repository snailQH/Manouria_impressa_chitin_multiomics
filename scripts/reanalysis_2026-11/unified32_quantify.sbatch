#!/bin/bash
#SBATCH --job-name=u32_quant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=$PROJECT/03.unified_analysis/slurm/u32_quant_%A_%a.out
#SBATCH --error=$PROJECT/03.unified_analysis/slurm/u32_quant_%A_%a.err
set -uo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/chitin_meta
BASE=$PROJECT/03.unified_analysis; OUT=$BASE/15_unified32; S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BASE/samples_all32.txt); T=20
[ -s $OUT/counts/${S}.counts.tsv ] && { echo "$S done already"; exit 0; }
echo "=== quantify $S $(date) node=$(hostname) ==="
bowtie2 -x $OUT/idx/nrgenes -1 $BASE/02_host_removal/$S/${S}_nonhost_R1.fq.gz -2 $BASE/02_host_removal/$S/${S}_nonhost_R2.fq.gz -p $T --no-unal --very-sensitive-local 2> $OUT/counts/${S}.bowtie2.log \
  | samtools view -@ 4 -bS - | samtools sort -@ 4 -m 2G -o $OUT/counts/${S}.bam -
samtools index -@ 4 $OUT/counts/${S}.bam
samtools idxstats $OUT/counts/${S}.bam | awk -F"\t" "\$1!=\"*\"{print \$1\"\t\"\$3}" > $OUT/counts/${S}.counts.tsv
grep "overall alignment rate" $OUT/counts/${S}.bowtie2.log; rm -f $OUT/counts/${S}.bam $OUT/counts/${S}.bam.bai; echo "=== $S done $(date) ==="
