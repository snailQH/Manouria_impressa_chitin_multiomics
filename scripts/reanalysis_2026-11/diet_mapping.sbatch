#!/bin/bash
#SBATCH --job-name=diet_map
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/diet_map_%A_%a.out
# Dietary-DNA quantification: map each sample's non-host read pairs to the three mushroom genomes
# (P. eryngii 183, H. marmoreus HM62, P. ostreatus PC9; one combined Bowtie2 index, contig names
# prefixed with the genome tag). Count properly paired, MAPQ>=10 pairs per genome.
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/chitin_meta
U=$PROJECT/03.unified_analysis
R=$PROJECT/05.reanalysis_2026-11/diet_mapping; mkdir -p $R/idx $R/counts; T=16
# build the combined index once (task 1 builds; others wait)
if [ ! -s $R/idx/index.done ]; then
  if [ "$SLURM_ARRAY_TASK_ID" = "1" ]; then
    D=$PROJECT/05.reanalysis_2026-11/diet_genomes
    for g in Pleurotus_eryngii_183 Hypsizygus_marmoreus_HM62 Pleurotus_ostreatus_PC9; do
      zcat $D/$g.fna.gz | awk -v g=$g '/^>/{print ">"g"|"substr($1,2); next}{print}'
    done > $R/idx/diet3.fna
    bowtie2-build --threads $T $R/idx/diet3.fna $R/idx/diet3 > $R/idx/build.log 2>&1
    touch $R/idx/index.done
  else
    while [ ! -e $R/idx/index.done ]; do sleep 60; done
  fi
fi
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $U/samples_all32.txt)
[ -s $R/counts/$S.tsv ] && { echo "$S done"; exit 0; }
echo "=== $S $(date) $(hostname)"
bowtie2 -x $R/idx/diet3 -1 $U/02_host_removal/$S/${S}_nonhost_R1.fq.gz -2 $U/02_host_removal/$S/${S}_nonhost_R2.fq.gz \
  -p $T --very-sensitive --no-unal 2> $R/counts/$S.bowtie2.log \
  | samtools view -@ 2 -b -f 2 -q 10 - | samtools sort -@ 2 -m 1G -o $R/counts/$S.bam -
samtools index $R/counts/$S.bam
# pairs per genome: strict = properly paired MAPQ>=10 (the BAM); lenient = any aligned mate from the full log.
tot=$(grep -m1 "reads; of these" $R/counts/$S.bowtie2.log | awk '{print $1}')
conc=$(grep -m1 "aligned concordantly exactly 1 time" $R/counts/$S.bowtie2.log | awk '{print $1}')
concm=$(grep -m1 "aligned concordantly >1 times" $R/counts/$S.bowtie2.log | awk '{print $1}')
{ echo -e "Sample\tGenome\tpairs_q10\tnonhost_pairs\tpairs_q10_per_M\tconcordant_pairs_total\tconcordant_multi_total"
  samtools view $R/counts/$S.bam | awk -F'\t' -v s=$S -v tot=$tot -v c=$conc -v cm=$concm \
    'and($2,64){split($3,a,"|"); n[a[1]]++} END{for(g in n) printf "%s\t%s\t%d\t%d\t%.4f\t%d\t%d\n", s,g,n[g],tot,n[g]/tot*1e6,c,cm; if(length(n)==0) printf "%s\tnone\t0\t%d\t0\t%d\t%d\n", s,tot,c,cm}'
} > $R/counts/$S.tsv
cat $R/counts/$S.tsv; echo "=== done $(date)"
