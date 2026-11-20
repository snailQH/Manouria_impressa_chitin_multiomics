#!/bin/bash
#SBATCH --job-name=u32_catalogue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=400G
#SBATCH --time=96:00:00
#SBATCH --output=$PROJECT/03.unified_analysis/slurm/u32_catalogue_%j.out
#SBATCH --error=$PROJECT/03.unified_analysis/slurm/u32_catalogue_%j.err
# ONE gene catalogue for all 32 samples (response to the 2026-11-10 review: per-cohort catalogues are separate reference spaces).
# Same tools/parameters as the per-cohort catalogues (cohort_catalogue.sbatch): CD-HIT 4.8.1 -c 0.95 -aS 0.9 on the Prodigal genes of
# all 32 samples; Bowtie2 --very-sensitive-local + samtools idxstats counts (array job unified32_quantify.sbatch, launched from here);
# eggNOG-mapper 2.1.12 (diamond, --dbmem) on the NR proteins.
set -uo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh
BASE=$PROJECT/03.unified_analysis; OUT=$BASE/15_unified32; LIST=$BASE/samples_all32.txt; THREADS=60
echo "###### UNIFIED-32 CATALOGUE START $(date) node=$(hostname) ######"
conda activate $CONDA_ROOT/envs/chitin_meta
SAMPLES=$(cat $LIST)
for S in $SAMPLES; do [ -s $BASE/04_gene_prediction/$S/${S}_proteins.faa ] || { echo "FATAL: $S has no gene prediction"; exit 1; }; done
if [ ! -s $OUT/catalog/all_proteins.faa ]; then
  # DS01 Prodigal headers carry no sample prefix (DS02/DS03 do) -> prefix every header with <sample>_ unless already present
  for S in $SAMPLES; do awk -v S="$S" '/^>/{ if (substr($0,2,length(S)+1)!=S"_") $0=">"S"_"substr($0,2) } {print}' $BASE/04_gene_prediction/$S/${S}_proteins.faa; done > $OUT/catalog/all_proteins.faa
  for S in $SAMPLES; do awk -v S="$S" '/^>/{ if (substr($0,2,length(S)+1)!=S"_") $0=">"S"_"substr($0,2) } {print}' $BASE/04_gene_prediction/$S/${S}_genes.fna;    done > $OUT/catalog/all_genes.fna
fi
DUP=$(grep "^>" $OUT/catalog/all_proteins.faa | awk "{print \$1}" | sort -S 50G --parallel=16 | uniq -d | wc -l)
echo "all proteins: $(grep -c "^>" $OUT/catalog/all_proteins.faa)  duplicate IDs: $DUP"; [ "$DUP" -ne 0 ] && { echo "FATAL: duplicate gene identifiers"; exit 1; }
[ -s $OUT/catalog/nr_proteins.faa ] || cd-hit     -i $OUT/catalog/all_proteins.faa -o $OUT/catalog/nr_proteins.faa -c 0.95 -aS 0.9 -M 0 -T $THREADS -d 0 > $OUT/catalog/cdhit_prot.log 2>&1
[ -s $OUT/catalog/nr_genes.fna ]    || cd-hit-est -i $OUT/catalog/all_genes.fna    -o $OUT/catalog/nr_genes.fna    -c 0.95 -aS 0.9 -M 0 -T $THREADS -d 0 > $OUT/catalog/cdhit_gene.log 2>&1
echo "NR proteins: $(grep -c "^>" $OUT/catalog/nr_proteins.faa)  NR genes: $(grep -c "^>" $OUT/catalog/nr_genes.fna)  $(date)"
if [ ! -s $OUT/idx/nrgenes.1.bt2l ] && [ ! -s $OUT/idx/nrgenes.1.bt2 ]; then
  bowtie2-build --threads $THREADS --large-index $OUT/catalog/nr_genes.fna $OUT/idx/nrgenes > $OUT/idx/build.log 2>&1; echo "index built $(date)"
fi
# quantification of the 32 samples as an array (3 nodes in parallel), independent of the annotation below
sbatch --array=1-32%6 $BASE/slurm/unified32_quantify.sbatch
# functional annotation of the NR proteins
conda activate $CONDA_ROOT/envs/eggnog
emapper.py --version 2>&1 | head -2
[ -s $OUT/kegg/unified32.emapper.annotations ] || emapper.py -i $OUT/catalog/nr_proteins.faa -o unified32 --output_dir $OUT/kegg --data_dir $REFDIR/eggnog --cpu $THREADS -m diamond --dbmem --override
echo "###### UNIFIED-32 CATALOGUE + eggNOG DONE $(date) ######"
