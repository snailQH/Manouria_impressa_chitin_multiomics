#!/bin/bash
#SBATCH --job-name=gh18_tax
#SBATCH --cpus-per-task=16
#SBATCH --mem=130G
#SBATCH --time=6:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/gh18_tax_%j.out
# Taxonomy of the contigs carrying the 643 dbCAN3 GH18 genes of the unified catalogue (all 32 samples), with
# Kraken2 on the same database as the read-level taxonomy (RefSeq complete genomes; confidence 0.2).
# Output: gh18_contigs.kraken (per contig: taxid) + gh18_contigs.report; joined to gene counts by gh18_contig_taxonomy.py.
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/taxonomy
U=$PROJECT/03.unified_analysis; R=$PROJECT/05.reanalysis_2026-11/gh18_taxonomy; mkdir -p $R; cd $R
cat ../dbcan/out/part_*/overview.tsv | awk -F'\t' '$1!="Gene ID" && $6>=2 && $7 ~ /(^|[^0-9A-Za-z_])GH18([^0-9]|$)/ {print $1}' | sort -u > gh18_genes.txt
# gene id = <sample>_<contig>_<n>  -> sample, contig
awk '{g=$1; n=split(g,a,"_"); s=a[1]; c=a[2]"_"a[3]; print s"\t"c"\t"g}' gh18_genes.txt > gh18_gene_contig.tsv
cut -f1,2 gh18_gene_contig.tsv | sort -u > contigs.tsv; echo "GH18 genes: $(wc -l < gh18_genes.txt) on $(wc -l < contigs.tsv) contigs"
: > gh18_contigs.fa
for s in $(cut -f1 contigs.tsv | sort -u); do awk -v s=$s -F'\t' '$1==s{print $2}' contigs.tsv > ids_$s.txt
  conda run -p $CONDA_ROOT/envs/chitin_meta samtools faidx $U/03_assembly/$s/final.contigs.fa $(cat ids_$s.txt) | sed "s/^>/>${s}_/" >> gh18_contigs.fa; rm ids_$s.txt; done
grep -c ">" gh18_contigs.fa
kraken2 --db $PROJECT/kraken2_db_local --threads 16 --confidence 0.2 --minimum-hit-groups 2 --use-names --report gh18_contigs.report --output gh18_contigs.kraken gh18_contigs.fa
echo "classified: $(awk '$1=="C"' gh18_contigs.kraken | wc -l) / $(wc -l < gh18_contigs.kraken)"
echo DONE
