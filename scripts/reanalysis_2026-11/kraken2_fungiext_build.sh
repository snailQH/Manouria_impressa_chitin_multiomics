#!/bin/bash
# Head node, background. Build an extended Kraken2 DB = existing full standard+protozoa+fungi
# libraries ($REFDIR/kraken2_db, RefSeq complete genomes) + ALL RefSeq fungal assemblies
# (any assembly level) + GenBank Neocallimastigomycota (anaerobic gut fungi) + the three dietary
# mushroom genomes. New sequences get kraken:taxid headers so no accession2taxid lookup is needed.
# Step 1 (this script): assemble library -> Step 2: kraken2_fungiext_build.sbatch on a 480 GB node.
set -u
DB=$REFDIR/kraken2_pluspf_fungiext_2026-11; OLD=$REFDIR/kraken2_db
R=$PROJECT/05.reanalysis_2026-11; LOG=$R/logs/kraken2_fungiext_prep.log
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/taxonomy
mkdir -p $DB/library $DB/taxonomy $DB/new_genomes; cd $DB
exec > >(tee -a $LOG) 2>&1; echo "=== prep start $(date)"
# reuse the old taxonomy and libraries via file symlinks (kraken2-build finds library/*/prelim_map.txt)
for f in $OLD/taxonomy/*; do ln -sf $f taxonomy/; done
for l in archaea bacteria viral human UniVec_Core protozoa fungi; do mkdir -p library/$l; for f in library.fna prelim_map.txt manifest.txt assembly_summary.txt; do [ -e $OLD/library/$l/$f ] && ln -sf $OLD/library/$l/$f library/$l/$f; done; done
# --- new fungal genomes
cd $DB/new_genomes
curl -sfL -o refseq_fungi_summary.txt https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt
curl -sfL -o genbank_fungi_summary.txt https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt
# RefSeq fungi not already in the old build (any assembly level, latest, full genome rep)
sed 's#.*/##; s#_[A-Za-z].*##' $OLD/library/fungi/manifest.txt | sed -E 's/^(GCF_[0-9]+\.[0-9]+).*/\1/' | sort -u > old_fungi_acc.txt
awk -F'\t' 'NR==FNR{a[$1];next} !/^#/ && !($1 in a) && $11=="latest" && $14=="Full"{print $1"\t"$6"\t"$8"\t"$20}' old_fungi_acc.txt refseq_fungi_summary.txt > refseq_new.tsv
# GenBank anaerobic gut fungi (Neocallimastigomycota genera)
GEN='Piromyces|Anaeromyces|Neocallimastix|Pecoramyces|Caecomyces|Orpinomyces|Feramyces|Aklioshbomyces|Khoyollomyces|Agriosomyces|Ghazallomyces|Joblinomyces|Liebetanzomyces|Oontomyces|Buwchfawromyces|Tahromyces|Capellomyces|Paucimyces|Aestipascuomyces|Cyllamyces'
awk -F'\t' -v g="^($GEN) " '!/^#/ && $8 ~ g && $11=="latest" && $14=="Full"{print $1"\t"$6"\t"$8"\t"$20}' genbank_fungi_summary.txt > gutfungi.tsv
# dietary mushrooms (GenBank) - taxids: P. eryngii 5323, H. marmoreus 39966, P. ostreatus 5322
printf 'GCA_001717165.1\t5323\tPleurotus eryngii 183\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/717/165/GCA_001717165.1_ASM171716v1\nGCA_013433165.1\t39966\tHypsizygus marmoreus HM62\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/433/165/GCA_013433165.1_ASM1343316v1\nGCA_014466165.1\t5322\tPleurotus ostreatus PC9\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/466/165/GCA_014466165.1_ASM1446616v1\n' > diet.tsv
cat refseq_new.tsv gutfungi.tsv diet.tsv | awk -F'\t' '!seen[$1]++' > to_add.tsv
echo "RefSeq fungi to add: $(wc -l < refseq_new.tsv); gut fungi: $(wc -l < gutfungi.tsv); diet: 3; total $(wc -l < to_add.tsv)"
mkdir -p fna; n=0
while IFS=$'\t' read acc taxid org url; do
  b=$(basename $url); out=fna/$acc.fna
  [ -s $out ] && continue
  curl -sfL "$url/${b}_genomic.fna.gz" | zcat | awk -v t=$taxid '/^>/{print ">"substr($1,2)"|kraken:taxid|"t; next}{print}' > $out.tmp && mv $out.tmp $out || { echo "FAILED $acc $org"; rm -f $out.tmp; }
  n=$((n+1)); [ $((n%50)) = 0 ] && echo "downloaded $n $(date)"
done < to_add.tsv
echo "downloaded files: $(ls fna | wc -l); failed: $(grep -c FAILED $LOG)"
cat fna/*.fna > fungi_ext.fna; echo "fungi_ext.fna $(grep -c '>' fungi_ext.fna) sequences, $(stat -c %s fungi_ext.fna) bytes"
cd $DB; kraken2-build --add-to-library new_genomes/fungi_ext.fna --db $DB --threads 16 --no-masking
echo "=== prep done $(date)"; ls library/
sbatch $R/kraken2_fungiext_build.sbatch
