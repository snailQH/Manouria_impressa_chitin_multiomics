#!/bin/bash
#SBATCH --job-name=cmp_gh18
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/cmp_gh18_%j.out
# Comparative GH18 locus analysis across Testudines.
#  For Mimp_v3 (control) and 8 RefSeq turtle genomes:
#   (a) hmmsearch of the dbCAN GH18 model against the annotated proteome  -> annotation-based GH18 gene list
#   (b) miniprot of the 12 M. impressa GH18-locus gene models + 6 flanking genes against the genome
#       -> annotation-independent GH18 loci and the orthologous interval (delimited by the flanking genes)
#  comparative_gh18_parse.py then counts loci per genome and inside the interval.
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh
R=$PROJECT/05.reanalysis_2026-11/comparative; REL=$PROJECT/01.Genome_Assembly/new-version/06_release
W=$R/gh18; mkdir -p $W/q $W/out; T=16
cd $W
# queries
conda activate $CONDA_ROOT/envs/binning
python3 - <<'PY'
ids=[f"MIMP3_g{n:05d}" for n in list(range(8465,8468))+list(range(8468,8480))+list(range(8480,8483))]
seq={}; cur=None
for l in open("$PROJECT/01.Genome_Assembly/new-version/06_release/Mimp_v3.proteins.faa"):
    if l.startswith(">"): cur=l[1:].split[0]; seq[cur]=[]
    else: seq[cur].append(l.strip)
out=open("q/mimp_locus_and_flanks.faa","w"); n=0
for k,v in seq.items:
    g=k.split(".")[0]
    if g in ids and k.endswith(".t1"): out.write(f">{g}\n"+"".join(v)+"\n"); n+=1
out.close; print("queries",n)
PY
hmmfetch $REFDIR/dbcan/dbCAN.hmm GH18.hmm > q/GH18.hmm
hmmstat q/GH18.hmm | tail -1
# genomes: control + 8 turtles
printf "Mimp_v3\t%s\t%s\t%s\n" $REL/Mimp_v3.fa $REL/Mimp_v3.proteins.faa $REL/Mimp_v3.genes.gff3 > genomes.tsv
while IFS=$'\t' read acc sp asm url; do d=$R/turtles/$acc; b=$(basename $url)
  for f in $d/${b}_genomic.fna.gz $d/${b}_protein.faa.gz $d/${b}_genomic.gff.gz; do [ -s ${f%.gz} ] || gunzip -k $f; done
  printf "%s\t%s\t%s\t%s\n" $acc $d/${b}_genomic.fna $d/${b}_protein.faa $d/${b}_genomic.gff >> genomes.tsv
done < $R/turtle_assemblies.tsv
cat genomes.tsv
while IFS=$'\t' read name fa faa gff; do
  echo "=== $name $(date)"
  conda activate $CONDA_ROOT/envs/binning
  [ -s out/$name.gh18.domtbl ] || hmmsearch --cpu 8 -E 1e-5 --domtblout out/$name.gh18.domtbl q/GH18.hmm $faa > /dev/null
  conda activate $CONDA_ROOT/envs/annot
  [ -s out/$name.miniprot.gff ] || miniprot -t $T --gff -N 30 --outs 0.3 $fa q/mimp_locus_and_flanks.faa > out/$name.miniprot.gff 2> out/$name.miniprot.log
  echo "  hmm hits: $(grep -vc '^#' out/$name.gh18.domtbl)  miniprot mRNAs: $(awk '$3=="mRNA"' out/$name.miniprot.gff | wc -l)"
done < genomes.tsv
conda activate $CONDA_ROOT/envs/binning
python3 $R/gh18/comparative_gh18_parse.py genomes.tsv out > comparative_gh18_summary.tsv
cat comparative_gh18_summary.tsv
echo "=== DONE $(date)"
