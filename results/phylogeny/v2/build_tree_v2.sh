#!/bin/bash
# Commands that produced gh18_v2_trim.faa and the IQ-TREE output in this directory.
# Run from this directory with mafft, trimal and iqtree on $PATH.
# (Originally run in a conda environment named "phylo" on the project compute server.)
set -uo pipefail
cd "$(dirname "$0")"
echo "###### TREE v2 START $(date) ######"
mafft --auto --thread 16 gh18_v2.faa > gh18_v2_aln.faa 2> mafft.log
trimal -in gh18_v2_aln.faa -out gh18_v2_trim.faa -automated1
iqtree -s gh18_v2_trim.faa -m MFP -B 1000 -alrt 1000 -T 16 --prefix gh18_v2 -redo > iqtree_stdout.log 2>&1
echo "###### TREE v2 DONE $(date) ######"
