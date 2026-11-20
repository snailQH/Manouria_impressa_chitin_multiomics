#!/bin/bash
# Pfam domain architecture of the twelve Mimp_v3 GH18-locus gene models (MIMP3_g08468-g08479) and the four human GH18 paralogs.
# hmmscan (HMMER 3.4) against Pfam-A (current release, $REFDIR/pfam) with the Pfam gathering thresholds (--cut_ga).
# Output: $PROJECT/05.reanalysis_2026-11/domain_architecture/{queries.faa,pfam.domtbl}
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/binning
W=$PROJECT/05.reanalysis_2026-11/domain_architecture; mkdir -p $W; cd $W
REL=$PROJECT/01.Genome_Assembly/new-version/06_release
python3 - <<'PY'
ids = {f"MIMP3_g{n:05d}" for n in range(8468, 8480)}; out = open('queries.faa', 'w'); keep = False
for l in open('$PROJECT/01.Genome_Assembly/new-version/06_release/Mimp_v3.proteins.faa'):
    if l.startswith('>'): keep = l[1:].split('.')[0].split[0] in ids
    if keep: out.write(l)
out.write(open('human_gh18.faa').read); out.close
PY
grep -c '>' queries.faa
hmmscan --cpu 8 --cut_ga --domtblout pfam.domtbl $REFDIR/pfam/Pfam-A.hmm queries.faa > pfam.hmmscan.txt
echo "=== DONE $(date): $(grep -vc '^#' pfam.domtbl) domain hits"
