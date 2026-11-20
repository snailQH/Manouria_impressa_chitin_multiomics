#!/bin/bash
#SBATCH --job-name=k2_build2
#SBATCH --cpus-per-task=60
#SBATCH --mem=400G
#SBATCH --time=3-00:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/kraken2_fungiext_build2_%j.out
# Needs a node with ~500 GB RAM: the Kraken2 index is held in memory.
# Third attempt. Build 1 (94933): hash write to NFS truncated at 1.2 GB. Build 2 (95357): hash correct but the added
# library's prelim_map had been deleted, so the 594 new genomes were skipped (seqid2taxid.map = base DB only). Build here with the DB
# directory on the node's local disk, verify hash.k2d against its own header, then copy to shared storage.
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/taxonomy; export LC_ALL=C
NFS=$REFDIR/kraken2_pluspf_fungiext_2026-11; LOC=/tmp/k2_fungiext_build; mkdir -p $LOC; cd $LOC
rm -f hash.k2d hash.k2d.tmp opts.k2d opts.k2d.tmp taxo.k2d.tmp
ln -sfn $NFS/library library; ln -sfn $NFS/taxonomy taxonomy; rm -f seqid2taxid.map taxo.k2d
# step 1 must rerun: the added library's prelim_map was lost before build 2 (its 594 genomes were silently skipped),
# so seqid2taxid.map / taxonomy/prelim_map.txt are regenerated from library/*/prelim_map.txt here
[ -s $NFS/library/added/prelim_map.txt ] || cat $NFS/library/added/prelim_map_*.txt > $NFS/library/added/prelim_map.txt
[ "$(wc -l < $NFS/library/added/prelim_map.txt)" -gt 300000 ] || { echo 'added prelim_map.txt too small'; exit 1; }
echo "=== build start $(date) $(hostname); local free: $(df -h /tmp | tail -1)"
kraken2-build --build --db $LOC --threads 60
echo "seqid2taxid.map lines: $(wc -l < seqid2taxid.map) (base 192,551 + added 378,443 expected)"; grep -c 'kraken:taxid|5323' seqid2taxid.map
python3 - <<'PY'
import struct, os, sys
h=open('hash.k2d','rb').read(32); cap,size,kb,vb=struct.unpack('<QQQQ',h); exp=32+cap*4; act=os.path.getsize('hash.k2d')
print(f'hash.k2d capacity {cap} size {size} expected bytes {exp} actual {act}')
sys.exit(0 if exp==act else 1)
PY
echo "=== hash verified $(date); copying to NFS"
cp -f hash.k2d $NFS/hash.k2d.copy && cp -f opts.k2d $NFS/opts.k2d && cp -f taxo.k2d $NFS/taxo.k2d
[ "$(stat -c %s $NFS/hash.k2d.copy)" = "$(stat -c %s hash.k2d)" ] && mv $NFS/hash.k2d.copy $NFS/hash.k2d || { echo "COPY SIZE MISMATCH"; exit 1; }
echo "=== copied $(date)"; ls -la $NFS/*.k2d
kraken2-inspect --db $NFS --skip-counts | head -5
cp -f seqid2taxid.map $NFS/seqid2taxid.map   # Bracken (kmer2read_distr) needs it in the DB directory — a job failed here for lack of it
rm -rf $LOC
bracken-build -d $NFS -t 60 -k 35 -l 150   # if this step alone fails, rerun kraken2_fungiext_bracken.sbatch (reuses database.kraken)
echo "=== bracken done $(date)"
kraken2-inspect --db $NFS --threads 16 | grep -iE "Pleurotus ostreatus|Pleurotus eryngii|Hypsizygus|Neocallimastigomycota|Piromyces" | head
echo "=== ALL DONE $(date)"
