#!/bin/bash
# Pfam-A (current release) to $REFDIR/pfam: the EBI FTP allows ~0.4 MB/s per connection from this site,
# so fetch 24 byte ranges in parallel (resumable), concatenate, verify the size against the server, gunzip, hmmpress.
set -u
D=$REFDIR/pfam; B=https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release; N=24
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/binning
mkdir -p $D/parts; cd $D
pkill -f "curl -sS -o Pfam-A.hmm.gz" 2>/dev/null; rm -f Pfam-A.hmm.gz
curl -sS -o Pfam.version.gz $B/Pfam.version.gz; curl -sS -o relnotes.txt $B/relnotes.txt; gunzip -f Pfam.version.gz
SIZE=$(curl -sIL $B/Pfam-A.hmm.gz | grep -i content-length | tail -1 | tr -dc 0-9); echo "=== start $(date) size $SIZE"
chunk=$(( (SIZE + N - 1) / N )); cd parts
for i in $(seq 0 $((N-1))); do
  s=$((i*chunk)); e=$((s+chunk-1)); [ $e -ge $SIZE ] && e=$((SIZE-1))
  ( f=part_$(printf %02d $i); want=$((e-s+1))
    for try in $(seq 1 10); do
      have=$(stat -c %s $f 2>/dev/null || echo 0); [ "$have" -ge "$want" ] && break
      curl -sL -r $((s+have))-$e "$B/Pfam-A.hmm.gz" >> $f || sleep 20
    done ) &
done
wait; cat part_* > ../Pfam-A.hmm.gz; cd ..; got=$(stat -c %s Pfam-A.hmm.gz); echo "size $got expected $SIZE"
if [ "$got" = "$SIZE" ]; then gunzip -f Pfam-A.hmm.gz && hmmpress -f Pfam-A.hmm > /dev/null && rm -rf parts && echo "PFAM_READY $(head -1 Pfam.version) $(grep -c '^NAME' Pfam-A.hmm) models $(date)"; else echo "SIZE MISMATCH"; fi
