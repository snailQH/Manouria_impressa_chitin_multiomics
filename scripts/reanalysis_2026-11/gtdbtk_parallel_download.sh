#!/bin/bash
# GTDB-Tk R226 package (141,442,235,198 bytes): the server allows ~0.4 MB/s per connection, so fetch
# 32 byte ranges in parallel (resumable), concatenate, verify the size, extract to $REFDIR/gtdbtk/release226.
set -u
D=$REFDIR/gtdbtk; U=https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
SIZE=141442235198; N=32; mkdir -p $D/parts; cd $D/parts
chunk=$(( (SIZE + N - 1) / N ))
echo "=== start $(date)"
for i in $(seq 0 $((N-1))); do
  s=$((i*chunk)); e=$((s+chunk-1)); [ $e -ge $SIZE ] && e=$((SIZE-1))
  ( f=part_$(printf %02d $i); want=$((e-s+1))
    for try in 1 2 3 4 5 6 7 8 9 10; do
      have=$(stat -c %s $f 2>/dev/null || echo 0); [ "$have" -ge "$want" ] && break
      curl -sL -r $((s+have))-$e "$U" >> $f || sleep 30
    done
    echo "$f $(stat -c %s $f) / $want" ) &
done
wait; echo "=== parts done $(date)"
cat part_* > ../gtdbtk_r226_data.tar.gz; got=$(stat -c %s ../gtdbtk_r226_data.tar.gz); echo "size $got expected $SIZE"
if [ "$got" = "$SIZE" ]; then cd $D && tar xzf gtdbtk_r226_data.tar.gz && rm -rf gtdbtk_r226_data.tar.gz parts && echo "=== extracted $(ls $D)"; else echo "SIZE MISMATCH"; fi
echo "=== end $(date)"
