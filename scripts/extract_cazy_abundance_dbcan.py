#!/usr/bin/env python3
"""Per-sample CAZy family abundance from dbCAN3 (run_dbcan 5.x) overview tables + Bowtie2 per-gene counts.
Same join and output format as extract_cazy_abundance.py (eggNOG), so the two can be compared directly.
A gene is assigned to the families in 'Recommend Results' (>= 2 of HMM / dbCAN-sub / DIAMOND agree; subfamily
suffixes such as GH5_2 / GH13_e123 are collapsed to the family), one read count per family per gene.
Usage: extract_cazy_abundance_dbcan.py --overview overview.tsv [overview2.tsv ...] --counts-dir DIR --out out.tsv"""
import argparse, os, re
from collections import defaultdict
ap = argparse.ArgumentParser; ap.add_argument('--overview', nargs='+', required=True); ap.add_argument('--counts-dir', required=True); ap.add_argument('--out', required=True)
a = ap.parse_args
gene2fam = {}
for f in a.overview:
    with open(f) as fh:
        hdr = fh.readline.rstrip('\n').split('\t'); i_rec = hdr.index('Recommend Results'); i_n = hdr.index('#ofTools')
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) <= i_rec or int(p[i_n]) < 2: continue
            fams = set(re.findall(r'\b(GH\d+|GT\d+|PL\d+|CE\d+|AA\d+|CBM\d+)', p[i_rec]))
            if fams: gene2fam[p[0]] = fams
samples, per, tot = [], {}, {}
for fn in sorted(os.listdir(a.counts_dir)):
    if not fn.endswith('.counts.tsv'): continue
    s = fn[:-len('.counts.tsv')]; samples.append(s); d = defaultdict(int); t = 0
    for l in open(os.path.join(a.counts_dir, fn)):
        g, c = l.rstrip('\n').split('\t')[:2]; c = int(c); t += c
        for f in gene2fam.get(g, ): d[f] += c
    per[s] = d; tot[s] = t
fams = sorted({f for d in per.values for f in d})
with open(a.out, 'w') as o:
    o.write('# Per-sample CAZy family abundance (Bowtie2 read counts)\n# CAZy families from dbCAN3 run_dbcan Recommend Results (>= 2 tools)\n')
    o.write('Family\t' + '\t'.join(samples) + '\tTotal\n')
    for f in fams:
        v = [per[s].get(f, 0) for s in samples]; o.write(f + '\t' + '\t'.join(map(str, v)) + f'\t{sum(v)}\n')
    o.write('TOTAL_MAPPED\t' + '\t'.join(str(tot[s]) for s in samples) + f'\t{sum(tot.values)}\n')
print(f'{len(gene2fam)} CAZyme genes, {len(fams)} families, {len(samples)} samples -> {a.out}')
