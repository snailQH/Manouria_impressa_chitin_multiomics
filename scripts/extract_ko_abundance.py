#!/usr/bin/env python3
"""Per-sample KEGG-ortholog abundance (read counts) for a cohort from eggNOG-mapper KEGG_ko
column + Bowtie2 per-gene counts. Same join as extract_cazy_abundance.py.
Usage: extract_ko_abundance.py --annotations X.emapper.annotations --counts-dir DIR --out out.tsv [--kos K01183,K12373,...]"""
import argparse, os, sys
from collections import defaultdict
ap = argparse.ArgumentParser; ap.add_argument('--annotations', required=True); ap.add_argument('--counts-dir', required=True)
ap.add_argument('--out', required=True); ap.add_argument('--kos', default='K01183,K12373,K01207,K01443,K02564,K00884,K01452,K03478,K03431,K04042,K02804')
a = ap.parse_args; want = set(a.kos.split(','))
gene2ko = {}; col = None
for line in open(a.annotations):
    if line.startswith('#'):
        if line.startswith('#query'): col = line.rstrip('\n').split('\t').index('KEGG_ko')
        continue
    p = line.rstrip('\n').split('\t')
    if len(p) <= (col or 11): continue
    kos = {k.replace('ko:', '') for k in p[col or 11].split(',') if k != '-'} & want
    if kos: gene2ko[p[0]] = kos
per = {}
for fn in sorted(os.listdir(a.counts_dir)):
    if fn.endswith('.counts.tsv'):
        d = {}
        for l in open(os.path.join(a.counts_dir, fn)):
            g, c = l.rstrip('\n').split('\t')[:2]
            if int(c): d[g] = int(c)
        per[fn[:-len('.counts.tsv')]] = d
samples = sorted(per); tab = defaultdict(lambda: defaultdict(int)); tot = {}
for s in samples:
    tot[s] = sum(per[s].values)
    for g, c in per[s].items:
        for k in gene2ko.get(g, ): tab[k][s] += c
with open(a.out, 'w') as o:
    o.write('KO\t' + '\t'.join(samples) + '\n')
    for k in sorted(want): o.write(k + '\t' + '\t'.join(str(tab[k].get(s, 0)) for s in samples) + '\n')
    o.write('TOTAL_MAPPED\t' + '\t'.join(str(tot[s]) for s in samples) + '\n')
print(f'wrote {a.out}; genes with wanted KOs: {len(gene2ko):,}', file=sys.stderr)
