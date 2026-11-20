#!/usr/bin/env python3
"""Server-side (head node): KEGG summaries of the unified 32-sample catalogue from the eggNOG-mapper annotations
and the per-sample Bowtie2 counts.
  unified32_kegg_level1_cpm.tsv  sample x KEGG level-1 category (reads per million mapped; a gene with several
                                 pathways is counted once per category it belongs to); categories from the map
                                 number (00010-00999 Metabolism; 01xxx Global/overview maps excluded; 02xxx
                                 Environmental information processing + Cellular processes (02010/02020/02024/02025/
                                 02026/02030/02040/02060 = EIP; 02xxx others cellular); 03xxx Genetic information
                                 processing; 04xxx Cellular processes / Organismal systems; 05xxx Human diseases)
  unified32_kegg_level2_cpm.tsv  sample x level-2 category by map-number range of the metabolism maps
  unified32_kegg_functional_summary.tsv  catalogue-level counts (annotated, KO, ko00520, chitin KOs)"""
import os, sys, collections
W = '$PROJECT/03.unified_analysis/15_unified32'
ANN = f'{W}/kegg/unified32.emapper.annotations'; CNT = f'{W}/counts'
EIP = {'02010','02020','02024','02025','02026','02030','02040','02060'}
L2 = [(10, 99, 'Carbohydrate metabolism'), (100, 199, 'Energy metabolism'), (200, 299, 'Lipid metabolism'), (300, 399, 'Nucleotide metabolism'),
      (400, 499, 'Amino acid metabolism'), (500, 599, 'Glycan biosynthesis and metabolism'), (600, 699, 'Metabolism of other amino acids / cofactors'),
      (700, 799, 'Metabolism of cofactors and vitamins'), (800, 899, 'Biosynthesis of other secondary metabolites'), (900, 999, 'Xenobiotics / terpenoids / polyketides')]
def level1(m):
    n = int(m[3:]) if m.startswith('map') else int(m)
    if n < 1000: return 'Metabolism'
    if 1000 <= n < 2000: return None
    if 2000 <= n < 3000: return 'Environmental information processing' if m[-5:] in EIP else 'Cellular processes'
    if 3000 <= n < 4000: return 'Genetic information processing'
    if 4000 <= n < 5000: return 'Cellular processes / Organismal systems'
    return 'Human diseases'
def level2(m):
    n = int(m[3:]) if m.startswith('map') else int(m)
    for a, b, lab in L2:
        if a <= n <= b: return lab
    return None
gene_l1, gene_l2 = {}, {}; n_ann = n_ko = n_520 = 0; chitin = collections.Counter
KOS = ['K01183','K12373','K01207','K01443','K02564','K00884','K01452','K03478','K03431','K04042','K02804','K00844','K02563','K00820']
for l in open(ANN):
    if l.startswith('#'): continue
    p = l.rstrip('\n').split('\t'); n_ann += 1
    if p[11] != '-':
        n_ko += 1
        for k in KOS:
            if ('ko:' + k) in p[11].split(','): chitin[k] += 1
    if p[12] != '-':
        maps = [m for m in p[12].split(',') if m.startswith('map')]
        if any(m == 'map00520' for m in maps): n_520 += 1
        l1 = {level1(m) for m in maps} - {None}; l2 = {level2(m) for m in maps} - {None}
        if l1: gene_l1[p[0]] = l1
        if l2: gene_l2[p[0]] = l2
samples = sorted(f[:-11] for f in os.listdir(CNT) if f.endswith('.counts.tsv'))
cats1 = sorted({c for v in gene_l1.values for c in v}); cats2 = [lab for _, _, lab in L2]
with open(f'{W}/unified32_kegg_level1_cpm.tsv', 'w') as o1, open(f'{W}/unified32_kegg_level2_cpm.tsv', 'w') as o2:
    o1.write('Sample\ttotal_mapped\t' + '\t'.join(cats1) + '\n'); o2.write('Sample\ttotal_mapped\t' + '\t'.join(cats2) + '\n')
    for s in samples:
        c1 = collections.Counter; c2 = collections.Counter; tot = 0
        for l in open(f'{CNT}/{s}.counts.tsv'):
            g, c = l.split('\t')[:2]; c = int(c); tot += c
            if c and g in gene_l1:
                for k in gene_l1[g]: c1[k] += c
            if c and g in gene_l2:
                for k in gene_l2[g]: c2[k] += c
        o1.write(s + f'\t{tot}\t' + '\t'.join(f'{c1[k] / tot * 1e6:.1f}' for k in cats1) + '\n')
        o2.write(s + f'\t{tot}\t' + '\t'.join(f'{c2[k] / tot * 1e6:.1f}' for k in cats2) + '\n')
        print(s, tot, flush=True)
with open(f'{W}/unified32_kegg_functional_summary.tsv', 'w') as o:
    o.write('metric\tvalue\n'); o.write(f'catalogue_proteins\t2850805\nannotated_eggNOG\t{n_ann}\nwith_KO\t{n_ko}\nin_ko00520\t{n_520}\n')
    for k in KOS: o.write(f'genes_{k}\t{chitin[k]}\n')
print('done')
