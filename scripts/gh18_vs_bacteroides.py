#!/usr/bin/env python3
"""Is the 1.7x GH18 difference between the H. tessellatus and P. ostreatus cohorts just a
Bacteroides bloom? Compares per-sample GH18 read counts (tables/01) with Bacteroides read
abundance and total genus-assigned reads (tables/11) for the two Novogene cohorts, and tests
whether GH18 normalised to Bacteroides still differs. Writes tables/34_gh18_vs_bacteroides.tsv."""
import os, pandas as pd, numpy as np
from scipy import stats
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
g = pd.read_csv(os.path.join(BASE, 'tables', '11_genus_abundance_cross_diet.tsv'), sep='\t', comment='#').set_index('Genus')
c = pd.read_csv(os.path.join(BASE, 'tables', '01_chitin_cazyme_cross_diet.tsv'), sep='\t', comment='#')
c = c[c.Dataset.isin(['H_tessellatus', 'P_ostreatus'])].copy
cols = {s: [x for x in g.columns if x.endswith('|' + s)][0] for s in c.Sample}
c['bacteroides_reads'] = [g.loc['g__Bacteroides', cols[s]] for s in c.Sample]
c['total_genus_reads'] = [g[cols[s]].sum for s in c.Sample]
c['bacteroides_frac'] = c.bacteroides_reads / c.total_genus_reads
c['GH18_per_total'] = c.GH18 / c.total_genus_reads * 1e6
c['GH18_per_bacteroides'] = c.GH18 / c.bacteroides_reads * 1e3
rows = []
for col, lab in [('GH18', 'GH18 read counts (as reported)'), ('bacteroides_frac', 'Bacteroides fraction of genus-assigned reads'),
                 ('GH18_per_total', 'GH18 per million genus-assigned reads'), ('GH18_per_bacteroides', 'GH18 per thousand Bacteroides reads')]:
    a = c[c.Dataset == 'H_tessellatus'][col]; b = c[c.Dataset == 'P_ostreatus'][col]
    U, p = stats.mannwhitneyu(a, b, alternative='two-sided')
    rows.append(dict(metric=lab, H_tessellatus_mean=a.mean, P_ostreatus_mean=b.mean, fold=a.mean / b.mean, U=U, P=p))
    print(f'{lab:48s} Ht={a.mean:10.3f} Po={b.mean:10.3f} fold={a.mean/b.mean:.2f} P={p:.3g}')
r, p = stats.spearmanr(c.GH18, c.bacteroides_reads); print(f'Spearman GH18 vs Bacteroides reads (n={len(c)}): rho={r:.2f} P={p:.3g}')
out = pd.DataFrame(rows)
with open(os.path.join(BASE, 'tables', '34_gh18_vs_bacteroides.tsv'), 'w') as fh:
    fh.write('# GH18 vs Bacteroides, Novogene cohorts only (DS02 n=13, DS03 n=7); vendor read-count tables 01 and 11.\n')
    fh.write(f'# Spearman GH18 vs Bacteroides reads across the 20 samples: rho={r:.2f}, P={p:.3g}\n')
    out.to_csv(fh, sep='\t', index=False)
c[['Dataset', 'Sample', 'GH18', 'bacteroides_reads', 'total_genus_reads', 'bacteroides_frac', 'GH18_per_total', 'GH18_per_bacteroides']].to_csv(
    os.path.join(BASE, 'tables', '34_gh18_vs_bacteroides_per_sample.tsv'), sep='\t', index=False)
