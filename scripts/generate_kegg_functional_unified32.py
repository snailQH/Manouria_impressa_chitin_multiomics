#!/usr/bin/env python3
"""Supp Fig S18 + Table S25 on the single 32-sample catalogue (replaces the provider-era DS02-vs-DS03 KEGG level-1
comparison and the per-cohort-catalogue KO summaries).
Inputs (results/unified32/, produced on the server by scripts/reanalysis_2026-11/unified32_kegg_summary.py):
  unified32_kegg_level1_cpm.tsv, unified32_kegg_level2_cpm.tsv (sample x category, reads per million mapped;
  categories from KEGG map numbers; global/overview maps excluded), unified32_kegg_functional_summary.tsv.
Outputs: tables/63_unified32_kegg_level1_level2_cpm.tsv (+ Kruskal-Wallis per category), tables/27_unified32_kegg_functional_summary.tsv,
         figures SuppFig_S18_KEGG_functional."""
import os, pandas as pd, numpy as np
from scipy import stats
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); U = os.path.join(BASE, 'results', 'unified32'); T = os.path.join(BASE, 'tables')
coh = lambda s: 'DS01' if s.startswith('ITF') else 'DS02' if s.startswith('HM') else 'DS03'
DIET = {'DS01': 'P. eryngii', 'DS02': 'H. tessellatus', 'DS03': 'P. ostreatus'}; COL = {'DS01': '#2E86AB', 'DS02': '#A23B72', 'DS03': '#F18F01'}
l1 = pd.read_csv(os.path.join(U, 'unified32_kegg_level1_cpm.tsv'), sep='\t').set_index('Sample'); l2 = pd.read_csv(os.path.join(U, 'unified32_kegg_level2_cpm.tsv'), sep='\t').set_index('Sample')
rows = []
for level, df in [('level 1', l1)]:   # level-2 by map-number range is not reliable (KEGG numbering is not range-coded) and is not reported
    for cat in [c for c in df.columns if c != 'total_mapped']:
        v = df[cat]; g = {c: v[[s for s in v.index if coh(s) == c]] for c in DIET}
        kw = stats.kruskal(*g.values).pvalue
        rows.append(dict(level=level, category=cat, **{f'{c}_mean_CPM': round(g[c].mean, 1) for c in DIET}, **{f'{c}_sd': round(g[c].std, 1) for c in DIET}, kruskal_P=round(kw, 4)))
res = pd.DataFrame(rows)
with open(os.path.join(T, '63_unified32_kegg_level1_cpm.tsv'), 'w') as fh:
    fh.write('# KEGG level-1 category abundance (reads per million mapped; a gene counted once per category) on the single 32-sample catalogue, eggNOG-mapper KEGG_Pathway column; categories by map-number range (00010-00999 metabolism; 02xxx environmental information processing / cellular processes; 03xxx genetic information processing; 04xxx cellular processes / organismal systems; 05xxx human diseases; 01xxx global/overview maps excluded). Kruskal-Wallis across the three cohorts.\n')
    res.to_csv(fh, sep='\t', index=False)
summ = pd.read_csv(os.path.join(U, 'unified32_kegg_functional_summary.tsv'), sep='\t'); summ.to_csv(os.path.join(T, '27_unified32_kegg_functional_summary.tsv'), sep='\t', index=False)
# figure: level-1 (left) and level-2 metabolism (right), cohort means with per-sample points
fig, ax = plt.subplots(figsize=(8.6, 5.6))
cats = [c for c in l1.columns if c != 'total_mapped']; y = np.arange(len(cats))[::-1]; w = 0.27
for i, c in enumerate(DIET):
    ss = [s for s in l1.index if coh(s) == c]; m = l1.loc[ss, cats] / 1e3
    ax.barh(y + (i - 1) * w, m.mean, w, color=COL[c], edgecolor='#333', linewidth=0.4, label=f'{DIET[c]} (n = {len(ss)})')
    for s in ss: ax.scatter(m.loc[s], y + (i - 1) * w + np.random.normal(0, 0.03, len(cats)), s=6, color='#333', alpha=0.5, zorder=3)
ax.set_yticks(y); ax.set_yticklabels(cats, fontsize=8.5)       # upright, unbroken category names
ax.set_xlabel('reads per thousand mapped'); ax.spines[['top', 'right']].set_visible(False)
ax.legend(frameon=False, fontsize=8, loc='upper right')
ax.set_title('Supplementary Figure S18. KEGG level-1 functional categories across the three diets\n'
             '(all 32 samples on the single catalogue; eggNOG KEGG pathways; global maps excluded)', fontweight='bold', fontsize=10)
fig.tight_layout
for _d in ('png', 'pdf'):  # a fresh clone has no figures/ yet
    os.makedirs(os.path.join(BASE, 'figures', _d), exist_ok=True)
for d in ('png', 'pdf'): fig.savefig(os.path.join(BASE, 'figures', d, f'SuppFig_S18_KEGG_functional.{d}'), dpi=300, bbox_inches='tight')
print(res[['level', 'category', 'DS01_mean_CPM', 'DS02_mean_CPM', 'DS03_mean_CPM', 'kruskal_P']].to_string(index=False))
