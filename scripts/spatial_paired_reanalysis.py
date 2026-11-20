#!/usr/bin/env python3
"""
Per-protein paired re-analysis of the segment-resolved proteomics (DS04).

Why this exists: the original spatial test (spatial_chitinase_analysis.py)
pooled every protein x replicate intensity as an independent observation
(~80 values per segment), which is pseudoreplication - proteins differ in
absolute abundance by orders of magnitude, so pooled tests mostly measure
which proteins happen to be in the set. Both simulated referees flagged it.

Here the PROTEIN is the unit of analysis:
  1. each protein's per-segment mean is centred on its own row mean
     (within-protein profile), removing between-protein abundance differences;
  2. Friedman test across the six segments on those paired profiles;
  3. Wilcoxon signed-rank, small intestine (mean of D/J/I) vs hindgut (C/R)
     and vs stomach, paired within protein;
  4. exact binomial test on where each protein peaks (3 of 6 segments are SI,
     so null P(peak in SI) = 0.5).

The GH18 set is also curated: entries the phylogeny places in the CTBS or
CHID1 clades are not chitinases and are dropped from the GH18 test.
All 18 runs come from ONE animal (three tissue pieces per segment), so even
the paired result describes one individual, not the species.

Outputs tables/31_spatial_paired_statistics.tsv
"""
import os, numpy as np, pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(BASE, 'tables')
SEG = ["Stomach", "Duodenum", "Jejunum", "Ileum", "Colon", "Rectum"]
SI, HG = ["Duodenum", "Jejunum", "Ileum"], ["Colon", "Rectum"]

# not chitinases: CTBS (di-N-acetylchitobiase) and CHID1-clade, per the GH18 tree
NOT_CHITINASE = {"A0A8C0FXB1", "A0A452J3X8"}

df = pd.read_csv(os.path.join(T, "02_proteomics_chitin_proteins.tsv"), sep='\t', comment='#')
df = df[df['Dataset'] == 'Proteomics_1st_GI'].copy
df['acc'] = df['Protein_ID'].str.split('|').str[1]

segmean = {}
for s in SEG:
    cols = [c for c in df.columns if (c.startswith(s + '_') or c.startswith(s + ' _'))
            and any(c.rstrip.endswith(f'_{i}') for i in range(1, 4))][:3]
    segmean[s] = df[cols].apply(pd.to_numeric, errors='coerce').mean(axis=1)
M = pd.DataFrame(segmean)
M['Category'] = df['Category'].values
M['acc'] = df['acc'].values

rows = []
for cat in ['Chitinase (GH18)', 'Hexosaminidase (GH20)', 'GlcNAc metabolism']:
    sub = M[M['Category'] == cat]
    dropped = []
    if cat.startswith('Chitinase'):
        dropped = sorted(set(sub['acc']) & NOT_CHITINASE)
        sub = sub[~sub['acc'].isin(NOT_CHITINASE)]
    sub = sub.dropna(subset=SEG)          # paired tests need complete profiles
    n = len(sub)
    prof = sub[SEG].sub(sub[SEG].mean(axis=1), axis=0)   # within-protein centring

    fr_stat, fr_p = stats.friedmanchisquare(*[prof[s].values for s in SEG])
    si = prof[SI].mean(axis=1); hg = prof[HG].mean(axis=1); st = prof['Stomach']
    w1 = stats.wilcoxon(si, hg); w2 = stats.wilcoxon(si, st)
    peaks = sub[SEG].idxmax(axis=1)
    k_si = int(peaks.isin(SI).sum)
    binom = stats.binomtest(k_si, n, 0.5, alternative='greater')
    rows += [
        dict(Category=cat, n_proteins=n, dropped=';'.join(dropped) or '-', Test='Friedman across 6 segments (paired by protein)',
             statistic=f'{fr_stat:.2f}', p_value=f'{fr_p:.3g}', effect=''),
        dict(Category=cat, n_proteins=n, dropped='', Test='Wilcoxon signed-rank SI vs hindgut (paired)',
             statistic=f'W={w1.statistic:.0f}', p_value=f'{w1.pvalue:.3g}',
             effect=f'median diff {np.median(si-hg):+.2f} log2'),
        dict(Category=cat, n_proteins=n, dropped='', Test='Wilcoxon signed-rank SI vs stomach (paired)',
             statistic=f'W={w2.statistic:.0f}', p_value=f'{w2.pvalue:.3g}',
             effect=f'median diff {np.median(si-st):+.2f} log2'),
        dict(Category=cat, n_proteins=n, dropped='', Test='Binomial: proteins peaking in SI (null 3/6)',
             statistic=f'{k_si}/{n}', p_value=f'{binom.pvalue:.3g}', effect=''),
    ]
    print(f'{cat}: n={n} dropped={dropped}  Friedman P={fr_p:.3g}  SI-vs-HG P={w1.pvalue:.3g} '
          f'(median {np.median(si-hg):+.2f})  SI-vs-St P={w2.pvalue:.3g}  peak-in-SI {k_si}/{n} P={binom.pvalue:.3g}')
    if cat.startswith('Chitinase'):
        print('  per-protein peak segment:', peaks.value_counts.to_dict)

out = pd.DataFrame(rows)
with open(os.path.join(T, '31_spatial_paired_statistics.tsv'), 'w') as fh:
    fh.write('# Per-protein PAIRED re-analysis of DS04 segment-resolved proteomics (one animal, 3 tissue pieces per segment).\n')
    fh.write('# Unit of analysis = protein (within-protein centred profile). Replaces the pooled protein x replicate tests in 21_spatial_statistics.tsv.\n')
    fh.write('# GH18 set curated: CTBS (A0A8C0FXB1) and CHID1-clade (A0A452J3X8) removed.\n')
    out.to_csv(fh, sep='\t', index=False)
print('wrote tables/31_spatial_paired_statistics.tsv')
