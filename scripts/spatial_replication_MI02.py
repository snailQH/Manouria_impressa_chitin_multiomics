#!/usr/bin/env python3
"""Does the proximal (small intestine > hindgut) host-chitinase pattern measured in MI01 (DS04, 3 pieces/segment) recur in the
second animal MI02 (DS05, 2 pieces/segment)? Same procedure as spatial_paired_reanalysis.py (tables/31): per-protein segment means,
within-protein centring, paired tests across the 6 shared GI segments; plus per-protein concordance of the SI−hindgut contrast
between the two animals. Output: tables/48_spatial_replication_MI02.tsv"""
import os, numpy as np, pandas as pd
from scipy import stats
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); T = os.path.join(BASE, 'tables')
SEG = ["Stomach", "Duodenum", "Jejunum", "Ileum", "Colon", "Rectum"]; SI, HG = ["Duodenum", "Jejunum", "Ileum"], ["Colon", "Rectum"]
NOT_CHITINASE = {"A0A8C0FXB1", "A0A452J3X8"}
df = pd.read_csv(os.path.join(T, "02_proteomics_chitin_proteins.tsv"), sep='\t', comment='#'); df['acc'] = df['Protein_ID'].str.split('|').str[1]
def profiles(dataset, nrep):
    d = df[df['Dataset'] == dataset].copy; seg = {}
    for s in SEG:
        cols = [c for c in d.columns if (c.startswith(s + '_') or c.startswith(s + ' _')) and any(c.rstrip.endswith(f'_{i}') for i in range(1, nrep + 1))]
        cols = [c for c in cols if d[c].notna.any][:nrep]
        seg[s] = d[cols].apply(pd.to_numeric, errors='coerce').mean(axis=1)
    M = pd.DataFrame(seg); M['Category'] = d['Category'].values; M['acc'] = d['acc'].values; return M
rows = []; per = {}
for animal, dataset, nrep in (('MI01', 'Proteomics_1st_GI', 3), ('MI02', 'Proteomics_2nd_MultiTissue', 2)):
    M = profiles(dataset, nrep)
    for cat in ['Chitinase (GH18)', 'Hexosaminidase (GH20)', 'GlcNAc metabolism']:
        sub = M[M['Category'] == cat]
        if cat.startswith('Chitinase'): sub = sub[~sub['acc'].isin(NOT_CHITINASE)]
        sub = sub.dropna(subset=SEG); n = len(sub)
        if n < 3: continue
        prof = sub[SEG].sub(sub[SEG].mean(axis=1), axis=0); si = prof[SI].mean(axis=1); hg = prof[HG].mean(axis=1)
        fr = stats.friedmanchisquare(*[prof[s].values for s in SEG]); w = stats.wilcoxon(si, hg)
        k = int((si > hg).sum); peaks = sub[SEG].idxmax(axis=1)
        rows.append(dict(animal=animal, dataset=dataset, pieces_per_segment=nrep, category=cat, n_proteins=n,
                         segment_means_centred_log2='; '.join(f'{s} {prof[s].mean:+.2f}' for s in SEG),
                         friedman_P=f'{fr.pvalue:.3g}', wilcoxon_SI_vs_hindgut_P=f'{w.pvalue:.3g}', median_SI_minus_hindgut_log2=f'{np.median(si - hg):+.2f}',
                         proteins_SI_gt_hindgut=f'{k}/{n}', peak_segment_counts=str(peaks.value_counts.to_dict)))
        if cat.startswith('Chitinase'): per[animal] = (si - hg).rename(animal).set_axis(sub['acc'].values)
# per-protein concordance of the SI-hindgut contrast between animals
both = pd.concat([per['MI01'], per['MI02']], axis=1).dropna
rho = stats.spearmanr(both['MI01'], both['MI02']); agree = int((np.sign(both['MI01']) == np.sign(both['MI02'])).sum)
rows.append(dict(animal='MI01 vs MI02', dataset='DS04 vs DS05', pieces_per_segment='3 vs 2', category='Chitinase (GH18) — per-protein SI−hindgut contrast', n_proteins=len(both),
                 segment_means_centred_log2='', friedman_P='', wilcoxon_SI_vs_hindgut_P=f'Spearman rho={rho.statistic:.2f}, P={rho.pvalue:.2g}', median_SI_minus_hindgut_log2='',
                 proteins_SI_gt_hindgut=f'same sign in {agree}/{len(both)}', peak_segment_counts=''))
out = pd.DataFrame(rows)
with open(os.path.join(T, '48_spatial_replication_MI02.tsv'), 'w') as fh:
    fh.write('# Replication check of the proximal host-chitinase pattern in the second animal (MI02, DS05, 2 tissue pieces per segment) with the\n# procedure of tables/31 (per-protein segment means, within-protein centring, paired tests over the six shared GI segments). GH18 curated (CTBS, CHID1-clade removed).\n')
    out.to_csv(fh, sep='\t', index=False)
print(out[['animal', 'category', 'n_proteins', 'friedman_P', 'wilcoxon_SI_vs_hindgut_P', 'median_SI_minus_hindgut_log2', 'proteins_SI_gt_hindgut']].to_string(index=False))
print(both.round(2).to_string)
