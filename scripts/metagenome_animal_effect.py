#!/usr/bin/env python3
"""Animal-aware re-test of the cross-diet metagenome results now that the faecal donors are known (collection series 1 = MI02,
series 2 = MI03; MI01 contributed no faeces). Inputs: tables/35 (GH18 CPM per sample), tables/38 (unified genus profiles).
Tests: (i) animal effect on GH18 within each cohort; (ii) diet effect on GH18 within each animal; (iii) KW on animal-centred GH18;
(iv) PERMANOVA (Bray–Curtis, genus) for diet with free permutations and with permutations restricted within animal (strata),
for animal (free and within diet), and for diet within each animal. Output: tables/49_metagenome_animal_effect.tsv"""
import os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, kruskal
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); T = os.path.join(BASE, 'tables'); rng = np.random.default_rng(20261110); NPERM = 9999
cz = pd.read_csv(os.path.join(T, '35_unified_cazy_per_sample.tsv'), sep='\t')
series = lambda s: (s.split('-')[1] if '-' in s else s.split('.')[1])
cz['animal'] = cz['Sample'].map(lambda s: 'MI02' if series(s) == '1' else 'MI03')
rows = []
def add(section, test, value, P='', note=''): rows.append(dict(section=section, test=test, value=value, P=P, note=note))
for (c, a), g in cz.groupby(['Cohort', 'animal']): add('design', f'{c} × {a}', f'n={len(g)}; GH18 CPM {g.GH18_CPM.mean:.0f} ± {g.GH18_CPM.std:.0f}')
for c, g in cz.groupby('Cohort'):
    a, b = g[g.animal == 'MI02'].GH18_CPM, g[g.animal == 'MI03'].GH18_CPM
    add('animal effect within cohort (GH18)', f'{c}: MI02 vs MI03', f'{a.mean:.0f} vs {b.mean:.0f} CPM (n={len(a)}, {len(b)})', f'{mannwhitneyu(a, b).pvalue:.3f}', 'Mann–Whitney U, two-sided')
add('animal effect (GH18)', 'MI02 vs MI03, all cohorts', f'{cz[cz.animal=="MI02"].GH18_CPM.mean:.0f} vs {cz[cz.animal=="MI03"].GH18_CPM.mean:.0f} CPM', f'{kruskal(*[v.values for _, v in cz.groupby("animal").GH18_CPM]).pvalue:.3f}', 'Kruskal–Wallis')
for an, g in cz.groupby('animal'):
    add('diet effect within animal (GH18)', f'{an}: 3 cohorts', '; '.join(f'{c} {v.mean:.0f} (n={len(v)})' for c, v in g.groupby('Cohort').GH18_CPM), f'{kruskal(*[v.values for _, v in g.groupby("Cohort").GH18_CPM]).pvalue:.3f}', 'Kruskal–Wallis')
    for A, B in [('DS01', 'DS03'), ('DS02', 'DS03'), ('DS01', 'DS02')]:
        a, b = g[g.Cohort == A].GH18_CPM, g[g.Cohort == B].GH18_CPM
        add('diet effect within animal (GH18)', f'{an}: {A} vs {B}', f'fold {a.mean / b.mean:.2f} (n={len(a)}, {len(b)})', f'{mannwhitneyu(a, b).pvalue:.3f}', 'Mann–Whitney U, two-sided')
cz['GH18_centred'] = cz.GH18_CPM - cz.groupby('animal').GH18_CPM.transform('mean')
add('diet effect adjusted for animal (GH18)', 'KW across cohorts on animal-centred GH18 CPM', '', f'{kruskal(*[v.values for _, v in cz.groupby("Cohort").GH18_centred]).pvalue:.4f}', 'each sample minus its animal mean')
# PERMANOVA
gen = pd.read_csv(os.path.join(T, '38_unified_genus.tsv'), sep='\t', index_col=0); gen = gen.loc[[i for i in gen.index if i != 'Others']]
X = gen.T; X = X.div(X.sum(axis=1), axis=0); samples = X.index.tolist; A = X.values; n = len(samples)
D = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n): D[i, j] = D[j, i] = np.abs(A[i] - A[j]).sum / (A[i] + A[j]).sum
meta = cz.set_index('Sample').loc[samples]
def pseudoF(D, lab):
    lab = np.asarray(lab); N = len(lab); SSW = 0
    for k in np.unique(lab):
        idx = np.where(lab == k)[0]; SSW += (D[np.ix_(idx, idx)] ** 2).sum / (2 * len(idx))
    SST = (D ** 2).sum / (2 * N); SSA = SST - SSW; a = len(np.unique(lab)); return (SSA / (a - 1)) / (SSW / (N - a)), SSA / SST
def permanova(D, lab, strata=None):
    F, R2 = pseudoF(D, lab); lab = np.asarray(lab); cnt = 0
    for _ in range(NPERM):
        p = lab.copy
        if strata is None: p = rng.permutation(p)
        else:
            st = np.asarray(strata)
            for s in np.unique(st): idx = np.where(st == s)[0]; p[idx] = rng.permutation(p[idx])
        if pseudoF(D, p)[0] >= F: cnt += 1
    return F, R2, (cnt + 1) / (NPERM + 1)
for name, lab, strata in [('diet, free permutations', meta.Cohort.values, None), ('diet, permutations within animal (strata = animal)', meta.Cohort.values, meta.animal.values),
                          ('animal, free permutations', meta.animal.values, None), ('animal, permutations within diet (strata = diet)', meta.animal.values, meta.Cohort.values)]:
    F, R2, P = permanova(D, lab, strata); add('PERMANOVA (Bray–Curtis, genus, 32 samples)', name, f'pseudo-F {F:.2f}, R² {R2:.3f}', f'{P:.4f}', f'{NPERM} permutations')
for an in ['MI02', 'MI03']:
    idx = np.where(meta.animal.values == an)[0]; F, R2, P = permanova(D[np.ix_(idx, idx)], meta.Cohort.values[idx])
    add('PERMANOVA within animal', f'diet within {an} (n={len(idx)})', f'pseudo-F {F:.2f}, R² {R2:.3f}', f'{P:.4f}', f'{NPERM} permutations')
out = pd.DataFrame(rows)
with open(os.path.join(T, '49_metagenome_animal_effect.tsv'), 'w') as fh:
    fh.write('# Animal-aware re-test of the cross-diet metagenome results (faecal donors: series 1 = MI02, series 2 = MI03; confirmed by the authors 2026-11-10).\n'); out.to_csv(fh, sep='\t', index=False)
print(out.to_string(index=False))
