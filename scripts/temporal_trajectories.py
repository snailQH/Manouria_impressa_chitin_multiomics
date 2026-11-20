#!/usr/bin/env python3
"""Temporal structure of the 32 faecal metagenomes (collection dates from the authors' registry, 2026-11-10).
The three diets were fed sequentially, so diet is confounded with time; this script asks whether the cross-diet
differences look like a step at the diet switches or like a continuous drift.
Inputs: tables/35 (GH18 CPM), tables/38 (genus profiles), tables/deposition/sample_registry_all.tsv (dates, animals).
Tests (all non-parametric, sample = unit, exploratory):
  (i)  GH18 CPM vs days since the first collection of the same cohort — Spearman, per cohort and per cohort × animal;
  (ii) GH18 CPM vs calendar time across all 32 samples (Spearman) versus the cohort effect (Kruskal–Wallis), and the
       within-cohort time trend after centring GH18 on the cohort mean;
  (iii) Bray–Curtis distance between consecutive collections of the same animal: within a diet versus across a diet switch
       (Mann–Whitney), with the time gap of each pair reported because the switch gaps are longer;
  (iv) within-diet drift: Spearman between Bray–Curtis distance and days apart for all same-animal, same-diet pairs.
Outputs: tables/69_temporal_trajectories.tsv, figures SuppFig_S4_temporal_trajectories.{png,pdf}."""
import pathlib, itertools, numpy as np, pandas as pd, matplotlib
from scipy.stats import spearmanr, kruskal, mannwhitneyu
from scipy.spatial.distance import pdist, squareform


def _registry(root):
    """The full registry is not published (it carries dissection and specimen notes); the public
    subset holds the animal, date and cohort these analyses need."""
    full = root / 'tables/deposition/sample_registry_all.tsv'
    return full if full.exists else root / 'tables/deposition/sample_registry_public.tsv'

matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = pathlib.Path(__file__).resolve.parents[1]
cz = pd.read_csv(R/'tables/35_unified_cazy_per_sample.tsv', sep='\t', comment='#')[['Cohort', 'Diet', 'Sample', 'GH18_CPM']]
reg = pd.read_csv(_registry(R), sep='\t', comment='#'); reg = reg[reg.sample_type.str.contains('faec', na=False)]
reg = reg.rename(columns={'sample_id': 'Sample', 'animal_id (unified)': 'animal'})[['Sample', 'animal', 'collection_date']]
d = cz.merge(reg, on='Sample'); assert len(d) == 32 and d.collection_date.notna.all
d['date'] = pd.to_datetime(d.collection_date); d['t'] = (d.date - d.date.min).dt.days
d['t_cohort'] = d.t - d.groupby('Cohort').t.transform('min'); d['gh18_c'] = d.GH18_CPM - d.groupby('Cohort').GH18_CPM.transform('mean')
G = pd.read_csv(R/'tables/38_unified_genus.tsv', sep='\t', index_col=0)[d.Sample].T
G = G.div(G.sum(axis=1), axis=0); BC = pd.DataFrame(squareform(pdist(G.values, 'braycurtis')), index=G.index, columns=G.index)
out = []
def add(sec, test, val, p='', note=''): out.append(dict(section=sec, test=test, value=val, P=p, note=note))
for c, g in d.groupby('Cohort'):
    r, p = spearmanr(g.t_cohort, g.GH18_CPM); add('GH18 vs time within cohort', c, f'rho {r:.2f} (n={len(g)}, span {g.t_cohort.max} d)', f'{p:.3f}')
    for a, ga in g.groupby('animal'):
        if len(ga) >= 4: r, p = spearmanr(ga.t_cohort, ga.GH18_CPM); add('GH18 vs time within cohort × animal', f'{c} {a}', f'rho {r:.2f} (n={len(ga)})', f'{p:.3f}')
r, p = spearmanr(d.t, d.GH18_CPM); add('GH18 vs calendar time, all 32', 'Spearman', f'rho {r:.2f}', f'{p:.4f}', 'time and diet are confounded')
k = kruskal(*[g.GH18_CPM for _, g in d.groupby('Cohort')]); add('GH18 vs cohort, all 32', 'Kruskal–Wallis', f'H {k.statistic:.2f}', f'{k.pvalue:.4f}')
r, p = spearmanr(d.t_cohort, d.gh18_c); add('GH18 (cohort-centred) vs days within cohort, all 32', 'Spearman', f'rho {r:.2f}', f'{p:.3f}', 'residual time trend after removing the cohort means')
# consecutive-pair distances per animal
pairs = []
for a, ga in d.sort_values('date').groupby('animal'):
    s = ga.Sample.tolist
    for x, y in zip(s[:-1], s[1:]):
        gx, gy = ga.set_index('Sample').loc[x], ga.set_index('Sample').loc[y]
        pairs.append(dict(animal=a, s1=x, s2=y, gap_days=int((gy.date - gx.date).days), switch=gx.Cohort != gy.Cohort, bc=BC.loc[x, y], cohorts=f'{gx.Cohort}>{gy.Cohort}'))
P = pd.DataFrame(pairs)
w, x = P[~P.switch].bc, P[P.switch].bc; u = mannwhitneyu(x, w, alternative='two-sided')
add('Bray–Curtis between consecutive collections (same animal)', 'within diet', f'median {w.median:.3f} (n={len(w)}; gap median {int(P[~P.switch].gap_days.median)} d, max {int(P[~P.switch].gap_days.max)})')
add('Bray–Curtis between consecutive collections (same animal)', 'across a diet switch', f'median {x.median:.3f} (n={len(x)}; gaps {", ".join(str(v) for v in P[P.switch].gap_days)} d)', f'{u.pvalue:.4f}', 'Mann–Whitney vs within-diet pairs; switch gaps are longer')
for _, q in P[P.switch].iterrows: add('Bray–Curtis across each switch', f'{q.animal} {q.cohorts} ({q.s1}→{q.s2})', f'{q.bc:.3f}', '', f'gap {q.gap_days} d')
# within-diet drift: BC vs days apart, all same-animal same-diet pairs
dr = []
for (c, a), g in d.groupby(['Cohort', 'animal']):
    for x, y in itertools.combinations(g.Sample, 2):
        dr.append(dict(Cohort=c, animal=a, days=abs(int((d.set_index('Sample').loc[x, 'date'] - d.set_index('Sample').loc[y, 'date']).days)), bc=BC.loc[x, y]))
DR = pd.DataFrame(dr); r, p = spearmanr(DR.days, DR.bc); add('Within-diet drift: Bray–Curtis vs days apart (same animal, same diet)', f'all pairs (n={len(DR)})', f'rho {r:.2f}', f'{p:.4f}')
for c, g in DR.groupby('Cohort'):
    r, p = spearmanr(g.days, g.bc); add('Within-diet drift: Bray–Curtis vs days apart (same animal, same diet)', f'{c} (n={len(g)})', f'rho {r:.2f}', f'{p:.3f}')
# within-diet pairs with long gaps vs switch pairs (gap-matched comparison)
long_w = DR[DR.days >= P[P.switch].gap_days.min]
if len(long_w): u2 = mannwhitneyu(P[P.switch].bc, long_w.bc, alternative='two-sided'); add('Gap-matched: switch pairs vs within-diet pairs with gap >= shortest switch gap', f'within-diet long-gap pairs n={len(long_w)}, median BC {long_w.bc.median:.3f}; switch median {P[P.switch].bc.median:.3f}', f'{u2.pvalue:.4f}', 'Mann–Whitney')
O = pd.DataFrame(out)
with open(R/'tables/69_temporal_trajectories.tsv', 'w') as fh:
    fh.write('# Temporal structure of the 32 faecal metagenomes (collection dates: authors\' registry 2026-11-10; diets fed sequentially). Exploratory, sample = unit. scripts/temporal_trajectories.py\n'); O.to_csv(fh, sep='\t', index=False)
print(O.to_string(index=False))
# ---- Supp Fig S5
col = {'DS01': '#2E86AB', 'DS02': '#A23B72', 'DS03': '#F18F01'}; lab = {'DS01': 'P. eryngii', 'DS02': 'H. tessellatus', 'DS03': 'P. ostreatus'}   # palette as generate_unified_figs.py
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={'width_ratios': [2, 1, 1]})
for a, m in (('MI02', 'o'), ('MI03', 's')):
    ga = d[d.animal == a].sort_values('date'); ax[0].plot(ga.date, ga.GH18_CPM, '-', color='0.75', lw=0.8, zorder=1)
    for c, g in ga.groupby('Cohort'): ax[0].scatter(g.date, g.GH18_CPM, marker=m, color=col[c], s=34, edgecolor='k', lw=0.4, zorder=3, label=f'{lab[c]} ({a})')
for c, g in d.groupby('Cohort'): ax[0].axvspan(g.date.min, g.date.max, color=col[c], alpha=0.07)
ax[0].set_ylabel('GH18 (CPM, dbCAN3)', fontsize=8); ax[0].set_title('A  GH18 abundance by collection date (circles MI02, squares MI03; shading = feeding period)', loc='left', fontsize=9)
ax[0].set_ylim(top=d.GH18_CPM.max * 1.18); ax[0].set_xlabel('collection date', fontsize=8); ax[0].tick_params(axis='x', labelsize=7)
from matplotlib.patches import Patch
h0, l0 = ax[0].get_legend_handles_labels
h0 += [Patch(color=col[c], alpha=0.25, label=f'{lab[c]} feeding period') for c in ('DS01', 'DS02', 'DS03')]
ax[0].legend(handles=h0, fontsize=6.2, ncol=3, frameon=False, loc='upper left',
             bbox_to_anchor=(0, -0.22))                                   # below the axes, never over the series
for (c, a), g in d.groupby(['Cohort', 'animal']):                       # keep the circle/square animal code of panel A
    ax[1].scatter(g.t_cohort, g.gh18_c, color=col[c], s=26, edgecolor='k', lw=0.3,
                  marker={'MI02': 'o', 'MI03': 's'}[a], label=lab[c] if a == 'MI02' else None)
ax[1].axhline(0, color='0.6', lw=0.8); ax[1].set_xlabel('days since the first collection of the cohort', fontsize=8); ax[1].set_ylabel('GH18 CPM − cohort mean', fontsize=8)
rr = O[O.section.str.startswith('GH18 (cohort-centred)')].iloc[0]; ax[1].set_title(f'B  Within-cohort time trend ({rr.value}, P = {rr.P})', loc='left', fontsize=9); ax[1].legend(fontsize=7, frameon=False)
ax[2].scatter(DR.days, DR.bc, color='0.5', s=12, label='same diet, same animal'); ax[2].scatter(P[P.switch].gap_days, P[P.switch].bc, color='k', marker='D', s=36, label='consecutive, across a diet switch')
ax[2].set_xlabel('days between collections', fontsize=8); ax[2].set_ylabel('Bray–Curtis (genus)', fontsize=8); ax[2].legend(fontsize=7, frameon=False, loc='lower right')
rr = O[O.test.str.startswith('all pairs')].iloc[0]; ax[2].set_title(f'C  Distance vs days apart (within-diet drift {rr.value})', loc='left', fontsize=9)
for a_ in ax:
    for s in ('top', 'right'): a_.spines[s].set_visible(False)
fig.suptitle('Supplementary Figure S4. Diet switches versus time: GH18 abundance and community distance along the two animals\' collection series', fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.95))
for _d in ('png', 'pdf'):  # a fresh clone has no figures/ yet
    (R/'figures'/_d).mkdir(parents=True, exist_ok=True)
for ext in ('png', 'pdf'): fig.savefig(R/f'figures/{ext}/SuppFig_S4_temporal_trajectories.{ext}', dpi=300)
print('written SuppFig_S4_temporal_trajectories')
