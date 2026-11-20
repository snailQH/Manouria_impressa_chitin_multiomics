#!/usr/bin/env python3
"""
Figures 2A-2D, 3A, 3B, 3D, 4A and Supp Figs S7, S19 from the UNIFIED pipeline tables
(tables/35-43, produced by scripts/unified_cross_cohort.py on the single 32-sample catalogue with dbCAN3 CAZy families, 2026-11-11). All 32 samples were
processed identically, so every panel is one pipeline, one unit, one axis.
Supersedes the vendor-table versions of these panels in generate_plots.py /
generate_diversity_plots.py and a superseded vendor-table script, neither of which is called.
"""
import os, sys, numpy as np, pandas as pd, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial.distance import pdist, squareform
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 10, 'figure.dpi': 300})
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); T = os.path.join(BASE, 'tables')
COH = ['DS01', 'DS02', 'DS03']
LAB = {'DS01': 'P. eryngii\n(King Oyster)', 'DS02': 'H. tessellatus\n(Beech)', 'DS03': 'P. ostreatus\n(Oyster)'}
COL = {'DS01': '#2E86AB', 'DS02': '#A23B72', 'DS03': '#F18F01'}
rng = np.random.default_rng(0)

def save(fig, name):
    from figure_ownership import may_write
    if not may_write(name, __file__): return
    for d in ('png', 'pdf'):
        os.makedirs(os.path.join(BASE, 'figures', d), exist_ok=True)
        fig.savefig(os.path.join(BASE, 'figures', d, f'{name}.{d}'), dpi=300, bbox_inches='tight')
    plt.close(fig); print(' ', name)

caz = pd.read_csv(os.path.join(T, '35_unified_cazy_per_sample.tsv'), sep='\t')
order = list(caz.Sample); coh_of = caz.set_index('Sample').Cohort
phy = pd.read_csv(os.path.join(T, '37_unified_phylum.tsv'), sep='\t', index_col=0)[order]
gen = pd.read_csv(os.path.join(T, '38_unified_genus.tsv'), sep='\t', index_col=0)[order]
alpha = pd.read_csv(os.path.join(T, '39_unified_alpha_diversity.tsv'), sep='\t', comment='#', nrows=32)
tests = pd.read_csv(os.path.join(T, '36_unified_cazy_tests.tsv'), sep='\t')
idx = pd.read_csv(os.path.join(T, '43_unified_cell_wall_indices.tsv'), sep='\t', comment='#', nrows=32)

def boxes(ax, values_by_cohort, ylabel, title=None):
    for i, c in enumerate(COH):
        v = np.asarray(values_by_cohort[c])
        bp = ax.boxplot([v], positions=[i], widths=0.55, patch_artist=True, showfliers=False)
        bp['boxes'][0].set_facecolor(COL[c]); bp['boxes'][0].set_alpha(0.55); bp['medians'][0].set_color('black')
        ax.scatter(i + rng.normal(0, 0.07, len(v)), v, c=COL[c], s=26, alpha=0.85, edgecolors='white', linewidth=0.5, zorder=3)
    ax.set_xticks(range(3)); ax.set_xticklabels([LAB[c] for c in COH], fontsize=8); ax.set_ylabel(ylabel)
    if title: ax.set_title(title, fontweight='bold', fontsize=10)

def bracket(ax, i, j, y, text):
    ax.plot([i, i, j, j], [y, y * 1.02, y * 1.02, y], color='k', lw=1); ax.text((i + j) / 2, y * 1.03, text, ha='center', fontsize=8)

# ---------- Fig 2A phylum ----------
rel = phy * 100; top = rel.mean(axis=1).nlargest(9).index; plot = rel.loc[top].copy; plot.loc['Other'] = rel.loc[~rel.index.isin(top)].sum
fig, ax = plt.subplots(figsize=(12, 5)); bottom = np.zeros(len(order)); cmap = plt.get_cmap('tab10')
for k, p in enumerate(plot.index):
    ax.bar(range(len(order)), plot.loc[p].values, bottom=bottom, color=cmap(k) if p != 'Other' else '#BBBBBB', label=p, width=0.85); bottom += plot.loc[p].values
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=7); ax.set_ylabel('Relative abundance (% of classified reads)')
x0 = 0
for c in COH:
    n = (coh_of[order] == c).sum; ax.axvspan(x0 - 0.5, x0 + n - 0.5, color=COL[c], alpha=0.06); ax.text(x0 + n / 2 - 0.5, 101, LAB[c].split('\n')[0], ha='center', fontsize=9, style='italic'); x0 += n
ax.set_ylim(0, 106); ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8, frameon=False)
ax.set_title('Phylum-level composition, all 32 samples, one pipeline (Kraken2 conf 0.1 + Bracken)', fontweight='bold', fontsize=10)
save(fig, 'Fig2A_phylum_barplot')

# ---------- Fig 2B PCoA ----------
D = squareform(pdist(gen.T.values, 'braycurtis')); n = len(order); J = np.eye(n) - np.ones((n, n)) / n; B = -0.5 * J @ (D ** 2) @ J
w, v = np.linalg.eigh(B); o = np.argsort(w)[::-1]; w, v = w[o], v[:, o]; pc = v[:, :2] * np.sqrt(np.clip(w[:2], 0, None)); ve = w[:2] / w[w > 0].sum * 100
fig, ax = plt.subplots(figsize=(6.2, 5.2))
for c in COH:
    m = (coh_of[order] == c).values; ax.scatter(pc[m, 0], pc[m, 1], c=COL[c], s=55, edgecolors='white', label=LAB[c].split('\n')[0], zorder=3)
    cx, cy = pc[m, 0].mean, pc[m, 1].mean
    for x, y in pc[m]: ax.plot([cx, x], [cy, y], color=COL[c], alpha=0.25, lw=0.8)
perm = pd.read_csv(os.path.join(T, '40_unified_permanova.tsv'), sep='\t'); r = perm.iloc[0]
ax.set_xlabel(f'PCo1 ({ve[0]:.1f}%)'); ax.set_ylabel(f'PCo2 ({ve[1]:.1f}%)'); ax.legend(frameon=False, fontsize=8)
ax.set_title(f'PCoA, Bray-Curtis on genus profiles\nPERMANOVA 3 cohorts: pseudo-F = {r.pseudo_F:.2f}, R² = {r.R2:.2f}, P = {r.P:.3f}', fontsize=9.5, fontweight='bold')
save(fig, 'Fig2B_PCoA')

# ---------- Fig 2C alpha ----------
fig, axes = plt.subplots(1, 3, figsize=(11, 4))
for ax, m, lab in zip(axes, ['Shannon', 'Simpson', 'Genera_observed'], ['Shannon index (genus)', 'Simpson index (genus)', 'Genera observed']):
    boxes(ax, {c: alpha[alpha.Cohort == c][m] for c in COH}, lab)
    kw = stats.kruskal(*[alpha[alpha.Cohort == c][m] for c in COH]); ax.set_title(f'Kruskal-Wallis P = {kw.pvalue:.3g}', fontsize=9)
axes[2].set_title('note: DS01 sequenced ~3x deeper — richness is depth-sensitive', fontsize=7, color='#555')
fig.suptitle('Alpha diversity by diet cohort (Bracken genus profiles, one pipeline)', fontweight='bold', fontsize=11); plt.tight_layout
save(fig, 'Fig2C_alpha_diversity')

# ---------- Fig 2D genus heatmap ----------
topg = gen.mean(axis=1).nlargest(20).index; H = np.log10(gen.loc[topg] * 100 + 0.01)
fig, ax = plt.subplots(figsize=(12, 7)); im = ax.imshow(H.values, aspect='auto', cmap='YlGnBu')
ax.set_yticks(range(len(topg))); ax.set_yticklabels(topg, fontsize=8, style='italic'); ax.set_xticks(range(n)); ax.set_xticklabels(order, rotation=90, fontsize=7)
for i, s in enumerate(order): ax.get_xticklabels[i].set_color(COL[coh_of[s]])
plt.colorbar(im, ax=ax, fraction=0.03, label='log10(% of classified reads)'); ax.set_title('Top 20 genera across all 32 samples (one pipeline)', fontweight='bold', fontsize=10)
save(fig, 'Fig2D_genus_heatmap')

# ---------- Fig 3A GH18 ----------
fig, ax = plt.subplots(figsize=(6, 4.6)); boxes(ax, {c: caz[caz.Cohort == c].GH18_CPM for c in COH}, 'GH18 reads per million mapped reads')
ymax = caz.GH18_CPM.max
for (a, b, i, j, h) in [('DS02', 'DS03', 1, 2, 1.05), ('DS01', 'DS03', 0, 2, 1.16), ('DS01', 'DS02', 0, 1, 1.27)]:
    t = tests[(tests.Family == 'GH18') & (tests.A == a) & (tests.B == b)].iloc[0]; bracket(ax, i, j, ymax * h, f'{t.fold_A_over_B:.2f}×, P = {t.P:.3g}')
ax.set_ylim(top=ymax * 1.36); ax.set_title('GH18 chitinase abundance across diets\n(all 32 samples, one catalogue: dbCAN3 CAZy + Bowtie2 counts)', fontweight='bold', fontsize=10)
save(fig, 'Fig3A_GH18_boxplot')

# ---------- Fig 3B CAZy heatmap ----------
fams = ['GH18', 'GH19', 'GH20', 'CBM5', 'CBM12', 'CBM50', 'AA10']; M = caz.set_index('Sample')[[f + '_CPM' for f in fams]].T; M.index = fams
fig, ax = plt.subplots(figsize=(12, 3.8)); im = ax.imshow(np.log10(M.values + 1), aspect='auto', cmap='viridis')
ax.set_yticks(range(len(fams))); ax.set_yticklabels(fams); ax.set_xticks(range(n)); ax.set_xticklabels(order, rotation=90, fontsize=7)
for i, s in enumerate(order): ax.get_xticklabels[i].set_color(COL[coh_of[s]])
plt.colorbar(im, ax=ax, fraction=0.02, label='log10(CPM + 1)'); ax.set_title('Chitin-related CAZy families, reads per million (one pipeline)', fontweight='bold', fontsize=10)
save(fig, 'Fig3B_CAZyme_heatmap')

# ---------- Fig 3D / 7C indices ----------
def index_fig(name, suptitle):
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    for ax, m, lab in zip(axes, ['CDI_pct', 'chitin_to_cellulose'], ['Chitin degradation index (% of GH reads)', 'Chitin : cellulose ratio (GH18 / GH5+6+9+48)']):
        boxes(ax, {c: idx[idx.Cohort == c][m] for c in COH}, lab); kw = stats.kruskal(*[idx[idx.Cohort == c][m] for c in COH]); ax.set_title(f'Kruskal-Wallis P = {kw.pvalue:.3g}', fontsize=9)
    fig.suptitle(suptitle, fontweight='bold', fontsize=10); plt.tight_layout; save(fig, name)
index_fig('Fig3D_digestibility_indices', 'Cell-wall degradation indices, all 32 samples on one annotation')
# Fig 7C: GH18 abundance versus Bacteroides — the diet signal is GH18-specific and only partly compositional (tables/41)
gb = pd.read_csv(os.path.join(T, '41_unified_gh18_vs_bacteroides.tsv'), sep='\t', comment='#', nrows=32)
for c_ in ['GH18_CPM', 'Bacteroides_frac', 'GH18_CPM_per_Bacteroides_frac']: gb[c_] = pd.to_numeric(gb[c_])
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
ax = axes[0]
for c in COH:
    g = gb[gb.Cohort == c]; ax.scatter(g.Bacteroides_frac * 100, g.GH18_CPM, color=COL[c], s=36, edgecolor='#333', linewidth=0.5, label=LAB[c].split('\n')[0])
rho = stats.spearmanr(gb.Bacteroides_frac, gb.GH18_CPM); ax.set_xlabel('Bacteroides (% of Bracken-assigned reads)'); ax.set_ylabel('GH18 (CPM, dbCAN3)')
ax.set_title(f'GH18 tracks Bacteroides: Spearman ρ = {rho.statistic:.2f}, P = {rho.pvalue:.1g}', fontsize=9); ax.legend(frameon=False, fontsize=8); ax.spines[['top', 'right']].set_visible(False)
boxes(axes[1], {c: gb[gb.Cohort == c].GH18_CPM_per_Bacteroides_frac for c in COH}, 'GH18 CPM per unit Bacteroides fraction')
kw = stats.kruskal(*[gb[gb.Cohort == c].GH18_CPM_per_Bacteroides_frac for c in COH]); axes[1].set_title(f'GH18 per unit Bacteroides: Kruskal-Wallis P = {kw.pvalue:.3g}', fontsize=9)
fig.suptitle('Diet signal in the GH18 pool: proportional to Bacteroides in the two Bacteroidota-led cohorts, in excess of it under P. eryngii', fontweight='bold', fontsize=10)
plt.tight_layout; save(fig, 'SuppFig_S19_diet_response')

# ---------- Supp S8 top genera per cohort ----------
fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), sharex=True)
xmax = max((gen[caz[caz.Cohort == c].Sample].mean(axis=1) * 100).max for c in COH)
for k, (ax, c) in enumerate(zip(axes, COH)):
    ss = caz[caz.Cohort == c].Sample; mm = (gen[ss].mean(axis=1) * 100).nlargest(10)[::-1]
    ax.barh(range(10), mm.values, color=COL[c]); ax.set_yticks(range(10)); ax.set_yticklabels(mm.index, fontsize=8, style='italic')
    ax.set_xlabel('% of classified reads (mean)'); ax.set_xlim(0, xmax * 1.08)   # one scale across panels
    ax.set_title(f'{"ABC"[k]}  ' + LAB[c].split('\n')[0], loc='left', fontsize=10)
fig.suptitle('Supplementary Figure S7. Top 10 genera per diet cohort (one pipeline; common x-axis)', fontweight='bold')
plt.tight_layout; save(fig, 'SuppFig_S7_genus_barplot')

# ---------- Fig 4A KEGG chitin-pathway KOs ----------
ko = pd.read_csv(os.path.join(T, '45_unified_ko_per_sample.tsv'), sep='\t'); kt = pd.read_csv(os.path.join(T, '46_unified_ko_tests.tsv'), sep='\t').set_index('KO')
panels = [('K01183', 'K01183 chitinase'), ('K12373', 'K12373 HEXA/B hexosaminidase'), ('K01207', 'K01207 nagZ hexosaminidase'), ('K01443', 'K01443 nagA deacetylase'),
          ('K02564', 'K02564 nagB deaminase'), ('K02804', 'K02804 nagE PTS transporter'), ('K00884', 'K00884 NAGK GlcNAc kinase'), ('K03431', 'K03431 glmM')]
fig, axes = plt.subplots(2, 4, figsize=(15, 6.2))   # keeps Figure 4 under 225 mm with its legend
for ax, (k, lab) in zip(axes.flat, panels):
    boxes(ax, {c: ko[ko.Cohort == c][k + '_CPM'] for c in COH}, 'reads per million')
    ax.set_title(f'{lab}\nKruskal-Wallis P = {kt.loc[k, "KW_P"]:.3g}', fontsize=8.5, fontweight='bold')
fig.suptitle('KEGG chitin → GlcNAc → glycolysis orthologs across diets (all 32 samples, eggNOG KO + Bowtie2 counts)', fontweight='bold', fontsize=11)
plt.tight_layout; save(fig, 'Fig4A_KEGG_chitin_KOs')
