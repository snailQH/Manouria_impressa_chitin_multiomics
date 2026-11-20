#!/usr/bin/env python3
"""
Generate microbiome diversity and composition plots (Figure 2 + supplementary).
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 11,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.spines.top': False, 'axes.spines.right': False,
})

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, 'results')

def _src(name):
    """A vendor-era table lives in results/ on the analysis server and in tables/ in the published
    repository. Prefer whichever exists, so the figure scripts run in both places."""
    p = os.path.join(RESULTS, name)
    return p if os.path.exists(p) else os.path.join(BASE, 'tables', name)

PNG_DIR = os.path.join(BASE, 'figures', 'png')
PDF_DIR = os.path.join(BASE, 'figures', 'pdf')

DIET_COLORS = {'P.eryngii': '#2E86AB', 'H.tessellatus': '#A23B72', 'P.ostreatus': '#F18F01'}
DIET_LABELS = {'P.eryngii': 'P. eryngii\n(King Oyster)', 'H.tessellatus': 'H. tessellatus\n(Beech)', 'P.ostreatus': 'P. ostreatus\n(Oyster)'}
PHYLUM_COLORS = sns.color_palette('Set3', 12) + sns.color_palette('Paired', 8)


def save_fig(fig, name):
    from figure_ownership import may_write
    if not may_write(name, __file__): return
    fig.savefig(os.path.join(PNG_DIR, f'{name}.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(PDF_DIR, f'{name}.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f'  {name}')


def get_diet(col_name):
    if 'P.eryngii' in col_name or 'ITF' in col_name: return 'P.eryngii'
    if 'H.tessellatus' in col_name or 'HM.' in col_name: return 'H.tessellatus'
    if 'P.ostreatus' in col_name or 'PM.' in col_name: return 'P.ostreatus'
    return 'Unknown'


def get_sample(col_name):
    return col_name.split('|')[-1] if '|' in col_name else col_name


# ============================================================
# Fig 2A: Phylum-level stacked barplot
# ============================================================

# ============================================================
# Phylum nomenclature harmonisation
# ============================================================
# The vendor tables use different phylum nomenclature per cohort: DS01 carries
# the pre-2021 names (Firmicutes, Bacteroidetes, ...) while DS02/DS03 use the
# current ICNP names (Bacillota, Bacteroidota, ...). These are the SAME taxa.
# Left unmerged, every DS01 row reads as zero against the DS02/DS03 names,
# which is why DS01 was previously dropped from the phylum figure entirely.
PHYLUM_SYNONYMS = {
    'p_Firmicutes': 'Bacillota',
    'p__Bacillota': 'Bacillota',
    'p_Bacteroidetes': 'Bacteroidota',
    'p__Bacteroidota': 'Bacteroidota',
    'p_Proteobacteria': 'Pseudomonadota',
    'p__Pseudomonadota': 'Pseudomonadota',
    'p_Actinobacteria': 'Actinomycetota',
    'p__Actinomycetota': 'Actinomycetota',
    'p_Synergistetes': 'Synergistota',
    'p__Synergistota': 'Synergistota',
    'p_Spirochaetes': 'Spirochaetota',
    'p__Spirochaetota': 'Spirochaetota',
    'p_Ascomycota': 'Ascomycota',
    'p__Ascomycota': 'Ascomycota',
}


def harmonise_genera(df):
    """Strip the differing rank prefixes (g_ vs g__) so cohorts share names."""
    df = df.copy
    df.index = [str(i).replace('g__', '').replace('g_', '') for i in df.index]
    return df.groupby(level=0).sum


def harmonise_phyla(df):
    """Merge synonymous phylum names so all cohorts share one nomenclature."""
    def canon(name):
        if name in PHYLUM_SYNONYMS:
            return PHYLUM_SYNONYMS[name]
        return name.replace('p__', '').replace('p_', '')
    df = df.copy
    df.index = [canon(i) for i in df.index]
    return df.groupby(level=0).sum


def plot_fig2a:
    df = pd.read_csv(_src('10_phylum_abundance_cross_diet.tsv'), sep='\t', index_col=0, comment='#')
    df = harmonise_phyla(df)

    # Normalize to relative abundance per sample
    df_rel = df.div(df.sum(axis=0), axis=1) * 100

    # Keep top 10 phyla, group rest as "Other"
    top_phyla = df_rel.mean(axis=1).nlargest(10).index.tolist
    other = df_rel.loc[~df_rel.index.isin(top_phyla)].sum
    df_plot = df_rel.loc[df_rel.index.isin(top_phyla)].copy
    df_plot.loc['Other'] = other

    # All three dietary cohorts, now that nomenclature is harmonised
    cols = [c for c in df_plot.columns
            if get_diet(c) in ['P.eryngii', 'H.tessellatus', 'P.ostreatus']]
    df_plot = df_plot[cols]

    # Order by diet then sample
    diet_rank = {'P.eryngii': 0, 'H.tessellatus': 1, 'P.ostreatus': 2}
    col_order = sorted(cols, key=lambda x: (diet_rank.get(get_diet(x), 9), get_sample(x)))
    df_plot = df_plot[col_order]

    fig, ax = plt.subplots(figsize=(12, 5))

    bottom = np.zeros(len(col_order))
    phyla_list = df_plot.index.tolist
    colors = PHYLUM_COLORS[:len(phyla_list)]

    for i, phylum in enumerate(phyla_list):
        vals = df_plot.loc[phylum].values
        ax.bar(range(len(col_order)), vals, bottom=bottom, label=phylum.replace('p__', ''),
               color=colors[i], edgecolor='white', linewidth=0.5)
        bottom += vals

    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels([get_sample(c) for c in col_order], rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('Relative Abundance (%)')
    ax.set_title('Phylum-Level Gut Microbiome Composition', fontweight='bold')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7, title='Phylum')

    # Diet group brackets, derived from the actual column order so this stays
    # correct however many cohorts are present.
    boundaries, labels = [], []
    prev, run_start = None, 0
    for i, c in enumerate(col_order + [None]):
        d = get_diet(c) if c is not None else '__end__'
        if d != prev:
            if prev is not None:
                labels.append((prev, run_start, i))
                if c is not None:
                    boundaries.append(i - 0.5)
            prev, run_start = d, i
    for b in boundaries:
        ax.axvline(x=b, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    pretty = {'P.eryngii': 'P. eryngii', 'H.tessellatus': 'H. tessellatus',
              'P.ostreatus': 'P. ostreatus'}
    for name, s, e in labels:
        ax.text((s + e - 1) / 2, 105, f'{pretty.get(name, name)} (n={e - s})',
                ha='center', fontsize=9, fontweight='bold', style='italic')

    save_fig(fig, 'Fig2A_phylum_barplot')


# ============================================================
# Fig 2B: PCoA / PCA ordination
# ============================================================
def plot_fig2b:
    df = pd.read_csv(_src('11_genus_abundance_cross_diet.tsv'), sep='\t', index_col=0, comment='#')

    df = harmonise_genera(df)
    ds_cols = [c for c in df.columns
               if get_diet(c) in ['P.eryngii', 'H.tessellatus', 'P.ostreatus']]
    df_sub = df[ds_cols].T

    # Normalize to relative abundance
    df_norm = df_sub.div(df_sub.sum(axis=1), axis=0)
    df_norm = df_norm.fillna(0)

    # Bray-Curtis distance -> PCoA (classical MDS)
    bc_dist = pdist(df_norm.values, metric='braycurtis')
    bc_sq = squareform(bc_dist)

    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, normalized_stress='auto')
    coords = mds.fit_transform(bc_sq)

    fig, ax = plt.subplots(figsize=(6, 5))

    diets = [get_diet(c) for c in ds_cols]
    for diet in ['P.eryngii', 'H.tessellatus', 'P.ostreatus']:
        mask = [d == diet for d in diets]
        x = coords[mask, 0]
        y = coords[mask, 1]
        ax.scatter(x, y, c=DIET_COLORS[diet], s=60, alpha=0.8, edgecolors='white',
                   linewidth=1, label=DIET_LABELS[diet].replace('\n', ' '), zorder=3)
        # Draw ellipse (convex hull approximation)
        if len(x) > 2:
            from matplotlib.patches import Ellipse
            cx, cy = np.mean(x), np.mean(y)
            sx, sy = np.std(x) * 1.5, np.std(y) * 1.5
            ell = Ellipse((cx, cy), sx*2, sy*2, alpha=0.15, color=DIET_COLORS[diet])
            ax.add_patch(ell)

    ax.set_xlabel('MDS1')
    ax.set_ylabel('MDS2')
    ax.set_title('PCoA (Bray-Curtis) — Genus-Level', fontweight='bold')
    ax.legend(loc='best', framealpha=0.9)

    save_fig(fig, 'Fig2B_PCoA')


# ============================================================
# Fig 2C: Alpha diversity boxplots
# ============================================================
def plot_fig2c:
    df = pd.read_csv(_src('12_alpha_diversity_cross_diet.tsv'), sep='\t', comment='#')

    fig, axes = plt.subplots(1, 3, figsize=(11, 4))

    metrics = ['Shannon', 'Simpson', 'Chao1']
    diets = ['P.eryngii', 'H.tessellatus', 'P.ostreatus']
    diet_labels_short = ['P. eryngii', 'H. tessellatus', 'P. ostreatus']

    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        data_groups = []
        for diet in diets:
            vals = df[df['Diet'] == diet][metric].dropna.values
            data_groups.append(vals)

        bp = ax.boxplot(data_groups, patch_artist=True, showfliers=False, widths=0.5)
        for i, (patch, diet) in enumerate(zip(bp['boxes'], diets)):
            patch.set_facecolor(DIET_COLORS[diet])
            patch.set_alpha(0.6)
            bp['medians'][i].set_color('black')
            bp['medians'][i].set_linewidth(2)
            # Strip
            vals = data_groups[i]
            jitter = np.random.normal(0, 0.06, len(vals))
            ax.scatter(np.full(len(vals), i+1) + jitter, vals,
                      c=DIET_COLORS[diet], s=20, alpha=0.7, edgecolors='white', linewidth=0.5, zorder=3)

        ax.set_xticklabels(diet_labels_short, fontsize=8, rotation=15, ha='right')
        ax.set_title(metric, fontweight='bold')

    # Note about Chao1
    axes[2].text(0.5, -0.2, '*Chao1 not directly comparable\nacross pipelines (Kraken2 vs MicroNR)',
                transform=axes[2].transAxes, fontsize=7, color='gray', ha='center', style='italic')

    fig.suptitle('Alpha Diversity Across Mushroom Diets', fontweight='bold', fontsize=12, y=1.02)
    plt.tight_layout
    save_fig(fig, 'Fig2C_alpha_diversity')


# ============================================================
# Fig 2D: Top genus heatmap
# ============================================================
def plot_fig2d:
    df = pd.read_csv(_src('11_genus_abundance_cross_diet.tsv'), sep='\t', index_col=0, comment='#')

    df = harmonise_genera(df)
    ds_cols = [c for c in df.columns
               if get_diet(c) in ['P.eryngii', 'H.tessellatus', 'P.ostreatus']]
    df_sub = df[ds_cols]

    # Normalize to relative abundance
    df_norm = df_sub.div(df_sub.sum(axis=0), axis=1) * 100

    # Top 15 genera by mean abundance (exclude "Others")
    means = df_norm.mean(axis=1)
    top_genera = means[means.index != 'Others'].nlargest(15).index.tolist
    df_plot = df_norm.loc[top_genera]

    # Order columns by diet
    col_order = sorted(ds_cols, key=lambda x: (get_diet(x), get_sample(x)))
    df_plot = df_plot[col_order]

    # Clean labels
    df_plot.index = [g.replace('g__', '') for g in df_plot.index]
    df_plot.columns = [get_sample(c) for c in df_plot.columns]

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.heatmap(df_plot, cmap='YlOrRd', ax=ax, linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'Relative Abundance (%)', 'shrink': 0.6},
                yticklabels=True, xticklabels=True)
    ax.set_title('Top 15 Genera Across Samples', fontweight='bold')
    ax.set_ylabel('')

    # Add vertical line between diets
    ht_count = sum(1 for c in col_order if get_diet(c) == 'H.tessellatus')
    ax.axvline(x=ht_count, color='black', linewidth=2)

    plt.setp(ax.get_yticklabels, style='italic', fontsize=9)
    plt.setp(ax.get_xticklabels, fontsize=7, rotation=45, ha='right')

    save_fig(fig, 'Fig2D_genus_heatmap')


# ============================================================
# Supp: KEGG functional comparison
# ============================================================
def plot_supp_kegg_functional:
    df = pd.read_csv(_src('13_kegg_functional_cross_diet.tsv'), sep='\t', comment='#')

    # Filter to Level1 categories
    level1 = df[df['Level'] == 'Level1'].copy
    if len(level1) == 0:
        print("  WARNING: No Level1 KEGG data found")
        return

    # Get DS02 and DS03 columns
    ht_cols = [c for c in level1.columns if 'H.tessellatus' in c]
    po_cols = [c for c in level1.columns if 'P.ostreatus' in c]

    if not ht_cols or not po_cols:
        print("  WARNING: Cannot find diet-specific columns")
        return

    level1['HT_mean'] = level1[ht_cols].astype(float).mean(axis=1)
    level1['PO_mean'] = level1[po_cols].astype(float).mean(axis=1)

    fig, ax = plt.subplots(figsize=(9, 5))
    y = np.arange(len(level1))
    height = 0.35

    ax.barh(y + height/2, level1['HT_mean'] / 1e6, height,
            label='H. tessellatus', color='#A23B72', alpha=0.8)
    ax.barh(y - height/2, level1['PO_mean'] / 1e6, height,
            label='P. ostreatus', color='#F18F01', alpha=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels(level1['Category'].values, fontsize=9)
    ax.set_xlabel('Mean Abundance (millions)')
    ax.set_title('KEGG Level 1 Functional Categories (DS02 vs DS03)', fontweight='bold')
    ax.legend
    ax.invert_yaxis

    save_fig(fig, 'SuppFig_S18_KEGG_functional')


# ============================================================
# Supp: Genus barplot (top 10 per diet)
# ============================================================
def plot_supp_genus_barplot:
    df = pd.read_csv(_src('11_genus_abundance_cross_diet.tsv'), sep='\t', index_col=0, comment='#')
    df = harmonise_genera(df)

    # DS02 + DS03 only
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for ax, diet, title in [(ax1, 'H.tessellatus', 'H. tessellatus (Beech)'),
                             (ax2, 'P.ostreatus', 'P. ostreatus (Oyster)')]:
        cols = [c for c in df.columns if get_diet(c) == diet]
        if not cols:
            continue
        means = df[cols].mean(axis=1)
        top10 = means[means.index != 'Others'].nlargest(10)
        top10.index = [g.replace('g__', '') for g in top10.index]

        color = DIET_COLORS[diet]
        bars = ax.barh(range(len(top10)), top10.values, color=color, alpha=0.8)
        ax.set_yticks(range(len(top10)))
        ax.set_yticklabels(top10.index, fontsize=9, style='italic')
        ax.set_xlabel('Mean Abundance')
        ax.set_title(title, fontweight='bold', fontsize=10, style='italic')
        ax.invert_yaxis

    fig.suptitle('Top 10 Genera Per Diet Group',
                fontweight='bold', fontsize=11, y=1.02)
    plt.tight_layout
    save_fig(fig, 'SuppFig_S7_genus_barplot')


if __name__ == '__main__':
    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(PDF_DIR, exist_ok=True)

    print('Generating microbiome diversity plots...')
    print('\n--- Figure 2 Panels ---')
    # Fig 2A-2D are now produced from the unified pipeline by generate_unified_figs.py;
    # the vendor-table versions below are retained for provenance but not run.
    # plot_fig2a; plot_fig2b; plot_fig2c; plot_fig2d

    print('\n--- Supplementary ---')
    # plot_supp_kegg_functional   # provider-era DS02-vs-DS03 version; Supp Fig S18 is produced by generate_kegg_functional_unified32.py — never re-enable
    # plot_supp_genus_barplot   # superseded by generate_unified_figs.py (Supp Fig S8)

    print(f'\nDone! Saved to {PNG_DIR}/ and {PDF_DIR}/')
