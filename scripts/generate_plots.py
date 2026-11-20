#!/usr/bin/env python3
"""
Generate publication-quality figures from chitin metabolism analysis data.
Outputs PNG (300 dpi) and PDF for each figure.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
import seaborn as sns
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Style settings
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
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

DIET_COLORS = {
    'P_eryngii': '#2E86AB',
    'H_tessellatus': '#A23B72',
    'P_ostreatus': '#F18F01',
}
DIET_LABELS = {
    'P_eryngii': 'P. eryngii\n(King Oyster)',
    'H_tessellatus': 'H. tessellatus\n(Beech)',
    'P_ostreatus': 'P. ostreatus\n(Oyster)',
}
TISSUE_ORDER = ['Stomach', 'Duodenum', 'Jejunum', 'Ileum', 'Colon', 'Rectum']


def save_fig(fig, name):
    from figure_ownership import may_write
    if not may_write(name, __file__): return
    fig.savefig(os.path.join(PNG_DIR, f'{name}.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(PDF_DIR, f'{name}.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f'  {name}.png / .pdf')


# ============================================================
# Fig 3A: GH18 abundance across diets (boxplot + strip)
# ============================================================
def plot_fig3a:
    df = pd.read_csv(_src('01_chitin_cazyme_cross_diet.tsv'), sep='\t', comment='#')

    # Split into two panels: DS01 (TPM) separate from DS02/DS03 (absolute counts)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4.5), gridspec_kw={'width_ratios': [1, 2.2]})

    # Panel left: DS02 vs DS03 (same units - directly comparable)
    for i, ds in enumerate(['H_tessellatus', 'P_ostreatus']):
        subset = df[df['Dataset'] == ds]['GH18'].values
        bp = ax1.boxplot([subset], positions=[i], widths=0.5,
                        patch_artist=True, showfliers=False)
        bp['boxes'][0].set_facecolor(DIET_COLORS[ds])
        bp['boxes'][0].set_alpha(0.6)
        bp['medians'][0].set_color('black')
        bp['medians'][0].set_linewidth(2)
        jitter = np.random.normal(0, 0.08, len(subset))
        ax1.scatter(np.full(len(subset), i) + jitter, subset,
                   c=DIET_COLORS[ds], s=30, alpha=0.8, edgecolors='white', linewidth=0.5, zorder=3)

    ax1.set_xticks([0, 1])
    ax1.set_xticklabels([DIET_LABELS['H_tessellatus'], DIET_LABELS['P_ostreatus']], fontsize=8)
    ax1.set_ylabel('GH18 Abundance (absolute gene counts)')
    ax1.set_title('DS02 vs DS03\n(same pipeline)', fontweight='bold', fontsize=10)

    # Significance bracket
    max_val = max(df[df['Dataset'] == 'H_tessellatus']['GH18'])
    ax1.annotate('', xy=(0, max_val * 1.05), xytext=(1, max_val * 1.05),
                arrowprops=dict(arrowstyle='-', lw=1.5))
    ax1.text(0.5, max_val * 1.08, '1.7x, P = 0.0012', ha='center', fontsize=9,
             fontweight='bold')
    # headroom so the significance annotation cannot collide with the subplot title
    ax1.set_ylim(top=max_val * 1.20)

    # Panel right: DS01 (TPM - separate axis)
    subset01 = df[df['Dataset'] == 'P_eryngii']['GH18'].values
    bp = ax2.boxplot([subset01], positions=[0], widths=0.5,
                    patch_artist=True, showfliers=False)
    bp['boxes'][0].set_facecolor(DIET_COLORS['P_eryngii'])
    bp['boxes'][0].set_alpha(0.6)
    bp['medians'][0].set_color('black')
    bp['medians'][0].set_linewidth(2)
    jitter = np.random.normal(0, 0.08, len(subset01))
    ax2.scatter(np.full(len(subset01), 0) + jitter, subset01,
               c=DIET_COLORS['P_eryngii'], s=30, alpha=0.8, edgecolors='white', linewidth=0.5, zorder=3)

    ax2.set_xticks([0])
    ax2.set_xticklabels([DIET_LABELS['P_eryngii']], fontsize=8)
    ax2.set_ylabel('GH18 Abundance (gene-level TPM)')
    ax2.set_title('DS01\n(different pipeline)', fontweight='bold', fontsize=10)

    fig.suptitle('GH18 Chitinase Abundance Across Mushroom Diets', fontweight='bold', fontsize=12, y=1.02)
    plt.tight_layout
    save_fig(fig, 'Fig3A_GH18_boxplot')


# ============================================================
# Fig 3B: Chitin CAZyme heatmap
# ============================================================
def plot_fig3b:
    df = pd.read_csv(_src('01_chitin_cazyme_cross_diet.tsv'), sep='\t', comment='#')

    families = ['GH18', 'GH19', 'GH20']
    # Pivot to matrix
    matrix_data = []
    for _, row in df.iterrows:
        for fam in families:
            matrix_data.append({
                'Sample': row['Sample'],
                'Dataset': row['Dataset'],
                'Family': fam,
                'Abundance': row[fam]
            })
    mdf = pd.DataFrame(matrix_data)
    pivot = mdf.pivot_table(index='Family', columns='Sample', values='Abundance')

    # Order columns by dataset
    order = (df.sort_values(['Dataset', 'Sample'])['Sample']).tolist
    pivot = pivot[[c for c in order if c in pivot.columns]]

    fig, ax = plt.subplots(figsize=(12, 3))
    sns.heatmap(pivot, cmap='YlOrRd', ax=ax, linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'Abundance', 'shrink': 0.6})
    ax.set_title('Chitin-Degrading CAZyme Families Across Samples', fontweight='bold')
    ax.set_ylabel('')

    # Add dataset color bar at bottom
    dataset_map = df.set_index('Sample')['Dataset'].to_dict
    cols = [c for c in order if c in pivot.columns]
    for i, col in enumerate(cols):
        ds = dataset_map.get(col, '')
        color = DIET_COLORS.get(ds, 'gray')
        ax.add_patch(plt.Rectangle((i, len(families)), 1, 0.3,
                                    facecolor=color, clip_on=False))

    save_fig(fig, 'Fig3B_CAZyme_heatmap')


# ============================================================
# Fig 3C: GH18 taxonomy (horizontal bar)
# ============================================================
def plot_fig3c:
    df = pd.read_csv(_src('04_gh18_taxonomy_distribution.tsv'), sep='\t', comment='#')
    df = df[df['Genus'] != 'Unclassified_genus'].head(12)

    fig, ax = plt.subplots(figsize=(7, 5))

    y = np.arange(len(df))
    height = 0.35

    ax.barh(y + height/2, df['DS02_GH18_gene_count'], height, label='H. tessellatus', color='#A23B72', alpha=0.8)
    ax.barh(y - height/2, df['DS03_GH18_gene_count'], height, label='P. ostreatus', color='#F18F01', alpha=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels([f'  {g}' for g in df['Genus']], fontsize=9, style='italic')
    ax.set_xlabel('GH18 Gene Count')
    ax.set_title('Taxonomic Origin of GH18 Chitinase Genes', fontweight='bold')
    ax.legend(loc='lower right', framealpha=0.9)
    ax.invert_yaxis

    save_fig(fig, 'Fig3C_GH18_taxonomy')


# ============================================================
# Fig 3D: Digestibility indices
# ============================================================
def plot_fig3d:
    df = pd.read_csv(_src('07_mushroom_digestibility_comparison.tsv'), sep='\t', comment='#')

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    metrics = [
        ('Chitin_Degradation_Index', 'Chitin Degradation Index\n(GH18+19+20)/total GH', '%'),
        ('Chitin_to_Cellulose_Ratio', 'Chitin-to-Cellulose Ratio\nGH18/cellulose GH', ''),
        ('Digestibility_Index', 'Overall Digestibility Index\ncell wall enzymes/total', '%'),
    ]

    datasets = ['DS01_P.eryngii', 'DS02_H.tessellatus', 'DS03_P.ostreatus']
    ds_colors = ['#2E86AB', '#A23B72', '#F18F01']
    ds_labels = ['P. eryngii', 'H. tessellatus', 'P. ostreatus']

    for idx, (col, title, unit) in enumerate(metrics):
        ax = axes[idx]
        for i, (ds, color, label) in enumerate(zip(datasets, ds_colors, ds_labels)):
            subset = df[df['Dataset'] == ds][col].dropna
            if len(subset) > 0:
                bp = ax.boxplot([subset.values], positions=[i], widths=0.5,
                                patch_artist=True, showfliers=False)
                bp['boxes'][0].set_facecolor(color)
                bp['boxes'][0].set_alpha(0.6)
                bp['medians'][0].set_color('black')
                jitter = np.random.normal(0, 0.08, len(subset))
                ax.scatter(np.full(len(subset), i) + jitter, subset.values,
                           c=color, s=20, alpha=0.7, edgecolors='white', linewidth=0.5, zorder=3)

        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(ds_labels, fontsize=8, rotation=15, ha='right')
        ax.set_title(title, fontsize=9, fontweight='bold')
        if unit == '%':
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x*100:.1f}%'))

    fig.suptitle('Mushroom Cell Wall Digestibility Comparison', fontweight='bold', fontsize=12, y=1.02)
    plt.tight_layout
    save_fig(fig, 'Fig3D_digestibility_indices')


# ============================================================
# Fig 4B: KEGG chitin pathway enzymes (cross-diet)
# ============================================================
def plot_fig4b:
    df = pd.read_csv(_src('05_kegg_chitin_pathway.tsv'), sep='\t', comment='#')

    # Filter to key genes
    genes = ['K01183', 'K12373', 'K01443', 'K02564']
    gene_names = {'K01183': 'Chitinase\n(K01183)', 'K12373': 'HEXA/B\n(K12373)',
                  'K01443': 'nagA\n(K01443)', 'K02564': 'nagB\n(K02564)'}

    ds02_cols = [c for c in df.columns if c.startswith('DS02_')]
    ds03_cols = [c for c in df.columns if c.startswith('DS03_')]

    fig, ax = plt.subplots(figsize=(8, 5))

    x = np.arange(len(genes))
    width = 0.35

    ds02_means = []
    ds03_means = []
    ds02_sds = []
    ds03_sds = []

    for ko in genes:
        row = df[df['KO_ID'] == ko]
        if len(row) > 0:
            vals2 = row[ds02_cols].values.flatten.astype(float)
            vals3 = row[ds03_cols].values.flatten.astype(float)
            ds02_means.append(np.mean(vals2))
            ds03_means.append(np.mean(vals3))
            ds02_sds.append(np.std(vals2))
            ds03_sds.append(np.std(vals3))
        else:
            ds02_means.append(0)
            ds03_means.append(0)
            ds02_sds.append(0)
            ds03_sds.append(0)

    bars1 = ax.bar(x - width/2, ds02_means, width, yerr=ds02_sds,
                   label='H. tessellatus', color='#A23B72', alpha=0.8, capsize=3)
    bars2 = ax.bar(x + width/2, ds03_means, width, yerr=ds03_sds,
                   label='P. ostreatus', color='#F18F01', alpha=0.8, capsize=3)

    ax.set_xticks(x)
    ax.set_xticklabels([gene_names[g] for g in genes])
    ax.set_ylabel('Mean Abundance')
    ax.set_title('Chitin Degradation Pathway Enzymes (KEGG)', fontweight='bold')
    ax.legend

    # Add fold change annotations
    for i in range(len(genes)):
        if ds03_means[i] > 0:
            fc = ds02_means[i] / ds03_means[i]
            max_h = max(ds02_means[i] + ds02_sds[i], ds03_means[i] + ds03_sds[i])
            ax.text(i, max_h * 1.05, f'{fc:.1f}×', ha='center', fontsize=8, fontweight='bold', color='#333')

    save_fig(fig, 'Fig4B_KEGG_chitin_enzymes')


# ============================================================
# Fig 5D: Host chitinase tissue expression (proteomics)
# ============================================================
def plot_fig5d:
    df = pd.read_csv(_src('02_proteomics_chitin_proteins.tsv'), sep='\t', comment='#')

    # Filter to chitinases from 1st proteomics batch
    chit = df[(df['Category'] == 'Chitinase (GH18)') & (df['Dataset'] == 'Proteomics_1st_GI')].copy

    if len(chit) == 0:
        print("  WARNING: No chitinase data found for tissue plot")
        return

    # Get tissue columns
    tissue_cols = {}
    for tissue in TISSUE_ORDER:
        cols = [c for c in chit.columns if tissue in c and c != 'Tissue_Pattern']
        if cols:
            tissue_cols[tissue] = cols

    if not tissue_cols:
        print("  WARNING: No tissue columns found")
        return

    # Calculate mean per tissue per protein
    plot_data = []
    for _, row in chit.iterrows:
        pid = row['Protein_ID'].split('|')[1] if '|' in str(row['Protein_ID']) else str(row['Protein_ID'])[:15]
        for tissue in TISSUE_ORDER:
            if tissue in tissue_cols:
                vals = []
                for c in tissue_cols[tissue]:
                    try:
                        v = float(row[c])
                        if not np.isnan(v):
                            vals.append(v)
                    except (ValueError, TypeError):
                        pass
                if vals:
                    plot_data.append({'Protein': pid, 'Tissue': tissue, 'Intensity': np.mean(vals)})

    if not plot_data:
        print("  WARNING: No plottable tissue data")
        return

    pdf = pd.DataFrame(plot_data)

    # Select top 6 proteins by max intensity
    top_prots = pdf.groupby('Protein')['Intensity'].max.nlargest(6).index.tolist
    pdf_top = pdf[pdf['Protein'].isin(top_prots)]

    fig, ax = plt.subplots(figsize=(8, 5))

    palette = sns.color_palette('Set2', len(top_prots))
    for i, prot in enumerate(top_prots):
        sub = pdf_top[pdf_top['Protein'] == prot]
        tissue_means = []
        for t in TISSUE_ORDER:
            t_data = sub[sub['Tissue'] == t]['Intensity']
            tissue_means.append(t_data.mean if len(t_data) > 0 else np.nan)
        ax.plot(range(len(TISSUE_ORDER)), tissue_means, 'o-', label=prot,
                color=palette[i], linewidth=2, markersize=6)

    ax.set_xticks(range(len(TISSUE_ORDER)))
    ax.set_xticklabels(TISSUE_ORDER, fontsize=9)
    ax.set_ylabel('log2 Intensity')
    ax.set_title('Host Chitinase Expression Along the GI Tract', fontweight='bold')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7, title='Protein ID')

    # Highlight small intestine region
    ax.axvspan(0.5, 3.5, alpha=0.1, color='#27AE60', zorder=0)
    ax.text(2, ax.get_ylim[0] + 0.3, 'Small Intestine', ha='center',
            fontsize=8, color='#27AE60', fontweight='bold')

    plt.tight_layout
    save_fig(fig, 'Fig7B_tissue_chitinase_expression')


# ============================================================
# Fig 5C: Chitin protein categories (bar chart)
# ============================================================
def plot_fig5c:
    categories = {
        'Chitinase\n(GH18)': 54,
        'Hexosaminidase\n(GH20)': 9,
        'GlcNAc\nmetabolism': 60,
        'Other': 1,
    }

    fig, ax = plt.subplots(figsize=(5, 4))
    colors = ['#2E86AB', '#F18F01', '#27AE60', '#95A5A6']
    bars = ax.bar(categories.keys, categories.values, color=colors, alpha=0.85, edgecolor='white', linewidth=1.5)

    for bar, val in zip(bars, categories.values):
        ax.text(bar.get_x + bar.get_width/2, bar.get_height + 1,
                str(val), ha='center', fontweight='bold', fontsize=11)

    ax.set_ylabel('Number of Proteins')
    ax.set_title('Chitin-Related Proteins Detected\n(DIA Proteomics, n=124 total)', fontweight='bold')
    ax.set_ylim(0, max(categories.values) * 1.15)

    save_fig(fig, 'Fig7A_chitin_protein_categories')


# ============================================================
# Fig 5A: Chr3 gene cluster diagram
# ============================================================
def plot_fig5a:
    """chr3 GH18 locus map: the three annotated gene models plus the five loci
    recovered by genome-level tBLASTn/miniprot (results/32_gh18_locus_search_record.tsv).
    Peptide support = proteotypic peptides of the species-specific DIA-NN search per
    Mimp_v3 model (tables/62_diann_peptides_per_locus.tsv), mapped onto the v1 loci via
    tables/45 (088402 = g08469 + g08470; 088393 = g08478 (+ g08477); L4-L6 share one group)."""
    loci = pd.read_csv(os.path.join(BASE, 'tables', '32_gh18_locus_search_record.tsv'), sep='\t', comment='#')
    # draw in Mimp_v3 chr4 coordinates (the reference of the paper); tables/45 maps every locus between versions
    v3 = pd.read_csv(os.path.join(BASE, 'tables', '45_gh18_locus_Mimp_v3_map.tsv'), sep='\t', comment='#').set_index('locus')
    loci['start'] = [int(v3.loc[l, 'Mimp_v3_gene_start']) for l in loci['locus']]
    loci['end'] = [int(v3.loc[l, 'Mimp_v3_gene_end']) for l in loci['locus']]
    loci['strand'] = ['-' if str(v3.loc[l, 'Mimp_v3_strand']) == '-' else '+' for l in loci['locus']]
    pep = pd.read_csv(os.path.join(BASE, 'tables', '62_diann_peptides_per_locus.tsv'), sep='\t', comment='#').set_index('model').peptides_in_own_group
    v1_to_v3 = {'088392': ['MIMP3_g08479'], '088393': ['MIMP3_g08478', 'MIMP3_g08477'], 'L3': ['MIMP3_g08476'], 'L4': ['MIMP3_g08475'], 'L5': ['MIMP3_g08474'],
                'L6': ['MIMP3_g08473'], 'L7': ['MIMP3_g08472'], '088402': ['MIMP3_g08469', 'MIMP3_g08470']}
    loci['unique_peptides'] = [sum(int(pep.get(m, 0)) for m in v1_to_v3[str(l)]) for l in loci['locus']]
    fig, ax = plt.subplots(figsize=(12, 3.6))
    chr_start, chr_end = 139.04, 139.27
    ax.plot([chr_start, chr_end], [0, 0], color='#CCC', linewidth=8, solid_capstyle='round', zorder=1)
    colors = {'annotated': '#2E86AB', 'unannotated': '#F39C12'}
    for i, r in loci.iterrows:
        s_, e_ = r['start'] / 1e6, r['end'] / 1e6
        y = 0.16 if r['strand'] == '+' else -0.16
        ax.add_patch(plt.Rectangle((s_, y - 0.12), e_ - s_, 0.24, facecolor=colors[r['status']],
                                   edgecolor='white', linewidth=1, zorder=2))
        ax.annotate('', xy=(e_ if r['strand'] == '+' else s_, y), xytext=(s_ if r['strand'] == '+' else e_, y),
                    arrowprops=dict(arrowstyle='->', color='white', lw=1.2), zorder=3)
        lab = f"{r['locus']}\n{r['clade']}"
        if r['catalytic_motif'] != 'intact':
            lab += f"\n({r['catalytic_motif']})"
        if r['unique_peptides'] > 0:
            lab += f"\n{int(r['unique_peptides'])} peptides" + (' (shared L4–L6)' if r['locus'] in ('L4', 'L5', 'L6') else '')
        if str(r['locus']) == '088402':
            lab = lab.replace(f"{r['clade']}", '088402a/b: two GH18 genes (g08469, g08470)')
        ax.text((s_ + e_) / 2, 0.42 if i % 2 == 0 else 0.95, lab, ha='center', va='bottom', fontsize=6.6)
    ax.set_xlim(chr_start - 0.005, chr_end + 0.005); ax.set_ylim(-0.6, 1.55)
    ax.set_xlabel('Mimp_v3 chromosome 4 position (Mb)')
    ax.set_title('GH18 chitinase locus on Mimp_v3 chromosome 4: 3 annotated gene models + 5 loci recovered by genome search; peptides = species-specific DIA-NN search',
                 fontweight='bold', fontsize=10)
    ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.plot([139.05, 139.07], [-0.5, -0.5], 'k-', linewidth=2); ax.text(139.06, -0.57, '20 kb', ha='center', fontsize=8)
    ax.text(139.268, -0.45, '+ strand above / − strand below the axis (locus on the minus strand)', ha='right', fontsize=7, color='#555')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=colors['annotated'], label='annotated in the published gene set'),
                       Patch(color=colors['unannotated'], label='recovered by genome search (tBLASTn + miniprot)')],
              loc='upper left', fontsize=7, frameon=False)
    save_fig(fig, 'Fig5A_gh18_locus_map')


# ============================================================
# RNA-Seq expression barplot
# ============================================================
def plot_rnaseq_expression:
    df = pd.read_csv(_src('09_rnaseq_chitinase_expression.tsv'), sep='\t', comment='#')

    # Filter to key genes
    key_genes = df[df['Category'].isin(['Chitinase_chr3', 'GlcNAc_pathway', 'Other_chitin'])].copy
    key_genes = key_genes[key_genes['IT3_TPM'] + key_genes['IT4_TPM'] > 1]  # filter very low

    fig, ax = plt.subplots(figsize=(10, 5))

    x = np.arange(len(key_genes))
    width = 0.35

    bars1 = ax.bar(x - width/2, key_genes['IT3_TPM'], width, label='IT3 (heart)', color='#2E86AB', alpha=0.8)
    bars2 = ax.bar(x + width/2, key_genes['IT4_TPM'], width, label='IT4 (lung)', color='#F18F01', alpha=0.8)

    ax.set_xticks(x)
    nice = {'Chitinase_chr3_1': '088392 (L1, CHIT1-clade)', 'Chitinase_chr3_2': '088393 (L2, CLP)', 'Chitinase_chr3_3': '088402 (the two GH18 genes)'}
    ax.set_xticklabels([nice.get(g, g) for g in key_genes['Gene_Name']], fontsize=8, rotation=30, ha='right')
    ax.set_ylabel('TPM')
    ax.set_title('Supplementary Figure S10. Host chitin-related gene expression (RNA-Seq, heart and lung)', fontweight='bold', fontsize=11)
    star = plt.Line2D([], [], ls='', marker='$*$', color='#E74C3C', ms=11, label='TPM > 100 in either tissue')
    ax.legend(handles=[bars1, bars2, star], fontsize=8, frameon=False)
    ax.spines[['top', 'right']].set_visible(False)

    # Highlight the highly expressed gene
    for i, (_, row) in enumerate(key_genes.iterrows):
        if row['IT3_TPM'] > 100 or row['IT4_TPM'] > 100:
            max_val = max(row['IT3_TPM'], row['IT4_TPM'])
            ax.annotate('*', xy=(i, max_val), xytext=(0, 2), textcoords='offset points', fontsize=14, ha='center', color='#E74C3C', fontweight='bold')

    plt.tight_layout
    save_fig(fig, 'SuppFig_S10_RNASeq_expression')


# ============================================================
# Fig 6: Integrated model (simplified schematic)
# ============================================================
def plot_fig7_model:
    """Figure 7A - the proposed pathway. Each enzyme annotation is joined by a leader line to the
    reaction arrow it acts on; the annotation column used to overrun the right edge of the canvas."""
    fig, ax = plt.subplots(figsize=(12, 5.8))   # trimmed so Figure 9 clears the 225 mm limit with its legend
    ax.set_xlim(0, 13); ax.set_ylim(0, 10); ax.axis('off')
    ax.set_title('Integrated chitin degradation model in $\\mathit{M.\\,impressa}$', fontweight='bold', fontsize=13, pad=16)

    X, W, H = 1.8, 4.4, 0.8
    boxes = [
        ('Fungal cell wall chitin\n(β-1,4-GlcNAc polymer)', '#8B4513'),
        ('Chitooligosaccharides', '#2E86AB'),
        ('Chitobiose $\\mathrm{(GlcNAc)_2}$', '#2E86AB'),
        ('N-acetylglucosamine (GlcNAc)', '#27AE60'),
        ('Fructose-6-P → glycolysis', '#F39C12'),
    ]
    ys = [8.7 - 1.75 * i for i in range(len(boxes))]
    for (txt, col), y in zip(boxes, ys):
        ax.add_patch(plt.Rectangle((X, y), W, H, facecolor=col, edgecolor='white', lw=2, alpha=0.88, zorder=2))
        ax.text(X + W / 2, y + H / 2, txt, ha='center', va='center', fontsize=9.5, color='white', fontweight='bold', zorder=3)
    mids = []
    for y0, y1 in zip(ys[:-1], ys[1:]):
        ax.annotate('', xy=(X + W / 2, y1 + H), xytext=(X + W / 2, y0), arrowprops=dict(arrowstyle='->', lw=2.2, color='#555'))
        mids.append((y1 + H + y0) / 2)

    BOTH, MICR = '#C0392B', '#2E86AB'
    enzymes = [
        ('Host GH18 (chr4 locus, 8 loci / 9 genes: CHIT1-, CHIA-array,\nand CHI3L-clade loci) + microbial GH18\n(643 GH18 genes in the catalogue;\nBacteroides the largest classified contributor)', BOTH),
        ('Microbial GH18\n(exochitinase)', MICR),
        ('GH20 hexosaminidase\n(K12373, most abundant)', MICR),
        ('nagA → nagB\n(microbial + host GNPDA1)', BOTH),
    ]
    for (txt, col), ym in zip(enzymes, mids):
        ax.annotate('', xy=(X + W / 2 + 0.2, ym), xytext=(6.75, ym),
                    arrowprops=dict(arrowstyle='-', color=col, lw=1.1, shrinkA=0, shrinkB=0))
        ax.text(6.85, ym, txt, fontsize=8, ha='left', va='center', color=col, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=col, alpha=0.95))

    # clear of the bottom metabolite box (x 1.8-6.2, y 1.7-2.5): "host + microbiome" used to run under it
    ax.text(0.2, 1.25, 'Enzyme source', fontsize=9, fontweight='bold', color='#333')
    for k, (col, lab) in enumerate([(MICR, 'microbiome'), (BOTH, 'host + microbiome')]):
        ax.plot([0.2, 0.75], [0.88 - 0.34 * k] * 2, color=col, lw=3.5)
        ax.text(0.85, 0.88 - 0.34 * k, lab, fontsize=8, color=col, va='center')
    ax.text(6.85, 1.5, 'Diet-associated (one catalogue, dbCAN3):\nGH18 pool lowest on P. ostreatus\n(2.1× below H. tessellatus, both animals)',
            fontsize=8, color='#555', style='italic', va='center',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#F5F5F5', edgecolor='#DDD'))
    save_fig(fig, 'Fig9A_integrated_model')


# ============================================================
# Supp Fig S6: QC and Assembly Summary Table
# ============================================================
def plot_supp_qc:
    """Supp Fig S6: per-sample QC / host removal / assembly / gene prediction for all 32
    samples under the unified pipeline (tables/42_unified_sample_summary.tsv)."""
    q = pd.read_csv(os.path.join(BASE, 'tables', '42_unified_sample_summary.tsv'), sep='\t', comment='#')
    qc = pd.read_csv(os.path.join(BASE, 'results', 'unified', 'qc_assembly_32.tsv'), sep='\t', comment='#')
    q = q.merge(qc.drop(columns=['Cohort']), on='Sample', how='left')   # per-sample QC / assembly statistics (unified pipeline)
    q = q.sort_values(['Cohort', 'Sample']).reset_index(drop=True)
    q['Host_aligned'] = q['Host_aligned'].fillna('n.d.')
    cols = ['Cohort', 'Sample', 'Raw_Gb', 'Clean_Gb', 'Q30_pct', 'Host_aligned', 'nonhost_read_pairs', 'Contigs',
            'Assembly_Mb', 'Contig_N50_bp', 'Predicted_genes', 'catalogue_mapping_rate']
    lab = ['Cohort', 'Sample', 'Raw (Gb)', 'Clean (Gb)', 'Q30 (%)', 'Host (%)', 'Non-host\npairs (M)', 'Contigs',
           'Assembly (Mb)', 'N50 (bp)', 'Genes', 'Mapped (%)']
    d = q[cols].copy
    # Table 1 / Supp Table S3 source: the merged per-sample overview (QC + assembly + mapping + Bracken genus-assigned fraction)
    t1 = q.copy; t1['bracken_genus_assigned_pct'] = (100 * t1['bracken_genus_assigned_reads'] / t1['nonhost_read_pairs']).round(2)
    with open(os.path.join(BASE, 'tables', '72_unified_sample_overview.tsv'), 'w') as fh:
        fh.write('# Per-sample overview of the unified pipeline, all 32 samples: results/unified/qc_assembly_32.tsv merged with tables/42_unified_sample_summary.tsv (Table 1 / Supplementary Table S3 / Supplementary Figure S6). scripts/generate_plots.py\n')
        t1.to_csv(fh, sep='\t', index=False)
    d['nonhost_read_pairs'] = (d['nonhost_read_pairs'] / 1e6).map('{:.1f}'.format)
    # one precision per column, and the % lives in the header, not in every cell
    for c in ('Raw_Gb', 'Clean_Gb'): d[c] = d[c].map('{:.2f}'.format)
    d['Q30_pct'] = d['Q30_pct'].map('{:.2f}'.format)
    for c in ('Host_aligned', 'catalogue_mapping_rate'): d[c] = d[c].astype(str).str.rstrip('%')
    d['Contigs'] = d['Contigs'].map('{:,}'.format); d['Predicted_genes'] = d['Predicted_genes'].map('{:,}'.format)
    d['Contig_N50_bp'] = d['Contig_N50_bp'].map('{:,}'.format)
    fig, ax = plt.subplots(figsize=(13, 10)); ax.axis('off')
    ax.set_title('Supplementary Figure S5. Sequencing QC, host removal, assembly and gene prediction — all 32 samples, one pipeline',
                 fontweight='bold', fontsize=11, pad=16)
    ax.text(0.5, -0.015, 'n.d. = host alignment not recorded for this sample; all percentages are of the column header quantity.',
            transform=ax.transAxes, ha='center', fontsize=7.5, color='#555')
    table = ax.table(cellText=d.values, colLabels=lab, cellLoc='center', loc='center')
    table.auto_set_font_size(False); table.set_fontsize(7.5); table.scale(1, 1.25)
    cohort_col = {'DS01': '#E3F0F7', 'DS02': '#F5E3EC', 'DS03': '#FDF0E0'}
    for (i, j), cell in table.get_celld.items:
        if i == 0:
            cell.set_facecolor('#1B3A5C'); cell.set_text_props(color='white', fontweight='bold')
        else:
            cell.set_facecolor(cohort_col[d.iloc[i - 1, 0]])
    save_fig(fig, 'SuppFig_S5_QC_assembly')


# ============================================================
# Supp Fig: Proteomics protein count overview
# ============================================================
def plot_supp_proteomics_overview:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Batch 1: proteins per tissue
    tissues_b1 = ['Stomach', 'Duodenum', 'Jejunum', 'Ileum', 'Colon', 'Rectum']
    counts_b1 = [7140, 7664, 7743, 7467, 7538, 6960]
    bars1 = ax1.bar(tissues_b1, counts_b1, color=sns.color_palette('Blues_d', 6), alpha=0.8, edgecolor='white')
    ax1.set_ylabel('Identified Proteins')
    ax1.set_title('Batch 1: Proteins per GI Segment\n(3 replicates each)', fontweight='bold', fontsize=10)
    ax1.set_xticklabels(tissues_b1, rotation=30, ha='right', fontsize=8)
    for bar, val in zip(bars1, counts_b1):
        ax1.text(bar.get_x + bar.get_width/2, bar.get_height + 50,
                f'{val:,}', ha='center', fontsize=7)

    # Batch 2: proteins per tissue
    tissues_b2 = ['Liver', 'Heart', 'Stomach', 'Duodenum', 'Jejunum', 'Ileum', 'Colon', 'Rectum', 'Muscle']
    counts_b2 = [8979, 8233, 7943, 6303, 8812, 8846, 7210, 8015, 6285]
    colors_b2 = ['#E74C3C', '#E74C3C'] + ['#3498DB']*6 + ['#E74C3C']
    bars2 = ax2.bar(tissues_b2, counts_b2, color=colors_b2, alpha=0.7, edgecolor='white')
    ax2.set_ylabel('Identified Proteins')
    ax2.set_title('Batch 2: Proteins per Tissue\n(2 replicates each)', fontweight='bold', fontsize=10)
    ax2.set_xticklabels(tissues_b2, rotation=30, ha='right', fontsize=8)
    for bar, val in zip(bars2, counts_b2):
        ax2.text(bar.get_x + bar.get_width/2, bar.get_height + 50,
                f'{val:,}', ha='center', fontsize=7)

    fig.suptitle('Proteomics Protein Identification Summary',
                fontweight='bold', fontsize=11, y=1.02)
    plt.tight_layout
    save_fig(fig, 'SuppFig_S14_proteomics_overview')


# ============================================================
# Supp Fig: Host genome overview
# ============================================================
def plot_supp_genome:
    fig, ax = plt.subplots(figsize=(13.5, 4.8))
    ax.axis('off')
    ax.set_title('Supplementary Figure S6. $\\mathit{M.\\,impressa}$ reference genome Mimp_v3 — summary', fontweight='bold', fontsize=11, pad=20)

    col_labels = ['Metric', 'Value']
    data = [
        ['Reference used here', 'Mimp_v3 = GWHHIVU00000000.2 (BioProject PRJCA038476)'],
        ['Assembly size', '2,417,402,780 bp (2.42 Gb); 361 sequences; 148 gaps'],
        ['Chromosomes', '26 (2n = 52), 94.5% of the sequence; scaffold N50 138.6 Mb, contig N50 86.9 Mb'],
        ['Quality', 'Merqury QV 67.9, k-mer completeness 95.1%; BUSCO vertebrata 99.3% complete'],
        ['Haplotypes', 'hap1 GWHKKMH00000000.1 (2.37 Gb, QV 69.1); hap2 GWHKKMI00000000.1 (2.38 Gb, QV 68.5)'],
        ['Protein-coding genes', '24,875 (29,927 transcripts; 97.6% complete BUSCOs in protein mode; 89.0% eggNOG-annotated)'],
        ['Repeats', '48.6% masked'],
        ['Sequencing', 'PacBio HiFi 28x (68.5 Gb) + Hi-C (9.8 Gb) + RNA-seq heart/lung; one female (MI01)'],
        ['GH18 chitinase locus', 'chr4:139.05-139.26 Mb (minus strand), 8 loci / 9 complete GH18 genes, MIMP3_g08469-g08479 (g08468 = TAF8)'],
    ]

    table = ax.table(cellText=data, colLabels=col_labels, cellLoc='left', loc='center',
                     colWidths=[0.3, 0.7])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    for (i, j), cell in table.get_celld.items:
        if i == 0:
            cell.set_facecolor('#1B3A5C')
            cell.set_text_props(color='white', fontweight='bold')
        elif i % 2 == 0:
            cell.set_facecolor('#F1F5F9')
        if j == 0 and i > 0:
            cell.set_text_props(fontweight='bold')
    # two-line values used to spill over the rules above and below: give those rows the height they need
    for r, row in enumerate(data, start=1):
        lines = max(v.count('\n') for v in row) + 1
        if lines > 1:
            for c in range(len(col_labels)):
                table[(r, c)].set_height(table[(r, c)].get_height * lines)

    save_fig(fig, 'SuppFig_S6_genome_summary')


# ============================================================
# Supp Fig: Host chitinase chr3 + RNA-Seq combined
# ============================================================
def plot_supp_host_combined:
    df = pd.read_csv(_src('09_rnaseq_chitinase_expression.tsv'), sep='\t', comment='#')
    chr3 = df[df['Category'] == 'Chitinase_chr3'].copy

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5), gridspec_kw={'width_ratios': [1.5, 1]})

    # Left: Chr3 cluster with expression
    # read from the tables rather than hardcoding: a pasted number cannot drift from its source
    _rna = pd.read_csv(os.path.join(BASE, 'tables', '09_rnaseq_chitinase_expression.tsv'), sep='\t', comment='#').set_index('Gene_ID')
    _bl = pd.read_csv(os.path.join(BASE, 'tables', '25_CHIA_CHIT1_blast_verification.tsv'), sep='\t', comment='#')
    _bl = _bl[_bl.Query.str.contains('CH3L1_HUMAN', na=False)].set_index('Subject')
    _MODELS = [('maker_gene_088392', 'GWHPHIVU010568.1', '088392 / L1', '571 aa'),
               ('maker_gene_088393', 'GWHPHIVU010569.1', '088393 / L2', '717 aa'),
               ('maker_gene_088402', 'GWHPHIVU010570.1', '088402', '1,160 aa; fused model')]
    genes = [{'name': f'{lab}\n({aa})', 'pos': i,
              'tpm_it3': float(_rna.IT3_TPM[gid]), 'tpm_it4': float(_rna.IT4_TPM[gid]),
              'ident': float(_bl.Pident[sub])} for i, (gid, sub, lab, aa) in enumerate(_MODELS)]
    x = [0, 1, 2]
    width = 0.35
    ax1.bar([p - width/2 for p in x], [g['tpm_it3'] for g in genes], width,
            label='IT3 (heart)', color='#2E86AB', alpha=0.8)
    ax1.bar([p + width/2 for p in x], [g['tpm_it4'] for g in genes], width,
            label='IT4 (lung)', color='#F18F01', alpha=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels([g['name'] for g in genes], fontsize=8)
    ax1.set_ylabel('TPM')
    ax1.set_title('A  Heart / lung RNA-Seq (log scale)', loc='left', fontweight='bold', fontsize=10)
    ax1.legend
    ax1.set_yscale('log')
    ax1.set_ylim(1, 500)

    # Right: Comparison table
    ax2.axis('off')
    col_labels = ['Gene', 'Identity to\nhuman CHI3L1', 'TPM\nIT3 (heart)', 'TPM\nIT4 (lung)']
    data = [[lab, f"{g['ident']:.1f}%", f"{g['tpm_it3']:.2f}", f"{g['tpm_it4']:.2f}"]
            for g, (_, _, lab, _) in zip(genes, _MODELS)]
    table = ax2.table(cellText=data, colLabels=col_labels, cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)
    for (i, j), cell in table.get_celld.items:
        if i == 0:
            cell.set_facecolor('#1B3A5C')
            cell.set_text_props(color='white', fontweight='bold')
        elif i == 2:  # Highlight the highly expressed gene
            cell.set_facecolor('#FFF3E0')

    ax2.set_title('B  Identity to human CHI3L1 and TPM', loc='left', fontweight='bold', fontsize=10, y=0.86)
    fig.suptitle('Supplementary Figure S9. Expression of the three annotated GH18-locus gene models',
                fontweight='bold', fontsize=11, y=1.02)
    plt.tight_layout
    save_fig(fig, 'SuppFig_S9_host_chitinase_detail')



def plot_fig4c:
    """Figure 4B - the cascade as ONE vertical column of metabolites, each step annotated to the
    right and joined to its own reaction arrow by a leader line. The earlier L-shaped layout left
    the enzyme boxes floating in white space with no visual link to the step they catalyse."""
    fig, ax = plt.subplots(figsize=(12, 7.6))
    ax.set_xlim(0, 12); ax.set_ylim(0, 10); ax.axis('off')
    ax.set_title('Chitin degradation cascade: enzymes, abundances and sources', fontweight='bold', fontsize=12, pad=12)

    # metabolite column — fills encode position in the chain, never enzyme source (red is reserved for the key)
    nodes = [
        ('Fungal cell wall\nCHITIN', '#8B4513'), ('Chitin\noligosaccharides', '#2E86AB'),
        ('Chitobiose $\\mathrm{(GlcNAc)_2}$', '#2E86AB'), ('GlcNAc', '#27AE60'),
        ('GlcNAc-6-P', '#27AE60'), ('GlcN-6-P', '#27AE60'), ('Fru-6-P\n(glycolysis)', '#F39C12'),
    ]
    X, W, H = 0.5, 2.8, 0.78
    ys = [9.05 - 1.33 * i for i in range(len(nodes))]
    for (txt, col), y in zip(nodes, ys):
        ax.add_patch(plt.Rectangle((X, y), W, H, facecolor=col, edgecolor='white', lw=2, alpha=0.9, zorder=2))
        ax.text(X + W / 2, y + H / 2, txt, ha='center', va='center', fontsize=9, color='white', fontweight='bold', zorder=3)
    mids = []
    for y0, y1 in zip(ys[:-1], ys[1:]):
        ax.annotate('', xy=(X + W / 2, y1 + H), xytext=(X + W / 2, y0), arrowprops=dict(arrowstyle='->', lw=2.4, color='#555'))
        mids.append((y1 + H + y0) / 2)

    ko = {}
    for l in open(os.path.join(BASE, 'tables', '46_unified_ko_tests.tsv')):
        if l.startswith('#') or l.startswith('KO'): continue
        f = l.rstrip('\n').split('\t'); ko[f[0]] = tuple(round(float(x)) for x in f[2:5])
    genes = {}
    for l in open(os.path.join(BASE, 'tables', '27_unified32_kegg_functional_summary.tsv')):
        f = l.rstrip('\n').split('\t')
        if len(f) >= 2 and f[0].startswith('genes_K'): genes[f[0][6:]] = int(f[1])
    cpm = lambda k: 'CPM %d / %d / %d' % ko[k]
    BOTH, MICR = '#C0392B', '#2E86AB'
    steps = [
        ('Endochitinase (GH18)\nEC 3.2.1.14 / K01183', BOTH,
         'Host GH18 locus (8 loci) + microbial GH18\n(643 dbCAN3 GH18 genes; Bacteroides largest share)\nK01183: %d genes, %s' % (genes['K01183'], cpm('K01183'))),
        ('Exochitinase (GH18)\nchitobiosidase', MICR,
         'Primarily microbial\n(Bacteroides, Butyrivibrio)'),
        ('β-N-acetylhexosaminidase (GH20)\nEC 3.2.1.52 / K12373', BOTH,
         'Most abundant chitin-pathway enzyme\nK12373: %d genes, %s\n+ host HEXD (chr15)' % (genes['K12373'], cpm('K12373'))),
        ('NAGK (K00884) / PTS nagE (K02804)', MICR,
         'GlcNAc kinase rare (K00884: %d genes, %s);\nphosphorylation mainly at import (nagE: %d genes)' % (genes['K00884'], cpm('K00884'), genes['K02804'])),
        ('nagA (K01443)\nGlcNAc-6-P deacetylase', BOTH,
         'K01443: %d genes, %s\n+ host NAGA (chr1)' % (genes['K01443'], cpm('K01443'))),
        ('nagB (K02564)\nGlcN-6-P deaminase', BOTH,
         'K02564: %d genes, %s\n+ host GNPDA1 (chr8, TPM = 395)' % (genes['K02564'], cpm('K02564'))),
    ]
    for (name, col, detail), ym in zip(steps, mids):
        ax.annotate('', xy=(X + W / 2 + 0.22, ym), xytext=(3.95, ym),      # leader to the arrow it catalyses
                    arrowprops=dict(arrowstyle='-', color=col, lw=1.1, shrinkA=0, shrinkB=0))
        ax.text(4.05, ym, name, fontsize=8, fontweight='bold', color=col, ha='left', va='center',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor=col, alpha=0.95))
        ax.text(7.45, ym, detail, fontsize=7, color='#555', ha='left', va='center', style='italic')
    ax.text(0.5, 0.28, 'CPM = cohort mean per million mapped reads, P. eryngii / H. tessellatus / P. ostreatus; '
            'gene counts from the single 32-sample catalogue (dbCAN3 / eggNOG-KEGG)', fontsize=7, color='#555', style='italic')
    ax.text(7.45, 9.85, 'Enzyme source', fontsize=8, fontweight='bold', color='#333')
    for k, (col, lab) in enumerate([(MICR, 'microbiome'), (BOTH, 'host + microbiome')]):
        ax.plot([7.45, 7.95], [9.55 - 0.3 * k] * 2, color=col, lw=4)
        ax.text(8.05, 9.55 - 0.3 * k, lab, fontsize=8, color=col, va='center', fontweight='bold')
    save_fig(fig, 'Fig4B_chitin_cascade')


# ============================================================
# Entry point
# ============================================================
if __name__ == '__main__':
    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(PDF_DIR, exist_ok=True)
    print('Generating publication-quality figures...')
    print('\n--- Main Figures ---')
    # Fig 3A, 3B, 3D are now produced from the unified pipeline by generate_unified_figs.py
    # plot_fig3a; plot_fig3b; plot_fig3d
    # plot_fig3c   # vendor two-cohort version — Fig 3C is produced by scripts/mag_gh18_analysis.py (genome-resolved, three diets); never re-enable
    # plot_fig4b   # vendor-table version superseded by Fig 4A in generate_unified_figs.py
    plot_fig5a
    # plot_fig5c; plot_fig5d   # proxy-database versions superseded by scripts/proteomics_diann_chitin.py (2026-11-11)
    plot_fig7_model
    plot_fig4c

    print('\n--- Supplementary Figures ---')
    plot_rnaseq_expression
    plot_supp_qc
    # plot_supp_proteomics_overview   # superseded by scripts/proteomics_diann_chitin.py (Supp S14)
    plot_supp_genome
    plot_supp_host_combined
    print(f'\nAll figures saved to:\n  {PNG_DIR}/\n  {PDF_DIR}/')
