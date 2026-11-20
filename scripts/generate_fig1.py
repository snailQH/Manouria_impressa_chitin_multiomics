#!/usr/bin/env python3
"""
Generate Figure 1 — Study design and multi-omics data overview.

Panel A: M. impressa and the three dietary mushroom species (photographic plate)
Panel B: Multi-omics experimental design and analysis workflow
Panel C: Data generation summary (sample counts and sequencing depth)

Source of panel A image: docs/initial_design/figure-x.jpg
  NOTE: provenance of this photographic plate must be confirmed before submission
  (see FIG1_NOTES.md). Swap PANEL_A_IMAGE below to replace it.

Outputs PNG (300 dpi) and PDF, matching scripts/generate_plots.py conventions.
"""

import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# Style settings — identical to scripts/generate_plots.py
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
PNG_DIR = os.path.join(BASE, 'figures', 'png')
PDF_DIR = os.path.join(BASE, 'figures', 'pdf')
PANEL_A_IMAGE = os.path.join(BASE, 'docs', 'initial_design', 'figure-x.jpg')

DIET_COLORS = {
    'P_eryngii': '#2E86AB',
    'H_tessellatus': '#A23B72',
    'P_ostreatus': '#F18F01',
}

# ---------------------------------------------------------------------------
# Data — every number below is traceable to a source file or server log.
# ---------------------------------------------------------------------------

# Metagenomics cohorts
METAGENOMICS = [
    # label, n, colour key, mean Gbp/sample (R1+R2), source
    ('$\\mathit{P.\\,eryngii}$\n(King Oyster)', 12, 'P_eryngii', 22.9,
     'unified pipeline fastp logs'),
    ('$\\mathit{H.\\,tessellatus}$\n(Beech)', 13, 'H_tessellatus', 7.6,
     'unified pipeline fastp logs'),
    ('$\\mathit{P.\\,ostreatus}$\n(Oyster)', 7, 'P_ostreatus', 7.9,
     'unified pipeline fastp logs'),
]

# Other omics layers
OTHER_LAYERS = [
    ('Proteomics I\n6 GI segments', 18, '#55707F'),
    ('Proteomics II\n9 tissues', 18, '#55707F'),
    ('RNA-Seq\nheart, lung', 2, '#55707F'),
    ('Reference\ngenome', 1, '#55707F'),
]


def save_fig(fig, name):
    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(PDF_DIR, exist_ok=True)
    fig.savefig(os.path.join(PNG_DIR, f'{name}.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(PDF_DIR, f'{name}.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f'  {name}.png / .pdf')


def _box(ax, x, y, w, h, text, facecolor, fontsize=8.5, fontweight='normal'):
    """Draw a rounded workflow box centred on (x, y)."""
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle='round,pad=0.012,rounding_size=0.02',
        facecolor=facecolor, edgecolor='#4A4A4A', linewidth=0.9,
    ))
    ax.text(x, y, text, ha='center', va='center',
            fontsize=fontsize, fontweight=fontweight, zorder=5)


def _arrow(ax, x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle='-|>', mutation_scale=11,
        color='#6B6B6B', linewidth=1.0,
        connectionstyle='arc3,rad=0.0', zorder=1,
    ))


def panel_a(ax):
    """Photographic plate: study species and the three dietary mushrooms."""
    if os.path.exists(PANEL_A_IMAGE):
        ax.imshow(mpimg.imread(PANEL_A_IMAGE))
    else:
        ax.text(0.5, 0.5, 'Panel A image not found:\n' + PANEL_A_IMAGE,
                ha='center', va='center', fontsize=9, color='#B00020')
    ax.axis('off')
    ax.set_title('Study species and dietary fungi', loc='left',
                 fontsize=11, fontweight='bold', pad=6)


def panel_b(ax):
    """Multi-omics experimental design and analysis workflow."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title('Multi-omics experimental design', loc='left',
                 fontsize=11, fontweight='bold', pad=6)

    host, meta, integ = '#DCE9F5', '#F7E3D3', '#E8E2F2'

    # Row 1 — data generation
    _box(ax, 0.16, 0.90, 0.28, 0.115,
         'Host genome Mimp_v3\n$\\it{Manouria\\ impressa}$\n2.42 Gb, 26 chr, QV 67.9', host)
    _box(ax, 0.47, 0.90, 0.20, 0.115, 'RNA-Seq\n(heart, lung)', host)
    _box(ax, 0.81, 0.90, 0.30, 0.115,
         'Gut metagenome\n3 fungal diets, n = 32 (2 animals)', meta)

    # Row 2 — processing
    _box(ax, 0.30, 0.66, 0.32, 0.105,
         'Host genome annotation\n24,875 genes, BUSCO 99.3%', host)
    _box(ax, 0.81, 0.66, 0.30, 0.105,
         'Microbial community\nand gene catalogue', meta)

    _arrow(ax, 0.16, 0.842, 0.24, 0.715)
    _arrow(ax, 0.47, 0.842, 0.36, 0.715)
    _arrow(ax, 0.81, 0.842, 0.81, 0.715)

    # Row 3 — integration
    _box(ax, 0.50, 0.40, 0.66, 0.135,
         'Chitin metabolism genes, enzymes and pathways\n'
         'host / microbe origin assignment', integ, fontsize=9)

    _arrow(ax, 0.30, 0.607, 0.38, 0.470)
    _arrow(ax, 0.81, 0.607, 0.64, 0.470)

    # Row 4 — validation
    _box(ax, 0.50, 0.14, 0.50, 0.125,
         'Validation by gut proteomics (DIA-NN vs Mimp_v3)\n'
         '2 animals, n = 36 (6 GI segments x 3; 9 tissues)', '#DDEEDD', fontsize=9)
    _arrow(ax, 0.50, 0.333, 0.50, 0.205)


def panel_c(ax):
    """Data generation summary: sample counts per omics layer."""
    labels, counts, colors = [], [], []
    for lab, n, key, _depth, _src in METAGENOMICS:
        labels.append(lab)
        counts.append(n)
        colors.append(DIET_COLORS[key])
    for lab, n, col in OTHER_LAYERS:
        labels.append(lab)
        counts.append(n)
        colors.append(col)

    bars = ax.bar(range(len(counts)), counts, color=colors,
                  edgecolor='#3A3A3A', linewidth=0.7, width=0.68)
    for bar, n in zip(bars, counts):
        ax.text(bar.get_x + bar.get_width / 2, n + 0.4, str(n),
                ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel('Samples (n)')
    ax.set_ylim(0, max(counts) * 1.25)
    ax.set_title('Data generation summary', loc='left',
                 fontsize=11, fontweight='bold', pad=6)

    # Separate metagenomics from the other omics layers
    ax.axvline(2.5, color='#9A9A9A', linestyle='--', linewidth=0.9)
    ax.text(1.0, max(counts) * 1.16, 'Metagenomics', ha='center',
            fontsize=8.5, style='italic', color='#4A4A4A')
    ax.text(4.5, max(counts) * 1.16, 'Host omics', ha='center',
            fontsize=8.5, style='italic', color='#4A4A4A')

    ax.text(0.0, -0.30,
            'Sequencing depth (unified pipeline): $\\it{P.\\,eryngii}$ 20.5–24.3 Gb/sample; '
            '$\\it{H.\\,tessellatus}$ and $\\it{P.\\,ostreatus}$ 5.8–9.5 Gb/sample; '
            'host reads 0.15–5.0% (one sample 11.4%); 32 samples from two animals (series 1 = MI02, series 2 = MI03).',
            transform=ax.transAxes, fontsize=7.5, color='#4A4A4A', va='top')


def build_figure1:
    fig = plt.figure(figsize=(11, 11))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.92, 1.05, 0.72], hspace=0.26)

    axes = [fig.add_subplot(gs[i]) for i in range(3)]
    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])

    # Panel letters
    for ax, letter in zip(axes, 'ABC'):
        pos = ax.get_position
        fig.text(pos.x0 - 0.045, pos.y1 + 0.006, letter,
                 fontsize=15, fontweight='bold', va='bottom')

    save_fig(fig, 'Fig1_Study_Design')


if __name__ == '__main__':
    print('Generating Figure 1...')
    build_figure1
    print('Done.')
