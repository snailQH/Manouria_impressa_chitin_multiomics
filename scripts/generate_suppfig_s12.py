#!/usr/bin/env python3
"""Supplementary Figure S12 - full-length human GH18 paralogs vs M. impressa GH18-locus proteins and
NCBI Testudines 'chitotriosidase-1-like' proteins: BLASTP identity and query coverage
(tables/25_CHIA_CHIT1_blast_verification.tsv). Replaces the earlier version built from
truncated queries, which depicted absence where full-length searches give strong alignments."""
import os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, pandas as pd, numpy as np
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
t = pd.read_csv(os.path.join(BASE, 'tables', '25_CHIA_CHIT1_blast_verification.tsv'), sep='\t', comment='#')
t['q'] = t.Query.str.split('|').str[2].str.replace('_HUMAN', '')
subj = {'GWHPHIVU010568.1': '$\\mathit{M.\\,impressa}$ 088392 / L1', 'GWHPHIVU010569.1': '$\\mathit{M.\\,impressa}$ 088393 / L2', 'GWHPHIVU010570.1': '$\\mathit{M.\\,impressa}$ 088402',
        'XP_007065194.2': '$\\mathit{C.\\,mydas}$ XP_007065194', 'XP_005308333.2': '$\\mathit{C.\\,picta}$ XP_005308333',
        'XP_032641725.1': '$\\mathit{C.\\,abingdonii}$ XP_032641725', 'XP_030417332.1': '$\\mathit{G.\\,evgoodei}$ XP_030417332'}
t['s'] = t.Subject.map(subj)
qs = ['CHIA', 'CHIT1', 'CH3L1', 'CH3L2']; ss = list(subj.values)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
for ax, col, title, vmin, vmax, cbl in ((axes[0], 'Pident', 'A  Percent identity', 40, 65, 'identity (%, scale clipped at 40-65)'),
                                        (axes[1], 'Qcov', 'B  Query coverage', 70, 100, 'query coverage (%)')):
    M = np.full((len(ss), len(qs)), np.nan)
    for i, s_ in enumerate(ss):
        for j, q in enumerate(qs):
            r = t[(t.s == s_) & (t.q == q)]
            if len(r): M[i, j] = r[col].max
    im = ax.imshow(M, cmap='YlGnBu', vmin=vmin, vmax=vmax, aspect='auto')
    for i in range(len(ss)):
        for j in range(len(qs)):
            if np.isnan(M[i, j]): continue
            # YlGnBu goes near-navy at the top of the scale: black annotations there are unreadable
            dark = (M[i, j] - vmin) / (vmax - vmin) > 0.55
            ax.text(j, i, f'{M[i, j]:.0f}', ha='center', va='center', fontsize=9,
                    color='white' if dark else '#222', fontweight='bold' if dark else 'normal')
    ax.set_xticks(range(len(qs))); ax.set_xticklabels(['CHIA', 'CHIT1', 'CHI3L1', 'CHI3L2'])
    ax.set_yticks(range(len(ss))); ax.set_yticklabels(ss, fontsize=8.5)
    ax.set_title(title, loc='left', fontsize=10, fontweight='bold')
    cb = plt.colorbar(im, ax=ax, fraction=0.046); cb.set_label(cbl, fontsize=8); cb.ax.tick_params(labelsize=7)
fig.suptitle('Supplementary Figure S12. Full-length human GH18 paralogs (queries) vs the annotated $\\mathit{M.\\,impressa}$ GH18-locus proteins\n'
             'and Testudines "chitotriosidase-1-like" proteins (BLASTP). '
             'Each turtle protein aligns comparably to all four paralogs: similarity cannot assign orthology; see the phylogeny (Fig. 5B)', fontsize=9.5)
plt.tight_layout
for d in ('png', 'pdf'):
    for _d in ('png', 'pdf'):  # a fresh clone has no figures/ yet
        os.makedirs(os.path.join(BASE, 'figures', _d), exist_ok=True)
    fig.savefig(os.path.join(BASE, 'figures', d, f'SuppFig_S12_CHIT1_turtle_comparison.{d}'), dpi=300, bbox_inches='tight')
print('SuppFig S12 regenerated')
