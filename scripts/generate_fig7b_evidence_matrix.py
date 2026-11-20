#!/usr/bin/env python3
"""Figure 7B - multi-omics evidence matrix, rebuilt from the current claims.
The previous PNG had no generating script and carried withdrawn statements
("AMCase host genome: Yes", "GH18 chitinase host genome: No")."""
import os, textwrap, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cols = ['Faecal\nmetagenome', 'Gut-wall\nproteome', 'Host genome', 'Heart/lung\nRNA-Seq']
rows = [
 ('Microbial GH18 chitinase',            ['Yes (643 genes, dbCAN3;\n317-665 CPM)', 'Low (5 gut-restricted prokaryotic\ngroups, mostly 1 peptide)', 'n/a', 'n/a']),
 ('Microbial GH20 hexosaminidase',       ['Yes (most abundant KO)', 'Trace', 'n/a', 'n/a']),
 ('Microbial nagA / nagB',               ['Yes', 'Trace', 'n/a', 'n/a']),
 ('Host CHIA-clade array (L4-L7)',       ['n/a', 'Yes (L4-L6 one group, 6 peptides; L7 4)', 'Yes (4 copies; 1-2 in\n8 other turtles)', 'No (0-9 reads)']),
 ('Host CHIA-related L3',                ['n/a', 'Yes (12 peptides; top abundance)', 'Yes (MIMP3_g08476)', 'No']),
 ('Host 088392 (CHIT1/CHI3L lineage)',   ['n/a', 'Yes (11 peptides)', 'Yes (MIMP3_g08479)', 'Low (TPM 4-32)']),
 ('Host 088393 (inactive CLP)',          ['n/a', 'Yes (17 peptides; flat across gut)', 'Yes (MIMP3_g08477-78)', 'Yes (TPM 179-324)']),
 ('Host 088402a/b (reptile GH18)',          ['n/a', 'Yes (2 of 3 models; 19 + 11 peptides)', 'Yes (MIMP3_g08469-70;\ng08468 = TAF8)', 'Low (TPM 2-3)']),
 ('Host GH20 + GlcNAc pathway',          ['n/a', 'Yes (HEXA/B/DC; 10 GlcNAc enzymes)', 'Yes', 'Yes']),
 ('Chitinase enzyme activity',           ['not measured', 'not measured', 'motif only', 'n/a']),
 ('Microbial GH18 per gut segment',      ['not measured\n(faeces only)', 'Yes, but no proximal-distal\ngradient (P = 0.94 / 0.10)', 'n/a', 'n/a']),
]
def level(t):
    t = t.lower
    if t.startswith('yes'): return 1.0
    if t.startswith(('partial', 'low', 'trace')): return 0.5
    if t.startswith(('no', 'not measured', 'motif')): return 0.0
    return np.nan
M = np.array([[level(c) for c in r[1]] for r in rows])
# Colour-blind-safe: the old green/red pair is indistinguishable under deuteranopia, so "observed"
# is blue — red vs blue is safe — and every cell also carries a glyph, so colour is never the only channel.
FILL = {1.0: '#1F5F8B', 0.5: '#F0A830', 0.0: '#A5352A'}
# Pin the font: matplotlib ships DejaVu Sans, which has these three glyphs in both weights. The
# shipped Figure 9B was built where the resolved font did not, so every icon rendered as a tofu box
# - in the TIFF the journal would have received. Do not remove this line.
plt.rcParams['font.family'] = 'DejaVu Sans'
GLYPH = {1.0: '\u2713', 0.5: '\u223c', 0.0: '\u2717'}
NA = '#ECECEC'

from matplotlib.font_manager import findfont, FontProperties
try:
    from fontTools.ttLib import TTFont as _TT
    _cm = set
    for _t in _TT(findfont(FontProperties(family=plt.rcParams['font.family'])), fontNumber=0)['cmap'].tables:
        _cm |= set(_t.cmap)
    _missing = [hex(ord(g)) for g in ('\u2713', '\u223c', '\u2717') if ord(g) not in _cm]
    assert not _missing, f'font lacks glyphs {_missing}: the evidence matrix would render tofu boxes'
except ImportError:
    pass   # fontTools absent: proceed, the pinned DejaVu Sans is the matplotlib default anyway
fig, ax = plt.subplots(figsize=(12.5, 6.6))   # trimmed so Figure 9 clears the 225 mm limit
ax.set_xlim(-0.5, len(cols) - 0.5); ax.set_ylim(len(rows) - 0.5, -0.5)
for i, (_, cells) in enumerate(rows):
    for j, c in enumerate(cells):
        lv = level(c); na = np.isnan(lv)
        ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, facecolor=NA if na else FILL[lv],
                                   edgecolor='white', lw=1.6))
        if na:
            ax.text(j, i, 'n/a', ha='center', va='center', fontsize=7.6, color='#888'); continue
        dark = lv in (0.0, 1.0)
        ax.text(j - 0.45, i, GLYPH[lv], ha='left', va='center', fontsize=10,
                color='white' if dark else '#5A3A00', fontweight='bold')
        # wrap so a long entry cannot run under the glyph or out of the cell
        w = '\n'.join('\n'.join(textwrap.wrap(ln, 30)) for ln in c.split('\n'))
        ax.text(j + 0.06, i, w, ha='center', va='center', fontsize=7.4,
                color='white' if dark else '#222', fontweight='bold' if dark else 'normal', linespacing=1.35)
ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, fontsize=9)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
for sp in ax.spines.values: sp.set_visible(False)
ax.tick_params(length=0)
ax.set_title('Multi-omics evidence matrix for chitin degradation in $\\mathit{M.\\,impressa}$',
             fontsize=11, fontweight='bold', pad=12)
key = [(1.0, 'observed'), (0.5, 'partial / low'), (0.0, 'not observed or not measured')]
handles = [plt.Line2D([], [], marker='s', ls='', ms=10, mfc=FILL[l], mec='none',
                      label=f'{GLYPH[l]}  {t}') for l, t in key]
handles.append(plt.Line2D([], [], marker='s', ls='', ms=10, mfc=NA, mec='none', label='   not applicable'))
ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.09), ncol=4,
          frameon=False, fontsize=8.5, handletextpad=0.4, columnspacing=1.6)
for d in ('png', 'pdf'):
    os.makedirs(os.path.join(BASE, 'figures', d), exist_ok=True)
    fig.savefig(os.path.join(BASE, 'figures', d, f'Fig9B_evidence_matrix.{d}'), dpi=300, bbox_inches='tight')
print('Fig7B regenerated')
