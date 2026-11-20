#!/usr/bin/env python3
"""
Generate Figure 5B — Maximum-likelihood phylogeny of vertebrate GH18 chitinases.

Input : results/phylogeny/v2/gh18_v2.contree  (IQ-TREE consensus tree, 162 sequences:
        154 public + 3 annotated M. impressa gene models + 5 unannotated loci)
Output: figures/png/SuppFig_S11_full_phylogeny.png (+ PDF)

The tree is rooted on the CTBS/CHID1 outgroup paralogs (as Methods 2.11 states),
tips are coloured by host lineage, all eight M. impressa GH18 loci are
highlighted, the mammalian CHIA / CHIT1 / CHI3L reference clades are labelled,
and UFBoot support is printed on EVERY node that the text relies on (not only
on nodes >= 90, which hid the contested 41% and 74% nodes in the first version).
"""

import os
import re

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo
from matplotlib.lines import Line2D

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREE = os.path.join(BASE, 'results', 'phylogeny', 'v2', 'gh18_v2.contree')
PNG_DIR = os.path.join(BASE, 'figures', 'png')
PDF_DIR = os.path.join(BASE, 'figures', 'pdf')

# Host lineage → colour. Order matters: first match wins.
LINEAGES = [
    ('Manouria impressa', '#D7263D', ['Manouria_impressa']),
    ('Testudines (turtles)', '#1B7F79', [
        'Chrysemys_picta', 'Pelodiscus_sinensis', 'Chelonia_mydas',
        'Chelonoidis_abingdonii', 'Gopherus_agassizii', 'Gopherus_evgoodei',
        'Terrapene', 'Trachemys', 'Mauremys', 'Malaclemys', 'Dermochelys',
    ]),
    ('Crocodilia', '#8E6C1F', ['Alligator', 'Crocodylus', 'Gavialis']),
    ('Squamata (lizards)', '#5FA05A', ['Anolis', 'Pogona', 'Python', 'Lacerta', 'Gekko']),
    ('Aves (birds)', '#E08A1E', ['Gallus', 'Meleagris', 'Taeniopygia', 'Coturnix', 'Anas']),
    ('Mammalia', '#2E5EAA', [
        'Homo_sapiens', 'Mus_musculus', 'Rattus', 'Bos_taurus', 'Sus_scrofa',
        'Canis', 'Felis', 'Ovis', 'Oryctolagus', 'Macaca', 'Pan_', 'Equus',
        'Capra', 'Cricetulus', 'Mesocricetus', 'Cavia', 'Pongo', 'Gorilla',
    ]),
    ('Amphibia', '#7D5BA6', ['Xenopus', 'Rana', 'Nanorana', 'Bufo']),
    ('Actinopterygii (fish)', '#4C8FB5', ['Danio_rerio', 'Oryzias', 'Salmo', 'Oncorhynchus', 'Takifugu']),
]
OTHER_COLOR = '#8A8A8A'

MIMP_LABELS = {
    'MimpChit1-active': '088392 / L1 (annotated; CHIT1-clade)',
    'MimpCLP-inactive': '088393 / L2 (annotated; inactive CLP)',
    'MimpChit2-dualactive': '088402b (annotated; GH18 domain of the 088402 pair)',
    'Mimp_L3_GT31type': 'locus L3 (unannotated; GT31-type)',
    'Mimp_L4_CHIA1': 'locus L4 (unannotated; CHIA-clade copy 1)',
    'Mimp_L5_CHIA2': 'locus L5 (unannotated; CHIA-clade copy 2)',
    'Mimp_L6_CHIA3': 'locus L6 (unannotated; CHIA-clade copy 3)',
    'Mimp_L7_CHIA4': 'locus L7 (unannotated; CHIA-clade copy 4)',
}
# Reference proteins used to label clades on the figure
CLADE_REFS = {
    'mammalian CHIA': ['CHIA_Homo_sapiens', 'Chia_Mus_musculus'],
    'mammalian CHIT1': ['CHIT1_Homo_sapiens', 'Chit1_Mus_musculus'],
    'mammalian CHI3L1/CHI3L2': ['CHI3L1_Homo_sapiens', 'CHI3L2_Homo_sapiens'],
}


def lineage_of(name):
    for label, color, keys in LINEAGES:
        for k in keys:
            if k in name:
                return label, color
    return 'Other / unassigned', OTHER_COLOR


def pretty(name):
    """Turn GENE_Genus_species_ACC into a compact readable label."""
    for key, desc in MIMP_LABELS.items:
        if name.startswith(key):
            return f'>> M. impressa  {desc}'
    parts = name.split('_')
    gene = parts[0]
    m = re.search(r'_([A-Z][a-z]+)_([a-z]+)_', name)
    if m:
        sp = f'{m.group(1)[0]}. {m.group(2)}'
    else:
        sp = ' '.join(parts[1:3])
    acc = parts[-1]
    return f'{gene}  {sp}  ({acc})'


def save_fig(fig, name):
    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(PDF_DIR, exist_ok=True)
    fig.savefig(os.path.join(PNG_DIR, f'{name}.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(PDF_DIR, f'{name}.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f'  {name}.png / .pdf')


def build:
    if not os.path.exists(TREE):
        raise SystemExit(f'Tree not found: {TREE}')

    tree = Phylo.read(TREE, 'newick')
    # Root on the distant CTBS/CHID1 paralogs, as the Methods state. IQ-TREE
    # writes an unrooted tree; the first version of this figure was drawn from
    # an arbitrary (ingroup) root, which a referee rightly objected to.
    outgroup = [t for t in tree.get_terminals if re.search(r'CTBS|CHID1', t.name or '')]
    if outgroup:
        tree.root_with_outgroup(*outgroup)
    tree.ladderize(reverse=True)

    terminals = tree.get_terminals
    n = len(terminals)
    print(f'  tips: {n}')

    # y position per tip; x position = root-to-node branch length
    ypos = {t: i for i, t in enumerate(terminals)}

    def depth(clade):
        return tree.distance(tree.root, clade)

    def yof(clade):
        if clade.is_terminal:
            return ypos[clade]
        kids = clade.clades
        return (yof(kids[0]) + yof(kids[-1])) / 2.0

    fig, ax = plt.subplots(figsize=(13, max(12, n * 0.13)))

    # Draw branches
    def draw(clade):
        x0 = depth(clade)
        y0 = yof(clade)
        for kid in clade.clades:
            x1 = depth(kid)
            y1 = yof(kid)
            ax.plot([x0, x0], [y0, y1], color='#4A4A4A', lw=0.7, zorder=1)
            ax.plot([x0, x1], [y1, y1], color='#4A4A4A', lw=0.7, zorder=1)
            draw(kid)

    draw(tree.root)

    # Tip labels
    xmax = max(depth(t) for t in terminals)
    used = set
    for t in terminals:
        label, color = lineage_of(t.name or '')
        used.add(label)
        is_mimp = 'Manouria_impressa' in (t.name or '')
        ax.text(depth(t) + xmax * 0.012, ypos[t], pretty(t.name or ''),
                va='center', ha='left',
                fontsize=6.4 if not is_mimp else 7.6,
                color=color,
                fontweight='bold' if is_mimp else 'normal', zorder=3)
        if is_mimp:
            ax.scatter([depth(t)], [ypos[t]], s=26, color=color,
                       zorder=4, edgecolor='white', linewidth=0.6)

    # Support values: every node that unites a M. impressa locus with a reference
    # paralog, plus every node with >= 5 tips. Low values are printed too.
    mimp_tips = [t for t in terminals if 'Manouria_impressa' in (t.name or '')]
    ref_tips = [t for t in terminals if any(k in (t.name or '') for ks in CLADE_REFS.values for k in ks)]
    key_nodes = set
    for m in mimp_tips:
        for r in ref_tips:
            key_nodes.add(tree.common_ancestor([m, r]))
        for anc in tree.get_path(m)[:-1]:
            if len(anc.get_terminals) <= 20:
                key_nodes.add(anc)
    for clade in tree.get_nonterminals:
        if clade.confidence is None:
            continue
        try:
            sup = float(clade.confidence)
        except (TypeError, ValueError):
            continue
        if clade in key_nodes or len(clade.get_terminals) >= 5:
            ax.text(depth(clade) - xmax * 0.004, yof(clade) + 0.55, f'{sup:.0f}',
                    fontsize=5.2, color='#B0322F' if sup >= 80 else '#B08A2F',
                    fontweight='bold' if clade in key_nodes else 'normal',
                    ha='right', va='bottom', zorder=3)

    # Clade labels for the mammalian reference paralogs
    for label, keys in CLADE_REFS.items:
        tips = [t for t in terminals if any(k in (t.name or '') for k in keys)]
        if len(tips) >= 2:
            mrca = tree.common_ancestor(tips)
            ys = [ypos[t] for t in mrca.get_terminals]
            ax.plot([xmax * 1.52] * 2, [min(ys) - 0.3, max(ys) + 0.3], color='#2E5EAA', lw=2.2)
            if max(ys) - min(ys) < 6:      # short clade: horizontal label so it cannot overlap a neighbour
                ax.text(xmax * 1.56, (min(ys) + max(ys)) / 2, label, va='center',
                        ha='left', fontsize=7, color='#2E5EAA')   # offset past the rotated labels of neighbouring clades
            else:
                ax.text(xmax * 1.535, (min(ys) + max(ys)) / 2, label, rotation=270, va='center',
                        ha='left', fontsize=7, color='#2E5EAA')

    ax.set_xlim(-xmax * 0.02, xmax * 1.70)
    ax.set_ylim(-1.5, n + 0.5)
    ax.invert_yaxis
    ax.axis('off')

    # Scale bar
    bar = round(xmax / 6, 1) or 0.1
    ax.plot([0, bar], [n + 0.0, n + 0.0], color='#333333', lw=1.4)
    ax.text(bar / 2, n + 0.9, f'{bar} substitutions/site',
            ha='center', va='top', fontsize=7.5, color='#333333')

    handles = [Line2D([0], [0], color=c, lw=3, label=l)
               for l, c, _ in LINEAGES if l in used]
    if 'Other / unassigned' in used:
        handles.append(Line2D([0], [0], color=OTHER_COLOR, lw=3, label='Other / unassigned'))
    ax.legend(handles=handles, loc='lower right', frameon=False,
              fontsize=8, title='Host lineage', title_fontsize=8.5)

    ax.set_title(
        'Maximum-likelihood phylogeny of vertebrate GH18 chitinases (162 sequences)\n'
        'IQ-TREE, ModelFinder best model, 1,000 UFBoot; rooted on CTBS/CHID1; '
        'support shown on all nodes ≥5 tips and on every node joining a M. impressa locus (bold)',
        fontsize=11, fontweight='bold', loc='left', pad=14)

    save_fig(fig, 'SuppFig_S11_full_phylogeny')


if __name__ == '__main__':
    print('Generating Figure 5B (GH18 phylogeny)...')
    build
    print('Done.')


# ---------------------------------------------------------------- condensed tree
def build_condensed:
    """Main-figure version: reference clades collapsed to wedges, every M. impressa
    locus shown individually. The full 162-tip tree stays as a supplementary figure —
    at 162 tips it is ~30 cm tall and illegible once scaled to a journal page."""
    tree = Phylo.read(TREE, 'newick')
    outgroup = [t for t in tree.get_terminals if re.search(r'CTBS|CHID1', t.name or '')]
    if outgroup:
        tree.root_with_outgroup(*outgroup)
    tree.ladderize(reverse=True)

    def has_mimp(c):
        return any('Manouria_impressa' in (t.name or '') for t in c.get_terminals)

    def paralog(name):
        # tip names are GENE_Genus_species_ACC with inconsistent case (CHIA / chia)
        g = (name or '').split('_')[0].upper
        if g.startswith('CHIA'): return 'CHIA'
        if g.startswith('CHIT1'): return 'CHIT1'
        if g.startswith('CHI3L'): return 'CHI3L1/2'
        if g.startswith('CTBS'): return 'CTBS'
        if g.startswith('CHID1'): return 'CHID1'
        if g.startswith('OVGP1'): return 'OVGP1'
        return 'GH18 (unnamed)'

    # a clade collapses if it holds no M. impressa locus and is paralog-uniform
    units = []                                    # (clade, kind) kind='tip'|'wedge'
    def walk(c):
        tips = c.get_terminals
        # collapse the LARGEST clade that contains no M. impressa locus: walking top-down,
        # the first such clade on each path becomes one wedge, which is what keeps the
        # display tree small (uniform-paralog collapsing left 86 units because the
        # paralogs interleave - which is itself one of the paper's points)
        if not has_mimp(c) and len(tips) >= 2:
            units.append((c, 'wedge')); return
        if c.is_terminal: units.append((c, 'tip')); return
        for k in c.clades: walk(k)
    walk(tree.root)

    ypos = {id(c): i for i, (c, _) in enumerate(units)}
    kind = {id(c): k for c, k in units}
    depth = lambda c: tree.distance(tree.root, c)

    def yof(c):
        if id(c) in ypos: return ypos[id(c)]
        kids = c.clades
        return (yof(kids[0]) + yof(kids[-1])) / 2.0

    fig, ax = plt.subplots(figsize=(11, max(5.0, len(units) * 0.165)))  # keeps Figure 5 under the 225 mm limit

    def draw(c):
        if id(c) in ypos: return
        x0, y0 = depth(c), yof(c)
        for kid in c.clades:
            x1, y1 = depth(kid), yof(kid)
            ax.plot([x0, x0], [y0, y1], color='#4A4A4A', lw=0.9, zorder=1)
            ax.plot([x0, x1], [y1, y1], color='#4A4A4A', lw=0.9, zorder=1)
            draw(kid)
    draw(tree.root)

    xmax = max(depth(t) for t in tree.get_terminals)
    for c, k in units:
        y = ypos[id(c)]; x = depth(c)
        if k == 'wedge':
            tips = c.get_terminals; dmax = max(depth(t) for t in tips)
            import collections
            pc = collections.Counter(paralog(t.name) for t in tips)
            top, ntop = pc.most_common(1)[0]
            lab, col = lineage_of(tips[0].name or '')
            ax.fill([x, dmax, dmax], [y, y - 0.34, y + 0.34], color=col, alpha=0.5,
                    edgecolor='#4A4A4A', lw=0.7, zorder=2)
            txt = f'{top}  ({len(tips)} sequences)' if ntop == len(tips) else \
                  f'{top} + {len(pc) - 1} other paralog{"s" if len(pc) > 2 else ""}  ({len(tips)} sequences)'
            ax.text(dmax + xmax * 0.012, y, txt, va='center', fontsize=8.5, color=col, zorder=3)
        else:
            nm = c.name or ''
            is_mimp = 'Manouria_impressa' in nm
            lab, col = lineage_of(nm)
            ax.plot([x, x], [y, y], marker='o', ms=4 if is_mimp else 0, color=col, zorder=4)
            ax.text(x + xmax * 0.012, y, pretty(nm), va='center', fontsize=8.5 if is_mimp else 7.5,
                    color=col, fontweight='bold' if is_mimp else 'normal', zorder=3)
    # UFBoot support on the nodes the text relies on: every ancestor of an M. impressa locus
    visible = {id(c) for c, _ in units}
    mimp_tips = [t for t in tree.get_terminals if 'Manouria_impressa' in (t.name or '')]
    shown = set
    for m in mimp_tips:
        for anc in tree.get_path(m)[:-1]:
            if id(anc) in shown or id(anc) in visible: continue
            try: sup = float(anc.confidence)
            except (TypeError, ValueError): continue
            shown.add(id(anc))
            ax.text(depth(anc) - xmax * 0.004, yof(anc) - 0.22, f'{sup:.0f}', fontsize=6,
                    color='#B0322F' if sup >= 80 else '#B08A2F', ha='right', va='bottom', zorder=3)
    handles = [Line2D([], [], marker='s', ls='', color=c, label=l)
               for l, c, _ in LINEAGES if l in {lineage_of(t.name or '')[0] for t in tree.get_terminals}]
    handles.append(Line2D([], [], marker='>', ls='', color='#8A8A8A', label='collapsed reference clade'))
    ax.legend(handles=handles, loc='lower right', frameon=False, fontsize=7.5, title='host lineage', title_fontsize=7.5)
    ax.set_ylim(-1, len(units)); ax.invert_yaxis; ax.set_xlim(-xmax * 0.02, xmax * 1.55)
    ax.set_yticks([]); ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.set_xlabel('substitutions per site', fontsize=9)
    ax.set_title('Maximum-likelihood phylogeny of vertebrate GH18 chitinases — reference clades collapsed\n'
                 f'(162 sequences, JTTDCMUT+G4; full tree in the supplementary figures)',
                 fontsize=10, fontweight='bold', loc='left')
    fig.tight_layout
    save_fig(fig, 'Fig5B_chitinase_phylogeny_condensed')
    print(f'  condensed to {len(units)} visible units')


if __name__ == '__main__':
    build_condensed
