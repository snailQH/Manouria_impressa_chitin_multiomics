#!/usr/bin/env python3
"""Figure 6 (promoted from Supplementary Figure S3 on 2026-11-17) — architecture of the GH18 chitinase locus across nine Testudines genomes.
Input: results/reanalysis_2026-11/comparative_gh18/<assembly>.loci.tsv (miniprot loci of the 12 Mimp_v3 locus models,
clustered, pseudogene-flagged; see scripts/reanalysis_2026-11/comparative_gh18_parse.py) and tables/56_comparative_gh18_locus_testudines.tsv.
(A) one row per species, loci drawn in locus order from the CHIT1-clade gene (088392) end so that all species read the same way
irrespective of strand; box colour = clade, hatched = pseudogene (frameshift / in-frame stop in the best alignment).
(B) copy counts of the four locus components per species."""
import csv, pathlib, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt; from matplotlib.patches import Patch
R = pathlib.Path(__file__).resolve.parents[1]; D = R/'results/reanalysis_2026-11/comparative_gh18'
rows = [r for r in csv.DictReader((l for l in open(R/'tables/56_comparative_gh18_locus_testudines.tsv') if not l.startswith('#')), delimiter='\t')]
COL = {'088392 CHIT1-clade':'#2f6db3','088393 CLP':'#8ab6e8','CLP fragment (g08477)':'#d6e4f5','L3 CHIA-related':'#f2a93b','CHIA-array':'#c8452c','088402a GH18+CBM14':'#3d7a4c','088402b GH18':'#9ccba6'}
def clade(l):
    c = l['clade'] if isinstance(l, dict) else l
    if isinstance(l, dict) and l['best_query'] == 'MIMP3_g08477': return 'CLP fragment (g08477)'
    return 'CHIA-array' if c.startswith('CHIA-array') else c
figA, ax = plt.subplots(figsize=(13, 6.2))
figB, bx = plt.subplots(figsize=(9, 4.6))
XMAX = []
yl = []
for i, r in enumerate(rows):
    loci = [l for l in csv.DictReader((x for x in open(D/f"{r['assembly']}.loci.tsv") if not x.startswith('#')), delimiter='\t') if l['in_interval'] == 'True']
    loci.sort(key=lambda l: int(l['start']))
    # orient every species so the CHIT1-clade gene (088392) lies left of the 088402 models, irrespective of strand
    s392 = [int(l['start']) for l in loci if l['clade'] == '088392 CHIT1-clade']; s402 = [int(l['start']) for l in loci if l['clade'].startswith('088402')]
    rev = bool(s392 and s402 and s392[0] > sum(s402) / len(s402))
    mx = max(int(l['end']) for l in loci); o = min(int(l['start']) for l in loci)
    segs = sorted([((mx - int(l['end'])) if rev else (int(l['start']) - o), int(l['end']) - int(l['start']), l) for l in loci], key=lambda t: t[0])
    # compress intergenic gaps > 60 kb to 60 kb and mark the break (Pelodiscus: 088392 lies 240 kb from the rest)
    y = len(rows) - 1 - i; shift = 0; prev_end = None
    for x0, w, l in segs:
        x0 -= shift
        if prev_end is not None and x0 - prev_end > 60000:
            shift += x0 - prev_end - 60000; x0 = prev_end + 60000
            ax.text((prev_end + 30000) / 1000, y, '//', ha='center', va='center', fontsize=9, color='0.4')
        prev_end = x0 + w
        x0, w = x0 / 1000, w / 1000
        ax.add_patch(plt.Rectangle((x0, y - 0.32), w, 0.64, facecolor=COL[clade(l)], edgecolor='k', lw=0.6,
                                   hatch='////' if l['pseudogene'] == 'True' else None))
        # array copies are labelled generically: cross-species assignment to a particular M. impressa paralog (L4-L7) is not reliable
        lab = {'088392 CHIT1-clade':'392','088393 CLP':'393','CLP fragment (g08477)':'393f','L3 CHIA-related':'L3','088402a GH18+CBM14':'402a','088402b GH18':'402b','CHIA-array':'CHIA'}[clade(l)]
        ax.text(x0 + w / 2, y + 0.42, lab + ('ψ' if l['pseudogene'] == 'True' else ''), ha='center', va='bottom', fontsize=6.5)
        XMAX.append(x0 + w)
    yl.append(f"$\\it{{{r['species'].replace(' ', '\\ ')}}}$\n{r['family']} · {r['diet']}")
ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels(yl, fontsize=7.5)
ax.set_xlim(-5, max(XMAX) * 1.04); ax.set_ylim(-0.7, len(rows) - 0.2); ax.set_xlabel('kb from the first locus (each species oriented with the CHIT1-clade gene on the left; // = gap > 60 kb compressed)', fontsize=8)
ax.legend(handles=[Patch(facecolor=c, edgecolor='k', label=k) for k, c in COL.items] + [Patch(facecolor='w', edgecolor='k', hatch='////', label='pseudogene (frameshift / stop)')],
          fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4, frameon=False)
ax.set_title('GH18 locus architecture in the orthologous interval of nine Testudines genomes\n(miniprot placements of the eleven GH18 $\\mathit{M.\\,impressa}$ models; TAF8 g08468 excluded)', fontsize=10, loc='left')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
comp = [('CHIA array (intact)', 'CHIA_array_intact', '#c8452c'), ('CHIA array (ψ)', 'CHIA_array_pseudo', '#e8a89c'),
        ('L3 CHIA-related', 'L3_CHIA_related', '#f2a93b'), ('088393 CLP', 'CLP_088393', '#8ab6e8'), ('088392 CHIT1-clade', 'CHIT1_088392', '#2f6db3'),
        ('088402a GH18+CBM14', 'GH18_088402a_CBM14', '#3d7a4c'), ('088402b GH18', 'GH18_088402b', '#9ccba6')]
import numpy as np
yy = np.arange(len(rows))[::-1]; left = np.zeros(len(rows))
for lab, k, c in comp:
    v = np.array([int(r[k]) for r in rows]); bx.barh(yy, v, left=left, color=c, edgecolor='k', lw=0.5, height=0.64, label=lab); left += v
bx.set_yticks(yy); bx.set_yticklabels([r['species'].split[0][0] + '. ' + ' '.join(r['species'].split[1:]) for r in rows], fontsize=7.5, style='italic')
bx.set_xlabel('GH18 loci in the interval (excluding the g08477 fragment)', fontsize=8); bx.legend(fontsize=7, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)
bx.set_xlim(0, left.max * 1.04)   # trim AFTER the bars exist, never before
bx.set_title('Copy number per locus component', fontsize=10, loc='left')
for s in ('top', 'right'): bx.spines[s].set_visible(False)
figA.tight_layout; figB.tight_layout
for ext in ('png', 'pdf'):
    for _d in ('png', 'pdf'):  # a fresh clone has no figures/ yet
        (R/'figures'/_d).mkdir(parents=True, exist_ok=True)
    figA.savefig(R/f'figures/{ext}/Fig6A_comparative_locus_architecture.{ext}', dpi=300, bbox_inches='tight')
    figB.savefig(R/f'figures/{ext}/Fig6B_comparative_locus_copynumber.{ext}', dpi=300, bbox_inches='tight')
print('written Fig6A_comparative_locus_architecture, Fig6B_comparative_locus_copynumber')
