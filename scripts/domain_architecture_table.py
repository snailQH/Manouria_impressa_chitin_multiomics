#!/usr/bin/env python3
"""Pfam domain architecture of the twelve Mimp_v3 gene models in the chr4 GH18 locus (MIMP3_g08468-g08479) and the four human
GH18 paralogs (CHIA, CHIT1, CHI3L1, CHI3L2). hmmscan (HMMER 3.4, --cut_ga) against Pfam-A 39.0 was run by
scripts/reanalysis_2026-11/domain_architecture_hmmscan.sh; inputs results/reanalysis_2026-11/domain_architecture/{queries.faa,pfam.domtbl}.
The catalytic motif is located by the GH18 consensus DxxDxDxE inside the Glyco_hydro_18 domain. No signal-peptide predictor was run.
Outputs: tables/70_gh18_locus_domain_architecture.tsv, figures SuppFig_S3_gh18_locus_domains.{png,pdf}."""
import re, pathlib, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt; from matplotlib.patches import Patch
R = pathlib.Path(__file__).resolve.parents[1]; D = R/'results/reanalysis_2026-11/domain_architecture'
seqs = {}; name = None
for l in open(D/'queries.faa'):
    if l.startswith('>'): name = l[1:].split[0]; seqs[name] = ''
    else: seqs[name] += l.strip
dom = {}
for l in open(D/'pfam.domtbl'):
    if l.startswith('#'): continue
    p = l.split; dom.setdefault(p[3], []).append((int(p[17]), int(p[18]), p[0], p[1].split('.')[0], float(p[12])))
LOCUS = {'MIMP3_g08468.t1': ('088402 (v1 N-terminal third)', 'TAF8'), 'MIMP3_g08469.t1': ('088402a (v1 domain 1)', 'reptile lineage'), 'MIMP3_g08470.t1': ('088402b (v1 domain 2)', 'reptile lineage'),
         'MIMP3_g08471.t1': ('fragment between 088402 and L7', '—'), 'MIMP3_g08472.t1': ('L7', 'CHIA array'), 'MIMP3_g08473.t1': ('L6', 'CHIA array'), 'MIMP3_g08474.t1': ('L5', 'CHIA array'),
         'MIMP3_g08475.t1': ('L4', 'CHIA array'), 'MIMP3_g08476.t1': ('L3', 'CHIA-related'), 'MIMP3_g08477.t1': ('088393 fragment (v1 C-terminal part)', 'CHI3L lineage'),
         'MIMP3_g08478.t1': ('088393 / L2', 'CHI3L lineage CLP'), 'MIMP3_g08479.t1': ('088392 / L1', 'CHIT1/CHI3L lineage')}
rows = []
for q in seqs:
    s = seqs[q]; ds = sorted(dom.get(q, [])); gh = [d for d in ds if d[2] == 'Glyco_hydro_18']
    motif = ''
    for a, b, *_ in gh:
        m = re.search(r'D..D.D.E', s[a - 1:b]); 
        if m: motif = f'{m.group(0)} ({a + m.start}-{a + m.end - 1})'; break
    if not motif and gh:
        for a, b, *_ in gh:
            m = re.search(r'D..[ANSG].[DEN].[EQ]|[DEN]..D.[DEN].[EQ]', s[a - 1:b])
            if m: motif = f'{m.group(0)} ({a + m.start}-{a + m.end - 1}) — degenerate'; break
    if not motif and gh: motif = 'no DxxDxDxE in the GH18 domain'
    locus, clade = LOCUS.get(q, ({'CH3L1': 'CHI3L1', 'CH3L2': 'CHI3L2'}.get(q.split('|')[2].replace('_HUMAN', ''), q.split('|')[2].replace('_HUMAN', '')) if '|' in q else q, 'human'))
    rows.append(dict(protein=q, locus=locus, clade=clade, length=len(s), domains='; '.join(f'{n} {acc} {a}-{b} (E={e:.1e})' for a, b, n, acc, e in ds) or 'none',
                     n_GH18_domains=len(gh), CBM14='yes' if any(d[2] == 'CBM_14' for d in ds) else 'no', catalytic_motif=motif))
with open(R/'tables/70_gh18_locus_domain_architecture.tsv', 'w') as fh:
    fh.write('# Pfam-A 39.0 (hmmscan, HMMER 3.4, gathering thresholds) domain architecture of the twelve Mimp_v3 models of the chr4 GH18 locus and the four human GH18 paralogs; catalytic motif = DxxDxDxE located\n'
             '# inside the Glyco_hydro_18 domain (position in the protein). g08468 is TAF8 (Bromo_TP + TAF8_C): the N-terminal third of the fused 1,160-aa model 088402 is a transcription-factor subunit, not a chitinase.\n'
             '# No signal-peptide prediction was performed. scripts/reanalysis_2026-11/domain_architecture_hmmscan.sh; scripts/domain_architecture_table.py\n')
    fh.write('\t'.join(rows[0].keys) + '\n'); [fh.write('\t'.join(map(str, r.values)) + '\n') for r in rows]
for r in rows: print(r['protein'], r['locus'], r['length'], r['n_GH18_domains'], r['CBM14'], r['catalytic_motif'], '|', r['domains'])
# ---- Supp Fig S4
COL = {'Glyco_hydro_18': '#c8452c', 'CBM_14': '#2f6db3', 'Bromo_TP': '#6D4C41', 'TAF8_C': '#D7A86E'}   # the two TAF8 domains were near-identical browns
top_to_bottom = sorted([q for q in seqs if q.startswith('MIMP3')], reverse=True) + [q for q in seqs if not q.startswith('MIMP3')]  # 088392 (g08479) first, human last
order = top_to_bottom[::-1]
fig, ax = plt.subplots(figsize=(11, 0.42 * len(order) + 2.2))
for i, q in enumerate(order):
    s = seqs[q]; L = len(s); ax.plot([0, L], [i, i], color='0.3', lw=1.2, zorder=1)
    for a, b, n, acc, e in dom.get(q, []):
        ax.add_patch(plt.Rectangle((a, i - 0.3), b - a + 1, 0.6, facecolor=COL.get(n, '#999'), edgecolor='k', lw=0.5, zorder=2))
    r = next(x for x in rows if x['protein'] == q)
    m = re.search(r'\((\d+)-(\d+)\)', r['catalytic_motif'])
    if m: ax.plot(int(m.group(1)), i + 0.42, marker='v', color='k', ms=5, zorder=3)
    hum = {'CH3L1': 'CHI3L1', 'CH3L2': 'CHI3L2'}
    gid = q.split('|')[-1].replace('.t1', '') if '|' not in q else q.split('|')[2].replace('_HUMAN', '')
    lab = f"{hum.get(r['locus'], r['locus'])}  ({hum.get(gid, gid)}; {L} aa)"
    ax.text(-15, i, lab, ha='right', va='center', fontsize=7.5)
    ax.text(L + 12, i, r['catalytic_motif'].split(' (')[0] if r['n_GH18_domains'] else '', va='center', fontsize=7)   # one font family throughout
XR = max(len(v) for v in seqs.values) * 1.30
ax.set_yticks([]); ax.set_xlim(-10, XR); ax.set_ylim(-0.8, len(order) - 0.2); ax.set_xlabel('amino-acid position', fontsize=8)
for sp in ('top', 'right', 'left'): ax.spines[sp].set_visible(False)
ax.legend(handles=[Patch(facecolor=c, edgecolor='k', label=k) for k, c in COL.items] + [plt.Line2D([], [], marker='v', color='k', ls='', label='DxxDxDxE catalytic motif')],
          fontsize=7, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=5)
ax.axhline(3.5, color='0.8', lw=0.8, ls='--')
ax.text(XR * 0.995, 2.95, 'human paralogs', ha='right', va='center', fontsize=7, color='0.4')
ax.text(XR * 0.995, 4.05, '$\\mathit{M.\\,impressa}$ chr4 locus (genomic order, 088392 end first)', ha='right', va='center', fontsize=7, color='0.4')
fig.suptitle('Supplementary Figure S3. Pfam domain architecture of the chr4 GH18-locus models (Mimp_v3) and the human GH18 paralogs', fontsize=9.5)
fig.tight_layout(rect=(0, 0, 1, 0.97))
for _d in ('png', 'pdf'):  # a fresh clone has no figures/ yet
    (R/'figures'/_d).mkdir(parents=True, exist_ok=True)
for ext in ('png', 'pdf'): fig.savefig(R/f'figures/{ext}/SuppFig_S3_gh18_locus_domains.{ext}', dpi=300)
print('written SuppFig_S3_gh18_locus_domains')
