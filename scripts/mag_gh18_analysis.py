#!/usr/bin/env python3
"""Genome-resolved view of the chitinolytic repertoire: 389 species-level MAGs (MetaBAT2 2.15 -> CheckM2 -> dRep 95 % ANI
-> GTDB-Tk R226 -> dbCAN3 >= 2 methods). Inputs in results/reanalysis_2026-11/mags/. Outputs:
  tables/64_mag_summary.tsv            one row per MAG: taxonomy, quality, size, CAZy counts (GH18/GH19/GH20/CBM50/GH5/GH9), chitin-route completeness
  tables/65_mag_gh18_by_taxon.tsv      GH18-carrying MAGs by phylum / genus, with per-cohort abundance (summed length-weighted depth of cluster members)
  tables/66_mag_gh18_abundance.tsv     per sample x GH18-carrying MAG: depth of the sample's member bin of that cluster
  figures Fig3C_GH18_taxonomy (genome-resolved genus share of GH18 abundance per cohort; the best-hit version is Supp S2) and SuppFig_S8_mag_overview."""
import os, re, pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); M = os.path.join(BASE, 'results', 'reanalysis_2026-11', 'mags'); T = os.path.join(BASE, 'tables')
coh = lambda s: 'DS01' if s.startswith('ITF') else 'DS02' if s.startswith('HM') else 'DS03'
DIET = {'DS01': 'P. eryngii', 'DS02': 'H. tessellatus', 'DS03': 'P. ostreatus'}; COL = {'DS01': '#2E86AB', 'DS02': '#A23B72', 'DS03': '#F18F01'}
gt = pd.concat([pd.read_csv(os.path.join(M, f), sep='\t') for f in ('gtdbtk.bac120.summary.tsv', 'gtdbtk.ar53.summary.tsv')]).set_index('user_genome')
tax = gt.classification.str.split(';', expand=True); tax.columns = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
tax = tax.apply(lambda c: c.str.replace(r'^[a-z]__', '', regex=True))
q = pd.read_csv(os.path.join(M, 'bin_quality.tsv'), sep='\t').set_index('Name')
cz = pd.read_csv(os.path.join(M, 'mag_cazymes_ge2.tsv'), sep='\t', header=None, names=['mag', 'gene', 'fams'])
def fams(s): return set(re.findall(r'(GH\d+|GT\d+|PL\d+|CE\d+|AA\d+|CBM\d+)', str(s)))
cz['famset'] = cz.fams.map(fams)
FAM = ['GH18', 'GH19', 'GH20', 'CBM50', 'CBM5', 'AA10', 'GH5', 'GH9', 'GH48']
cnt = pd.DataFrame({f: cz.groupby('mag').famset.apply(lambda x: sum(f in s for s in x)) for f in FAM}).reindex(gt.index).fillna(0).astype(int)
cnt['CAZymes'] = cz.groupby('mag').size.reindex(gt.index).fillna(0).astype(int)
# chitin route: GH18 (chitinase) + GH20 (hexosaminidase) in the same genome = chitin -> GlcNAc; nagA/nagB not in dbCAN (KEGG) -> flagged from GH only
summ = tax.join(q[['Completeness', 'Contamination', 'Genome_Size', 'quality']]).join(cnt)
summ['chitin_to_GlcNAc'] = (summ.GH18 > 0) & (summ.GH20 > 0)
summ.index.name = 'MAG'; summ.to_csv(os.path.join(T, '64_mag_summary.tsv'), sep='\t')
# cluster membership -> per-sample abundance of each species cluster (sum of depth of member bins in that sample)
cdb = pd.read_csv(os.path.join(M, 'Cdb.csv')); wdb = pd.read_csv(os.path.join(M, 'Wdb.csv'))
rep_of = dict(zip(wdb.cluster, wdb.genome.str.replace('.fa', '', regex=False))); cdb['rep'] = cdb.secondary_cluster.map(rep_of); cdb['bin'] = cdb.genome.str.replace('.fa', '', regex=False)
dep = pd.read_csv(os.path.join(M, 'bin_depth_per_sample.tsv'), sep='\t'); dep = dep.merge(cdb[['bin', 'rep']], on='bin', how='inner')
ab = dep.groupby(['rep', 'sample']).mean_depth.sum.unstack(fill_value=0.0)
# normalise per sample to relative depth (fraction of all bins' depth in that sample)
rel = ab / ab.sum(axis=0)
gh = summ[summ.GH18 > 0]
print(f'MAGs: {len(summ)}; with GH18: {len(gh)}; with GH18+GH20: {int(summ.chitin_to_GlcNAc.sum)}')
print('GH18 MAGs by phylum:', gh.phylum.value_counts.to_dict)
print('GH18 MAGs by genus (top):', gh.genus.replace('', 'unassigned').value_counts.head(10).to_dict)
# per cohort: share of GH18-MAG relative abundance by genus
rows = []
for c in DIET:
    ss = [s for s in rel.columns if coh(s) == c]; r = rel[ss].reindex(gh.index).fillna(0).mean(axis=1)
    tot = r.sum
    rg = r * gh.GH18.reindex(r.index); totg = rg.sum   # GH18-copy-weighted (the share quoted in the text and drawn in Fig 3C)
    for mag, v in r.items: rows.append(dict(Cohort=c, MAG=mag, phylum=gh.phylum[mag], genus=gh.genus[mag] or 'unassigned', species=gh.species[mag] or '', GH18_genes=int(gh.GH18[mag]), GH20_genes=int(gh.GH20[mag]), quality=gh.quality[mag], mean_rel_depth=round(v, 5), share_of_GH18_MAG_abundance=round(v / tot, 4) if tot else 0, share_GH18_gene_weighted=round(rg[mag] / totg, 4) if totg else 0))
t65 = pd.DataFrame(rows); t65.to_csv(os.path.join(T, '65_mag_gh18_by_taxon.tsv'), sep='\t', index=False)
rel.reindex(gh.index).to_csv(os.path.join(T, '66_mag_gh18_abundance.tsv'), sep='\t')
gshare = t65.groupby(['Cohort', 'genus']).share_of_GH18_MAG_abundance.sum.unstack(0).fillna(0)
top = gshare.mean(axis=1).sort_values(ascending=False).head(12).index.tolist
print(gshare.loc[top].round(3))
# GH18-gene-weighted share: each MAG's abundance x its GH18 gene count
rows2 = []
for c in DIET:
    ss = [s for s in rel.columns if coh(s) == c]; r = (rel[ss].reindex(gh.index).fillna(0).mean(axis=1) * gh.GH18); r = r / r.sum
    rows2.append(r.groupby(gh.genus.replace('', 'unassigned')).sum.rename(c))
gw = pd.concat(rows2, axis=1).fillna(0); topw = gw.mean(axis=1).sort_values(ascending=False).head(12).index.tolist
fig, ax = plt.subplots(figsize=(10, 5.2)); x = np.arange(len(topw)); w = 0.27
for i, c in enumerate(DIET): ax.bar(x + (i - 1) * w, gw.loc[topw, c] * 100, w, color=COL[c], edgecolor='#333', linewidth=0.5, label=DIET[c])
ax.set_xticks(x); ax.set_xticklabels([(f'$\\it{{{g.split("_")[0].replace("-", chr(0x2011))}}}$' + (f' ({g.split("_",1)[1]})' if '_' in g else '')) if g != 'unassigned' else 'no genus (GTDB)' for g in topw], rotation=40, ha='right', fontsize=9)
ax.set_ylabel('share of GH18 gene abundance carried by MAGs (%)'); ax.legend(frameon=False, fontsize=9); ax.spines[['top', 'right']].set_visible(False)
ax.set_title(f'Genome-resolved origin of GH18 chitinase genes: {len(gh)} of {len(summ)} species-level MAGs carry GH18\n(MAG relative depth × GH18 gene count, per diet; GTDB-Tk R226 genus)', fontweight='bold', fontsize=10)
fig.tight_layout
for _d in ('png', 'pdf'):  # a fresh clone has no figures/ yet
    os.makedirs(os.path.join(BASE, 'figures', _d), exist_ok=True)
for d in ('png', 'pdf'): fig.savefig(os.path.join(BASE, 'figures', d, f'Fig3C_GH18_taxonomy.{d}'), dpi=300)
# Supp S9: MAG overview
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axes[0]; ph = summ.phylum.value_counts.head(10); ax.barh(range(len(ph)), ph.values, color='#5A7A58'); ax.set_yticks(range(len(ph))); ax.set_yticklabels(ph.index, fontsize=8); ax.invert_yaxis; ax.set_xlabel('MAGs'); ax.set_title(f'A  {len(summ)} species-level MAGs by phylum', loc='left', fontsize=10, fontweight='bold')
ax = axes[1]
for lab, col, keep in (('carries GH18', '#B5541C', summ.GH18 > 0), ('no GH18', '#999999', summ.GH18 == 0)):
    ax.scatter(summ.Completeness[keep], summ.Contamination[keep], s=12, c=col, alpha=0.7, label=lab)
ax.set_xlabel('CheckM2 completeness (%)'); ax.set_ylabel('contamination (%)')
ax.legend(frameon=False, fontsize=8, markerscale=1.6); ax.set_title('B  MAG quality', loc='left', fontsize=10, fontweight='bold')
ax = axes[2]; lv = pd.Series({'genus': (summ.genus != '').sum, 'species': (summ.species != '').sum}); ax.bar(['assigned to a\nknown genus', 'assigned to a\nknown species'], lv.values, color='#2E6F9E'); ax.set_ylabel('MAGs'); ax.set_title('C  GTDB-Tk novelty', loc='left', fontsize=10, fontweight='bold')
for i, v in enumerate(lv.values): ax.text(i, v + 3, f'{v} / {len(summ)}', ha='center', fontsize=9)
for a in axes: a.spines[['top', 'right']].set_visible(False)
fig.suptitle('Supplementary Figure S8. Genome-resolved overview of the MAG set', fontweight='bold', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
for d in ('png', 'pdf'): fig.savefig(os.path.join(BASE, 'figures', d, f'SuppFig_S8_mag_overview.{d}'), dpi=300)
print('route complete (GH18+GH20) MAGs by genus:', summ[summ.chitin_to_GlcNAc].genus.replace('', 'unassigned').value_counts.head(8).to_dict)
print('GH18 genes in MAGs total:', int(summ.GH18.sum), '; MAG-GH18 share of catalogue GH18 (643):', round(summ.GH18.sum / 643, 2))
