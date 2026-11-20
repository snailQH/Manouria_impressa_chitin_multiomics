#!/usr/bin/env python3
"""Taxonomic origin of the GH18 chitinase genes of the unified 32-sample catalogue, all three cohorts
(replaces the providers' DIAMOND-NR gene taxonomy, which existed for DS02/DS03 only).
Primary assignment = taxonomy of each gene's eggNOG-mapper seed ortholog (closest homolog in eggNOG 5; NCBI lineage
via the taxdump), which is the same kind of evidence as the providers' DIAMOND best hit and exists for 592/643 genes.
Cross-check = Kraken2 (confidence 0.2, same database as the read taxonomy) of the contig carrying each gene; only
70/624 contigs are classifiable, mirroring the ~90 % unclassified reads, so it is reported but not used for the figure.
Inputs (05.reanalysis_2026-11/gh18_taxonomy/): gh18_seed_taxonomy.tsv, gh18_contigs.kraken/.report, gh18_gene_contig.tsv,
and the per-gene counts (results/unified32/gh18_gene_counts.tsv). Abundance = the gene's read count (CPM). Outputs: tables/57_gh18_gene_taxonomy_unified32.tsv (per gene) and
tables/04_gh18_taxonomy_by_cohort.tsv (genus share of GH18 reads per cohort; replaces the provider-era Table 4),
figure SuppFig_S2_GH18_taxonomy_besthit (three cohorts; the genome-resolved MAG version is Figure 3C, scripts/mag_gh18_analysis.py)."""
import os, re, pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(BASE, 'results', 'reanalysis_2026-11', 'gh18_taxonomy')
coh = lambda s: 'DS01' if s.startswith('ITF') else 'DS02' if s.startswith('HM') else 'DS03'
DIET = {'DS01': 'P. eryngii', 'DS02': 'H. tessellatus', 'DS03': 'P. ostreatus'}; COL = {'DS01': '#2E86AB', 'DS02': '#A23B72', 'DS03': '#F18F01'}
# kraken output: C/U, seqid, "name (taxid N)", length, lca string
kr = {}
for l in open(os.path.join(R, 'gh18_contigs.kraken')):
    p = l.rstrip('\n').split('\t'); kr[p[1]] = p[2]
gc = pd.read_csv(os.path.join(R, 'gh18_gene_contig.tsv'), sep='\t', header=None, names=['sample', 'contig', 'gene'])
gc['contig_id'] = gc['sample'] + '_' + gc['contig']; gc['kraken'] = gc.contig_id.map(kr).fillna('unclassified (taxid 0)')
# rank from the Kraken report (name -> rank/lineage)
rep = pd.read_csv(os.path.join(R, 'gh18_contigs.report'), sep='\t', header=None, names=['pct', 'clade', 'direct', 'rank', 'taxid', 'name'])
rep['name'] = rep['name'].str.strip; rank = dict(zip(rep.taxid, rep['rank']))
# lineage walk: build parent map from the indented report (indentation = depth)
lines = open(os.path.join(R, 'gh18_contigs.report')).read.split('\n'); stack = []; lineage = {}
for l in lines:
    if not l.strip: continue
    p = l.split('\t'); name = p[5]; depth = (len(name) - len(name.lstrip(' '))) // 2; name = name.strip; taxid = int(p[4]); rk = p[3]
    stack = stack[:depth]; stack.append((taxid, rk, name)); lineage[taxid] = list(stack)
def genus_phylum(s):
    m = re.search(r'\(taxid (\d+)\)', s); t = int(m.group(1)) if m else 0
    lin = lineage.get(t, []); g = next((n for ti, rk, n in lin if rk == 'G'), None); ph = next((n for ti, rk, n in lin if rk == 'P'), None)
    dom = next((n for ti, rk, n in lin if rk == 'D'), None)
    if t == 0: return 'unclassified', 'unclassified', 'unclassified'
    return (g or 'unassigned at genus'), (ph or 'unassigned at phylum'), (dom or 'root')
gc[['kraken_genus', 'kraken_phylum', 'kraken_domain']] = gc.kraken.map(genus_phylum).apply(pd.Series)
seed = pd.read_csv(os.path.join(R, 'gh18_seed_taxonomy.tsv'), sep='\t').set_index('gene')
gc = gc.set_index('gene').join(seed[['seed_taxid', 'genus', 'family', 'phylum', 'domain']]).reset_index
gc[['genus', 'family', 'phylum', 'domain']] = gc[['genus', 'family', 'phylum', 'domain']].fillna('unassigned')
cnt = pd.read_csv(os.path.join(BASE, 'results', 'unified32', 'gh18_gene_counts.tsv'), sep='\t', index_col=0)   # gene x sample read counts
gc = gc.set_index('gene'); gc['reads_total'] = cnt.sum(axis=1).reindex(gc.index).fillna(0).astype(int)
gc.reset_index.to_csv(os.path.join(BASE, 'tables', '57_gh18_gene_taxonomy_unified32.tsv'), sep='\t', index=False)
# per-cohort genus share of GH18 reads (CPM-weighted: each sample's GH18 reads normalised by its mapped total)
tot = pd.read_csv(os.path.join(BASE, 'results', 'unified32', 'unified32_mapping.tsv'), sep='\t', header=None, names=['s', 'pairs', 'rate'])
mapped = pd.read_csv(os.path.join(BASE, 'results', 'unified32', 'unified32_cazy_abundance_dbcan.tsv'), sep='\t', comment='#').set_index('Family').loc['TOTAL_MAPPED']
rows = []
for c in ['DS01', 'DS02', 'DS03']:
    ss = [s for s in cnt.columns if coh(s) == c]
    cpm = (cnt[ss] / mapped[ss] * 1e6).mean(axis=1)                      # mean CPM per gene within cohort
    df = pd.DataFrame({'cpm': cpm}).join(gc[['genus', 'phylum']]); tot_cpm = df.cpm.sum
    g = df.groupby('genus').cpm.sum.sort_values(ascending=False)
    for genus, v in g.items: rows.append(dict(Cohort=c, Diet=DIET[c], genus=genus, phylum=df[df.genus == genus].phylum.iloc[0], GH18_CPM=round(v, 2), share=round(v / tot_cpm, 4), n_genes=int((df.genus == genus).sum)))
t4 = pd.DataFrame(rows); t4.to_csv(os.path.join(BASE, 'tables', '04_gh18_taxonomy_by_cohort.tsv'), sep='\t', index=False)
# figure: top genera (by mean share across cohorts) + unassigned/unclassified, three cohorts
top = t4[~t4.genus.isin(['unassigned'])].groupby('genus').share.mean.sort_values(ascending=False).head(12).index.tolist
fig, ax = plt.subplots(figsize=(10, 5.2)); cats = top + ['unassigned']; x = np.arange(len(cats)); w = 0.27
for i, c in enumerate(['DS01', 'DS02', 'DS03']):
    sub = t4[t4.Cohort == c].set_index('genus').share.reindex(cats).fillna(0) * 100
    ax.bar(x + (i - 1) * w, sub.values, w, color=COL[c], edgecolor='#333', linewidth=0.5, label=f'{DIET[c]} (n = {sum(1 for s in cnt.columns if coh(s) == c)})')
ax.set_xticks(x); ax.set_xticklabels([f'$\\it{{{g}}}$' if g != 'unassigned' else 'no homolog /\nunassigned' for g in cats], rotation=40, ha='right', fontsize=9)
ax.set_ylabel('Share of GH18 reads (%)'); ax.legend(frameon=False, fontsize=9)
ax.set_title('Best-hit taxonomy of GH18 chitinase genes, all three diets (cross-check of Figure 3C)\n(643 dbCAN3 GH18 genes of the unified catalogue; closest eggNOG homolog, read-count weighted)', fontweight='bold', fontsize=10)
ax.spines[['top', 'right']].set_visible(False); fig.tight_layout
for _d in ('png', 'pdf'):  # a fresh clone has no figures/ yet
    os.makedirs(os.path.join(BASE, 'figures', _d), exist_ok=True)
for d in ('png', 'pdf'): fig.savefig(os.path.join(BASE, 'figures', d, f'SuppFig_S2_GH18_taxonomy_besthit.{d}'), dpi=300)
print(t4[t4.genus.isin(cats)].pivot(index='genus', columns='Cohort', values='share').reindex(cats).round(3))
print('seed-assigned genes:', (gc.genus != 'unassigned').sum, '/', len(gc), '; Kraken2-classified contigs:', (gc.kraken_domain != 'unclassified').sum)
print(t4.groupby('Cohort').apply(lambda d: d.groupby('phylum').share.sum.sort_values(ascending=False).head(4).round(3).to_dict))
