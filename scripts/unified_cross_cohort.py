#!/usr/bin/env python3
"""
Cross-cohort tables from the UNIFIED pipeline (all 32 samples processed identically:
fastp -> Bowtie2 host removal -> MEGAHIT -> Prodigal -> CD-HIT -> eggNOG-mapper /
Bowtie2 read counts -> Kraken2+Bracken conf 0.1). Inputs are the per-cohort outputs
collected in results/unified/. Replaces the vendor-derived tables 01/10/11/12/16 for
every cross-diet statement in the manuscript.

Outputs (tables/):
  35_unified_cazy_per_sample.tsv     CAZy family read counts + counts per million mapped (CPM)
  36_unified_cazy_tests.tsv          Mann-Whitney between every cohort pair, on CPM, for GH18/GH19/GH20/CBM50/AA10
  37_unified_phylum.tsv, 38_unified_genus.tsv   Bracken relative abundance (fraction of classified reads)
  39_unified_alpha_diversity.tsv     Shannon/Simpson/richness at genus level + Kruskal-Wallis / pairwise MW
  40_unified_permanova.tsv           Bray-Curtis PERMANOVA (999 perms) on genus profiles, 3-group + pairwise
  41_unified_gh18_vs_bacteroides.tsv GH18 CPM vs Bacteroides fraction, all 32 samples
  42_unified_sample_summary.tsv      reads, host-removed, mapping rate, classified fraction per sample
"""
import os, re, itertools, numpy as np, pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# UNIFIED_DIR / TABLES_DIR override: used by scripts/reanalysis_2026-11/build_unified32_tables.sh to run the
# identical statistics on the single 32-sample catalogue (results/unified32 -> tables/unified32/).
U = os.environ.get('UNIFIED_DIR', os.path.join(BASE, 'results', 'unified')); T = os.environ.get('TABLES_DIR', os.path.join(BASE, 'tables'))
# The manuscript's cross-diet numbers (tables/35-46: GH18 CPM 471/665/317, dbCAN3) come from
# results/unified32_dbcan. The defaults above point at results/unified — the OLDER per-cohort run —
# so a bare `python3 scripts/unified_cross_cohort.py` used to overwrite tables/35-46 with per-cohort
# eggNOG numbers (390/355/291) and no warning. Writing tables/ now requires the dbCAN input explicitly.
if os.path.realpath(T) == os.path.realpath(os.path.join(BASE, 'tables')) and not os.path.basename(U.rstrip('/')).endswith('unified32_dbcan'):
    raise SystemExit(
        'refusing to write the manuscript tables from ' + os.path.relpath(U, BASE) + '.\n'
        'tables/35-46 are the dbCAN3 single-catalogue set. Rebuild them with:\n'
        '    CAZY_SOURCE=dbcan bash scripts/reanalysis_2026-11/build_unified32_tables.sh\n'
        'To run these statistics on another input, send the output elsewhere:\n'
        '    UNIFIED_DIR=results/unified32 TABLES_DIR=tables/unified32 python3 scripts/unified_cross_cohort.py')
os.makedirs(T, exist_ok=True)
COH = {'DS01': 'P. eryngii', 'DS02': 'H. tessellatus', 'DS03': 'P. ostreatus'}
FAM = ['GH18', 'GH19', 'GH20', 'CBM5', 'CBM12', 'CBM50', 'AA10']
rng = np.random.default_rng(1)

# ---- CAZy per sample ----
rows = []
for c in COH:
    t = pd.read_csv(os.path.join(U, f'{c}_cazy_abundance.tsv'), sep='\t', comment='#').set_index('Family')
    samples = [s for s in t.columns if s != 'Total']
    tot = t.loc['TOTAL_MAPPED', samples]
    for s in samples:
        r = dict(Cohort=c, Diet=COH[c], Sample=s, mapped_reads=int(tot[s]))
        for f in FAM:
            v = int(t.loc[f, s]) if f in t.index else 0
            r[f] = v; r[f + '_CPM'] = v / tot[s] * 1e6
        rows.append(r)
caz = pd.DataFrame(rows)
caz.to_csv(os.path.join(T, '35_unified_cazy_per_sample.tsv'), sep='\t', index=False)

tests = []
for f in ['GH18', 'GH19', 'GH20', 'CBM50', 'AA10']:
    for a, b in itertools.combinations(COH, 2):
        x = caz[caz.Cohort == a][f + '_CPM']; y = caz[caz.Cohort == b][f + '_CPM']
        Uv, p = stats.mannwhitneyu(x, y, alternative='two-sided')
        tests.append(dict(Family=f, A=a, B=b, A_mean_CPM=x.mean, A_sd=x.std(ddof=1), B_mean_CPM=y.mean, B_sd=y.std(ddof=1),
                          fold_A_over_B=x.mean / y.mean, U=Uv, P=p, rank_biserial=1 - 2 * Uv / (len(x) * len(y))))
    kw = stats.kruskal(*[caz[caz.Cohort == c][f + '_CPM'] for c in COH])
    tests.append(dict(Family=f, A='all three', B='Kruskal-Wallis', U=kw.statistic, P=kw.pvalue))
pd.DataFrame(tests).to_csv(os.path.join(T, '36_unified_cazy_tests.tsv'), sep='\t', index=False)

# ---- taxonomy ----
def load_bracken(level):
    frames = []
    for c in COH:
        b = pd.read_csv(os.path.join(U, f'{c}_bracken_{level}.tsv'), sep='\t', header=None, names=['Sample', 'Taxon', 'reads', 'frac'])
        b['Cohort'] = c; frames.append(b)
    b = pd.concat(frames)
    wide = b.pivot_table(index='Taxon', columns='Sample', values='frac', fill_value=0.0)
    return b, wide
_, phy = load_bracken('P'); braw, gen = load_bracken('G')
order = list(caz.Sample)
phy = phy.reindex(columns=order, fill_value=0.0); gen = gen.reindex(columns=order, fill_value=0.0)
phy.loc[phy.mean(axis=1).sort_values(ascending=False).index].to_csv(os.path.join(T, '37_unified_phylum.tsv'), sep='\t')
gen.loc[gen.mean(axis=1).sort_values(ascending=False).index].to_csv(os.path.join(T, '38_unified_genus.tsv'), sep='\t')

# ---- alpha diversity (genus, Bracken fractions) ----
alpha = []
for s in order:
    p = gen[s].values; p = p[p > 0]; p = p / p.sum
    alpha.append(dict(Sample=s, Cohort=caz.set_index('Sample').loc[s, 'Cohort'], Shannon=-(p * np.log(p)).sum,
                      Simpson=1 - (p ** 2).sum, Genera_observed=int((gen[s] > 0).sum)))
alpha = pd.DataFrame(alpha)
arows = []
for m in ['Shannon', 'Simpson', 'Genera_observed']:
    kw = stats.kruskal(*[alpha[alpha.Cohort == c][m] for c in COH])
    arows.append(dict(Metric=m, Comparison='Kruskal-Wallis 3 cohorts', statistic=kw.statistic, P=kw.pvalue, means='; '.join(f'{c} {alpha[alpha.Cohort==c][m].mean:.2f}±{alpha[alpha.Cohort==c][m].std(ddof=1):.2f}' for c in COH)))
    for a, b in itertools.combinations(COH, 2):
        Uv, p = stats.mannwhitneyu(alpha[alpha.Cohort == a][m], alpha[alpha.Cohort == b][m], alternative='two-sided')
        arows.append(dict(Metric=m, Comparison=f'{a} vs {b} Mann-Whitney', statistic=Uv, P=p, means=''))
alpha.to_csv(os.path.join(T, '39_unified_alpha_diversity.tsv'), sep='\t', index=False)
with open(os.path.join(T, '39_unified_alpha_diversity.tsv'), 'a') as fh:
    fh.write('\n# Tests\n'); pd.DataFrame(arows).to_csv(fh, sep='\t', index=False)

# ---- PERMANOVA ----
def permanova(D, groups, nperm=999):
    n = len(groups); g = np.asarray(groups); labels = np.unique(g)
    def F_of(gv):
        sst = (D ** 2).sum / (2 * n); ssw = 0
        for lab in labels:
            idx = np.where(gv == lab)[0]; sub = D[np.ix_(idx, idx)]; ssw += (sub ** 2).sum / (2 * len(idx))
        a = len(labels); return ((sst - ssw) / (a - 1)) / (ssw / (n - a)), 1 - ssw / sst
    F, R2 = F_of(g); cnt = sum(F_of(rng.permutation(g))[0] >= F for _ in range(nperm))
    return F, R2, (cnt + 1) / (nperm + 1)
X = gen.T.loc[order].values; D = squareform(pdist(X, 'braycurtis'))
grp = caz.set_index('Sample').loc[order, 'Cohort'].values
prow = []
F, R2, p = permanova(D, grp); prow.append(dict(Comparison='3 cohorts', pseudo_F=F, R2=R2, P=p, n=len(order)))
for a, b in itertools.combinations(COH, 2):
    idx = np.where(np.isin(grp, [a, b]))[0]; F, R2, p = permanova(D[np.ix_(idx, idx)], grp[idx])
    prow.append(dict(Comparison=f'{a} vs {b}', pseudo_F=F, R2=R2, P=p, n=len(idx)))
# dispersion check: mean distance to group centroid (PCoA-free approximation: mean within-group distance)
for c in COH:
    idx = np.where(grp == c)[0]; sub = D[np.ix_(idx, idx)]
    prow.append(dict(Comparison=f'{c} mean within-group Bray-Curtis', pseudo_F=np.nan, R2=np.nan, P=np.nan, n=len(idx), note=sub[np.triu_indices(len(idx), 1)].mean))
pd.DataFrame(prow).to_csv(os.path.join(T, '40_unified_permanova.tsv'), sep='\t', index=False)

# ---- GH18 vs Bacteroides ----
bact = gen.loc['Bacteroides'] if 'Bacteroides' in gen.index else pd.Series(0, index=order)
g2 = caz.set_index('Sample').loc[order, ['Cohort', 'GH18_CPM', 'mapped_reads']].copy
g2['Bacteroides_frac'] = bact.loc[order].values
g2['GH18_CPM_per_Bacteroides_frac'] = g2.GH18_CPM / g2.Bacteroides_frac.replace(0, np.nan)
rho, pr = stats.spearmanr(g2.GH18_CPM, g2.Bacteroides_frac)
brow = [dict(metric='Spearman GH18_CPM vs Bacteroides fraction (n=32)', value=rho, P=pr)]
for a, b in itertools.combinations(COH, 2):
    for m in ['GH18_CPM', 'Bacteroides_frac', 'GH18_CPM_per_Bacteroides_frac']:
        x = g2[g2.Cohort == a][m].dropna; y = g2[g2.Cohort == b][m].dropna
        Uv, p = stats.mannwhitneyu(x, y, alternative='two-sided')
        brow.append(dict(metric=f'{m}: {a} vs {b}', value=x.mean / y.mean, P=p, A_mean=x.mean, B_mean=y.mean))
g2.reset_index.to_csv(os.path.join(T, '41_unified_gh18_vs_bacteroides.tsv'), sep='\t', index=False)
with open(os.path.join(T, '41_unified_gh18_vs_bacteroides.tsv'), 'a') as fh:
    fh.write('\n# Tests (value = fold A/B or rho)\n'); pd.DataFrame(brow).to_csv(fh, sep='\t', index=False)

# ---- sample summary ----
ms = pd.read_csv(os.path.join(U, 'mapping_summary.tsv'), sep='\t', header=None, names=['Cohort', 'Sample', 'nonhost_read_pairs', 'catalogue_mapping_rate'])
cl = braw.groupby('Sample').reads.sum   # Bracken genus-assigned reads
ms['bracken_genus_assigned_reads'] = ms.Sample.map(cl)
ms.to_csv(os.path.join(T, '42_unified_sample_summary.tsv'), sep='\t', index=False)

# ---- console summary ----
print('GH18 CPM by cohort:'); print(caz.groupby('Cohort').GH18_CPM.agg(['mean', 'std', 'min', 'max', 'count']).round(1))
print(pd.DataFrame(tests).query("Family=='GH18'")[['A', 'B', 'fold_A_over_B', 'P']].to_string(index=False))
print('\nAlpha:'); print(pd.DataFrame(arows)[['Metric', 'Comparison', 'P', 'means']].to_string(index=False))
print('\nPERMANOVA:'); print(pd.DataFrame(prow).to_string(index=False))
print('\nBacteroides:'); print(pd.DataFrame(brow).to_string(index=False))
print('\nTop phyla mean per cohort:')
for c in COH:
    ss = caz[caz.Cohort == c].Sample; print(c, phy[ss].mean(axis=1).sort_values(ascending=False).head(5).round(3).to_dict)
print('\nTop genera mean per cohort:')
for c in COH:
    ss = caz[caz.Cohort == c].Sample; print(c, gen[ss].mean(axis=1).sort_values(ascending=False).head(6).round(3).to_dict)

# ---- composition-matched indices on the unified annotation (CDI, chitin:cellulose) ----
irows = []
for c in COH:
    t = pd.read_csv(os.path.join(U, f'{c}_cazy_abundance.tsv'), sep='\t', comment='#').set_index('Family')
    samples = [s for s in t.columns if s != 'Total']
    gh = t.loc[[i for i in t.index if re.fullmatch(r'GH\d+', i)], samples].sum
    chit = t.loc[[f for f in ['GH18', 'GH19', 'GH20'] if f in t.index], samples].sum
    cell = t.loc[[f for f in ['GH5', 'GH6', 'GH9', 'GH48'] if f in t.index], samples].sum
    for s in samples:
        irows.append(dict(Cohort=c, Sample=s, total_GH_reads=int(gh[s]), CDI_pct=100 * chit[s] / gh[s],
                          chitin_to_cellulose=t.loc['GH18', s] / cell[s]))
idx = pd.DataFrame(irows)
itests = []
for m in ['CDI_pct', 'chitin_to_cellulose']:
    kw = stats.kruskal(*[idx[idx.Cohort == c][m] for c in COH])
    itests.append(dict(Index=m, Comparison='Kruskal-Wallis', statistic=kw.statistic, P=kw.pvalue,
                       means='; '.join(f'{c} {idx[idx.Cohort==c][m].mean:.3f}' for c in COH)))
    for a, b in itertools.combinations(COH, 2):
        Uv, p = stats.mannwhitneyu(idx[idx.Cohort == a][m], idx[idx.Cohort == b][m], alternative='two-sided')
        itests.append(dict(Index=m, Comparison=f'{a} vs {b}', statistic=Uv, P=p, means=''))
idx.to_csv(os.path.join(T, '43_unified_cell_wall_indices.tsv'), sep='\t', index=False)
with open(os.path.join(T, '43_unified_cell_wall_indices.tsv'), 'a') as fh:
    fh.write('\n# CDI = (GH18+GH19+GH20)/all GH reads; chitin:cellulose = GH18/(GH5+GH6+GH9+GH48). Tests:\n')
    pd.DataFrame(itests).to_csv(fh, sep='\t', index=False)
print('\nIndices:'); print(pd.DataFrame(itests).to_string(index=False))

# ---- KEGG chitin-pathway KOs per sample (CPM) ----
KO_NAMES = {'K01183': 'chitinase', 'K12373': 'HEXA_B hexosaminidase', 'K01207': 'nagZ hexosaminidase', 'K01443': 'nagA', 'K02564': 'nagB',
            'K00884': 'NAGK GlcNAc kinase', 'K01452': 'chitin deacetylase', 'K03478': 'chbG deacetylase', 'K03431': 'glmM', 'K04042': 'glmU', 'K02804': 'nagE PTS transporter'}
krows = []
for c in COH:
    t = pd.read_csv(os.path.join(U, f'{c}_ko_abundance.tsv'), sep='\t').set_index('KO'); samples = list(t.columns); tot = t.loc['TOTAL_MAPPED']
    for s_ in samples:
        r = dict(Cohort=c, Sample=s_)
        for k in KO_NAMES: r[k + '_CPM'] = t.loc[k, s_] / tot[s_] * 1e6
        krows.append(r)
ko = pd.DataFrame(krows); ko.to_csv(os.path.join(T, '45_unified_ko_per_sample.tsv'), sep='\t', index=False)
kt = []
for k, nm in KO_NAMES.items:
    kw = stats.kruskal(*[ko[ko.Cohort == c][k + '_CPM'] for c in COH])
    kt.append(dict(KO=k, name=nm, DS01_mean=ko[ko.Cohort=='DS01'][k+'_CPM'].mean, DS02_mean=ko[ko.Cohort=='DS02'][k+'_CPM'].mean, DS03_mean=ko[ko.Cohort=='DS03'][k+'_CPM'].mean, KW_P=kw.pvalue,
                   **{f'P_{a}_vs_{b}': stats.mannwhitneyu(ko[ko.Cohort==a][k+'_CPM'], ko[ko.Cohort==b][k+'_CPM'], alternative='two-sided').pvalue for a, b in itertools.combinations(COH, 2)}))
pd.DataFrame(kt).to_csv(os.path.join(T, '46_unified_ko_tests.tsv'), sep='\t', index=False)
print('\nKO CPM means:'); print(pd.DataFrame(kt).round(4).to_string(index=False))
