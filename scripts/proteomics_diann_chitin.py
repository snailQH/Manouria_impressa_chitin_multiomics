#!/usr/bin/env python3
"""Chitin-related proteins in the species-specific DIA-NN search (36 runs, two animals), replacing the proxy-database
analyses (tables 02/07/18/19/22/31/47, Figures 5C/5D/6, Supp S13, S17).
Inputs: results/reanalysis_2026-11/diann/report.{pg_matrix,pr_matrix}.tsv (DIA-NN 1.8.1, MaxLFQ, 1 % FDR, MBR),
        results/reanalysis_2026-11/mimp_v3_annotation/{Mimp_v3.eggnog.annotations.tsv.gz, gene_id_map.tsv, gh18_locus_proteins.faa},
        results/reanalysis_2026-11/diann/meta_detected_annotations.tsv (eggNOG rows of the detected metagenome proteins).
Categories: GH18 locus (chr4 models MIMP3_g08468-g08479; catalytic vs inactive CLP), other host GH18 (CHID1, CTBS, OVGP1),
GH20 hexosaminidase, GlcNAc metabolism (NAGK, GNPDA, GNPNAT1, GFPT, PGM3, UAP1, AMDHD2, RENBP), metagenome (any detected
catalogue protein, labelled by its CAZy/KO).
Outputs: tables/58_diann_chitin_related_proteins.tsv, 59_diann_spatial_by_segment.tsv, 60_diann_peak_segment.tsv,
         61_diann_spatial_tests.tsv, 62_diann_peptides_per_locus.tsv; figures Fig5C, Fig5D, Fig6, SuppFig S13–S17."""
import os, sys, re, gzip, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy import stats
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); D = os.path.join(BASE, 'results', 'reanalysis_2026-11', 'diann')
A = os.path.join(BASE, 'results', 'reanalysis_2026-11', 'mimp_v3_annotation'); T = os.path.join(BASE, 'tables'); FIG = os.path.join(BASE, 'figures')
def save(fig, name):
    from figure_ownership import may_write
    if not may_write(name, __file__): return
    for d in ('png', 'pdf'):
        os.makedirs(os.path.join(FIG, d), exist_ok=True)   # a fresh clone has no figures/ yet
        fig.savefig(os.path.join(FIG, d, f'{name}.{d}'), dpi=300, bbox_inches='tight')
    plt.close(fig); print('  ', name)
# ---------------- runs
run_meta = pd.read_csv(os.path.join(T, '50_diann_mimpv3_run_summary.tsv'), sep='\t', comment='#').set_index('run')
SEG = ['Stomach', 'Duodenum', 'Jejunum', 'Ileum', 'Colon', 'Rectum']; SI, HG = ['Duodenum', 'Jejunum', 'Ileum'], ['Colon', 'Rectum']
REGION = {'Stomach': ['Stomach'], 'Small intestine': SI, 'Hindgut': HG}
runof = lambda c: re.search(r'(25NP\d+)', c).group(1)
# ---------------- host annotation
idmap = pd.read_csv(os.path.join(A, 'gene_id_map.tsv'), sep='\t'); aug2rel = dict(zip(idmap.augustus_id, idmap.release_id))
ann = pd.read_csv(gzip.open(os.path.join(A, 'Mimp_v3.eggnog.annotations.tsv.gz'), 'rt'), sep='\t', comment='#', header=None,
                  names=['query', 'seed', 'evalue', 'score', 'OGs', 'lvl', 'COG', 'Description', 'Preferred_name', 'GOs', 'EC', 'KO', 'Pathway', 'Module', 'Reaction', 'rclass', 'BRITE', 'TC', 'CAZy', 'BiGG', 'PFAMs'])
ann['gene'] = ann['query'].str.replace(r'\.t\d+$', '', regex=True).map(aug2rel); ann = ann.dropna(subset=['gene']).drop_duplicates('gene').set_index('gene')
# a silently failed id join returns '-' for every annotation lookup and quietly empties the
# non-locus categories, so require that most augustus ids actually mapped to a release id
_matched = ann['query'].str.replace(r'\.t\d+$', '', regex=True).isin(aug2rel).mean if len(ann) else 0
assert len(ann) > 10000, f'host annotation join matched only {len(ann)} genes'
LOCUS = {'MIMP3_g08468': ('g08468 TAF8 (v1 088402 N-term; not GH18)', 'other'), 'MIMP3_g08469': ('088402a GH18+CBM14', 'catalytic'), 'MIMP3_g08470': ('088402b GH18', 'catalytic'), 'MIMP3_g08471': ('g08471 (GH18 fragment)', 'catalytic'),
         'MIMP3_g08472': ('L7 CHIA array', 'catalytic'), 'MIMP3_g08473': ('L6 CHIA array', 'catalytic'), 'MIMP3_g08474': ('L5 CHIA array', 'catalytic'), 'MIMP3_g08475': ('L4 CHIA array', 'catalytic'),
         'MIMP3_g08476': ('L3 CHIA-related', 'catalytic'), 'MIMP3_g08477': ('088393 fragment (g08477)', 'inactive CLP'), 'MIMP3_g08478': ('L2 088393 CLP', 'inactive CLP'), 'MIMP3_g08479': ('L1 088392 CHIT1-clade', 'catalytic')}
GLCNAC = {'NAGK', 'GNPDA1', 'GNPDA2', 'GNPNAT1', 'GFPT1', 'GFPT2', 'PGM3', 'UAP1', 'AMDHD2', 'RENBP'}
def host_category(genes):
    locus = [g for g in genes if g in LOCUS and LOCUS[g][1] != 'other']   # g08468 (TAF8) is in the interval but is not a chitinase
    if locus: return 'GH18 locus (' + LOCUS[locus[0]][1] + ')', ' + '.join(LOCUS[g][0] for g in locus)
    names = {str(ann.Preferred_name.get(g, '-')) for g in genes}; caz = {str(ann.CAZy.get(g, '-')) for g in genes}; pf = ' '.join(str(ann.PFAMs.get(g, '')) for g in genes)
    if 'GH18' in ' '.join(caz) or 'Glyco_hydro_18' in pf: return 'other host GH18', '/'.join(sorted(names - {'-'})) or genes[0]
    if 'GH20' in ' '.join(caz) or names & {'HEXA', 'HEXB', 'HEXDC', 'HEXDCL'}: return 'GH20 hexosaminidase', '/'.join(sorted(names - {'-'})) or genes[0]
    if names & GLCNAC: return 'GlcNAc metabolism', '/'.join(sorted(names & GLCNAC))
    return None, None
# ---------------- metagenome annotation
# 21 eggNOG columns. `ann.columns[:-1]` is 20 names (PFAMs is dropped once 'gene' becomes the index), and with
# header=None pandas then took column 1 as the index and shifted every name one place left — CAZy held BiGG,
# and the index was the seed ortholog, so every lookup by protein id missed and the whole metagenome half of
# this analysis read '-'. Name all 21 columns explicitly.
EGGNOG_COLS = ['query', 'seed', 'evalue', 'score', 'OGs', 'lvl', 'COG', 'Description', 'Preferred_name', 'GOs', 'EC',
               'KO', 'Pathway', 'Module', 'Reaction', 'rclass', 'BRITE', 'TC', 'CAZy', 'BiGG', 'PFAMs']
mann = pd.read_csv(os.path.join(D, 'meta_detected_annotations.tsv'), sep='\t', header=None, names=EGGNOG_COLS)
assert mann['query'].str.contains('k141').all, 'meta annotation column alignment'
mann = mann.set_index('query')
def meta_label(ids):
    cz = {str(mann.CAZy.get(i, '-')) for i in ids} - {'-'}; ko = {str(mann.KO.get(i, '-')) for i in ids} - {'-'}; nm = {str(mann.Preferred_name.get(i, '-')) for i in ids} - {'-'}
    return ','.join(sorted(cz)) or ','.join(sorted(k.replace('ko:', '') for k in ko)) or '-', '/'.join(sorted(nm)) or '-'
# ---------------- protein groups
pg = pd.read_csv(os.path.join(D, 'report.pg_matrix.tsv'), sep='\t'); cols = [c for c in pg.columns if c.endswith('.mzML')]
M = pg.set_index('Protein.Group')[cols]; M.columns = [runof(c) for c in cols]; M = np.log2(M.astype(float))
rows = []
for pgid in M.index:
    src = pgid.split('|')[0]; ids = [x.split('|')[1] for x in pgid.split(';')]
    if src == 'host':
        genes = [re.sub(r'\.t\d+$', '', i) for i in ids]; cat, lab = host_category(genes)
        if cat is None: continue
        fam = str(ann.CAZy.get(genes[0], '-')); ko = str(ann.KO.get(genes[0], '-')).replace('ko:', '')
    elif src == 'meta':
        cat = 'metagenome'; fam, lab = meta_label(ids); ko = ','.join(sorted({str(mann.KO.get(i, '-')).replace('ko:', '') for i in ids} - {'-'}))
    else: continue
    r = dict(protein_group=pgid, source=src, category=cat, label=lab, CAZy=fam, KO=ko, n_runs_detected=int(M.loc[pgid].notna.sum))
    for run in M.columns: r[f'{run}|{run_meta.animal[run]}|{run_meta.tissue[run]}_{run_meta.piece[run]}'] = round(M.loc[pgid, run], 3) if pd.notna(M.loc[pgid, run]) else np.nan
    rows.append(r)
chit = pd.DataFrame(rows)
# peptides per protein group (precursor matrix: proteotypic stripped sequences)
pr = pd.read_csv(os.path.join(D, 'report.pr_matrix.tsv'), sep='\t', usecols=['Protein.Group', 'Proteotypic', 'Stripped.Sequence'])
pep = pr.groupby('Protein.Group').agg(peptides=('Stripped.Sequence', 'nunique'), proteotypic_peptides=('Proteotypic', lambda x: int(((x == 1)).sum)))
chit = chit.merge(pep, left_on='protein_group', right_index=True, how='left')
order = ['GH18 locus (catalytic)', 'GH18 locus (inactive CLP)', 'other host GH18', 'GH20 hexosaminidase', 'GlcNAc metabolism', 'metagenome']
chit['category'] = pd.Categorical(chit.category, order); chit = chit.sort_values(['category', 'label']).reset_index(drop=True)
with open(os.path.join(T, '58_diann_chitin_related_proteins.tsv'), 'w') as fh:
    fh.write('# Chitin-related protein groups in the species-specific DIA-NN search (36 runs; log2 MaxLFQ per run, columns run|animal|tissue_piece). Host categories from the Mimp_v3 eggNOG annotation; metagenome groups labelled by CAZy family / KO of the catalogue proteins. peptides = distinct stripped sequences; proteotypic = unique to the group.\n')
    chit.to_csv(fh, sep='\t', index=False)
print('chitin-related groups:', chit.category.value_counts.to_dict)
# ---------------- per-animal segment means, peaks, paired tests
runcols = [c for c in chit.columns if '|' in c]
def animal_matrix(animal):
    cs = [c for c in runcols if c.split('|')[1] == animal]; X = chit.set_index('protein_group')[cs]
    tissue = pd.Series({c: c.split('|')[2].rsplit('_', 1)[0] for c in cs})
    seg = pd.DataFrame({t: X[[c for c in cs if tissue[c] == t]].mean(axis=1) for t in sorted(set(tissue))})
    return seg
seg_rows, peak_rows, test_rows = [], [], []
for animal in ['MI01', 'MI02']:
    seg = animal_matrix(animal); segs = [s for s in SEG if s in seg.columns]; nongut = [c for c in seg.columns if c not in SEG]
    seg = seg.join(chit.set_index('protein_group')[['category', 'label']])
    for cat in order:
        sub = seg[seg.category == cat][segs].dropna(thresh=4)
        if len(sub) == 0: continue
        for s in segs: seg_rows.append(dict(animal=animal, category=cat, segment=s, n_proteins=int(sub[s].notna.sum), mean_log2=round(sub[s].mean, 3), sd=round(sub[s].std, 3)))
        cen = sub.sub(sub.mean(axis=1), axis=0)
        si = cen[[s for s in SI if s in cen]].mean(axis=1); hg = cen[[s for s in HG if s in cen]].mean(axis=1); st = cen['Stomach'] if 'Stomach' in cen else None
        ok = si.notna & hg.notna; d = (si - hg)[ok]
        if len(d) >= 3:
            w = stats.wilcoxon(d, alternative='two-sided'); full = cen.dropna; fr = stats.friedmanchisquare(*[full[s] for s in segs]) if len(full) >= 3 else None
            test_rows.append(dict(animal=animal, category=cat, n_proteins=int(ok.sum), SI_minus_hindgut_median_log2=round(d.median, 3), n_higher_in_SI=int((d > 0).sum), wilcoxon_P=round(w.pvalue, 4),
                                  friedman_chi2=round(fr.statistic, 2) if fr else np.nan, friedman_P=round(fr.pvalue, 4) if fr else np.nan,
                                  SI_minus_stomach_median_log2=round((si - st)[ok & st.notna].median, 3) if st is not None else np.nan,
                                  wilcoxon_P_vs_stomach=round(stats.wilcoxon((si - st)[ok & st.notna]).pvalue, 4) if st is not None and (ok & st.notna).sum >= 3 else np.nan,
                                  gut_minus_nongut_median_log2=round((sub.mean(axis=1) - seg.loc[sub.index, nongut].mean(axis=1)).median, 3) if nongut else np.nan))
        for pgid, r in sub.iterrows:
            peak_rows.append(dict(animal=animal, category=cat, label=seg.label[pgid], protein_group=pgid, peak_segment=r.idxmax, peak_region='Stomach' if r.idxmax == 'Stomach' else 'Small intestine' if r.idxmax in SI else 'Hindgut',
                                  SI_minus_hindgut_log2=round(float(si.get(pgid, np.nan) - hg.get(pgid, np.nan)), 3), **{f'{s}_log2': round(r[s], 2) for s in segs}))
pd.DataFrame(seg_rows).to_csv(os.path.join(T, '59_diann_spatial_by_segment.tsv'), sep='\t', index=False)
pd.DataFrame(peak_rows).to_csv(os.path.join(T, '60_diann_peak_segment.tsv'), sep='\t', index=False)
tests = pd.DataFrame(test_rows)
with open(os.path.join(T, '61_diann_spatial_tests.tsv'), 'w') as fh:
    fh.write('# Paired per-protein spatial statistics on the species-specific DIA-NN quantities: each protein group centred on its own segment mean (MI01: mean of 3 pieces per segment; MI02: mean of 2); small intestine = duodenum+jejunum+ileum, hindgut = colon+rectum; Wilcoxon signed-rank two-sided; Friedman across the six segments on proteins quantified in all; gut vs non-gut (MI02: liver, heart, muscle).\n')
    tests.to_csv(fh, sep='\t', index=False)
print(tests.to_string(index=False))
# peptides per locus model
locus_rows = []
pr2 = pd.read_csv(os.path.join(D, 'report.pr_matrix.tsv'), sep='\t', usecols=['Protein.Group', 'Protein.Ids', 'Proteotypic', 'Stripped.Sequence'])
for g, (lab, kind) in LOCUS.items:
    sub = pr2[pr2['Protein.Ids'].str.contains(g)]; grp = pr2[pr2['Protein.Group'].str.contains(g)]
    locus_rows.append(dict(model=g, locus=lab, clade=kind, peptides_any_group=int(sub['Stripped.Sequence'].nunique), peptides_in_own_group=int(grp['Stripped.Sequence'].nunique), proteotypic=int(grp[grp.Proteotypic == 1]['Stripped.Sequence'].nunique), protein_group=';'.join(grp['Protein.Group'].unique[:1])))
loc = pd.DataFrame(locus_rows); loc.to_csv(os.path.join(T, '62_diann_peptides_per_locus.tsv'), sep='\t', index=False); print(loc.to_string(index=False))

# ---------------- Figures
COLC = {'GH18 locus (catalytic)': '#B5541C', 'GH18 locus (inactive CLP)': '#C9A227', 'other host GH18': '#8A6D3B', 'GH20 hexosaminidase': '#2E6F9E', 'GlcNAc metabolism': '#5A7A58', 'metagenome': '#7A7A7A'}
# Fig 5C: categories. The 213 metagenome groups are now split by annotation — lumping them in one
# grey bar hid the fact that some of them are chitin-active, and the bar swamped the host categories.
mg = chit[chit.category == 'metagenome']
mg_gh18 = int(mg.KO.astype(str).str.contains('K01183').sum); mg_gh20 = int(mg.KO.astype(str).str.contains('K12373').sum)
bars = [(c, int((chit.category == c).sum), COLC[c]) for c in order if c != 'metagenome']
bars += [('metagenome GH18 (K01183)', mg_gh18, '#9C6B4F'), ('metagenome GH20 (K12373)', mg_gh20, '#4F7FA3'),
         ('metagenome, other CAZy/KO', len(mg) - mg_gh18 - mg_gh20, '#9A9A9A')]
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.barh(range(len(bars)), [b[1] for b in bars], color=[b[2] for b in bars], edgecolor='#333', linewidth=0.5)
for i, b in enumerate(bars): ax.text(b[1] + 2, i, f'{b[1]}', va='center', fontsize=9)
ax.set_yticks(range(len(bars))); ax.set_yticklabels([b[0] for b in bars], fontsize=9); ax.invert_yaxis
ax.set_xlabel('protein groups detected (36 runs, two animals)'); ax.set_xscale('symlog', linthresh=10)
ax.set_xticks([0, 2, 5, 10, 50, 200]); ax.set_xticklabels(['0', '2', '5', '10', '50', '200'])
ax.axhline(len(order) - 1.5, color='#BBB', lw=0.8, ls='--')
ax.text(ax.get_xlim[1], len(order) - 1.7, 'host  ', ha='right', va='bottom', fontsize=7.5, color='#777')
ax.text(ax.get_xlim[1], len(order) - 1.3, 'microbial  ', ha='right', va='top', fontsize=7.5, color='#777')
ax.set_title('Chitin-related protein groups, species-specific DIA-NN search\n(Mimp_v3 proteome + metagenome CAZy/ko00520 proteins; 1 % FDR; symlog x)', fontweight='bold', fontsize=10)
ax.spines[['top', 'right']].set_visible(False)
save(fig, 'Fig7A_chitin_protein_categories')
# Fig 5D: GH18 locus groups by tissue, both animals
locus_pgs = chit[chit.category.astype(str).str.startswith('GH18 locus')].set_index('protein_group')
fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={'width_ratios': [6, 9]})
for ax, animal, tissues in [(axes[0], 'MI01', SEG), (axes[1], 'MI02', SEG + ['Liver', 'Heart', 'Muscle'])]:
    seg = animal_matrix(animal); x = np.arange(len(tissues)); n = len(locus_pgs); w = 0.8 / n
    for i, (pgid, r) in enumerate(locus_pgs.iterrows):
        vals = [seg.loc[pgid, t] if t in seg.columns else np.nan for t in tissues]
        # dots, not bars: a bar encodes magnitude from zero, and log2 MaxLFQ has no meaningful zero —
        # with the axis starting at 10 the bars also hid ~46 % of their own length
        ax.plot(x + (i - n / 2 + 0.5) * w, vals, marker='o', ms=5, ls='', label=r.label,
                color=plt.cm.tab10(i % 10), mec='#333', mew=0.4)
    ax.set_xticks(x); ax.set_xticklabels(tissues, rotation=30, ha='right', fontsize=9); ax.set_ylabel('log$_2$ MaxLFQ (mean of pieces)'); ax.grid(axis='y', color='#EEE', lw=0.6); ax.set_axisbelow(True)
    ax.set_title(f'{animal} — GH18-locus protein groups' + (' (6 GI segments × 3 pieces)' if animal == 'MI01' else ' (6 GI segments + 3 non-gut tissues × 2 pieces)'), fontsize=10, fontweight='bold'); ax.spines[['top', 'right']].set_visible(False)
axes[1].legend(fontsize=7, ncol=4, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.28)); fig.tight_layout; save(fig, 'Fig7B_tissue_chitinase_expression')
# ---------------- microbial chitin-active proteins along the gut (tables/73)
# Asked 2026-11-16: the search database carried 65,450 metagenome proteins, so the microbial side of the
# relay is testable from these runs. It had read as "search background" only because the annotation join
# above was broken. Non-gut tissue (liver/heart/muscle, MI02) is the negative control: a microbial protein
# group detected there is background, not gut microbiota.
NONGUT = ['Liver', 'Heart', 'Muscle']
meta_rows = chit[chit.source == 'meta'].set_index('protein_group')
runs_of = lambda an, segs: [r for r in M.columns if run_meta.animal[r] == an and run_meta.tissue[r] in segs]
micro = {}
for fam, ko in (('GH18', 'K01183'), ('GH20', 'K12373')):
    sel = meta_rows[meta_rows.KO.astype(str).str.contains(ko, na=False)]
    nongut_runs = runs_of('MI02', NONGUT)
    keep, bg, euk = [], [], []
    for pgid in sel.index:
        ids = [x.split('|')[1] for x in pgid.split(';')]
        ogs = ' '.join(str(mann.OGs.get(i, '')) for i in ids)
        if '|Eukaryota' in ogs and '|Bacteria' not in ogs and '|Archaea' not in ogs:
            euk.append(pgid)                      # e.g. Hymenoptera CHIA: dietary/environmental, not gut bacteria
        elif M.loc[pgid, nongut_runs].notna.any:
            bg.append(pgid)                       # present in liver/heart/muscle: search background
        else:
            keep.append(pgid)
    micro[fam] = dict(groups=keep, background=bg, eukaryotic=euk)
mrows = []
for fam, d in micro.items:
    for an in ('MI01', 'MI02'):
        for reg, segs in REGION.items:
            rs = runs_of(an, segs)
            if not rs or not d['groups']: continue
            det = int(M.loc[d['groups'], rs].notna.sum.sum); tot = len(rs) * len(d['groups'])
            mrows.append(dict(family=fam, animal=an, region=reg, n_groups=len(d['groups']), n_runs=len(rs),
                              detections=det, group_runs=tot, detection_pct=round(100 * det / tot, 1)))
micro_det = pd.DataFrame(mrows)
# is detection frequency different between regions? (counts, so a chi-square on detected/not)
mtests = []
for fam in micro:
    for an in ('MI01', 'MI02'):
        q = micro_det[(micro_det.family == fam) & (micro_det.animal == an)]
        if len(q) < 2: continue
        tab = np.array([[r.detections, r.group_runs - r.detections] for r in q.itertuples])
        chi = stats.chi2_contingency(tab)
        mtests.append(dict(family=fam, animal=an, test='detection frequency across the three regions',
                           statistic=round(chi.statistic, 2), P=f'{chi.pvalue:.3g}'))
pept = pd.read_csv(os.path.join(D, 'report.pr_matrix.tsv'), sep='\t', usecols=['Protein.Group', 'Stripped.Sequence', 'Proteotypic'])
# narrowest eggNOG orthologous-group clade of each detected protein: ties the gut-wall proteome to the
# taxonomy of the faecal metagenome, where Bacteroides carries the largest GH18 share
def og_clade(q):
    parts = [x.split('|')[-1] for x in str(mann.OGs.get(q, '')).split(',') if '|' in x]
    return parts[-1] if parts else '-'
clade_rows = []
for fam, d in micro.items:
    for kind, lst in (('gut-restricted', d['groups']), ('background (also non-gut)', d['background'])):
        for pgid in lst:
            for pid in [x.split('|')[1] for x in pgid.split(';')]:
                clade_rows.append(dict(family=fam, status=kind, protein=pid, clade=og_clade(pid),
                                       KO=str(mann.KO.get(pid, '-')).replace('ko:', ''), CAZy=str(mann.CAZy.get(pid, '-'))))
clade_df = pd.DataFrame(clade_rows)
# can any detected protein be traced to an individual MAG? (contig ids are <sample>_<contig>_<gene>)
MAGD = os.path.join(BASE, 'results', 'reanalysis_2026-11', 'mags')
clade_df['MAG'] = '-'; clade_df['MAG_GTDB'] = '-'
if os.path.exists(os.path.join(MAGD, 'mag_cazymes_ge2.tsv')):
    cz_mag = pd.read_csv(os.path.join(MAGD, 'mag_cazymes_ge2.tsv'), sep='\t', header=None, names=['mag', 'gene', 'fams'])
    cz_mag['key'] = cz_mag.mag.str.replace(r'\.bin\.\d+$', '', regex=True) + '_' + cz_mag.gene
    k2m = dict(zip(cz_mag.key, cz_mag.mag))
    gtd = pd.read_csv(os.path.join(MAGD, 'gtdbtk.bac120.summary.tsv'), sep='\t').set_index('user_genome')
    clade_df['MAG'] = clade_df.protein.map(k2m).fillna('-')
    clade_df['MAG_GTDB'] = clade_df.MAG.map(lambda m: gtd.classification.get(m, '-') if m != '-' else '-')
    n_binned = int((clade_df.MAG != '-').sum)
    print(f'  {n_binned} of {len(clade_df)} detected GH18/GH20 proteins lie in a binned MAG:')
    for r in clade_df[clade_df.MAG != '-'].itertuples:
        print(f'    {r.protein} ({r.family}) -> {r.MAG}  {r.MAG_GTDB.split(";g__")[-1] if ";g__" in r.MAG_GTDB else r.MAG_GTDB}')
print('  detected microbial GH18/GH20 proteins by clade:')
print('   ', clade_df[clade_df.family == 'GH18'].clade.value_counts.to_dict)
with open(os.path.join(T, '73_diann_microbial_chitin_spatial.tsv'), 'w') as fh:
    fh.write('# Microbial (metagenome-database) chitin-active protein groups in the 36 gut-wall runs, species-specific DIA-NN search.\n'
             '# Negative control: groups also detected in MI02 liver/heart/muscle are listed as background and excluded from the regional counts.\n'
             '# detection_pct = detected group-runs / (groups x runs) in that region. Most groups rest on a single peptide: these are\n'
             '# presence calls at the detection limit, NOT abundance measurements. scripts/proteomics_diann_chitin.py\n')
    micro_det.to_csv(fh, sep='\t', index=False)
    fh.write('#\n# per-group evidence\n')
    ev = []
    for fam, d in micro.items:
        for kind, lst in (('gut-restricted prokaryotic', d['groups']), ('background (also non-gut)', d['background']),
                          ('eukaryotic (excluded)', d['eukaryotic'])):
            for pgid in lst:
                pp = pept[pept['Protein.Group'] == pgid]
                ev.append(dict(family=fam, status=kind, protein_group=pgid,
                               label=meta_rows.label.get(pgid, '-'), CAZy=meta_rows.CAZy.get(pgid, '-'),
                               peptides=pp['Stripped.Sequence'].nunique,
                               proteotypic_peptides=pp[pp.Proteotypic == 1]['Stripped.Sequence'].nunique,
                               n_runs_detected=int(M.loc[pgid].notna.sum)))
    ev = pd.DataFrame(ev); ev.to_csv(fh, sep='\t', index=False)
    fh.write('#\n# regional tests\n'); pd.DataFrame(mtests).to_csv(fh, sep='\t', index=False)
    fh.write('#\n# narrowest eggNOG clade of each detected protein (coarse: the clade of the seed ortholog group, not a genome assignment)\n')
    clade_df.to_csv(fh, sep='\t', index=False)
print('  tables/73_diann_microbial_chitin_spatial.tsv')
for t_ in mtests: print('   ', t_)
print(micro_det.to_string(index=False))

# Fig 6: MI01 spatial panels
seg1 = animal_matrix('MI01').join(chit.set_index('protein_group')[['category', 'label']])
seg2 = animal_matrix('MI02').join(chit.set_index('protein_group')[['category', 'label']]); cats_plot = ['GH18 locus (catalytic)', 'GH18 locus (inactive CLP)', 'GH20 hexosaminidase', 'GlcNAc metabolism']
fig = plt.figure(figsize=(16.5, 10)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], width_ratios=[1, 1, 1.5])
ax = fig.add_subplot(gs[0, 0])
for cat in cats_plot:
    sub = seg1[seg1.category == cat][SEG].dropna(thresh=4); cen = sub.sub(sub.mean(axis=1), axis=0)
    ax.errorbar(range(6), cen.mean, yerr=cen.sem, marker='o', color=COLC[cat], label=f'{cat} (n = {len(sub)})', capsize=3)
ax.axhline(0, color='#999', lw=0.8); ax.set_xticks(range(6)); ax.set_xticklabels(SEG, rotation=30, ha='right'); ax.set_ylabel('within-protein centred log$_2$ MaxLFQ'); ax.legend(fontsize=8, frameon=False); ax.set_title('A  Segment profile by category (MI01, mean ± s.e.m.)', loc='left', fontweight='bold', fontsize=10); ax.spines[['top', 'right']].set_visible(False)
ax = fig.add_subplot(gs[0, 1:])
hm = seg1[seg1.category.isin(cats_plot)][SEG].dropna(thresh=4); z = hm.sub(hm.mean(axis=1), axis=0).div(hm.std(axis=1).replace(0, np.nan), axis=0); z = z.loc[z.idxmax(axis=1).map({s: i for i, s in enumerate(SEG)}).sort_values.index]
im = ax.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2); ax.set_xticks(range(6)); ax.set_xticklabels(SEG, rotation=30, ha='right'); ax.set_yticks(range(len(z)))
ax.set_yticklabels([f"{seg1.label[i]} [{str(seg1.category[i]).split(' (')[0]}]"[:48] for i in z.index], fontsize=6); plt.colorbar(im, ax=ax, fraction=0.03, label='row z-score')
ax.set_title(f'B  All chitin-related host proteins (n = {len(z)}), row z-scores, ordered by peak segment', loc='left', fontweight='bold', fontsize=10)
ax = fig.add_subplot(gs[1, 0]); pk = pd.DataFrame(peak_rows); pk = pk[(pk.animal == 'MI01') & pk.category.isin(cats_plot)]
ct = pd.crosstab(pk.category, pk.peak_region).reindex(index=cats_plot, columns=['Stomach', 'Small intestine', 'Hindgut']).fillna(0)
ct.plot.bar(stacked=True, ax=ax, color=['#C9A227', '#B5541C', '#2E6F9E'], edgecolor='#333', linewidth=0.4); ax.set_xticklabels(['GH18\ncatalytic', 'GH18\nCLP', 'GH20', 'GlcNAc\nmetab.'], rotation=0, fontsize=7.5); ax.set_xlabel(''); ax.set_ylabel('proteins'); ax.legend(fontsize=8, frameon=False, title='peak region')
ax.set_title('C  Peak segment per protein (MI01)', loc='left', fontweight='bold', fontsize=10); ax.spines[['top', 'right']].set_visible(False)
ax = fig.add_subplot(gs[1, 1]); cat_sub = seg1[seg1.category == 'GH18 locus (catalytic)'][SEG]; cen = cat_sub.sub(cat_sub.mean(axis=1), axis=0)
ax.boxplot([cen[s].dropna for s in SEG], widths=0.6, patch_artist=True, boxprops=dict(facecolor='#F0D9C7'), medianprops=dict(color='#B5541C'))
for i, s in enumerate(SEG): ax.scatter(np.random.normal(i + 1, 0.05, len(cen[s].dropna)), cen[s].dropna, s=14, color='#B5541C', zorder=3)
ax.set_xticks(range(1, 7)); ax.set_xticklabels(SEG, rotation=30, ha='right'); ax.set_ylabel('centred log$_2$ MaxLFQ'); t1 = tests[(tests.animal == 'MI01') & (tests.category == 'GH18 locus (catalytic)')].iloc[0]
ax.set_title(f'D  Catalytic GH18 loci (n = {int(t1.n_proteins)})\nSI − hindgut median {t1.SI_minus_hindgut_median_log2:+.2f} log₂, Wilcoxon P = {t1.wilcoxon_P}', loc='left', fontweight='bold', fontsize=10); ax.spines[['top', 'right']].set_visible(False)
ax = fig.add_subplot(gs[1, 2])
# ---- E: the two systems side by side as one heatmap.
# M. impressa columns are measured here (per-region mean of the within-protein-centred log2 MaxLFQ of the
# catalytic GH18-locus groups, one column per animal); the pangolin columns are a qualitative reading of the
# published profile (Cheng et al. 2023), scored -1 / 0 / +1, and the tortoise microbial column is not measured.
def region_profile(seg_df):
    sub = seg_df[seg_df.category == 'GH18 locus (catalytic)'][SEG].dropna(thresh=4)
    cen = sub.sub(sub.mean(axis=1), axis=0)                       # centre each protein on its own mean
    out = {}
    for r, cols_ in REGION.items:
        block = cen[cols_]
        out[r] = float(block.mean(axis=1).mean) if len(cols_) > 1 else float(block.iloc[:, 0].mean)
    return out
mi01, mi02 = region_profile(seg1), region_profile(seg2)
PANGOLIN_HOST = {'Stomach': 1.0, 'Small intestine': 0.0, 'Hindgut': -1.0}     # host AMCase: gastric peak
PANGOLIN_MICR = {'Stomach': -1.0, 'Small intestine': 0.0, 'Hindgut': 1.0}     # microbial chitinase: hindgut peak
# measured microbial side: detection frequency of the gut-restricted microbial GH18 groups (tables/73)
micro_col = {an: {r.region: r.detection_pct for r in micro_det[(micro_det.family == 'GH18') & (micro_det.animal == an)].itertuples}
             for an in ('MI01', 'MI02')}
micro_P = {t['animal']: t['P'] for t in mtests if t['family'] == 'GH18'}
cols = [('host', PANGOLIN_HOST, 'q'), ('microbial', PANGOLIN_MICR, 'q'),
        ('host\nMI01', mi01, 'm'), ('host\nMI02', mi02, 'm'),
        ('microbial\nMI01', micro_col['MI01'], 'd'), ('microbial\nMI02', micro_col['MI02'], 'd')]
rows = list(REGION)
HM = np.full((len(rows), len(cols)), np.nan)
for j, (_, d, kind) in enumerate(cols):
    if d is None: continue
    v = np.array([d[r] for r in rows], dtype=float)
    HM[:, j] = v / (np.abs(v).max or 1)                          # scale each column to [-1, 1] so the shapes compare
im = ax.imshow(np.ma.masked_invalid(HM), cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
for j, (_, d, kind) in enumerate(cols):
    for i, r in enumerate(rows):
        txt = {'q': lambda v: {1.0: 'high', 0.0: 'mid', -1.0: 'low'}[v], 'm': lambda v: f'{v:+.2f}', 'd': lambda v: f'{v:.0f} %'}[kind](d[r])
        ax.text(j, i, txt, ha='center', va='center', fontsize=7.5,
                color='white' if abs(HM[i, j]) > 0.6 else '#222', fontweight='bold' if kind in 'md' else 'normal')
ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[0] for c in cols], fontsize=7)
def group_bracket(ax, x0, x1, y, label, h=0.045):
    """Square bracket spanning x0..x1 in AXES coordinates, opening downward, label underneath."""
    ax.plot([x0, x0, x1, x1], [y - h, y, y, y - h], color='#333', lw=1.1, clip_on=False, transform=ax.transAxes,
            solid_joinstyle='miter')
    ax.text((x0 + x1) / 2, y - h - 0.04, label, ha='center', va='top', fontsize=7.8, fontweight='bold',
            clip_on=False, transform=ax.transAxes)
f = lambda xd: (xd + 0.5) / len(cols)                              # data column -> axes fraction
group_bracket(ax, f(-0.45), f(1.45), -0.17, 'Pangolin\n(Cheng et al. 2023)')
group_bracket(ax, f(1.55), f(5.45), -0.17, '$\\mathit{M.\\,impressa}$\n(this study)')
ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=8.5)
for k in range(len(cols) - 1): ax.axvline(k + .5, color='white', lw=2)
for k in range(len(rows) - 1): ax.axhline(k + .5, color='white', lw=2)
ax.axvline(1.5, color='#333', lw=1.6)                              # separates the two studies
ax.axvline(3.5, color='#999', lw=1.0, ls=':')                      # host | microbial within this study
ax.set_title('E  Where chitinase sits along the gut: pangolin vs $\\mathit{M.\\,impressa}$', loc='left', fontweight='bold', fontsize=10)
fig.suptitle('Spatial distribution of chitinase-related host proteins along the GI tract — species-specific DIA-NN search, one animal (MI01); protein abundance, not activity', fontsize=11, fontweight='bold'); fig.tight_layout; save(fig, 'Fig8_spatial_chitinase_GI')
# Supp S14: protein groups per run by source; S15/S16: host/meta/contaminant fractions; S17: PCA
rm = run_meta.copy; rm['label'] = rm.animal + ' ' + rm.tissue + '_' + rm.piece.astype(str)
TIS_ORDER = SEG + ['Liver', 'Heart', 'Muscle']
rm = rm.sort_values(['animal', 'tissue', 'piece'],                       # replicates ascend with the segments
                    key=lambda c: c.map(TIS_ORDER.index) if c.name == 'tissue' else c)
x = np.arange(len(rm))
fig, ax = plt.subplots(figsize=(14, 4.5))
ax.bar(x, rm.host, color='#2E6F9E', label='host (Mimp_v3)')
ax.bar(x, rm.metagenome, bottom=rm.host, color='#B5541C', label='metagenome')
ax.bar(x, rm.contaminant, bottom=rm.host + rm.metagenome, color='#999', label='contaminant')
ax.set_xticks(x); ax.set_xticklabels(rm.label, rotation=90, fontsize=7); ax.set_ylabel('protein groups')
ax.legend(frameon=False, fontsize=8, loc='upper left', bbox_to_anchor=(1.0, 1.0))
ax.set_title('Supplementary Figure S14. Protein groups per run by database source (DIA-NN, 1 % FDR); precursors 32,741–82,361 per run',
             fontweight='bold', fontsize=10)
ax.text(0, -0.42, 'the metagenome and contaminant layers are 1.4 % and 0.3 % of identifications and are barely visible at this scale — '
        'Supplementary Figure S15 gives the totals, S16 the metagenome layer alone',
        transform=ax.transAxes, fontsize=7.5, color='#555', va='top')
ax.spines[['top', 'right']].set_visible(False)
save(fig, 'SuppFig_S14_proteomics_overview')
# S15: a pie of three numbers hid the 0.3 % wedge entirely; one 100 % stacked bar shows the split and keeps the counts legible
fig, ax = plt.subplots(figsize=(10, 2.3)); tot = rm[['host', 'metagenome', 'contaminant']].sum; pct = tot / tot.sum * 100
left = 0.0; small = 0
for k, col in (('host', '#2E6F9E'), ('metagenome', '#B5541C'), ('contaminant', '#999999')):
    ax.barh(0, pct[k], left=left, color=col, edgecolor='white', height=0.5)
    lab = f'{k}\n{int(tot[k]):,} ({pct[k]:.1f} %)'
    if pct[k] > 10: ax.text(left + pct[k] / 2, 0, lab, ha='center', va='center', color='white', fontweight='bold', fontsize=9)
    else:
        dx, ty, ha = ((-5, 0.70, 'right'), (2, 1.05, 'left'))[small]; small += 1   # fan the two thin wedges apart
        ax.annotate(lab, xy=(left + pct[k] / 2, 0.26), xytext=(left + pct[k] / 2 + dx, ty), ha=ha, fontsize=8,
                    arrowprops=dict(arrowstyle='-', color=col, lw=1))
    left += pct[k]
ax.set_xlim(0, 100); ax.set_ylim(-0.45, 1.5); ax.set_yticks([]); ax.set_xlabel('% of protein-group identifications')
for sp in ('top', 'right', 'left'): ax.spines[sp].set_visible(False)
ax.set_title(f'Supplementary Figure S15. Protein-group identifications by database source, all 36 runs '
             f'({int(tot.sum):,} identifications)', fontsize=9.5, fontweight='bold', loc='left')
fig.tight_layout; save(fig, 'SuppFig_S15_proteomics_remapping')
# S16: the metagenome layer alone, coloured by animal, with the non-gut tissues marked
fig, ax = plt.subplots(figsize=(14, 4.0))
acol = {'MI01': '#2E6F9E', 'MI02': '#7FB3D5'}
bars = ax.bar(x, rm.metagenome, color=[acol[a] for a in rm.animal], edgecolor='#333', linewidth=0.4)
for b, t in zip(bars, rm.tissue):
    if t not in SEG: b.set_hatch('///')
ax.set_xticks(x); ax.set_xticklabels(rm.label, rotation=90, fontsize=7); ax.set_ylabel('metagenome protein groups')
ax.set_xlim(-0.8, len(rm) - 0.2)
handles = [plt.Rectangle((0, 0), 1, 1, color=acol[a], label=a) for a in ('MI01', 'MI02')]
handles.append(plt.Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='#333', hatch='///', label='non-gut tissue'))
ax.legend(handles=handles, frameon=False, fontsize=8, loc='upper left', bbox_to_anchor=(1.0, 1.0))
ax.set_title('Supplementary Figure S16. Metagenome-origin protein groups per run', fontweight='bold', fontsize=10, loc='left')
ax.text(0, -0.52, 'The count is 50–112 in every run and essentially constant across gut segments and in non-gut tissue: '
        'a database-search background, not a spatial measurement.', transform=ax.transAxes, fontsize=8, color='#555', va='top')
ax.spines[['top', 'right']].set_visible(False)
save(fig, 'SuppFig_S16_per_sample_origin')
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
NONGUT = ['Liver', 'Heart', 'Muscle']
ORDER = SEG + NONGUT
# the gut series is ordinal, so it gets a sequential ramp; tab10 placed Ileum red beside Jejunum green
TIS_COLOR = {t: plt.cm.viridis(v) for t, v in zip(SEG, np.linspace(0.05, 0.92, len(SEG)))}
TIS_COLOR.update({t: c for t, c in zip(NONGUT, ['#B5541C', '#8C6BB1', '#777777'])})
TIS_MARKER = {**{t: 'o' for t in SEG}, **{t: '^' for t in NONGUT}}
seen = []
for k, (ax, animal) in enumerate(zip(axes, ['MI01', 'MI02'])):
    cs = [c for c in M.columns if run_meta.animal[c] == animal]; X = M[cs].dropna; X = X.sub(X.mean(axis=1), axis=0).T.values
    u, s_, vt = np.linalg.svd(X, full_matrices=False); pc = u[:, :2] * s_[:2]; var = s_ ** 2 / (s_ ** 2).sum
    tis = [run_meta.tissue[c] for c in cs]; uniq = sorted(set(tis), key=ORDER.index)
    for t in uniq:
        idx = [j for j, tt in enumerate(tis) if tt == t]
        ax.scatter(pc[idx, 0], pc[idx, 1], label=t, s=60, color=TIS_COLOR[t], marker=TIS_MARKER[t], edgecolor='#333')
        if t not in seen: seen.append(t)
    ax.set_xlabel(f'PC1 ({var[0] * 100:.0f} %)'); ax.set_ylabel(f'PC2 ({var[1] * 100:.0f} %)')
    ax.set_title(f'{"AB"[k]}  {animal}: {X.shape[1]:,} protein groups quantified in every run', loc='left', fontsize=10, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)
handles = [plt.Line2D([], [], ls='', marker=TIS_MARKER[t], color=TIS_COLOR[t], mec='#333', ms=8, label=t)
           for t in sorted(seen, key=ORDER.index)]
fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.995, 0.5), frameon=False, fontsize=8)
fig.suptitle('Supplementary Figure S17. Proteomics PCA by tissue (log₂ MaxLFQ, species-specific search)', fontweight='bold')
fig.tight_layout(rect=[0, 0, 0.93, 1]); save(fig, 'SuppFig_S17_proteomics_PCA')
# Supp S13: peptide coverage of the locus models
seqs = {}; cur = None
for l in open(os.path.join(A, 'gh18_locus_proteins.faa')):
    if l.startswith('>'): cur = l[1:].strip.replace('host|', ''); seqs[cur] = ''
    else: seqs[cur] += l.strip
# one shared residue axis: twelve independent axes made peptide positions incomparable between tracks
LMAX = max(len(seqs.get(g + '.t1', '')) for g in LOCUS) * 1.02
fig, axes = plt.subplots(len(LOCUS), 1, figsize=(12.5, 0.55 * len(LOCUS) + 1.4), sharex=True)
for ax, (g, (lab, kind)) in zip(axes, LOCUS.items):
    s = seqs.get(g + '.t1', ''); L = max(len(s), 1); ax.set_xlim(0, LMAX); ax.set_ylim(0, 1); ax.set_yticks([])
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.add_patch(plt.Rectangle((0, 0.3), L, 0.4, color='#E8E2F2' if kind == 'inactive CLP' else '#DCE9F5'))
    peps = pr2[pr2['Protein.Ids'].str.contains(g)]['Stripped.Sequence'].unique; n = 0
    for p in peps:
        i = s.find(p)
        if i >= 0: ax.add_patch(plt.Rectangle((i, 0.3), len(p), 0.4, color='#B5541C', alpha=0.7)); n += 1
    m = re.search(r'D.D.E', s)
    if m: ax.axvline(m.start, color='#2E6F9E', lw=1.2); ax.text(m.start + 3, 0.78, 'DxDxE', fontsize=6, color='#2E6F9E')
    ax.text(LMAX * 1.005, 0.5, f'{g} {lab}: {n} peptide{"" if n == 1 else "s"}, {len(s)} aa', fontsize=7, va='center')
axes[-1].set_xlabel('residue (one shared scale)')
key = [plt.Rectangle((0, 0), 1, 1, color='#DCE9F5', label='gene model at the locus (protein length)'),
       plt.Rectangle((0, 0), 1, 1, color='#E8E2F2', label='inactive CLP'),
       plt.Rectangle((0, 0), 1, 1, color='#B5541C', alpha=0.7, label='detected peptide'),
       plt.Line2D([], [], color='#2E6F9E', lw=1.2, label='DxDxE catalytic motif')]
fig.legend(handles=key, loc='lower center', ncol=4, frameon=False, fontsize=7.5, bbox_to_anchor=(0.45, -0.02))
fig.suptitle('Supplementary Figure S13. Peptide coverage of the twelve GH18-locus gene models '
             '(DIA-NN precursors mapped to Mimp_v3 proteins)', fontsize=10, fontweight='bold')
fig.tight_layout(rect=[0, 0.03, 1, 1]); save(fig, 'SuppFig_S13_gh18_locus_peptide_coverage')
print('done')
