#!/usr/bin/env python3
"""Species-specific DIA-NN re-search (2026-11-10): summarise the 36-run search of the Orbitrap Astral data against the Mimp_v3 proteome
(29,927) + metagenome CAZy/ko00520 proteins (65,450) + cRAP (116). Inputs: results/reanalysis_2026-11/diann/report.{stats,pg_matrix}.tsv
(DIA-NN 1.8.1, 1 % precursor and protein-group FDR, MBR). Outputs:
  tables/50_diann_mimpv3_run_summary.tsv   per run: animal, tissue, precursors, protein groups by source (host / metagenome / contaminant)
  tables/51_diann_gh18_locus_quant.tsv     GH18 locus protein groups (MIMP3_g08468-g08479) x run, log2 MaxLFQ
  tables/52_diann_gh18_spatial.tsv         per animal: within-locus centred 6-segment profiles, SI-vs-hindgut paired test, gut vs non-gut"""
import os, re, numpy as np, pandas as pd
from scipy import stats
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); D = os.path.join(BASE, 'results', 'reanalysis_2026-11', 'diann'); T = os.path.join(BASE, 'tables')
RUN = {**{r: ('MI01', 'DS04', t, k) for r, t, k in [('25NP14550003', 'Stomach', 1), ('25NP14550002', 'Stomach', 2), ('25NP14550001', 'Stomach', 3), ('25NP14550015', 'Duodenum', 1), ('25NP14550014', 'Duodenum', 2), ('25NP14550013', 'Duodenum', 3), ('25NP14550027', 'Jejunum', 1), ('25NP14550026', 'Jejunum', 2), ('25NP14550025', 'Jejunum', 3), ('25NP14550039', 'Ileum', 1), ('25NP14550038', 'Ileum', 2), ('25NP14550037', 'Ileum', 3), ('25NP14550051', 'Colon', 1), ('25NP14550050', 'Colon', 2), ('25NP14550049', 'Colon', 3), ('25NP14550063', 'Rectum', 1), ('25NP14550062', 'Rectum', 2), ('25NP14550061', 'Rectum', 3)]},
       **{r: ('MI02', 'DS05', t, k) for r, t, k in [('25NP24100001', 'Liver', 1), ('25NP24100002', 'Liver', 2), ('25NP24100003', 'Heart', 1), ('25NP24100004', 'Heart', 2), ('25NP24100005', 'Stomach', 1), ('25NP24100006', 'Stomach', 2), ('25NP24100007', 'Duodenum', 1), ('25NP24100008', 'Duodenum', 2), ('25NP24100009', 'Jejunum', 1), ('25NP24100010', 'Jejunum', 2), ('25NP24100011', 'Ileum', 1), ('25NP24100012', 'Ileum', 2), ('25NP24100013', 'Colon', 1), ('25NP24100014', 'Colon', 2), ('25NP24100015', 'Rectum', 1), ('25NP24100016', 'Rectum', 2), ('25NP24100017', 'Muscle', 1), ('25NP24100018', 'Muscle', 2)]}}
LOCUS = {'MIMP3_g08468': 'TAF8 (v1 088402 N-terminal third; not GH18)', 'MIMP3_g08469': '088402a (GH18+CBM14)', 'MIMP3_g08470': '088402b (GH18)', 'MIMP3_g08471': 'GH18 fragment', 'MIMP3_g08472': 'L7 (CHIA array)', 'MIMP3_g08473': 'L6 (CHIA array)', 'MIMP3_g08474': 'L5 (CHIA array)', 'MIMP3_g08475': 'L4 (CHIA array)', 'MIMP3_g08476': 'L3 (CHIA-related)', 'MIMP3_g08477': '088393 fragment (g08477)', 'MIMP3_g08478': 'L2 (CHI3L-clade inactive CLP, 088393)', 'MIMP3_g08479': 'L1 (CHIT1-clade, 088392)'}
SEG = ['Stomach', 'Duodenum', 'Jejunum', 'Ileum', 'Colon', 'Rectum']; SI, HG = ['Duodenum', 'Jejunum', 'Ileum'], ['Colon', 'Rectum']
runof = lambda c: re.search(r'(25NP\d+)', c).group(1)
st = pd.read_csv(os.path.join(D, 'report.stats.tsv'), sep='\t'); st['run'] = st['File.Name'].map(runof)
pg = pd.read_csv(os.path.join(D, 'report.pg_matrix.tsv'), sep='\t'); cols = [c for c in pg.columns if c.endswith('.mzML')]
pg['source'] = pg['Protein.Group'].str.split('|').str[0]
rows = []
for c in cols:
    r = runof(c); a, b, t, k = RUN[r]; det = pg[c].notna
    rows.append(dict(run=r, animal=a, batch=b, tissue=t, piece=k, precursors=int(st.loc[st.run == r, 'Precursors.Identified'].iloc[0]),
                     protein_groups=int(det.sum), host=int((det & (pg.source == 'host')).sum), metagenome=int((det & (pg.source == 'meta')).sum), contaminant=int((det & (pg.source == 'crap')).sum)))
summ = pd.DataFrame(rows).sort_values(['batch', 'run'])
with open(os.path.join(T, '50_diann_mimpv3_run_summary.tsv'), 'w') as fh:
    fh.write(f"# DIA-NN 1.8.1 library-free, 1 % precursor + protein-group FDR, MBR; database = Mimp_v3 proteome (29,927) + metagenome CAZy/ko00520 proteins (65,450) + cRAP (116). Protein groups in total: {len(pg)} (host {int((pg.source=='host').sum)}, metagenome {int((pg.source=='meta').sum)}, contaminant {int((pg.source=='crap').sum)}).\n")
    summ.to_csv(fh, sep='\t', index=False)
# ---- GH18 loci
gh = pg[pg['Protein.Group'].str.contains(r'MIMP3_g0846[89]|MIMP3_g0847\d')].copy
M = gh.set_index('Protein.Group')[cols]; M.columns = [runof(c) for c in cols]; M = np.log2(M.astype(float))
M.index = [';'.join(sorted(set(re.findall(r'MIMP3_g\d+', i)))) for i in M.index]
def label(idx): return ' + '.join(LOCUS.get(g, g) for g in idx.split(';'))
out = M.copy; out.insert(0, 'locus', [label(i) for i in M.index])
out.columns = ['locus'] + [f'{r}|{RUN[r][0]}|{RUN[r][2]}_{RUN[r][3]}' for r in M.columns]
with open(os.path.join(T, '51_diann_gh18_locus_quant.tsv'), 'w') as fh:
    fh.write('# log2 MaxLFQ protein-group quantities of the chr4 GH18-locus gene models (MIMP3_g08468-g08479) per run; columns run|animal|tissue_piece. Shared peptides merge the near-identical tandem CHIA copies g08473/g08474/g08475 into one group.\n'); out.to_csv(fh, sep='\t')
# ---- spatial tests per animal
res = []
for animal in ('MI01', 'MI02'):
    runs = {r: v for r, v in RUN.items if v[0] == animal}
    seg = pd.DataFrame({s: M[[r for r, v in runs.items if v[2] == s]].mean(axis=1) for s in SEG}).dropna
    prof = seg.sub(seg.mean(axis=1), axis=0); si = prof[SI].mean(axis=1); hg = prof[HG].mean(axis=1)
    for g in prof.index:
        res.append(dict(animal=animal, group=g, locus=label(g), **{f'{s}_centred_log2': round(prof.loc[g, s], 3) for s in SEG}, SI_minus_hindgut_log2=round(si[g] - hg[g], 3)))
    w = stats.wilcoxon(si, hg) if len(prof) >= 5 else None; fr = stats.friedmanchisquare(*[prof[s].values for s in SEG]) if len(prof) >= 5 else None
    res.append(dict(animal=animal, group='SUMMARY', locus=f'{len(prof)} GH18 protein groups with complete 6-segment profiles',
                    **{f'{s}_centred_log2': round(prof[s].mean, 3) for s in SEG}, SI_minus_hindgut_log2=round(float(np.median(si - hg)), 3),
                    note=f"SI>hindgut in {int((si > hg).sum)}/{len(prof)}; Wilcoxon P={w.pvalue:.3g}; Friedman P={fr.pvalue:.3g}" if w else 'too few groups'))
    if animal == 'MI02':
        gut = seg.mean(axis=1)
        for t in ('Liver', 'Heart', 'Muscle'):
            v = M[[r for r, vv in runs.items if vv[2] == t]].mean(axis=1)
            for g in M.index: res.append(dict(animal=animal, group=g, locus=label(g), note=f'{t}: log2 {v[g]:.2f} vs gut mean {gut.get(g, np.nan):.2f}' if pd.notna(v[g]) else f'{t}: not detected (gut mean {gut.get(g, np.nan):.2f})'))
sp = pd.DataFrame(res)
with open(os.path.join(T, '52_diann_gh18_spatial.tsv'), 'w') as fh:
    fh.write('# Within-locus centred 6-segment profiles (log2 MaxLFQ, mean of tissue pieces) of the GH18-locus protein groups from the species-specific DIA-NN search; paired SI (duodenum, jejunum, ileum) vs hindgut (colon, rectum) per animal; MI02 gut vs non-gut tissues.\n'); sp.to_csv(fh, sep='\t', index=False)
print(summ.groupby('batch')[['precursors', 'protein_groups', 'host', 'metagenome', 'contaminant']].mean.round(0).to_string)
print(sp[sp.group == 'SUMMARY'][['animal', 'locus', 'SI_minus_hindgut_log2', 'note']].to_string(index=False))
