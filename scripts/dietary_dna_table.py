#!/usr/bin/env python3
"""Table 53 — dietary mushroom DNA in the 32 faecal metagenomes.
Input: results/reanalysis_2026-11/diet_mapping/<sample>.tsv (server job diet_mapping.sbatch: non-host pairs ->
Bowtie2 --very-sensitive against a combined index of the three mushroom genomes; properly paired, MAPQ >= 10).
Output: tables/53_dietary_dna_mapping.tsv — pairs and pairs per million non-host pairs per genome, split into the
species that cohort was fed vs the other two (the two Pleurotus share conserved regions, so 'other' is non-zero
for the Pleurotus cohorts but ~0 for the H. marmoreus cohort)."""
import glob, os, pandas as pd
from scipy import stats
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fs = sorted(glob.glob(os.path.join(BASE, 'results/reanalysis_2026-11/diet_mapping/*.tsv')))
t = pd.concat([pd.read_csv(f, sep='\t') for f in fs])
coh = lambda s: 'DS01' if s.startswith('ITF') else 'DS02' if s.startswith('HM') else 'DS03'
diet = {'DS01': 'Pleurotus_eryngii_183', 'DS02': 'Hypsizygus_marmoreus_HM62', 'DS03': 'Pleurotus_ostreatus_PC9'}
short = {'Pleurotus_eryngii_183': 'P_eryngii', 'Hypsizygus_marmoreus_HM62': 'H_marmoreus', 'Pleurotus_ostreatus_PC9': 'P_ostreatus'}
rows = []
for s in sorted({os.path.basename(f)[:-4] for f in fs}):
    sub = t[t.Sample == s]; tot = int(sub.nonhost_pairs.iloc[0])
    r = dict(Cohort=coh(s), Sample=s, fed_species=short[diet[coh(s)]], nonhost_pairs=tot, concordant_pairs_any_mapq=int(sub.concordant_pairs_total.iloc[0]))
    for g, sh in short.items:
        v = int(sub[sub.Genome == g].pairs_q10.sum); r[f'{sh}_pairs_q10'] = v; r[f'{sh}_per_M'] = round(v / tot * 1e6, 3)
    fed = diet[coh(s)]
    r['fed_species_per_M'] = round(sub[sub.Genome == fed].pairs_q10.sum / tot * 1e6, 3)
    r['other_species_per_M'] = round(sub[sub.Genome != fed].pairs_q10.sum / tot * 1e6, 3)
    rows.append(r)
d = pd.DataFrame(rows)
out = os.path.join(BASE, 'tables', '53_dietary_dna_mapping.tsv')
with open(out, 'w') as fh:
    d.to_csv(fh, sep='\t', index=False)
    w = stats.wilcoxon(d.fed_species_per_M, d.other_species_per_M, alternative='greater')
    fh.write(f'\n# fed species > other two species, paired Wilcoxon one-sided, n = 32: W = {w.statistic:.0f}, P = {w.pvalue:.2e}\n')
    fh.write(f'# fed species detected (>= 1 pair) in {(d.fed_species_per_M > 0).sum}/32 samples; median per cohort (per M non-host pairs): '
             + '; '.join(f'{c} {v:.1f}' for c, v in d.groupby("Cohort").fed_species_per_M.median.items) + '\n')
    fh.write(f'# maximum dietary-DNA fraction of non-host reads: {d.fed_species_per_M.max / 1e6:.2e} (HM.1.5)\n')
print(d.groupby('Cohort')[['fed_species_per_M', 'other_species_per_M']].median.round(2)); print('P =', f'{w.pvalue:.1e}')
