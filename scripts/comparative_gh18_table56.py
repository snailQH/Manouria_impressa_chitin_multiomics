#!/usr/bin/env python3
"""Table 56 — the GH18 locus across nine Testudines genomes, from the per-genome loci files written by
scripts/reanalysis_2026-11/comparative_gh18_parse.py (results/reanalysis_2026-11/comparative_gh18/<assembly>.loci.tsv +
comparative_gh18_summary.tsv). Locus components (best-matching Mimp_v3 model): 088392 CHIT1-clade (g08479); 088393 CLP (g08478) and its
243-aa GH18 fragment model (g08477, 'CLP fragment'); L3 CHIA-related (g08476); CHIA array L4-L7 (g08472-75); 088402a GH18+CBM14 (g08469) and
088402b GH18-only (g08470) — the two chitinase genes inside the fused 1,160-aa model 088402, whose N-terminal third is TAF8 (g08468; Pfam
Bromo_TP + TAF8_C) and is therefore not counted. pseudo = frameshift or in-frame stop in the best miniprot alignment."""
import csv, pathlib
R = pathlib.Path(__file__).resolve.parents[1]; D = R/'results/reanalysis_2026-11/comparative_gh18'
META = {'Mimp_v3': ('Manouria impressa', 'Testudinidae', 'mycophagous'), 'GCF_003597395.2': ('Chelonoidis abingdonii', 'Testudinidae', 'herbivorous'),
        'GCF_007399415.2': ('Gopherus evgoodei', 'Testudinidae', 'herbivorous'), 'GCF_009764565.3': ('Dermochelys coriacea', 'Dermochelyidae', 'gelatinous zooplankton'),
        'GCF_011386835.1': ('Chrysemys picta bellii', 'Emydidae', 'omnivorous'), 'GCF_013100865.1': ('Trachemys scripta elegans', 'Emydidae', 'omnivorous'),
        'GCF_015237465.2': ('Chelonia mydas', 'Cheloniidae', 'herbivorous (adult)'), 'GCF_016161935.1': ('Mauremys reevesii', 'Geoemydidae', 'omnivorous'),
        'GCF_052426305.1': ('Pelodiscus sinensis', 'Trionychidae', 'carnivorous')}
SHORT = {'088392 CHIT1-clade': '088392', '088393 CLP': '088393', 'L3 CHIA-related': 'L3', '088402a GH18+CBM14': '088402a', '088402b GH18': '088402b'}
summ = {r['genome']: r for r in csv.DictReader(open(D/'comparative_gh18_summary.tsv'), delimiter='\t')}
hdr = ['assembly', 'species', 'family', 'diet', 'interval_chrom', 'interval_kb', 'flanks_found', 'GH18_genes_annotated_genome', 'GH18_genes_annotated_interval',
       'loci_miniprot_interval', 'CHIA_array_copies', 'CHIA_array_intact', 'CHIA_array_pseudo', 'L3_CHIA_related', 'L3_intact', 'CLP_088393', 'CLP_fragment_g08477',
       'CHIT1_088392', 'GH18_088402a_CBM14', 'GH18_088402b', 'complete_GH18_genes_intact', 'order']
rows = []
for asm, (sp, fam, diet) in META.items:
    loci = [l for l in csv.DictReader((x for x in open(D/f'{asm}.loci.tsv') if not x.startswith('#')), delimiter='\t') if l['in_interval'] == 'True']
    loci.sort(key=lambda l: int(l['start'])); s = summ[asm]
    def lab(l):
        if l['best_query'] == 'MIMP3_g08477': return '393frag'
        c = l['clade']; return SHORT.get(c, c.split[-1])
    def n(pred): return sum(1 for l in loci if pred(l))
    arr = lambda l: l['clade'].startswith('CHIA-array'); ps = lambda l: l['pseudogene'] == 'True'
    rows.append([asm, sp, fam, diet, s['interval_chrom'], s['interval_kb'], f"{s['flank_up_found']}+{s['flank_dn_found']}", s['GH18_genes_hmm'], s['GH18_genes_hmm_in_interval'],
                 len(loci), n(arr), n(lambda l: arr(l) and not ps(l)), n(lambda l: arr(l) and ps(l)),
                 n(lambda l: l['clade'] == 'L3 CHIA-related'), n(lambda l: l['clade'] == 'L3 CHIA-related' and not ps(l)),
                 n(lambda l: l['best_query'] == 'MIMP3_g08478'), n(lambda l: l['best_query'] == 'MIMP3_g08477'),
                 n(lambda l: l['clade'] == '088392 CHIT1-clade'), n(lambda l: l['clade'] == '088402a GH18+CBM14'), n(lambda l: l['clade'] == '088402b GH18'),
                 n(lambda l: l['best_query'] != 'MIMP3_g08477' and not ps(l)),
                 ' - '.join(lab(l) + ('ψ' if ps(l) else '') for l in loci)])
with open(R/'tables/56_comparative_gh18_locus_testudines.tsv', 'w') as fh:
    fh.write('# GH18 locus across Testudines: hmmsearch (dbCAN GH18 model, E<=1e-10, >=35 % model coverage) on RefSeq proteomes; miniprot of the Mimp_v3 locus models\n'
             '# (identity >= 0.75, <=40 kb) against each genome; orthologous interval = span of the 3 upstream + 3 downstream flanking genes (MIMP3_g08465-67, g08480-82).\n'
             '# Loci labelled by best-matching M. impressa model; pseudo = frameshift or in-frame stop in the best alignment. g08468 (TAF8, Pfam Bromo_TP+TAF8_C) and the 191-aa fragment g08471 are\n'
             '# not counted; g08477 (243-aa GH18 fragment next to the CLP g08478) is listed separately. complete_GH18_genes_intact = loci other than the g08477 fragment without frameshift/stop.\n'
             '# Scripts: scripts/reanalysis_2026-11/comparative_gh18.sbatch, comparative_gh18_parse.py, scripts/comparative_gh18_table56.py (2026-11-15; Pfam 39.0 domain check)\n')
    fh.write('\t'.join(hdr) + '\n'); [fh.write('\t'.join(map(str, r)) + '\n') for r in rows]
for r in rows: print(r[1], *r[9:])
