#!/usr/bin/env python3
"""Fungal / eukaryotic profiling of the 32 non-host read sets against the extended Kraken2 database
($REFDIR/kraken2_pluspf_fungiext_2026-11: base libraries + 594 added fungal genomes incl. 50 Neocallimastigomycota
and the three diet mushrooms; built 2026-11-15, a job; profiling job scripts/reanalysis_2026-11/fungal_profiling.sbatch).
Input (copied from the server): results/reanalysis_2026-11/fungal_profiling/{kraken/*.report,*.log, bracken/*.S.bracken}
Output: tables/67_fungal_profiling_per_sample.tsv  (per sample: reads, classified %, Fungi / Neocallimastigomycota / diet-mushroom
        clade reads and reads per million, top fungal genus), tables/68_fungal_taxa_prevalence.tsv (fungal genera + AGF genera:
        prevalence, per-cohort median ppm), figures SuppFig_S1_fungal_profiling.{png,pdf}.
Usage: fungal_profiling_tables.py [indir]  (default results/reanalysis_2026-11/fungal_profiling)"""
import sys, re, glob, pathlib, collections, pandas as pd, numpy as np, matplotlib


def _registry(root):
    """The full registry is not published (it carries dissection and specimen notes); the public
    subset holds the animal, date and cohort these analyses need."""
    full = root / 'tables/deposition/sample_registry_all.tsv'
    return full if full.exists else root / 'tables/deposition/sample_registry_public.tsv'

matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = pathlib.Path(__file__).resolve.parents[1]
IN = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else R/'results/reanalysis_2026-11/fungal_profiling'
reg = pd.read_csv(_registry(R), sep='\t', comment='#')
reg = reg[reg.sample_type.str.contains('faec', na=False)].set_index('sample_id')
COH = {'DS01': 'P. eryngii', 'DS02': 'H. tessellatus', 'DS03': 'P. ostreatus'}
DIET = {'DS01': 'Pleurotus eryngii', 'DS02': 'Hypsizygus marmoreus', 'DS03': 'Pleurotus ostreatus'}
AGF = 'Neocallimastigomycetes'  # NCBI taxonomy node (class, taxid 451455) carrying the anaerobic gut fungi in the Kraken2 report
CLADES = ['Fungi', AGF, 'Pleurotus eryngii', 'Pleurotus ostreatus', 'Hypsizygus marmoreus', 'Metazoa', 'Viridiplantae']
# The database holds archaea, bacteria, viral, human, protozoa and fungi: there is no plant library and
# the only metazoan is human. Viridiplantae is therefore structurally zero and Metazoa means human reads;
# both are reported only to make the database's blind spots explicit, never as evidence of absence.
def read_report(f):
    """Kraken2 report -> dict name->(clade_reads, rank, depth) and total; names de-indented."""
    d = {}; rows = []
    for l in open(f):
        p = l.rstrip('\n').split('\t')
        if len(p) < 6: continue
        name = p[5]; depth = (len(name) - len(name.lstrip(' '))) // 2; name = name.strip
        rows.append((name, int(p[1]), int(p[2]), p[3], int(p[4]), depth))
        d.setdefault(name, (int(p[1]), p[3], depth, int(p[2])))
    return d, rows
def fungal_genera(rows):
    """genus-level clade reads for genera under Fungi (uses the indentation tree)."""
    out = {}; in_f = False; fdepth = None
    for name, clade, direct, rank, taxid, depth in rows:
        if name == 'Fungi': in_f, fdepth = True, depth; continue
        if in_f and depth <= fdepth: in_f = False
        if in_f and rank == 'G': out[name] = clade
    return out
recs = []; gen = {}
for f in sorted(glob.glob(str(IN/'kraken/*.report'))):
    s = pathlib.Path(f).stem; d, rows = read_report(f)
    log = open(IN/f'kraken/{s}.log').read
    tot = int(re.search(r'(\d+) sequences \(', log).group(1)); cls = int(re.search(r'(\d+) sequences classified', log).group(1))
    coh = reg.loc[s, 'dataset'] if s in reg.index else {'ITF': 'DS01', 'HM': 'DS02', 'PM': 'DS03'}[re.match(r'[A-Z]+', s).group]
    r = dict(Cohort=coh, Diet=COH[coh], Sample=s, animal=reg.loc[s, 'animal_id (unified)'] if s in reg.index else '', collection_date=reg.loc[s, 'collection_date'] if s in reg.index else '',
             read_pairs=tot, classified_pct=round(100 * cls / tot, 2))
    for c in CLADES:
        n = d.get(c, (0,))[0]; r[c.replace(' ', '_') + '_reads'] = n; r[c.replace(' ', '_') + '_ppm'] = round(1e6 * n / tot, 3)
    # reads left at the Fungi node itself (LCA = Fungi: k-mers shared across many fungal genomes, e.g. rDNA) versus reads resolved to a fungal genus
    r['Fungi_unresolved_at_node_reads'] = d.get('Fungi', (0, '', 0, 0))[3]
    r['fed_species_reads'] = d.get(DIET[coh], (0,))[0]; r['fed_species_ppm'] = round(1e6 * r['fed_species_reads'] / tot, 3)
    g = fungal_genera(rows); gen[s] = g
    top = sorted(g.items, key=lambda kv: -kv[1])[:3]
    r['Fungi_genus_resolved_reads'] = sum(g.values); r['Fungi_genus_resolved_ppm'] = round(1e6 * sum(g.values) / tot, 3)
    r['fungal_genera_detected'] = sum(1 for v in g.values if v >= 10)
    r['top_fungal_genera'] = '; '.join(f'{k} {v}' for k, v in top)
    recs.append(r)
df = pd.DataFrame(recs).sort_values(['Cohort', 'Sample'])
# AGF genera (Neocallimastigomycota) present in the reports
agf_names = set
for f in glob.glob(str(IN/'kraken/*.report')):
    _, rows = read_report(f); in_n = False; nd = None
    for name, clade, direct, rank, taxid, depth in rows:
        if name == AGF: in_n, nd = True, depth; continue
        if in_n and depth <= nd: in_n = False
        if in_n and rank == 'G': agf_names.add(name)
G = pd.DataFrame(gen).fillna(0).astype(int)  # genera x samples
ppm = G.div(df.set_index('Sample').read_pairs.reindex(G.columns), axis=1) * 1e6
prev = pd.DataFrame({'genus': G.index, 'AGF': G.index.isin(agf_names), 'samples_ge10_reads': (G >= 10).sum(axis=1).values, 'total_reads': G.sum(axis=1).values})
for c in ('DS01', 'DS02', 'DS03'):
    cols = df[df.Cohort == c].Sample; prev[f'median_ppm_{c}'] = ppm[cols].median(axis=1).round(3).values; prev[f'max_ppm_{c}'] = ppm[cols].max(axis=1).round(3).values
prev = prev.sort_values('total_reads', ascending=False)
hdr = (f"# Kraken2 (confidence 0.2, min-hit-groups 2) + Bracken against the extended database (base libraries + 594 fungal genomes incl. 50 Neocallimastigomycota "
       f"and the three diet mushrooms); {len(df)} samples; reads = non-host read pairs. Fungi clade reads: median {df.Fungi_ppm.median:.1f} ppm "
       f"(range {df.Fungi_ppm.min:.1f}-{df.Fungi_ppm.max:.1f}), of which a median {100*df.Fungi_unresolved_at_node_reads.sum/df.Fungi_reads.sum:.0f}% sit unresolved at the Fungi node and "
       f"{df.Fungi_genus_resolved_ppm.median:.1f} ppm (median) reach a genus; Neocallimastigomycetes (anaerobic gut fungi): {int((df.Neocallimastigomycetes_reads>0).sum)}/{len(df)} samples with >=1 read, "
       f"median {df.Neocallimastigomycetes_ppm.median:.2f} ppm, max {df.Neocallimastigomycetes_ppm.max:.2f}; fed mushroom species: {int((df.fed_species_reads>0).sum)}/{len(df)} samples, "
       f"median {df.fed_species_ppm.median:.2f} ppm.\n"
       f"# Database scope: the library set is archaea, bacteria, viral, human, protozoa and fungi. There is NO plant library, so Viridiplantae_reads is "
       f"structurally 0 in every sample and is not evidence that no plant DNA is present; the only metazoan genome is human, so Metazoa_reads means human reads. "
       f"scripts/fungal_profiling_tables.py\n")
with open(R/'tables/67_fungal_profiling_per_sample.tsv', 'w') as fh: fh.write(hdr); df.to_csv(fh, sep='\t', index=False)
with open(R/'tables/68_fungal_taxa_prevalence.tsv', 'w') as fh:
    fh.write(f"# Genus-level clade reads under Fungi in the Kraken2 reports (same run as tables/67); AGF = Neocallimastigomycota genus; {len(agf_names)} AGF genera in the database tree. scripts/fungal_profiling_tables.py\n"); prev.to_csv(fh, sep='\t', index=False)
print(hdr); print(prev.head(15).to_string(index=False))
# ---- Supplementary Figure S1
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={'width_ratios': [1, 1, 1.6]})
col = {'DS01': '#2E86AB', 'DS02': '#A23B72', 'DS03': '#F18F01'}   # same cohort palette as generate_unified_figs.py
for i, (k, lab) in enumerate([('Fungi_genus_resolved_ppm', 'Fungi, genus-resolved — reads per million'), ('Neocallimastigomycetes_ppm', 'Neocallimastigomycetes (AGF) — reads per million')]):
    for j, c in enumerate(('DS01', 'DS02', 'DS03')):
        v = df[df.Cohort == c][k]; x = j + np.random.default_rng(1).uniform(-0.18, 0.18, len(v))
        ax[i].scatter(x, v, color=col[c], s=22, zorder=3, edgecolor='k', lw=0.3); ax[i].hlines(v.median, j - 0.3, j + 0.3, color='k', lw=1.2)
    ax[i].set_xticks(range(3)); ax[i].set_xticklabels([COH[c] for c in ('DS01', 'DS02', 'DS03')], style='italic', fontsize=8); ax[i].set_ylabel(lab, fontsize=8)
    ax[i].set_title('AB'[i] + '  ' + lab.split(' — ')[0], loc='left', fontsize=9); ax[i].set_yscale('symlog', linthresh=1); ax[i].set_ylim(bottom=0)
    ax[i].set_ylabel(lab + '\n(symlog: linear below 1 ppm, so true zeros are shown)', fontsize=8)
    ax[i].plot([], [], color='k', lw=1.2, label='cohort median'); ax[i].legend(frameon=False, fontsize=7, loc='upper right')
    for s in ('top', 'right'): ax[i].spines[s].set_visible(False)
topg = prev[~prev.AGF].head(8).genus.tolist; agf = prev[prev.AGF].head(4).genus.tolist
M = ppm.reindex(topg + agf).fillna(0)[df.Sample]
im = ax[2].imshow(np.log10(M.values + 0.1), aspect='auto', cmap='YlOrBr')
ax[2].set_yticks(range(len(M.index))); ax[2].set_yticklabels([g + (' (AGF)' if g in agf else '') for g in M.index], fontsize=7, style='italic')
ax[2].set_xticks(range(len(M.columns))); ax[2].set_xticklabels(M.columns, rotation=90, fontsize=6)
for j, s in enumerate(M.columns): ax[2].get_xticklabels[j].set_color(col[df.set_index('Sample').loc[s, 'Cohort']])
plt.colorbar(im, ax=ax[2], fraction=0.03, pad=0.01).set_label('log10 (reads per million + 0.1)', fontsize=7)
ax[2].set_title('C  Most abundant fungal genera and the AGF genera, per sample', loc='left', fontsize=9)
# the sample labels are coloured by cohort — say so, rather than leaving the reader to infer it
ax[2].legend(handles=[plt.Line2D([], [], marker='s', ls='', color=col[c], label=COH[c]) for c in ('DS01', 'DS02', 'DS03')],
             title='sample label colour', loc='upper left', bbox_to_anchor=(1.12, 1.0), frameon=False, fontsize=7, title_fontsize=7)
fig.suptitle('Supplementary Figure S1. Fungal and anaerobic-gut-fungal (Neocallimastigomycetes) signal in the 32 faecal metagenomes (extended Kraken2 database)', fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.95))
for _d in ('png', 'pdf'):  # a fresh clone has no figures/ yet
    (R/'figures'/_d).mkdir(parents=True, exist_ok=True)
for ext in ('png', 'pdf'): fig.savefig(R/f'figures/{ext}/SuppFig_S1_fungal_profiling.{ext}', dpi=300)
print('written SuppFig_S1_fungal_profiling')
