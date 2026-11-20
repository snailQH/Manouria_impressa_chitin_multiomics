#!/usr/bin/env python3
"""Summarise the comparative GH18 locus analysis (see comparative_gh18.sbatch).
Usage: comparative_gh18_parse.py genomes.tsv outdir > summary.tsv ; also writes outdir/<name>.loci.tsv
Per genome: GH18 proteins by hmmsearch (E<=1e-10, >=35% of the HMM covered), their genes and positions (from the GFF);
miniprot placements of the 12 M. impressa locus models (identity >= 0.75, alignment span <= 40 kb — single M. impressa models span <= 13 kb)
clustered into loci by reciprocal overlap (>= 50 % of the shorter alignment); a locus is flagged pseudogene when its best alignment
carries a frameshift or in-frame stop; the orthologous interval =
span of the flanking-gene placements (g08465-67 upstream, g08480-82 downstream) on the chromosome carrying most GH18 loci."""
import sys, re, gzip, collections
genomes, outdir = sys.argv[1], sys.argv[2]
# Pfam 39.0 (2026-11-15): g08468 = TAF8 (Bromo_TP + TAF8_C), not GH18 — the N-terminal segment of the fused v1 model 088402; g08469 = GH18 + CBM14
# (v1 088402 domain 1); g08470 = GH18 only (domain 2); g08471 = 191-aa GH18 fragment. Non-GH18 / fragment models are kept as syntenic markers only.
LOCUS = {'MIMP3_g08468':'non-GH18','MIMP3_g08469':'088402a GH18+CBM14','MIMP3_g08470':'088402b GH18','MIMP3_g08471':'non-GH18',
         'MIMP3_g08472':'CHIA-array L7','MIMP3_g08473':'CHIA-array L6','MIMP3_g08474':'CHIA-array L5','MIMP3_g08475':'CHIA-array L4',
         'MIMP3_g08476':'L3 CHIA-related','MIMP3_g08477':'088393 CLP','MIMP3_g08478':'088393 CLP','MIMP3_g08479':'088392 CHIT1-clade'}
FLANK_UP = {'MIMP3_g08465','MIMP3_g08466','MIMP3_g08467'}; FLANK_DN = {'MIMP3_g08480','MIMP3_g08481','MIMP3_g08482'}
def opn(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)
def gff_proteins(gff):
    """protein_id -> (chrom, start, end, gene_name) from CDS lines (RefSeq) or mRNA/CDS (Mimp_v3)."""
    prot = {}; mrna_gene = {}
    for l in opn(gff):
        if l.startswith('#'): continue
        p = l.rstrip('\n').split('\t')
        if len(p) < 9: continue
        a = dict(kv.split('=',1) for kv in p[8].split(';') if '=' in kv)
        if p[2] == 'mRNA': mrna_gene[a.get('ID')] = a.get('gene') or a.get('Parent','').replace('gene-','')
        if p[2] == 'CDS':
            pid = a.get('protein_id') or (a.get('Parent','').replace('mRNA-','') if a.get('Parent','').startswith('MIMP3') else a.get('Parent'))
            if not pid: continue
            g = a.get('gene') or mrna_gene.get(a.get('Parent'), a.get('Parent'))
            c, s, e = p[0], int(p[3]), int(p[4])
            if pid in prot: prot[pid] = (c, min(prot[pid][1], s), max(prot[pid][2], e), g)
            else: prot[pid] = (c, s, e, g)
            PRODUCT[g] = a.get('product', '')
    return prot
PRODUCT = {}
def hmm_hits(domtbl):
    hits = {}
    for l in open(domtbl):
        if l.startswith('#'): continue
        p = l.split; qlen = int(p[5]); ev = float(p[6]); h1, h2 = int(p[15]), int(p[16])
        if ev <= 1e-10 and (h2 - h1 + 1) / qlen >= 0.35: hits[p[0]] = min(ev, hits.get(p[0], 1))
    return hits
def miniprot(gff):
    rows = []
    for l in open(gff):
        if l.startswith('#'): continue
        p = l.rstrip('\n').split('\t')
        if p[2] != 'mRNA': continue
        a = dict(kv.split('=',1) for kv in p[8].split(';') if '=' in kv)
        q = a['Target'].split[0]; pseudo = ('Frameshift' in a) or ('StopCodon' in a)
        rows.append((p[0], int(p[3]), int(p[4]), p[6], q, float(a['Identity']), pseudo))
    return rows
def merge(rows):
    rows = sorted(rows); loci = []
    for r in rows:
        placed = False
        for L in loci:
            if L['chrom'] != r[0]: continue
            ov = min(L['end'], r[2]) - max(L['start'], r[1]) + 1
            if ov >= 0.5 * min(r[2] - r[1] + 1, L['end'] - L['start'] + 1):
                L['hits'].append(r); L['start'] = min(L['start'], r[1]); L['end'] = max(L['end'], r[2]); placed = True; break
        if not placed: loci.append(dict(chrom=r[0], start=r[1], end=r[2], hits=[r]))
    for L in loci:
        best = max(L['hits'], key=lambda r: r[5]); L['best_query'] = best[4]; L['identity'] = best[5]; L['pseudo'] = best[6]
        L['clade'] = LOCUS.get(best[4], best[4]); L['n_queries'] = len({r[4] for r in L['hits']})
    return sorted(loci, key=lambda L: (L['chrom'], L['start']))
print('\t'.join(['genome','GH18_proteins_hmm','GH18_genes_hmm','GH18_loci_miniprot_genome','interval_chrom','interval_start','interval_end','interval_kb','flank_up_found','flank_dn_found','GH18_genes_hmm_in_interval','annotated_products_in_interval','GH18_loci_miniprot_in_interval','CHIA_array_loci_in_interval','CHIA_array_intact','interval_loci_clades']))
for l in open(genomes):
    name, fa, faa, gff = l.rstrip('\n').split('\t')
    PRODUCT.clear; prot = gff_proteins(gff); prot_product = dict(PRODUCT); hits = hmm_hits(f'{outdir}/{name}.gh18.domtbl')
    genes = {}
    for pid, ev in hits.items:
        if pid in prot: c, s, e, g = prot[pid]; genes.setdefault(g, (c, s, e, ev))
    mp = miniprot(f'{outdir}/{name}.miniprot.gff')
    gh = [r for r in mp if r[4] in LOCUS and LOCUS[r[4]] != 'non-GH18' and r[5] >= 0.75 and (r[2] - r[1]) <= 40000]
    loci = merge(gh)
    fl = [r for r in mp if (r[4] in FLANK_UP or r[4] in FLANK_DN) and r[5] >= 0.40]
    # chromosome carrying most GH18 loci (ties -> the one with flank hits)
    cnt = collections.Counter(L['chrom'] for L in loci); chrom = None
    if cnt:
        top = max(cnt.values); cands = [c for c, n in cnt.items if n == top]
        chrom = max(cands, key=lambda c: sum(1 for r in fl if r[0] == c))
    up = [r for r in fl if r[0] == chrom and r[4] in FLANK_UP]; dn = [r for r in fl if r[0] == chrom and r[4] in FLANK_DN]
    cl = [L for L in loci if L['chrom'] == chrom]
    if up and dn: s0 = min(min(r[1] for r in up), min(r[1] for r in dn)); e0 = max(max(r[2] for r in up), max(r[2] for r in dn))
    elif cl: s0 = min(L['start'] for L in cl) - 300000; e0 = max(L['end'] for L in cl) + 300000
    else: s0 = e0 = 0
    inl = [L for L in cl if L['start'] >= s0 and L['end'] <= e0]
    ing = [g for g, (c, s, e, ev) in genes.items if c == chrom and s >= s0 and e <= e0]
    chia = sum(1 for L in inl if L['clade'].startswith('CHIA-array')); chia_ok = sum(1 for L in inl if L['clade'].startswith('CHIA-array') and not L['pseudo'])
    prods = collections.Counter(re.sub(r'[- ]like.*|,.*| isoform.*', '', prot_product.get(g, '?')) for g in ing)
    print('\t'.join(map(str, [name, len(hits), len(genes), len(loci), chrom, s0, e0, round((e0 - s0) / 1000, 1), len({r[4] for r in up}), len({r[4] for r in dn}), len(ing), '; '.join(f'{k} x{v}' for k, v in prods.most_common), len(inl), chia, chia_ok, ';'.join(f"{L['clade']}@{L['start']}-{L['end']}({L['identity']:.2f}{',pseudo' if L['pseudo'] else ''})" for L in inl)])))
    with open(f'{outdir}/{name}.loci.tsv', 'w') as fh:
        fh.write('chrom\tstart\tend\tbest_query\tclade\tidentity\tpseudogene\tn_queries\tin_interval\thmm_genes_overlapping\tproducts\n')
        for L in loci:
            ov = [g for g, (c, s, e, ev) in genes.items if c == L['chrom'] and s <= L['end'] and e >= L['start']]
            fh.write(f"{L['chrom']}\t{L['start']}\t{L['end']}\t{L['best_query']}\t{L['clade']}\t{L['identity']:.3f}\t{L['pseudo']}\t{L['n_queries']}\t{L in inl}\t{','.join(ov)}\t{' | '.join(prot_product.get(g, '') for g in ov)}\n")
        fh.write('# annotated GH18 genes in the interval not overlapped by a miniprot locus:\n')
        for g in ing:
            if not any(genes[g][1] <= L['end'] and genes[g][2] >= L['start'] for L in inl): fh.write(f"# {g}\t{genes[g][0]}:{genes[g][1]}-{genes[g][2]}\t{prot_product.get(g, '')}\n")
