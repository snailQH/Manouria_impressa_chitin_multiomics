#!/usr/bin/env python3
"""
Summarise an eggNOG-mapper annotation file into the compact TSV consumed by
scripts/compare_harmonised_cohorts.py.

Run once per cohort on the server, e.g.

    python3 summarise_emapper.py 07_kegg/nr.emapper.annotations DS01_kegg_summary.tsv 2052794
    python3 summarise_emapper.py 09_ds03_harmonised/kegg/ds03.emapper.annotations DS03_kegg_summary.tsv <n>

The third argument is the number of proteins in the input catalogue, which the
annotation file itself does not record (it lists only annotated proteins).
"""
import sys, collections

F = sys.argv[1]
OUT = sys.argv[2]
CATALOG_N = int(sys.argv[3]) if len(sys.argv) > 3 else 0

CHITIN_KO = {
    "K01183": "chitinase (EC 3.2.1.14)",
    "K12373": "HEXA_B, beta-N-acetylhexosaminidase",
    "K01443": "nagA, N-acetylglucosamine-6-P deacetylase",
    "K02564": "nagB, glucosamine-6-P deaminase",
    "K01452": "chitin deacetylase (EC 3.5.1.41)",
    "K03478": "chbG, chitin disaccharide deacetylase (EC 3.5.1.105)",
    "K00884": "nagK, N-acetylglucosamine kinase",
    "K01207": "nagZ, beta-N-acetylhexosaminidase",
    "K03431": "glmM, phosphoglucosamine mutase",
    "K04042": "glmU, GlcNAc-1-P uridyltransferase",
}
hdr = None
n_total = n_ko = n_path = n_cazy = n_cog = 0
ko_counts = collections.Counter
cog_counts = collections.Counter
ko00520 = 0
cazy_fam = collections.Counter

with open(F) as fh:
    for line in fh:
        if line.startswith("#"):
            if line.startswith("#query"):
                hdr = line.rstrip("\n").split("\t")
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 21:
            continue
        n_total += 1
        ko = p[11]
        if ko and ko != "-":
            n_ko += 1
            for k in ko.split(","):
                k = k.replace("ko:", "").strip
                if k:
                    ko_counts[k] += 1
        pw = p[12]
        if pw and pw != "-":
            n_path += 1
            if "ko00520" in pw or "map00520" in pw:
                ko00520 += 1
        cog = p[6]
        if cog and cog != "-":
            n_cog += 1
            for c in cog.strip:
                cog_counts[c] += 1
        cz = p[18]
        if cz and cz != "-":
            n_cazy += 1
            for f in cz.split(","):
                f = f.strip
                if f:
                    cazy_fam[f.split("(")[0]] += 1

with open(OUT, "w") as o:
    o.write("# eggNOG-mapper functional annotation of the DS01 non-redundant gene catalogue\n")
    o.write("# Source: 07_kegg/nr.emapper.annotations (eggNOG-mapper 2.1.12, DIAMOND, eggNOG 5.0.2)\n")
    o.write("Metric\tValue\n")
    o.write(f"Proteins_in_catalogue\t{CATALOG_N}\n")
    o.write(f"Proteins_annotated\t{n_total}\n")
    o.write(f"With_KEGG_KO\t{n_ko}\n")
    o.write(f"With_KEGG_Pathway\t{n_path}\n")
    o.write(f"With_COG_category\t{n_cog}\n")
    o.write(f"With_CAZy_annotation\t{n_cazy}\n")
    o.write(f"In_ko00520_amino_sugar_pathway\t{ko00520}\n")
    o.write("#\n# Chitin-metabolism KEGG orthologs\nKO\tDescription\tGene_count\n")
    for k, d in CHITIN_KO.items:
        o.write(f"{k}\t{d}\t{ko_counts.get(k,0)}\n")
    o.write("#\n# Top COG functional categories\nCOG\tGene_count\n")
    for c, v in cog_counts.most_common(12):
        o.write(f"{c}\t{v}\n")
    o.write("#\n# Top CAZy families (eggNOG-assigned)\nFamily\tGene_count\n")
    for f, v in cazy_fam.most_common(15):
        o.write(f"{f}\t{v}\n")

print(f"annotated: {n_total:,} / {CATALOG_N:,}  ({100*n_total/CATALOG_N:.1f}%)" if CATALOG_N else f"annotated: {n_total:,}")
print(f"  KEGG KO      : {n_ko:,}")
print(f"  KEGG pathway : {n_path:,}")
print(f"  COG category : {n_cog:,}")
print(f"  CAZy         : {n_cazy:,}")
print(f"  ko00520      : {ko00520:,}")
print("\nchitin KOs:")
for k, d in CHITIN_KO.items:
    print(f"  {k}  {ko_counts.get(k,0):>8,}  {d}")
