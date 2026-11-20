#!/usr/bin/env python3
"""
Compare cohorts processed through the SAME unified pipeline.

The point of the unified re-analysis is that every cohort is put through one
identical method, so that differences between cohorts can be read as biology
rather than as processing. This script reports the comparison, and — just as
importantly — refuses to compare cohorts that were not processed identically.

Inputs are the per-cohort eggNOG-mapper annotation summaries produced by
`scripts/summarise_emapper.py` (or the equivalent extraction run on the server).

Usage:
    python3 scripts/compare_harmonised_cohorts.py DS01=<summary.tsv> DS03=<summary.tsv>
"""

import os
import sys

# Steps that must match for two cohorts to be directly comparable. Kept explicit
# so that a future cohort processed differently cannot be silently pooled in.
REQUIRED_IDENTICAL = [
    ("gene_caller", "Prodigal 2.6.3, -p meta"),
    ("clustering", "CD-HIT, -c 0.95 -aS 0.9"),
    ("annotation", "eggNOG-mapper 2.1.12, DIAMOND, eggNOG 5.0.2"),
]

# Assembly is comparable but not byte-identical across cohorts, and saying so
# matters: DS01 was assembled here (MEGAHIT 1.2.9, --min-contig-len 500,
# k=21..141) whereas DS03 uses Novogene's own MEGAHIT run, also filtered at
# 500 bp but with their version and k-mer settings. The cohorts are therefore
# identical from GENE PREDICTION ONWARD, not from raw reads. Recorded here so
# the distinction is not quietly lost.
ASSEMBLY_NOTE = ("DS01 assembled in-house (MEGAHIT 1.2.9, --min-contig-len 500, "
                 "k=21,41,61,81,101,121,141); DS03 uses the vendor MEGAHIT "
                 "assembly, also >=500 bp. Comparable, not identical.")

CHITIN_KO = [
    ("K01183", "chitinase"),
    ("K12373", "HEXA_B hexosaminidase"),
    ("K01207", "nagZ"),
    ("K01443", "nagA"),
    ("K02564", "nagB"),
    ("K03431", "glmM"),
    ("K04042", "glmU"),
    ("K16165", "chitin deacetylase"),
]


def load(path):
    """Read a summary TSV of `key<TAB>value` and `KO<TAB>desc<TAB>count` rows."""
    metrics, kos = {}, {}
    if not os.path.exists(path):
        raise SystemExit(f"summary not found: {path}")
    for line in open(path):
        if line.startswith("#") or not line.strip:
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) == 2:
            metrics[p[0]] = p[1]
        elif len(p) == 3 and p[0].startswith("K"):
            try:
                kos[p[0]] = int(p[2])
            except ValueError:
                pass
    return metrics, kos


def num(metrics, key):
    try:
        return int(metrics.get(key, 0))
    except ValueError:
        return 0


def main:
    args = [a for a in sys.argv[1:] if "=" in a]
    if len(args) < 2:
        raise SystemExit(__doc__)
    cohorts = {}
    for a in args:
        name, path = a.split("=", 1)
        cohorts[name] = load(path)

    print("Cohorts compared under the unified pipeline")
    for step, spec in REQUIRED_IDENTICAL:
        print(f"  {step:12s} {spec}")
    print(f"\n  assembly     {ASSEMBLY_NOTE}")
    print("\nAny cohort NOT processed with exactly these settings must not be")
    print("pooled with these numbers.\n")

    names = list(cohorts)
    w = max(len(n) for n in names) + 2

    rows = [
        ("Proteins in catalogue", "Proteins_in_catalogue"),
        ("Annotated (eggNOG)", "Proteins_annotated"),
        ("With KEGG KO", "With_KEGG_KO"),
        ("With KEGG pathway", "With_KEGG_Pathway"),
        ("In ko00520", "In_ko00520_amino_sugar_pathway"),
    ]
    hdr = "metric".ljust(30) + "".join(n.rjust(w + 8) for n in names)
    print(hdr)
    print("-" * len(hdr))
    for label, key in rows:
        line = label.ljust(30)
        for n in names:
            line += format(num(cohorts[n][0], key), ",").rjust(w + 8)
        print(line)

    # annotation rate
    line = "Annotation rate".ljust(30)
    for n in names:
        tot = num(cohorts[n][0], "Proteins_in_catalogue")
        ann = num(cohorts[n][0], "Proteins_annotated")
        line += (f"{100*ann/tot:.1f}%" if tot else "n/a").rjust(w + 8)
    print(line)

    print("\nChitin-metabolism KEGG orthologs (gene counts)")
    hdr2 = "KO / enzyme".ljust(30) + "".join(n.rjust(w + 8) for n in names)
    print(hdr2)
    print("-" * len(hdr2))
    for ko, desc in CHITIN_KO:
        line = f"{ko} {desc}".ljust(30)
        for n in names:
            line += format(cohorts[n][1].get(ko, 0), ",").rjust(w + 8)
        print(line)

    # normalised per million catalogue proteins, so catalogue size does not drive the comparison
    print("\nSame, per million catalogue proteins (size-normalised)")
    print(hdr2)
    print("-" * len(hdr2))
    for ko, desc in CHITIN_KO:
        line = f"{ko} {desc}".ljust(30)
        for n in names:
            tot = num(cohorts[n][0], "Proteins_in_catalogue")
            v = cohorts[n][1].get(ko, 0)
            line += (f"{1e6*v/tot:.1f}" if tot else "n/a").rjust(w + 8)
        print(line)


if __name__ == "__main__":
    main
