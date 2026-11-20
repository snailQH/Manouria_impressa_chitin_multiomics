#!/usr/bin/env python3
"""
Build per-sample CAZy family abundance for a cohort.

Joins two things, both produced by the unified pipeline:
  1. gene -> CAZy family, from the eggNOG-mapper CAZy column
  2. gene -> per-sample read count, from Bowtie2 mapping onto the gene catalogue

The result is per-sample abundance for each CAZy family, in read counts, which
is the same unit the Novogene cohorts report — so cohorts can be plotted on one
axis without any normalisation caveat.

CAZy families come from eggNOG rather than dbCAN deliberately: the dbCAN
database is no longer served upstream, whereas eggNOG annotation exists for
every cohort under identical settings. Using one source for all cohorts matters
more here than which source it is.

Usage:
    python3 scripts/extract_cazy_abundance.py \
        --annotations <cohort>.emapper.annotations \
        --counts-dir  <dir of *.counts.tsv> \
        --out         <cohort>_cazy_abundance.tsv
"""

import argparse
import os
import re
import sys
from collections import defaultdict

# Families reported in the manuscript. Prefix match, so GH18 also catches
# subfamily notation such as GH18_1 if eggNOG emits it.
FAMILIES_OF_INTEREST = ["GH18", "GH19", "GH20", "CBM5", "CBM12", "CBM50", "AA10"]


def parse_cazy_field(value):
    """eggNOG CAZy column looks like 'GH18' or 'GH18,CBM5' or 'GH13_28'."""
    out = set
    if not value or value == "-":
        return out
    for part in value.split(","):
        part = part.strip
        if not part:
            continue
        # normalise GH13_28 -> GH13 for family-level counting, keep both
        out.add(part)
        m = re.match(r"([A-Z]+\d+)", part)
        if m:
            out.add(m.group(1))
    return out


def load_gene_families(path):
    """gene id -> set of CAZy families, from the eggNOG annotation table."""
    gene2fam = {}
    cazy_col = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("#query"):
                    hdr = line.rstrip("\n").split("\t")
                    if "CAZy" in hdr:
                        cazy_col = hdr.index("CAZy")
                continue
            p = line.rstrip("\n").split("\t")
            col = cazy_col if cazy_col is not None else 18
            if len(p) <= col:
                continue
            fams = parse_cazy_field(p[col])
            if fams:
                gene2fam[p[0]] = fams
    return gene2fam


def load_counts(counts_dir):
    """sample -> {gene id: read count}, from samtools idxstats output."""
    per_sample = {}
    for fn in sorted(os.listdir(counts_dir)):
        if not fn.endswith(".counts.tsv"):
            continue
        sample = fn[: -len(".counts.tsv")]
        d = {}
        with open(os.path.join(counts_dir, fn)) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                try:
                    c = int(p[1])
                except ValueError:
                    continue
                if c:
                    d[p[0]] = c
        per_sample[sample] = d
    return per_sample


def main:
    ap = argparse.ArgumentParser
    ap.add_argument("--annotations", required=True)
    ap.add_argument("--counts-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args

    gene2fam = load_gene_families(a.annotations)
    print(f"genes with a CAZy assignment: {len(gene2fam):,}", file=sys.stderr)

    per_sample = load_counts(a.counts_dir)
    if not per_sample:
        raise SystemExit(f"no *.counts.tsv found in {a.counts_dir}")
    print(f"samples with counts: {len(per_sample)}", file=sys.stderr)

    # sanity check: do the identifier spaces actually intersect?
    any_sample = next(iter(per_sample.values))
    overlap = len(set(list(gene2fam)[:50000]) & set(any_sample))
    if overlap == 0:
        raise SystemExit(
            "FATAL: no gene identifiers shared between the annotation file and the\n"
            "count files. These must come from the SAME catalogue build — check that\n"
            "both were produced from the rebuilt, sample-prefixed catalogue."
        )

    samples = sorted(per_sample)
    table = defaultdict(lambda: defaultdict(int))
    total = defaultdict(int)
    for s in samples:
        for gene, c in per_sample[s].items:
            total[s] += c
            for fam in gene2fam.get(gene, ):
                table[fam][s] += c

    fams = [f for f in FAMILIES_OF_INTEREST if f in table]
    extra = sorted(f for f in table if f not in fams and table[f])
    with open(a.out, "w") as o:
        o.write("# Per-sample CAZy family abundance (Bowtie2 read counts)\n")
        o.write("# CAZy families from the eggNOG-mapper CAZy column\n")
        o.write("Family\t" + "\t".join(samples) + "\tTotal\n")
        for fam in fams + extra:
            row = [str(table[fam].get(s, 0)) for s in samples]
            o.write(fam + "\t" + "\t".join(row) + "\t" + str(sum(table[fam].values)) + "\n")
        o.write("TOTAL_MAPPED\t" + "\t".join(str(total[s]) for s in samples)
                + "\t" + str(sum(total.values)) + "\n")

    print(f"wrote {a.out}", file=sys.stderr)
    for fam in fams:
        vals = [table[fam].get(s, 0) for s in samples]
        print(f"  {fam:7s} total={sum(vals):>12,}  per-sample mean={sum(vals)/len(vals):>12,.0f}",
              file=sys.stderr)


if __name__ == "__main__":
    main
