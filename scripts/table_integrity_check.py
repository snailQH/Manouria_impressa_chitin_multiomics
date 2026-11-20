#!/usr/bin/env python3
"""Catch tables that silently misstate their data.

The figure equivalent (scripts/figure_integrity_check.py) looks for data drawn outside the view.
For tables the failure modes are different, so this checks:

  1. every tables/*.tsv parses, with one consistent column count
  2. no duplicated identifier in a key-looking first column
  3. no column that is 100 % a single missing-like value (the signature of a failed join:
     exactly how the metagenome annotation bug read '-' for all 213 groups)
  4. no numeric column that is entirely NaN
  5. derived columns agree with their inputs (X_CPM == X / mapped_reads * 1e6)
  6. every tables/... path cited in the manuscript exists on disk

Findings are printed for a human to judge; some are legitimate (a family genuinely absent).
    python3 scripts/table_integrity_check.py
"""
import os, re, sys, glob
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MISSING = {'-', '', 'nan', 'NaN', 'None', 'n/a', 'NA', 'unknown', 'unassigned'}
# checked and expected — kept here so a NEW instance of either pattern stands out
EXPECTED = {
    # multi-block TSVs: a data block plus a '# Tests' block. The workbook builder reads line by
    # line, so each block ships with its own header; only pandas-style parsing sees them as ragged.
    ('41_unified_gh18_vs_bacteroides.tsv', 'RAGGED'),
    ('73_diann_microbial_chitin_spatial.tsv', 'RAGGED'),
    ('27_DS01_kegg_functional_summary.tsv', 'RAGGED'),
    ('28_DS03_kegg_functional_summary.tsv', 'RAGGED'),
    ('44_DS02_kegg_functional_summary.tsv', 'RAGGED'),
    # genuinely absent, not a failed join: dbCAN3 calls no CBM12 in any cohort
    ('35_unified_cazy_per_sample.tsv', 'ALL ONE VALUE'),
    # structural: the Kraken2 database has no plant library, so Viridiplantae cannot be non-zero
    ('67_fungal_profiling_per_sample.tsv', 'ALL ONE VALUE'),
}
out = []


def rows_of(path):
    return [l for l in open(path, encoding='utf-8', errors='replace') if not l.startswith('#') and l.strip]


def main:
    files = sorted(glob.glob(os.path.join(BASE, 'tables', '**', '*.tsv'), recursive=True))
    files = [f for f in files if 'archive' not in f]
    print(f'checking {len(files)} tables\n')
    for f in files:
        rel = os.path.relpath(f, BASE)
        body = rows_of(f)
        if len(body) < 2:
            continue
        widths = {l.count('\t') + 1 for l in body}
        if len(widths) > 1:
            out.append((rel, 'RAGGED', f'inconsistent column counts {sorted(widths)} — will not parse as one table'))
            continue
        try:
            t = pd.read_csv(f, sep='\t', comment='#', low_memory=False)
        except Exception as e:
            out.append((rel, 'UNPARSEABLE', str(e)[:90])); continue
        if t.empty:
            continue
        key = t.columns[0]
        if t[key].dtype == object and t[key].nunique < len(t) and len(t) > 3:
            dup = t[key][t[key].duplicated].unique[:3]
            if t[key].nunique > len(t) * 0.5:
                out.append((rel, 'DUP KEY', f'{key}: {len(t) - t[key].nunique} repeated, e.g. {list(dup)}'))
        for c in t.columns:
            v = t[c]
            if v.isna.all:
                out.append((rel, 'ALL EMPTY', f'{c}: every value NaN ({len(t)} rows)')); continue
            u = set(v.dropna.astype(str).unique)
            if len(u) == 1 and (u & MISSING or u <= {'0', '0.0'}):
                out.append((rel, 'ALL ONE VALUE', f"{c}: every row is '{list(u)[0]}' ({len(t)} rows)"))
        # derived-column arithmetic
        denom = next((d for d in ('mapped_reads', 'nonhost_read_pairs', 'read_pairs') if d in t.columns), None)
        if denom:
            for c in t.columns:
                if not c.endswith('_CPM'):
                    continue
                src = c[:-4]
                if src in t.columns and pd.api.types.is_numeric_dtype(t[src]) and pd.api.types.is_numeric_dtype(t[denom]):
                    want = t[src] / t[denom] * 1e6
                    bad = (want - t[c]).abs > np.maximum(1e-6, want.abs * 1e-6)
                    if bad.any:
                        out.append((rel, 'DERIVED MISMATCH', f'{c} != {src}/{denom}*1e6 in {int(bad.sum)} row(s)'))
    # citations -> files
    _man_path = os.path.join(BASE, 'docs', 'MANUSCRIPT_DRAFT.md')
    if os.path.exists(_man_path):
        man = open(_man_path, encoding='utf-8').read
    else:
        # the published repository ships the analysis, not the manuscript. Checks 1-5 still apply,
        # so fall through to the normal report and skip only the citation resolution.
        print('note: docs/MANUSCRIPT_DRAFT.md is not part of this tree; skipping the citation check')
        man = ''
    cited = sorted({m for m in re.findall(r'`(tables/[A-Za-z0-9_./-]+)`', man)})
    for c in cited:
        if os.path.exists(os.path.join(BASE, c)):
            continue
        # `tables/06` is shorthand; resolve it, and flag only a dead or ambiguous reference
        stem = os.path.basename(c)
        hits = sorted(os.path.basename(f) for f in glob.glob(os.path.join(BASE, 'tables', stem + '_*')))
        if len(hits) == 1:
            continue
        kind = 'CITED BUT MISSING' if not hits else 'AMBIGUOUS CITATION'
        out.append((c, kind, 'no such table' if not hits else f'{len(hits)} files share this number: {hits}'))
    new = [f for f in out if (os.path.basename(f[0]), f[1]) not in EXPECTED]
    known = len(out) - len(new)
    print(f'{len(new)} unexpected finding(s); {known} known and documented:\n')
    for rel, kind, msg in new:
        print(f'  [{kind:17s}] {os.path.basename(rel):52s} {msg}')
    print(f'\n{len(cited)} table paths cited in the manuscript; '
          f'{sum(1 for c in cited if os.path.exists(os.path.join(BASE, c)))} resolve')
    return 0


if __name__ == '__main__':
    sys.exit(main)
