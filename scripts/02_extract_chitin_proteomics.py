#!/usr/bin/env python3
"""
02_extract_chitin_proteomics.py

Extract chitin-related proteins from proteomics vendor results (Datasets 04 & 05).
Searches annot_file.xls for chitin/chitinase/hexosaminidase/GlcNAc keywords.

Outputs:
  - results/02_proteomics_chitin_proteins.tsv
  - results/02_proteomics_chitin_report.txt

Usage:
  ssh <user>@<analysis-server> "python3 /path/to/this/script.py"
"""

import os
import sys
import csv
import re
from pathlib import Path

BASE = "$PROJECT/02.metagenomics-proteomics"

DATASETS = {
    "Proteomics_1st_GI": {
        "annot_file": os.path.join(BASE,
            "04_Manouria_impressa_Proteomics_First_Sequence_<order-id>",
            "<provider-delivery-dir>/report/BGI_result",
            "3.Diff_protein_statistics/annot_file.xls"),
        "description": "1st proteomics: 6 GI segments × 3 replicates",
    },
    "Proteomics_2nd_MultiTissue": {
        "annot_file": os.path.join(BASE,
            "05_Manouria_impressa_Proteomics_Second_Sequence_20251209",
            "<provider-delivery-dir>/<order-id>/report/BGI_result",
            "3.Diff_protein_statistics/annot_file.xls"),
        "description": "2nd proteomics: 9 tissues × 2 replicates",
    },
}

CHITIN_KEYWORDS = [
    r'chitin', r'chitinase', r'chitobiase', r'hexosaminidase',
    r'GlcNAc', r'N-acetylglucosamine', r'GH18', r'GH19', r'GH20',
    r'KOG2806',  # Chitinase KOG category
    r'K01183',   # KEGG chitinase
    r'K12310',   # Di-N-acetylchitobiase
    r'3\.2\.1\.14',  # EC chitinase
    r'3\.2\.1\.52',  # EC hexosaminidase
]

OUTPUT_DIR = "$PROJECT/results"


def search_annot_file(filepath, dataset_name):
    """Search annot_file.xls for chitin-related proteins."""
    results = []
    pattern = re.compile('|'.join(CHITIN_KEYWORDS), re.IGNORECASE)

    if not os.path.exists(filepath):
        # Try to find it
        parent = os.path.dirname(os.path.dirname(filepath))
        for root, dirs, files in os.walk(parent):
            for f in files:
                if f == 'annot_file.xls':
                    filepath = os.path.join(root, f)
                    break

    if not os.path.exists(filepath):
        print(f"  WARNING: annot_file.xls not found for {dataset_name}")
        return [], []

    headers = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if i == 0:
                headers = row
                continue
            line_text = '\t'.join(row)
            if pattern.search(line_text):
                results.append(row)

    return headers, results


def main:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("CHITIN-RELATED PROTEOMICS EXTRACTION REPORT")
    report_lines.append("=" * 80)

    all_proteins = []

    for ds_name, info in DATASETS.items:
        report_lines.append(f"\n{'=' * 60}")
        report_lines.append(f"Dataset: {ds_name}")
        report_lines.append(f"Description: {info['description']}")
        report_lines.append(f"{'=' * 60}")

        headers, hits = search_annot_file(info["annot_file"], ds_name)

        if not hits:
            report_lines.append(f"  No chitin-related proteins found (or file not accessible)")
            continue

        report_lines.append(f"  Found {len(hits)} chitin-related protein entries")
        report_lines.append(f"  Columns: {len(headers)}")

        # Identify key column indices
        protein_col = 0  # Usually first column
        # Find intensity columns (contain sample names like Stomach_1, Duodenum_1 etc)
        intensity_cols = []
        annotation_cols = []
        for idx, h in enumerate(headers):
            h_lower = h.lower
            if any(tissue in h_lower for tissue in ['stomach', 'duodenum', 'jejunum', 'ileum',
                                                      'colon', 'rectum', 'liver', 'heart', 'muscle']):
                intensity_cols.append(idx)
            elif idx > 0 and not any(tissue in h_lower for tissue in ['stomach', 'duodenum', 'jejunum',
                                                                       'ileum', 'colon', 'rectum',
                                                                       'liver', 'heart', 'muscle']):
                annotation_cols.append(idx)

        report_lines.append(f"  Intensity columns: {len(intensity_cols)}")
        report_lines.append(f"  Sample names: {[headers[i] for i in intensity_cols]}")

        # Process each hit
        for row in hits:
            protein_id = row[protein_col] if len(row) > protein_col else "unknown"

            # Get intensities
            intensities = {}
            for idx in intensity_cols:
                if idx < len(row):
                    sample_name = headers[idx]
                    try:
                        intensities[sample_name] = float(row[idx])
                    except (ValueError, TypeError):
                        intensities[sample_name] = None

            # Get annotations (join remaining columns)
            annot_text = ""
            for idx in annotation_cols:
                if idx < len(row) and row[idx].strip:
                    annot_text += f" | {headers[idx]}: {row[idx]}"

            # Classify the protein
            row_text = '\t'.join(row).lower
            if re.search(r'chitinase|3\.2\.1\.14|K01183|acidic mammalian', row_text, re.IGNORECASE):
                category = "Chitinase (GH18)"
            elif re.search(r'hexosaminidase|3\.2\.1\.52', row_text, re.IGNORECASE):
                category = "Hexosaminidase (GH20)"
            elif re.search(r'chitobiase|K12310', row_text, re.IGNORECASE):
                category = "Di-N-acetylchitobiase"
            elif re.search(r'CHID1|chitinase domain', row_text, re.IGNORECASE):
                category = "Chitinase-like (CHID1)"
            elif re.search(r'GlcNAc|N-acetylglucosamine|O-GlcNAc', row_text, re.IGNORECASE):
                category = "GlcNAc metabolism"
            else:
                category = "Other chitin-related"

            # Calculate tissue pattern
            if intensities:
                valid_vals = {k: v for k, v in intensities.items if v is not None}
                if valid_vals:
                    max_tissue = max(valid_vals, key=valid_vals.get)
                    min_tissue = min(valid_vals, key=valid_vals.get)
                    avg_val = sum(valid_vals.values) / len(valid_vals)
                    pattern = f"Highest: {max_tissue} ({valid_vals[max_tissue]:.1f}), Avg: {avg_val:.1f}"
                else:
                    pattern = "No valid intensities"
            else:
                pattern = "No intensity data"

            report_lines.append(f"\n  Protein: {protein_id}")
            report_lines.append(f"    Category: {category}")
            report_lines.append(f"    Pattern: {pattern}")

            # Store for output
            entry = {
                "Dataset": ds_name,
                "Protein_ID": protein_id,
                "Category": category,
                "Tissue_Pattern": pattern,
            }
            for sample, val in intensities.items:
                entry[sample] = val if val is not None else ""
            entry["Annotations"] = annot_text.strip(" |")

            all_proteins.append(entry)

    # Write combined output
    if all_proteins:
        # Get all possible column names
        all_cols = ["Dataset", "Protein_ID", "Category", "Tissue_Pattern"]
        sample_cols = set
        for p in all_proteins:
            for k in p.keys:
                if k not in all_cols and k != "Annotations":
                    sample_cols.add(k)
        all_cols.extend(sorted(sample_cols))
        all_cols.append("Annotations")

        outfile = os.path.join(OUTPUT_DIR, "02_proteomics_chitin_proteins.tsv")
        with open(outfile, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=all_cols, delimiter='\t',
                                    extrasaction='ignore')
            writer.writeheader
            writer.writerows(all_proteins)

        report_lines.append(f"\n\n{'=' * 60}")
        report_lines.append(f"SUMMARY")
        report_lines.append(f"{'=' * 60}")
        report_lines.append(f"Total chitin-related proteins: {len(all_proteins)}")

        # Count by category
        from collections import Counter
        cats = Counter(p["Category"] for p in all_proteins)
        for cat, count in cats.most_common:
            report_lines.append(f"  {cat}: {count}")

        report_lines.append(f"\nOutput: {outfile}")

    # Write report
    report_text = "\n".join(report_lines)
    report_file = os.path.join(OUTPUT_DIR, "02_proteomics_chitin_report.txt")
    with open(report_file, 'w') as f:
        f.write(report_text)

    print(report_text)
    print(f"\nReport saved: {report_file}")


if __name__ == "__main__":
    main
