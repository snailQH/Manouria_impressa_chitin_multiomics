#!/usr/bin/env python3
"""
01_extract_chitin_cazyme_data.py

Extract chitin-related CAZyme data from all 3 metagenomics vendor results.
Focuses on GH18, GH19, GH20, CBM5, CBM12, AA10 families.

Outputs:
  - results/01_cazyme_chitin_summary.tsv  (combined chitin CAZyme table)
  - results/01_cazyme_all_GH_summary.tsv  (all GH family abundances)
  - results/01_chitin_cazyme_report.txt   (text summary)

Usage:
  ssh <user>@<analysis-server> "python3 /path/to/this/script.py"
  OR run locally after scp-ing the annotation files
"""

import os
import sys
import glob
import csv
from collections import defaultdict
from pathlib import Path

# Server paths
BASE = "$PROJECT/02.metagenomics-proteomics"
DATASETS = {
    "P_eryngii": {
        "path": os.path.join(BASE, "01_Pleurotus_eryngii_12samples-metagenomics/report/upload/9_annotation/CAZy"),
        "diet": "Pleurotus eryngii (King Oyster Mushroom)",
        "samples_prefix": "ITF",
    },
    "H_tessellatus": {
        "path": os.path.join(BASE, "02_Hypsizygus_tessellatus_metagenomics/<provider-delivery-dir>/04.Annotation/CAZy"),
        "diet": "Hypsizygus tessellatus (Brown Beech Mushroom)",
        "samples_prefix": "HM",
    },
    "P_ostreatus": {
        "path": os.path.join(BASE, "03_Pleurotus_ostreatus_2024-12-09-Manouria_impressa_metagenomics/04.Annotation/CAZy"),
        "diet": "Pleurotus ostreatus (Oyster Mushroom)",
        "samples_prefix": "PM",
    },
}

# Chitin-related CAZyme families
CHITIN_FAMILIES = {
    "GH18": "Chitinase (endo-type, primary chitin degradation)",
    "GH19": "Chitinase (plant-type)",
    "GH20": "Beta-N-acetylhexosaminidase (GlcNAc release)",
    "GH23": "Lysozyme/chitinase (peptidoglycan + chitin)",
    "CBM5": "Chitin-binding module (bacterial)",
    "CBM12": "Chitin-binding module",
    "CBM14": "Chitin-binding module (invertebrate)",
    "AA10": "LPMO (lytic polysaccharide monooxygenase, chitin oxidation)",
}

OUTPUT_DIR = "$PROJECT/results"


def find_cazyme_files(cazy_dir):
    """Find CAZyme abundance files in vendor output directory."""
    candidates = []
    for root, dirs, files in os.walk(cazy_dir):
        for f in files:
            fp = os.path.join(root, f)
            fl = f.lower
            if any(kw in fl for kw in ['absolute', 'abundance', 'genenums', 'anno']):
                if fl.endswith(('.xls', '.tsv', '.txt', '.csv')):
                    candidates.append(fp)
    return candidates


def parse_abundance_file(filepath):
    """
    Parse a CAZyme abundance file. Expected format:
    First column: CAZyme family name (e.g., GH18)
    Remaining columns: sample abundances
    Returns: (headers, data_dict) where data_dict[family] = {sample: count}
    """
    data = {}
    headers = []

    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            reader = csv.reader(f, delimiter='\t')
            for i, row in enumerate(reader):
                if not row or len(row) < 2:
                    continue
                if i == 0:
                    headers = row[1:]  # sample names
                    continue
                family = row[0].strip
                if not family:
                    continue
                values = {}
                for j, val in enumerate(row[1:]):
                    if j < len(headers):
                        try:
                            values[headers[j]] = float(val)
                        except (ValueError, TypeError):
                            values[headers[j]] = 0.0
                data[family] = values
    except Exception as e:
        print(f"  Warning: Could not parse {filepath}: {e}", file=sys.stderr)

    return headers, data


def main:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("CHITIN-RELATED CAZyme EXTRACTION REPORT")
    report_lines.append(f"Generated from vendor metagenomics results")
    report_lines.append("=" * 80)

    all_chitin_data = []  # For combined output
    all_gh_data = []  # All GH families

    for dataset_name, info in DATASETS.items:
        cazy_dir = info["path"]
        diet = info["diet"]

        report_lines.append(f"\n{'=' * 60}")
        report_lines.append(f"Dataset: {dataset_name}")
        report_lines.append(f"Diet: {diet}")
        report_lines.append(f"CAZy path: {cazy_dir}")
        report_lines.append(f"{'=' * 60}")

        if not os.path.exists(cazy_dir):
            report_lines.append(f"  ERROR: Directory not found!")
            # Try alternative paths
            alt_paths = glob.glob(os.path.join(os.path.dirname(cazy_dir), "**/CAZy*"), recursive=True)
            if alt_paths:
                report_lines.append(f"  Found alternatives: {alt_paths[:3]}")
            continue

        # Find abundance files
        files = find_cazyme_files(cazy_dir)
        report_lines.append(f"  Found {len(files)} potential CAZyme files:")
        for fp in files:
            report_lines.append(f"    - {os.path.basename(fp)} ({os.path.getsize(fp)} bytes)")

        # Try to parse each file
        best_headers = []
        best_data = {}
        for fp in files:
            headers, data = parse_abundance_file(fp)
            if len(data) > len(best_data):
                best_headers = headers
                best_data = data
                best_file = fp

        if not best_data:
            report_lines.append(f"  WARNING: No parseable CAZyme abundance data found!")
            continue

        report_lines.append(f"\n  Best file: {os.path.basename(best_file)}")
        report_lines.append(f"  Total CAZyme families: {len(best_data)}")
        report_lines.append(f"  Samples: {best_headers}")

        # Extract chitin-related families
        report_lines.append(f"\n  --- Chitin-Related Families ---")
        found_chitin = False
        for family, description in CHITIN_FAMILIES.items:
            if family in best_data:
                found_chitin = True
                values = best_data[family]
                avg = sum(values.values) / len(values) if values else 0
                total = sum(values.values)
                report_lines.append(f"  {family}: avg={avg:.1f}, total={total:.1f} -- {description}")

                for sample, count in values.items:
                    all_chitin_data.append({
                        "Dataset": dataset_name,
                        "Diet": diet,
                        "Sample": sample,
                        "CAZy_Family": family,
                        "Description": description,
                        "Abundance": count,
                    })

        if not found_chitin:
            report_lines.append(f"  No chitin families found by exact name match.")
            # Try partial matching
            report_lines.append(f"  Trying partial match for 'GH1' or 'GH2' families...")
            for family in sorted(best_data.keys):
                if family.startswith("GH"):
                    values = best_data[family]
                    avg = sum(values.values) / len(values) if values else 0
                    report_lines.append(f"    {family}: avg={avg:.1f}")

        # Extract ALL GH families
        report_lines.append(f"\n  --- All GH (Glycoside Hydrolase) Families ---")
        gh_count = 0
        for family in sorted(best_data.keys):
            if family.startswith("GH"):
                gh_count += 1
                values = best_data[family]
                avg = sum(values.values) / len(values) if values else 0
                total = sum(values.values)
                report_lines.append(f"  {family}: avg={avg:.1f}, total={total:.1f}")

                for sample, count in values.items:
                    all_gh_data.append({
                        "Dataset": dataset_name,
                        "Diet": diet,
                        "Sample": sample,
                        "CAZy_Family": family,
                        "Abundance": count,
                    })

        report_lines.append(f"  Total GH families found: {gh_count}")

        # Total CAZyme summary
        total_cazyme = sum(sum(v.values) for v in best_data.values)
        gh18_total = sum(best_data.get("GH18", {}).values)
        if total_cazyme > 0:
            ratio = gh18_total / total_cazyme * 100
            report_lines.append(f"\n  Chitinolytic Index: GH18/total = {gh18_total:.0f}/{total_cazyme:.0f} = {ratio:.2f}%")

    # Write combined chitin CAZyme table
    if all_chitin_data:
        outfile = os.path.join(OUTPUT_DIR, "01_cazyme_chitin_summary.tsv")
        with open(outfile, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=["Dataset", "Diet", "Sample", "CAZy_Family", "Description", "Abundance"],
                                    delimiter='\t')
            writer.writeheader
            writer.writerows(all_chitin_data)
        report_lines.append(f"\n\nOutput: {outfile} ({len(all_chitin_data)} rows)")
    else:
        report_lines.append("\n\nWARNING: No chitin CAZyme data extracted!")

    # Write all GH family table
    if all_gh_data:
        outfile = os.path.join(OUTPUT_DIR, "01_cazyme_all_GH_summary.tsv")
        with open(outfile, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=["Dataset", "Diet", "Sample", "CAZy_Family", "Abundance"],
                                    delimiter='\t')
            writer.writeheader
            writer.writerows(all_gh_data)
        report_lines.append(f"Output: {outfile} ({len(all_gh_data)} rows)")

    # Write report
    report_text = "\n".join(report_lines)
    report_file = os.path.join(OUTPUT_DIR, "01_chitin_cazyme_report.txt")
    with open(report_file, 'w') as f:
        f.write(report_text)

    print(report_text)
    print(f"\nReport saved: {report_file}")


if __name__ == "__main__":
    main
