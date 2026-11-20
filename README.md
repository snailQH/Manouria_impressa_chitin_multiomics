# Chitin metabolism in the impressed tortoise (*Manouria impressa*) — analysis code and result tables

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24389230.svg)](https://doi.org/10.5281/zenodo.24389230)

Code and intermediate results for **"Multi-omics study indicates cooperative host–microbiome chitin degradation
in the mycophagous tortoise *Manouria impressa*"** (Wang, Yang & Li).

Every figure in the manuscript has exactly one generating script here, and every number it reports comes from a table
in `tables/`. The tables were produced from the raw data by the job scripts in `scripts/reanalysis_2026-11/`, which ran
on a compute cluster against data too large to publish — so those tables are shipped as the starting point, and the
figure scripts run from them on a laptop. A few tables (the evidence tiers, the vendor-pipeline comparison) are
hand-assembled summaries of other tables and say so in their own header.

## Data (not in this repository)

| Data | Accession | Link |
|---|---|---|
| 32 faecal metagenomes | NGDC GSA **CRA054275** | https://ngdc.cncb.ac.cn/gsa/browse/CRA054275 |
| ↳ BioProject | **PRJCA080007** | https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA080007 |
| ↳ BioSamples | SAMC9104754–SAMC9104785 | https://ngdc.cncb.ac.cn/search/?dbId=biosample&q=SAMC9104754 |
| ↳ Runs | CRR3904968–CRR3904999 | https://ngdc.cncb.ac.cn/gsa/browse/CRA054275/CRR3904968 |
| DIA proteomics (36 runs) | PRIDE/ProteomeXchange **PXD087475** | https://www.ebi.ac.uk/pride/archive/projects/PXD087475 — doi:[10.6019/PXD087475](https://doi.org/10.6019/PXD087475) |
| Host genome Mimp_v3 | NGDC GWH **GWHHIVU00000000.2** | https://ngdc.cncb.ac.cn/gwh/Assembly/GWHHIVU00000000.2 |
| ↳ Haplotype 1 | GWHKKMH00000000.1 | https://ngdc.cncb.ac.cn/gwh/Assembly/GWHKKMH00000000.1 |
| ↳ Haplotype 2 | GWHKKMI00000000.1 | https://ngdc.cncb.ac.cn/gwh/Assembly/GWHKKMI00000000.1 |
| ↳ BioProject | PRJCA038476 | https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA038476 |
| ↳ Mitogenome | GenBase **C_AA510199.1** | https://ngdc.cncb.ac.cn/genbase/browse/C_AA510199.1 |
| Genome assembly pipeline | — | https://github.com/snailQH/Manouria_impressa_genome_v3 — doi:[10.5281/zenodo.24056325](https://doi.org/10.5281/zenodo.24056325) |

Every link above was checked and resolves.

## Layout

```
scripts/                     one script per figure/table group; the header of each says what it writes
  reanalysis_2026-11/        the unified 32-sample re-analysis (Slurm): catalogue, dbCAN3, MAGs, DIA-NN, fungal profiling
tables/                      every result table cited in the manuscript (TSV)
  unified32/                 the eggNOG-annotated robustness set — NOT the manuscript's numbers (see its README)
  archive_2026-11-11…/       superseded per-cohort tables, kept for provenance; every file says so in its header
results/unified/             assembly and QC summaries read directly by the figure scripts
results/reanalysis_2026-11/  small intermediates (MAG taxonomy, comparative loci, Mimp_v3 gene map)
results/phylogeny/v2/        vertebrate GH18 alignment and ML tree (Additional file 3 of the paper)
docs/databases/              which reference database was used, its version, and where it came from
```

## Reproducing a figure

Each script states its inputs and outputs in its docstring:

```bash
python3 scripts/generate_unified_figs.py      # Figures 2A-D, 3A/3B/3D, 4A, Supp Figs S7/S19
python3 scripts/mag_gh18_analysis.py          # Figure 3C, Supp Fig S8, tables 64-66
python3 scripts/generate_plots.py             # Figures 4B, 5A, 9A, Supp Figs S5/S6/S9/S10
python3 scripts/generate_fig5b_phylogeny.py   # Figure 5B (condensed) and Supp Fig S11 (full 162-tip tree)
python3 scripts/generate_suppfig_S3_comparative_locus.py   # Figure 6, the locus across nine Testudines
python3 scripts/proteomics_diann_chitin.py    # Figures 7A/7B and 8, Supp Figs S13-S17, tables 58-62 and 73
python3 scripts/temporal_trajectories.py      # Supp Fig S4, table 69
```

Of the commands above, all but one run from what this repository ships. `proteomics_diann_chitin.py`
additionally needs the DIA-NN search output, which is deposited as PRIDE **PXD087475** rather than
duplicated here. Put the report matrices in `results/reanalysis_2026-11/diann/` and the script runs;
the tables it produces (`tables/58`–`62`, `73`) are shipped, so its results are available either way.

Nine scripts are the upstream pipeline steps rather than figure code — `unified_cross_cohort.py`,
`extract_cazy_abundance*.py`, `extract_ko_abundance.py`, `summarise_emapper.py`,
`compare_harmonised_cohorts.py`, `dietary_dna_table.py`, `proteomics_mimpv3_reanalysis.py` — and
they read the gene catalogue and annotation output, which are too large to publish. They are here as
the record of how the shipped tables were made, not as something a reader can run unchanged.

Anything that starts from the reads themselves — the catalogue, the annotation, the MAGs — ran as
batch jobs on a compute cluster (`scripts/reanalysis_2026-11/*.sbatch`, kept as the record of what
was run) against the metagenomes in GSA **CRA054275**. Paths in those job scripts are placeholders
(`$PROJECT`, `$REFDIR`, `$CONDA_ROOT`, `$TOOLS`, `$NFS`) — substitute your own.

Not every table here is a current result. Vendor-era and superseded tables are kept for provenance
and each one opens with a header saying what replaced it; `tables/07` and `tables/21` carry results
that were **withdrawn**, and say so. Cite `tables/35`–`46` for every cross-diet number.

## Checks

Three scripts exist to catch the kinds of error that do not raise an exception:

```bash
python3 scripts/figure_integrity_check.py scripts/generate_*.py   # data drawn outside the visible axes
python3 scripts/table_integrity_check.py                          # ragged tables, failed joins, derived-column drift
python3 scripts/verify_permanova.py                               # recomputes the published PERMANOVA independently
```

`scripts/figure_ownership.py` records which script owns which figure file, so two scripts cannot
silently write the same figure and leave the result to run order.

## Citing

Please cite the paper. If you use the code or tables directly, the archived release is:

> Wang J, Yang D, Li Q. Analysis code and result tables for chitin metabolism in *Manouria impressa*. Zenodo. doi:10.5281/zenodo.24389230 (concept DOI, resolves to the latest release)

## Licence

**Code** (`scripts/`): **PolyForm Noncommercial License 1.0.0** — free to use, modify and share for any
noncommercial purpose, including academic research and teaching. Commercial use needs a separate licence
from the authors; please contact the corresponding author.

**Result tables and data** (`tables/`, `results/`): **CC BY 4.0**, unchanged — the numbers behind the
published figures stay freely reusable, with attribution.

Full terms in [`LICENSE`](LICENSE).
