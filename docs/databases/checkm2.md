# CheckM2 database — `$REFDIR/checkm2`

| Item | Value |
|---|---|
| Purpose | Completeness/contamination of MetaBAT2 bins (MAG re-analysis, Nov 2026) |
| Downloaded | 2026-11-11 with `checkm2 database --download --path $REFDIR/checkm2` (CheckM2 from bioconda, env `$CONDA_ROOT/envs/binning`; log `05.reanalysis_2026-11/logs/checkm2_db.log`) |
| Source | Zenodo record 5571251 (CheckM2 DIAMOND database, `uniref100.KO.1.dmnd`, DB version **v2** — 2.7 GB) via the CheckM2 downloader |
| Files | `CheckM2_database/uniref100.KO.1.dmnd` (+ `.md5`, `version`) |
| Usage | `checkm2 predict --threads 32 --input bins/ -x fa --output-directory checkm2_out --database_path $REFDIR/checkm2/CheckM2_database/uniref100.KO.1.dmnd` |

Confirm the exact version string in `CheckM2_database/version` (or `checkm2 database --current`) after download and record it here.
