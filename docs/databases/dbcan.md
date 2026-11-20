# dbCAN3 database — `$REFDIR/dbcan`

| Item | Value |
|---|---|
| Purpose | CAZyme annotation of the unified 32-sample protein catalogue and MAGs (chitin paper re-analysis, Nov 2026) |
| Downloaded | 2026-11-11 00:49–00:59 UTC |
| How | `conda activate $CONDA_ROOT/envs/dbcan; run_dbcan database --db_dir $REFDIR/dbcan` (run_dbcan **5.2.9**, bioconda `dbcan=5.2.9`; log `05.reanalysis_2026-11/logs/dbcan_db.log`) |
| Source | dbCAN3 server, `https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/<file>` = the **"db_current"** release served on the download date (the dbCAN team does not version-stamp the files; record the download date) |
| Files | `dbCAN.hmm` (875 CAZy family/subfamily HMMs, HMMER3 format) · `CAZy.dmnd` (2.18 GB DIAMOND index of CAZyDB) · `dbCAN-sub.hmm` (5.13 GB subfamily HMMs) · `fam-substrate-mapping.tsv` · `TF.hmm`/`TF.dmnd` · `TCDB.dmnd` · `PUL.dmnd`, `dbCAN-PUL/`, `dbCAN-PUL.xlsx` · `peptidase_db.dmnd` · `sulfatlas_db.dmnd` |
| Missing | `STP.hmm` — server returned HTTP 403 on every attempt. Only needed for the CGC/STP (signal-transduction) step, which we do not run; CAZyme annotation (`--methods diamond,hmm,dbCANsub`) is unaffected |
| Tools in the env | diamond 2.2.8, pyhmmer 0.12.3 (used instead of hmmer binaries) |

## Usage (as run)

```bash
run_dbcan CAZyme_annotation --input_raw_data proteins.faa --mode protein \
  --db_dir $REFDIR/dbcan --methods diamond,hmm,dbCANsub --threads 16 --output_dir out/
```
Job: `05.reanalysis_2026-11/dbcan_unified32.sbatch` (8 chunks of `15_unified32/catalog/nr_proteins.faa`).
Earlier project runs (tables 05/06, 2026-04) used dbCAN3 web/standalone v4 — the two are not identical; this directory is the reference for every Nov-2026 result.
