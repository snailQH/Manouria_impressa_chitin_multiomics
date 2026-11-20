# GTDB-Tk reference data — `$REFDIR/gtdbtk`

| Item | Value |
|---|---|
| Purpose | Taxonomic classification of dereplicated MAGs (MAG re-analysis, Nov 2026) |
| Release | **GTDB R226** (April 2025), GTDB-Tk package `gtdbtk_r226_data.tar.gz` (~110 GB extracted → `release226/`) |
| Downloaded | 2026-11-11 by `mags/mag_setup_and_launch.sh`; log `05.reanalysis_2026-11/logs/gtdbtk_db.log` |
| Source | `https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz` (mirror: `https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/`) |
| Tool | GTDB-Tk **2.4.0** (bioconda, env `$CONDA_ROOT/envs/binning`) — R226 requires GTDB-Tk ≥ 2.4.0 |
| Usage | `export GTDBTK_DATA_PATH=$REFDIR/gtdbtk/release226; gtdbtk classify_wf --genome_dir mags/ -x fa --out_dir gtdbtk_out --cpus 32 --skip_ani_screen` (or `--mash_db` for ANI screen) |

**Env note (2026-11-11):** GTDB-Tk 2.4.0 in `$CONDA_ROOT/envs/binning` ships with numpy 2.x, where `ndarray.tostring` no longer exists; `gtdbtk/markers.py` was patched in place (`.tostring` → `.tobytes`). Rebuilding the env needs either this patch or `numpy<2`. CheckM2 1.1.0 in the same env requires AVX (TensorFlow) — it needs a host whose CPU provides AVX.
