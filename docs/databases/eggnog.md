# eggNOG-mapper database — `$REFDIR/eggnog`

| Item | Value |
|---|---|
| Database version | **eggNOG 5.0.2** (`eggnog.db` table `version` = 5.0.2; files dated 2021-03-02 on the eggNOG server) |
| Mapper version | emapper **2.1.12** (`$CONDA_ROOT/envs/eggnog`), DIAMOND 2.2.5, MMseqs2 18.8cc5c |
| Downloaded | 2026-09-02 (folder mtime) with `download_eggnog_data.py --data_dir $REFDIR/eggnog -y` |
| Source | `http://eggnog5.embl.de/download/emapperdb-5.0.2/` (`eggnog.db.gz`, `eggnog_proteins.dmnd.gz`, `eggnog.taxa.tar.gz`) |
| Files | `eggnog.db` 41.4 GB (annotation SQLite) · `eggnog_proteins.dmnd` 9.3 GB (DIAMOND index of eggNOG 5 proteins) · `eggnog.taxa.db` + `.traverse.pkl` |
| Used by | all unified-pipeline functional annotation (`03.unified_analysis/*/kegg/*.emapper.annotations`, incl. `15_unified32/kegg/unified32.emapper.annotations`): `emapper.py -i nr_proteins.faa -m diamond --dbmem --override --data_dir $REFDIR/eggnog` (emapper 2.1.12 defaults: e-value 0.001, score 60, `--tax_scope auto`, `--target_orthologs all`, `--go_evidence non-electronic`) |

KO, CAZy, COG and pathway assignments in the chitin manuscript derive from this release (KEGG mapping as frozen in eggNOG 5.0.2).
