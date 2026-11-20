# Reference databases used by the project

Each database folder on the server carries a copy of its file here as `README_DATABASE.md`
(what it is, where and how it was downloaded, version, size, how it is used). Keep both in sync.

| Server folder | README | Version |
|---|---|---|
| `$REFDIR/kraken2_db` and `…/Chitin_Metabolism_Analysis/kraken2_db_local` | [kraken2_db.md](kraken2_db.md) | RefSeq complete genomes, built 2025-10-21 (standard + protozoa + fungi, 100 GB) |
| `$REFDIR/kraken2_pluspf_fungiext_2026-11` | [kraken2_pluspf_fungiext.md](kraken2_pluspf_fungiext.md) | base libraries + all RefSeq fungi + Neocallimastigomycota + 3 diet genomes (build started 2026-11-11) |
| `$REFDIR/eggnog` | [eggnog.md](eggnog.md) | eggNOG 5.0.2 / emapper 2.1.12 |
| `$REFDIR/dbcan` | [dbcan.md](dbcan.md) | dbCAN3 `db_current` downloaded 2026-11-11, run_dbcan 5.2.9 |
| `$REFDIR/checkm2` | [checkm2.md](checkm2.md) | CheckM2 DIAMOND DB (downloaded 2026-11-11) |
| `$REFDIR/gtdbtk` | [gtdbtk.md](gtdbtk.md) | GTDB R226, GTDB-Tk 2.4.0 |
| `05.reanalysis_2026-11/db` | [diann_searchdb.md](diann_searchdb.md) | Mimp_v3 proteome + 65,450 metagenome proteins + cRAP (95,493 seqs) |
| `05.reanalysis_2026-11/comparative/turtles` | [comparative_turtles.md](comparative_turtles.md) | 8 RefSeq turtle genomes (2026-11-11) |
| `05.reanalysis_2026-11/diet_genomes` | [diet_genomes.md](diet_genomes.md) | 3 GenBank mushroom genomes (2026-11-11) |
