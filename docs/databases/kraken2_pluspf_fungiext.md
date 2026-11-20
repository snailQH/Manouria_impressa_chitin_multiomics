# Extended Kraken2 database — `$REFDIR/kraken2_pluspf_fungiext_2026-11`

| Item | Value |
|---|---|
| Purpose | Fungal profiling of the 32 faecal metagenomes: dietary mushrooms and anaerobic gut fungi are absent from the standard build (`kraken2_db`, RefSeq complete genomes only) |
| Built | started 2026-11-11 (`05.reanalysis_2026-11/kraken2_fungiext_build.sh` → `kraken2_fungiext_build.sbatch`; logs `logs/kraken2_fungiext_prep.log`, `logs/kraken2_fungiext_build_<job>.out`) |
| Base libraries | file-symlinks to `$REFDIR/kraken2_db/library/{archaea,bacteria,viral,human,UniVec_Core,protozoa,fungi}` and its `taxonomy/` (NCBI taxdump 2025-10-21) — i.e. the same sequence set as the published database |
| Added (`library/added`, source list `new_genomes/to_add.tsv`) | (1) every RefSeq fungal assembly not already present — any assembly level, `version_status=latest`, `genome_rep=Full` (`refseq_new.tsv`, from `https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt` on the build date); (2) all GenBank Neocallimastigomycota assemblies (20 genera incl. *Piromyces*, *Anaeromyces*, *Neocallimastix*, *Pecoramyces*, *Caecomyces*, *Orpinomyces*; `gutfungi.tsv`, from `genomes/genbank/fungi/assembly_summary.txt`); (3) the three dietary genomes GCA_001717165.1 *P. eryngii* (taxid 5323), GCA_013433165.1 *H. marmoreus* (39966), GCA_014466165.1 *P. ostreatus* (5322) |
| Header scheme | added sequences carry `|kraken:taxid|<taxid>` (taxid from assembly_summary column 6), so no accession2taxid lookup is required; added without masking (`--no-masking`) |
| Build | `kraken2-build --build --threads 60` (k = 35, l = 31 defaults, no size cap) then `bracken-build -k 35 -l 150` on a 480 GB node |
| Tool | Kraken2 / Bracken from `$CONDA_ROOT/envs/taxonomy` |

Runs on the 32 non-host read sets use the same parameters as the published run (`--confidence 0.2 --minimum-hit-groups 2`, Bracken `-r 150`), so bacterial/archaeal results are directly comparable and only the fungal fraction changes.

**Build history (2026-11):** build 1 — hash write truncated on NFS (1.2 of 136 GB); build 2 — hash correct but the 594 added genomes were silently skipped (the added library's prelim_map had been deleted; seqid2taxid.map = base DB only, 192,551 entries); build 3 (local disk on a high-memory compute node, 30 h) — hash verified byte-for-byte (23.8 G minimizers, 135,850,032,304 bytes), seqid2taxid.map 570,994 entries, but the Bracken step failed because seqid2taxid.map had not been copied to NFS; finished by a later job (`kraken2_fungiext_bracken.sbatch`: map regenerated from library/*/prelim_map.txt, Bracken k 35 / read length 150 on the existing database.kraken). Profiling a job guards against a map < 500 k entries.
