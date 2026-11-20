# DIA-NN search database — `05.reanalysis_2026-11/db`

| Item | Value |
|---|---|
| File | `searchdb_Mimp_v3_meta_crap.faa` — 95,493 sequences; md5 `cc6d1592aca3addc4cb907e70b482ce9` |
| Built | 2026-11-10 07:52 UTC (`diann_search.sbatch` header documents the composition) |
| Components | (1) `host_Mimp_v3.faa` — all **Mimp_v3** predicted protein isoforms (29,927 sequences; GWHHIVU00000000.2 annotation, 24,875 genes incl. the eight chr4 GH18 loci MIMP3_g08469–g08479); (2) `meta_prefixed.faa` — 65,450 metagenome proteins from the unified catalogue restricted to CAZy/KO-annotated genes (`meta_subset_ids.txt`, IDs prefixed by sample); (3) `crap_prefixed.faa` — cRAP contaminants, 116 sequences (`crap.fasta` from `ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta`, downloaded 2026-11-10) |
| Used by | DIA-NN 1.8.1 library-free search of all 36 .raw files (DS04 18 runs MI01, DS05 18 runs MI02): `--fasta-search --predictor --met-excision --cut K*,R* --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --min-pr-charge 2 --max-pr-charge 4 --min-pr-mz 300 --max-pr-mz 1800 --min-fr-mz 200 --max-fr-mz 1800 --qvalue 0.01 --matrices` (full command in `diann_search.sbatch`) |
| Output | `05.reanalysis_2026-11/diann/report.*` → `tables/50–52` |
