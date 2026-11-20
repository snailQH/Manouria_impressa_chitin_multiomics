# Pfam-A — `$REFDIR/pfam/`

| | |
|---|---|
| Database | Pfam-A protein family HMMs (EMBL-EBI / InterPro) |
| Release | **39.0** (2026-09; based on UniProtKB 2026_02; 53,275 families) — `Pfam.version`, `relnotes.txt` |
| Source | https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz (658,389,433 bytes gz) |
| Downloaded | 2026-11-15, 24 parallel byte ranges (`scripts/reanalysis_2026-11/pfam_parallel_download.sh`; size verified against the server), then `hmmpress` (HMMER 3.4, env `binning`) |
| Files | `Pfam-A.hmm` + `.h3{f,i,m,p}` |
| Used for | domain architecture of the Mimp_v3 GH18-locus models and human GH18 paralogs (`scripts/reanalysis_2026-11/domain_architecture_hmmscan.sh` → `tables/70`, Supp Fig S18): identified MIMP3_g08468 as TAF8 |
| Citation | Mistry et al. (2021) *Nucleic Acids Research* 49:D412–D419, doi:10.1093/nar/gkaa913 |

Usage: `hmmscan --cut_ga --domtblout out.domtbl $REFDIR/pfam/Pfam-A.hmm query.faa`
