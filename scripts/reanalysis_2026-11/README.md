# Re-analyses started 2026-11-10

| job | script | purpose | output (server) |
|---|---|---|---|
| 94802 (+ array) | `unified32_catalogue.sbatch` → `unified32_quantify.sbatch` | ONE non-redundant gene catalogue for all 32 samples (CD-HIT 95 %/90 %, same as per-cohort), Bowtie2 counts per sample, eggNOG-mapper annotation | `03.unified_analysis/15_unified32/{catalog,idx,counts,kegg}` |
| 94788 (array 36) → 94789 | `convert_raw.sbatch` → `diann_search.sbatch` | Thermo .raw → mzML (ThermoRawFileParser 2.0; DIA-NN 1.8.1 Linux cannot read .raw) → DIA-NN library-free re-search of all 36 runs against Mimp_v3 proteome (29,927) + metagenome CAZy/ko00520 protein subset (65,450, from the three current catalogues) + cRAP (116); 1 % FDR, MBR, Met-ox variable | `05.reanalysis_2026-11/{mzml,diann,db}` |

| 94830 (array 32) | `unified32_quantify.sbatch` | per-sample Bowtie2 counts against the unified catalogue (19/32 done 2026-11-11 01:00 UTC; DS01 samples ~4–5 h each) | `15_unified32/counts/*.counts.tsv` |
| 94856 (array 8) | `dbcan_unified32.sbatch` | dbCAN3 (run_dbcan 5.2.9: HMMER dbCAN-HMM + DIAMOND CAZy + dbCAN-sub) on the 2.85 M unified NR proteins in 8 chunks; DB `$REFDIR/dbcan` (downloaded 2026-11-11, `db_current`; STP.hmm 403 — not needed) | `05.reanalysis_2026-11/dbcan/out/part_*` |
| (array 32, `mags/`) | `mag_bin_sample.sbatch` ← `mag_setup_and_launch.sh` | per-sample MAG binning: Bowtie2 reads → own MEGAHIT assembly, `jgi_summarize_bam_contig_depths`, MetaBAT2 **2.15** (≥1,500 bp; env `metabat215` — the 2.18 bioconda build needs AVX2 and crashes with *Illegal instruction* on the Sandy Bridge compute nodes); CheckM2 1.1.0 / dRep 4.0.3 / GTDB-Tk 2.4.0 in env `binning`; CheckM2 DB + GTDB-Tk R226 downloaded to `$REFDIR/{checkm2,gtdbtk}` | `05.reanalysis_2026-11/mags/<S>/bins/` |
| 94861 (array 32) | `diet_mapping.sbatch` | dietary-DNA: non-host pairs → combined index of *P. eryngii* 183, *H. marmoreus* HM62, *P. ostreatus* PC9 (GenBank; `diet_genomes/`); properly paired MAPQ ≥ 10 pairs per genome, per million non-host pairs | `05.reanalysis_2026-11/diet_mapping/counts/<S>.tsv` |
| (head node → sbatch) | `kraken2_fungiext_build.sh` → `kraken2_fungiext_build.sbatch` | extended Kraken2 DB `$REFDIR/kraken2_pluspf_fungiext_2026-11`: the published DB's libraries (RefSeq complete genomes; fungi = 133 genomes, no Agaricales, no Neocallimastigomycota) + 541 further RefSeq fungal assemblies + 50 GenBank Neocallimastigomycota + the 3 diet genomes; then Bracken 150-mer. Fungal profiling of the 32 samples follows with the published parameters (conf 0.2, min-hit-groups 2) | `05.reanalysis_2026-11/logs/kraken2_fungiext_*` |
| — | `comparative/download.sh` | 8 RefSeq turtle genomes (+proteins, GFF) for the GH18-locus comparison; see `docs/databases/comparative_turtles.md` | `05.reanalysis_2026-11/comparative/turtles/` |

Every database folder carries a `README_DATABASE.md` (source URL, download date, version, size); the repo copies are in `docs/databases/`.

Search database: `05.reanalysis_2026-11/db/searchdb_Mimp_v3_meta_crap.faa` (95,493 sequences; headers prefixed host|, meta|, crap|; md5 cc6d1592aca3addc4cb907e70b482ce9). The metagenome subset = catalogue proteins with any eggNOG CAZy family or a ko00520 chitin/GlcNAc KO (K01183 chitinase, K01207/K12373 β-hexosaminidase, K01452/K03478 chitin deacetylase, K00884 GlcNAc kinase, K00844 hexokinase, K01443 NagA, K02564 NagB, K02563/K00820 GlmS, K13748, K18676, K22298, K00847).

Follow-ups once finished: regenerate tables 35–44 and Figures 2–4 from the unified counts (new script); peptide → Mimp_v3 gene table from the DIA-NN report, GH18 locus quantification per run (replaces the proxy-database tables 02/18/19/31/47/48).

