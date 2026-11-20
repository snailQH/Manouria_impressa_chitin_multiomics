# Kraken2 database — `$REFDIR/kraken2_db` (and the copy `kraken2_db_local`)

| Item | Value |
|---|---|
| Type | Locally built Kraken2 **full-size** database (no `--max-db-size` cap) + Bracken 150-mer distribution |
| Build date | 2025-10-21 → 2025-10-22 (file timestamps; built by uid 1002, not the chitin project) |
| Built with | `kraken2-build --download-taxonomy`, `--download-library <lib>` for each library below, `--build`; `bracken-build -l 150` |
| Libraries | archaea, bacteria, viral, human (GRCh38), UniVec_Core (= Kraken2 "standard") **+ protozoa + fungi** → equivalent in scope to the Kraken2 **PlusPF** collection |
| Source | NCBI RefSeq via `kraken2-build --download-library` (RefSeq **Complete Genome / Chromosome** assemblies only — the Kraken2 default; scaffold/contig-level assemblies are *not* included) |
| Genomes in the build | archaea 679 · bacteria 59,927 (141,909 sequences, 243 GB fasta) · fungi **133** (67 chromosome + 66 complete; 69 genera) · protozoa 53 · viral 15,001 · human 1 (`library/*/manifest.txt`) |
| Taxonomy | NCBI taxdump + `nucl_gb` / `nucl_wgs.accession2taxid` downloaded 2025-10-21 (`taxonomy/`) |
| Size | `hash.k2d` 107,186,287,208 bytes (100 GB); `taxo.k2d` 4.4 MB; `seqid2taxid.map` 7.5 MB; `database150mers.kmer_distrib` (Bracken, read length 150) |
| `opts.k2d` md5 | `88858ab7fa79a869f6045ca0fc9a3b7c` |

## `kraken2_db_local`

`$PROJECT/kraken2_db_local` (copied 2026-09-02) is a **byte-identical copy of the built index**
(`opts.k2d` md5 identical, `taxo.k2d` identical, `hash.k2d` prefix identical) without the `library/` and `taxonomy/`
source folders. It is the database used for all published Kraken2/Bracken results of the chitin project
(`03.unified_analysis/08_taxonomy_conf01`, `run_taxonomy.sh`: `--confidence 0.2 --minimum-hit-groups 2`, Bracken `-r 150 -l S`).

## What this means for fungal profiling

* Fungi **are** in the database, so the fungal fractions in the published reports (e.g. 0.03 % of ITF-1-1 reads) are real classifications, not a database omission.
* But only RefSeq complete/chromosome fungal genomes are present: **no *Pleurotus*, *Hypsizygus* or any Agaricales** (the dietary mushrooms), and **no Neocallimastigomycota** (anaerobic gut fungi). Reads from those taxa are unclassified or forced to the nearest Basidiomycota/Fungi node. Dietary and gut-fungal signals therefore need the extended database `kraken2_pluspf_fungiext_2026-11` (see its README) and/or direct read mapping to the mushroom genomes (`05.reanalysis_2026-11/diet_genomes`).

## How to rebuild

```bash
conda activate $CONDA_ROOT/envs/taxonomy          # kraken2 2.1.x, bracken
kraken2-build --download-taxonomy --db DB
for l in archaea bacteria viral human UniVec_Core protozoa fungi; do kraken2-build --download-library $l --db DB; done
kraken2-build --build --db DB --threads 60              # run on a 480 GB Slurm node
bracken-build -d DB -t 60 -k 35 -l 150
```
