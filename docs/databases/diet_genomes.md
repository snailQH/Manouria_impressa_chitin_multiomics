# Dietary mushroom genomes — `05.reanalysis_2026-11/diet_genomes`

Downloaded 2026-11-11 from NCBI GenBank (`comparative/download.sh`), `_genomic.fna.gz` only.
Purpose: quantify dietary fungal DNA in the 32 non-host faecal read sets (Bowtie2 mapping) and extend the Kraken2 fungal library.

| File | Species (diet cohort) | GenBank assembly |
|---|---|---|
| `Pleurotus_eryngii_183.fna.gz` | *Pleurotus eryngii* (DS01, <provider-delivery-dir> | GCA_001717165.1 ASM171716v1 (strain 183) |
| `Hypsizygus_marmoreus_HM62.fna.gz` | *Hypsizygus marmoreus* (DS02, <provider-delivery-dir> manuscript currently "H. tessellatus") | GCA_013433165.1 ASM1343316v1 (HM62) |
| `Pleurotus_ostreatus_PC9.fna.gz` | *Pleurotus ostreatus* (DS03, <provider-delivery-dir> | GCA_014466165.1 ASM1446616v1 (PC9; RefSeq twin GCF_014466165.1) |

None of these taxa is in the standard Kraken2 fungi library (RefSeq complete/chromosome genomes only), which is why direct mapping is used.
