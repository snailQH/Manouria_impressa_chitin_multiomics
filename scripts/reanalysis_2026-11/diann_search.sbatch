#!/bin/bash
#SBATCH --job-name=diann_v3
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=96:00:00
#SBATCH --output=$PROJECT/05.reanalysis_2026-11/logs/diann_%j.out
# DIA-NN 1.8.1 library-free re-search of all 36 Orbitrap Astral runs against the Mimp_v3 proteome (29,927) + metagenome CAZyme/ko00520
# subset (65,450) + cRAP (116). Settings follow the earlier run_diann.sh (vendor-like): trypsin K*,R*, 1 missed cleavage, 7-30 aa,
# charge 2-4, carbamidomethyl fixed, Met oxidation variable, 1% FDR, MBR, matrices.
set -uo pipefail
R=$PROJECT/05.reanalysis_2026-11; D=$TOOLS/diann-1.8.1; export LD_LIBRARY_PATH=$D:${LD_LIBRARY_PATH:-} LC_ALL=C LANG=C
n=$(ls $R/mzml/*.mzML | wc -l); [ "$n" -eq 36 ] || { echo "only $n mzML files — abort"; exit 1; }
F=""; for f in $R/mzml/*.mzML; do F="$F --f $f"; done
###### DIA-NN START $(date) $(hostname) ######"
$D/diann-1.8.1 $F --lib "" --threads 64 --verbose 1 --out $R/diann/report.tsv --out-lib $R/diann/report-lib.tsv --qvalue 0.01 --matrices \
  --temp $R/diann/temp --fasta $R/db/searchdb_Mimp_v3_meta_crap.faa --fasta-search --predictor --met-excision --cut K*,R* --missed-cleavages 1 \
  --min-pep-len 7 --max-pep-len 30 --min-pr-charge 2 --max-pr-charge 4 --min-pr-mz 300 --max-pr-mz 1800 --min-fr-mz 200 --max-fr-mz 1800 \
  --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --reanalyse --smart-profiling --pg-level 2 --mass-acc 15 --mass-acc-ms1 15 --relaxed-prot-inf
###### DIA-NN DONE $(date) ######"; ls -la $R/diann | head
