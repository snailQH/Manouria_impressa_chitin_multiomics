#!/bin/bash
# Rebuild the cross-cohort tables (35-43, 45-46) on the SINGLE 32-sample gene catalogue.
#  1. server: eggNOG CAZy / KO x Bowtie2 counts -> unified32_{cazy,ko}_abundance.tsv (extract_*_abundance.py, unchanged)
#  2. local : split the 32 columns by cohort into results/unified32/<cohort>_*.tsv, reuse the Bracken files
#             (taxonomy was not re-run) and the new mapping rates, then run unified_cross_cohort.py with
#             UNIFIED_DIR/TABLES_DIR -> tables/unified32/ for side-by-side comparison with tables/35-46.
set -euo pipefail
cd "$(dirname "$0")/../.."
# CAZY_SOURCE=eggnog (default): eggNOG CAZy column -> results/unified32, tables/unified32 (robustness set)
# CAZY_SOURCE=dbcan           : dbCAN3 Recommend Results (>= 2 tools) -> results/unified32_dbcan, tables/ (MANUSCRIPT set)
CAZY_SOURCE=${CAZY_SOURCE:-eggnog}
S=<user>@<slurm-head>; P=$PROJECT/03.unified_analysis; W=$P/15_unified32
scp -q scripts/extract_cazy_abundance.py scripts/extract_ko_abundance.py $S:$W/
ssh $S "cd $W && n=\$(ls counts/*.counts.tsv | wc -l) && [ \$n = 32 ] || { echo only \$n/32 count files; exit 1; }
  python3 extract_cazy_abundance.py --annotations kegg/unified32.emapper.annotations --counts-dir counts --out unified32_cazy_abundance.tsv
  python3 extract_ko_abundance.py --annotations kegg/unified32.emapper.annotations --counts-dir counts --out unified32_ko_abundance.tsv
  for f in counts/*.bowtie2.log; do s=\$(basename \$f .bowtie2.log); echo -e \"\$s\t\$(grep -m1 'reads; of these' \$f | cut -d' ' -f1)\t\$(grep overall \$f | cut -d' ' -f1)\"; done > unified32_mapping.tsv"
mkdir -p results/unified32 tables/unified32 results/unified32_dbcan
scp -q $S:$W/unified32_cazy_abundance.tsv $S:$W/unified32_ko_abundance.tsv $S:$W/unified32_mapping.tsv results/unified32/
[ -s results/unified32/unified32_cazy_abundance_dbcan.tsv ] || scp -q $S:$W/unified32_cazy_abundance_dbcan.tsv results/unified32/
CAZY_SOURCE=$CAZY_SOURCE python3 - <<'PY'
import pandas as pd, shutil, os
src=os.environ['CAZY_SOURCE']; R='results/unified32_dbcan' if src=='dbcan' else 'results/unified32'; U='results/unified'; E='results/unified32'
os.makedirs(R,exist_ok=True)
coh=lambda s: 'DS01' if s.startswith('ITF') else 'DS02' if s.startswith('HM') else 'DS03'
for kind,idx in [('cazy','Family'),('ko','KO')]:
    f=f'{E}/unified32_cazy_abundance_dbcan.tsv' if (kind=='cazy' and src=='dbcan') else f'{E}/unified32_{kind}_abundance.tsv'
    t=pd.read_csv(f,sep='\t',comment='#')
    samples=[c for c in t.columns if c not in (idx,'Total')]
    for c in ['DS01','DS02','DS03']:
        cols=[s for s in samples if coh(s)==c]; sub=t[[idx]+cols].copy
        if kind=='cazy': sub['Total']=sub[cols].sum(axis=1)
        sub.to_csv(f'{R}/{c}_{kind}_abundance.tsv',sep='\t',index=False)
m=pd.read_csv(f'{E}/unified32_mapping.tsv',sep='\t',header=None,names=['Sample','pairs','rate'])
m.insert(0,'Cohort',m.Sample.map(coh)); m.to_csv(f'{R}/mapping_summary.tsv',sep='\t',index=False,header=False)
for f in os.listdir(U):
    if 'bracken' in f or f=='qc_assembly_32.tsv': shutil.copy(f'{U}/{f}',f'{R}/{f}')
print('inputs ready:', sorted(os.listdir(R)))
PY
if [ "$CAZY_SOURCE" = dbcan ]; then
  UNIFIED_DIR=$PWD/results/unified32_dbcan TABLES_DIR=$PWD/tables python3 scripts/unified_cross_cohort.py
  python3 scripts/metagenome_animal_effect.py
else
  UNIFIED_DIR=$PWD/results/unified32 TABLES_DIR=$PWD/tables/unified32 python3 scripts/unified_cross_cohort.py
fi
echo; echo "=== GH18 CPM by cohort: tables/ (manuscript set) vs tables/unified32/ (eggNOG robustness set)"
python3 - <<'PY'
import pandas as pd
for lab,d in [('tables (manuscript)','tables'),('unified32 eggNOG','tables/unified32')]:
    t=pd.read_csv(f'{d}/35_unified_cazy_per_sample.tsv',sep='\t')
    g=t.groupby('Cohort')[['GH18_CPM','GH19_CPM','GH20_CPM']].mean.round(1)
    print(lab); print(g.to_string); print
    s=pd.read_csv(f'{d}/36_unified_cazy_tests.tsv',sep='\t'); print(s[s.Family=='GH18'][['A','B','fold_A_over_B','P']].round(4).to_string(index=False)); print
PY
