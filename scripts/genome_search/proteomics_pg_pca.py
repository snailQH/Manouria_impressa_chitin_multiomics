#!/usr/bin/env python3
"""Protein-group PCA per proteomics batch from the DIA-NN report (PG.MaxLFQ, log2, protein groups
quantified in every run of the batch). Runs are labelled by tissue from the vendor's sample order
(DS04) or sample_relation.txt (DS05). Writes tables/47_proteomics_pg_pca.tsv. Run on the server."""
import sys, numpy as np, pandas as pd
out = sys.argv[1]; rows = []
DS04 = ['Stomach_3','Stomach_2','Stomach_1','Duodenum_3','Duodenum_2','Duodenum_1','Jejunum_3','Jejunum_2','Jejunum_1','Ileum_3','Ileum_2','Ileum_1','Colon_3','Colon_2','Colon_1','Rectum_3','Rectum_2','Rectum_1']
for batch, f, rel in (('DS04', sys.argv[2], None), ('DS05', sys.argv[3], sys.argv[4])):
    d = pd.read_csv(f, sep='\t', usecols=['Run', 'Protein.Group', 'PG.MaxLFQ'])
    d['Run'] = d.Run.str.replace(r'.*/', '', regex=True).str.replace('.mzML', '', regex=False)
    m = d.drop_duplicates(['Run', 'Protein.Group']).pivot(index='Protein.Group', columns='Run', values='PG.MaxLFQ')
    m = m.dropna; m = m[(m > 0).all(axis=1)]; X = np.log2(m.values); X = X - X.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(X.T, full_matrices=False); pcs = U[:, :3] * S[:3]; var = S ** 2 / (S ** 2).sum
    runs = list(m.columns)
    if rel: tissue = dict(l.split for l in open(rel) if l.strip)
    else: tissue = dict(zip(sorted(runs), DS04))
    for r, p in zip(runs, pcs): rows.append(dict(Batch=batch, Run=r, Tissue=tissue[r], PC1=p[0], PC2=p[1], PC3=p[2], var1=var[0], var2=var[1], var3=var[2], n_protein_groups=len(m)))
    print(batch, len(runs), 'runs', len(m), 'complete protein groups; var', var[:3].round(3))
pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
