#!/usr/bin/env python3
"""Report, for each M. impressa GH18 locus in the v2 tree, its sister/clade context and the
support of the nodes joining it to the mammalian reference paralogs. Rooted on CTBS/CHID1."""
import re
from Bio import Phylo
t = Phylo.read('results/phylogeny/v2/gh18_v2.contree', 'newick')
og = [x for x in t.get_terminals if re.search(r'CTBS|CHID1', x.name)]
t.root_with_outgroup(*og)
terms = t.get_terminals
def find(s):
    r = [x for x in terms if s in x.name]; return r[0] if r else None
refs = {'CHIA': find('CHIA_Homo'), 'CHIT1': find('CHIT1_Homo'), 'CHI3L1': find('CHI3L1_Homo'), 'CHI3L2': find('CHI3L2_Homo')}
mimp = [x for x in terms if 'Manouria' in x.name]
def sz(c): return len(c.get_terminals)
for m in sorted(mimp, key=lambda x: x.name):
    print(f"\n=== {m.name} ===")
    path = t.get_path(m)[:-1]
    for anc in reversed(path):
        if sz(anc) > 30: break
        others = sorted(x.name.split('_Manouria')[0][:40] if 'Manouria' in x.name else '_'.join(x.name.split('_')[:3])[:40] for x in anc.get_terminals)
        print(f"   node n={sz(anc):2d} support={anc.confidence}: {others[:8]}{' ...' if len(others)>8 else ''}")
    for k, r in refs.items:
        ca = t.common_ancestor([m, r]); print(f"   MRCA with human {k}: {sz(ca)} tips, support {ca.confidence}")
print("\n=== reference clades ===")
for k, keys in {'CHIA': ['CHIA_Homo', 'Chia_Mus', 'Chia_Rattus', 'CHIA_Bos'], 'CHIT1': ['CHIT1_Homo', 'Chit1_Mus'], 'CHI3L': ['CHI3L1_Homo', 'CHI3L2_Homo']}.items:
    tips = [find(x) for x in keys]; ca = t.common_ancestor(tips)
    print(f"{k}: mammalian core n={sz(ca)} support={ca.confidence}")
    path = t.get_path(ca)
    node = ca
    # walk outward listing additions
    prev = set(x.name for x in ca.get_terminals)
    for anc in reversed(t.get_path(ca)[:-1]):
        cur = set(x.name for x in anc.get_terminals)
        if len(cur) > 45: break
        new = sorted(n.split('_')[0] + '_' + '_'.join(n.split('_')[1:3]) for n in cur - prev)
        print(f"    +{len(cur - prev):2d} at support {anc.confidence} (n={len(cur)}): {[n[:38] for n in new][:6]}")
        prev = cur
