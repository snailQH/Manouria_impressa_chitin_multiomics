#!/usr/bin/env python3
"""
Independently re-derive the PERMANOVA the manuscript publishes (Results 3.2).

WHAT THIS CHECKED BEFORE, AND WHY THAT WAS WRONG (2026-11-16 audit)
  It recomputed the statistic from tables/11_genus_abundance_cross_diet.tsv — a VENDOR table — and
  compared it with results/19_permanova_results.txt (pseudo-F 6.0957, R^2 0.2530), then printed
  "reproduction of published values: OK". Neither is what the manuscript reports. The published value
  is pseudo-F 4.61, R^2 0.241 on the unified 32-sample genus profiles (tables/38 -> tables/49), so the
  script was giving confident assurance about a superseded number computed from a forbidden source.

It now recomputes from tables/38_unified_genus.tsv with an independent implementation and compares
against tables/49_metagenome_animal_effect.tsv, which is what the manuscript cites.

adonis2-equivalent one-way PERMANOVA:
  SS_total  = sum of squared distances / n
  SS_within = sum over groups of (within-group squared distances / n_g)
  pseudo-F  = ((SS_total - SS_within)/(a-1)) / (SS_within/(N-a))
  R^2       = (SS_total - SS_within) / SS_total
"""
import numpy as np
import pandas as pd
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def bray_curtis(X):
    n = X.shape[0]; D = np.zeros((n, n))
    for i in range(n):
        num = np.abs(X[i] - X).sum(axis=1); den = (X[i] + X).sum(axis=1)
        D[i] = np.where(den > 0, num / den, 0.0)
    return D


def permanova(D, lab):
    N = len(lab); groups = sorted(set(lab)); a = len(groups)
    ss_t = (D ** 2).sum / 2 / N
    ss_w = 0.0
    for g in groups:
        idx = [i for i, l in enumerate(lab) if l == g]
        sub = D[np.ix_(idx, idx)]; ss_w += (sub ** 2).sum / 2 / len(idx)
    F = ((ss_t - ss_w) / (a - 1)) / (ss_w / (N - a))
    return F, (ss_t - ss_w) / ss_t


def main:
    gen = pd.read_csv(os.path.join(BASE, 'tables', '38_unified_genus.tsv'), sep='\t', index_col=0)
    gen = gen.loc[[i for i in gen.index if i != 'Others']]
    X = gen.T.values.astype(float)
    coh = ['DS01' if s.startswith('ITF') else 'DS02' if s.startswith('HM') else 'DS03' for s in gen.columns]
    F, R2 = permanova(bray_curtis(X), coh)
    print(f'recomputed from tables/38 ({X.shape[0]} samples, {X.shape[1]} genera): pseudo-F = {F:.2f}, R^2 = {R2:.3f}')
    t49 = open(os.path.join(BASE, 'tables', '49_metagenome_animal_effect.tsv')).read
    pub = [l for l in t49.split('\n') if 'PERMANOVA' in l and 'free permutations' in l]
    print('published  (tables/49):', pub[0].split('\t')[2] if pub else 'NOT FOUND')
    ok = pub and f'{F:.2f}' in pub[0] and f'{R2:.3f}' in pub[0]
    print('verification:', 'OK — matches the published value' if ok else 'MISMATCH — investigate before citing')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main)
