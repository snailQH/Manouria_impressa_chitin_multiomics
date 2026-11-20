# Robustness set — eggNOG CAZy annotation

**These are not the manuscript's numbers.** Everything here is the same 32-sample unified gene
catalogue as the published tables, annotated with the **eggNOG** CAZy column instead of **dbCAN3**.
It exists to show which conclusions survive a change of annotation method, and nothing else.

| | manuscript (`tables/35`, `36`) | this directory |
|---|---|---|
| CAZy source | dbCAN3 (run_dbcan 5.2.9, ≥ 2 methods) | eggNOG CAZy column |
| GH18 CPM, DS01 / DS02 / DS03 | **471 / 665 / 317** | 391 / 339 / 305 |
| GH18 genes | **643** | 195 shared with dbCAN3 |

dbCAN3 is the manuscript's source because eggNOG under-calls GH18: the two agree on only 195 of the
643 dbCAN3 GH18 genes, and eggNOG misses 372 genes that Pfam independently supports.

**What this set establishes.** The *P. ostreatus* (DS03) GH18 deficit is robust — it holds under
dbCAN3, eggNOG and Pfam alike. The DS01-versus-DS02 ordering is **not**: it reverses between dbCAN3
and eggNOG, and the manuscript therefore never asserts it.

Cite `tables/35`–`46` for every cross-diet number. Cite this directory only when discussing
annotation-method sensitivity.
