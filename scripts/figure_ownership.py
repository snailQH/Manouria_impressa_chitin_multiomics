#!/usr/bin/env python3
"""Which script owns which figure file.

Fourteen figure names were written by two scripts each, so the file on disk was decided by
run order rather than by intent. On 2026-11-16 that put the superseded vendor-table Figure 7C
(withdrawn digestibility index, "three vendor pipelines" footnote) into the submission bundle,
because a superseded vendor-table script ran after scripts/generate_unified_figs.py.

Every figure-writing script consults this table through `may_write`. A script that is not the
owner skips the file and says so, instead of silently overwriting the manuscript's version.
Set FIGURE_OWNERSHIP_OFF=1 to bypass (provenance rebuilds only).
"""
import os

OWNER = {
    # unified 32-sample catalogue — the manuscript's cross-diet set
    'Fig2A_phylum_barplot': 'generate_unified_figs.py',
    'Fig2B_PCoA': 'generate_unified_figs.py',
    'Fig2C_alpha_diversity': 'generate_unified_figs.py',
    'Fig2D_genus_heatmap': 'generate_unified_figs.py',
    'Fig3A_GH18_boxplot': 'generate_unified_figs.py',
    'Fig3B_CAZyme_heatmap': 'generate_unified_figs.py',
    'Fig3D_digestibility_indices': 'generate_unified_figs.py',
    'SuppFig_S19_diet_response': 'generate_unified_figs.py',
    # the pre-restructure name of the same panel: generate_fig7c.py (vendor tables) must not claim it
    'Fig7C_diet_response': 'generate_unified_figs.py',
    'SuppFig_S7_genus_barplot': 'generate_unified_figs.py',
    # genome-resolved (MAG) version of the GH18 taxonomy
    'Fig3C_GH18_taxonomy': 'mag_gh18_analysis.py',
    # species-specific DIA-NN search, not the provider proxy-database search
    'Fig7A_chitin_protein_categories': 'proteomics_diann_chitin.py',
    'Fig7B_tissue_chitinase_expression': 'proteomics_diann_chitin.py',
    'SuppFig_S14_proteomics_overview': 'proteomics_diann_chitin.py',
    # KEGG profile on the single catalogue
    'SuppFig_S18_KEGG_functional': 'generate_kegg_functional_unified32.py',
}


def may_write(name, caller):
    """True if `caller` (__file__ of the calling script) may write figure `name`."""
    owner = OWNER.get(name)
    if owner is None or os.environ.get('FIGURE_OWNERSHIP_OFF'):
        return True
    if os.path.basename(caller) == owner:
        return True
    print(f'  skip {name}: owned by scripts/{owner} (see scripts/figure_ownership.py)')
    return False
