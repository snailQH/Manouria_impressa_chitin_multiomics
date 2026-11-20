#!/usr/bin/env python3
"""Catch figures that silently misrepresent their data.

Motivated by Figure 6B (2026-11-17): an axis limit set before the bars were drawn pinned the view
to 0-1 and disabled autoscaling, so every stacked bar was clipped at 1.0 and only the first
component showed. The numbers were right; the axis was wrong, and nothing errored.

Wraps Figure.savefig and, for every axes, compares what was drawn (dataLim) with what is visible
(viewLim). Anything drawn outside the view is reported with the fraction lost. Deliberate zooms
exist, so this prints findings for a human to judge rather than failing the build.

    python3 scripts/figure_integrity_check.py [generator.py ...]
"""
import os, sys, runpy, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

FINDINGS = []
_orig = matplotlib.figure.Figure.savefig


def _check(fig, fname):
    name = os.path.basename(str(fname))
    if not name.endswith('.png'):
        return
    for i, ax in enumerate(fig.axes):
        try:
            d, v = ax.dataLim, ax.viewLim
            if not np.all(np.isfinite(d.get_points)):
                continue
            if d.width <= 0 and d.height <= 0:
                continue
            for axis, lo_d, hi_d, lo_v, hi_v, span in (
                    ('x', d.x0, d.x1, v.x0, v.x1, max(d.width, 1e-12)),
                    ('y', d.y0, d.y1, v.y0, v.y1, max(d.height, 1e-12))):
                lo_v, hi_v = min(lo_v, hi_v), max(lo_v, hi_v)       # inverted axes are fine
                out = max(0.0, lo_v - lo_d) + max(0.0, hi_d - hi_v)
                if out / span > 0.02:
                    FINDINGS.append((name, i, axis, round(100 * out / span, 1),
                                     (round(lo_d, 3), round(hi_d, 3)), (round(lo_v, 3), round(hi_v, 3))))
        except Exception:
            pass


def patched(self, fname, *a, **k):
    _check(self, fname)
    return _orig(self, fname, *a, **k)


matplotlib.figure.Figure.savefig = patched

if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    for script in sys.argv[1:]:
        sys.argv = [script]
        try:
            runpy.run_path(script, run_name='__main__')
        except SystemExit:
            pass
        except Exception as e:
            print(f'  !! {os.path.basename(script)} raised {type(e).__name__}: {e}')
        plt.close('all')
    print(f'\n{"=" * 92}\ndata drawn outside the visible axes ({len(FINDINGS)} finding(s)):')
    for n, i, axis, pct, dl, vl in sorted(FINDINGS, key=lambda f: -f[3]):
        print(f'  {n:52s} axes{i} {axis}: {pct:5.1f}% outside   data {dl} view {vl}')
