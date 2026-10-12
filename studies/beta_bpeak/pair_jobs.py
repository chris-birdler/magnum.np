#!/usr/bin/env python
"""
Paired field images for the A effect (decision Chris 2026-10-12, option 3): the corner K1 (A 20 pJ/m at d = 1 µm, i.e.
d/l_ex 212; L_eff 30 l_ex = K_eff 1.0 kJ/m³; K_σ 1.5 K_eff; φ 0.62) again with d/l_ex 300 (A 10 pJ/m), the same seeds
2011 … 2013 and all other arguments equal (snapshots at 9 / 50 / 150 mT, accumulators). In this model A enters only
as d/l_ex, so the pair differs only in A. The Herzer-cell grid scales with the box, so the random axes are not the
same cells: the pair shares seed, stress axes, packing and history, not the cell-by-cell microstructure.

    python pair_jobs.py > jobs/pair_K1A10.txt
"""
import shlex

from rerun_jobs import job_lines, set_arg

src = dict(job_lines("jobs/sobol_corners.txt"))
for k in range(3):
    tok = set_arg(shlex.split(src["K1_%d" % k]), "d_lex", ["300.0"])
    print("K1A10_%d %s" % (k, " ".join(tok)))
