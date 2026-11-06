#!/usr/bin/env python
"""
Counts for decision D0 of PLAN 4.21 (the tie of the cut cells), written to results/plan421_d0_counts.json, at K5 (the
run geometry of K5f12 of gate 1) and at K3 at q 16 (geometry.fcc_box of jobs/plan421.txt; a few minutes):
  n_cut          cut cells (particle cells with f below 1)
  n_thin         thin cut cells (f below l_ex/dx): the cells that V1 (thin_face) ties
  n_v2           cut cells with the centre outside the sphere: the cells that V2 (ghost_mirror) and option (a) tie
  nearest_cells  option (a): the distance from the mirror point of each V2 cell (the arithmetic of tie_ghost_mirror)
                 to the nearest cell of the same particle with the centre inside (field_terms_extra._nearest_cell), in
                 cells: median and largest
    python count_d0_421.py runs/plan421_gate1
"""
import json
import math
import pathlib
import shlex
import sys

import numpy as np

import scan421
from field_terms_extra import _nearest_cell
from geometry import FCC_SITES, fcc_box

HERE = pathlib.Path(__file__).resolve().parent


def counts(ids, frac, N, d, dx, dx_lex):
    a, R = N * dx, 0.5 * d
    cut = (ids >= 0) & (frac < 1.0)
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rr = np.full((N, N, N), np.inf)
    sites = FCC_SITES * a
    for p, s in enumerate(sites):
        sel = ids == p
        dd = [C[sel] - s[i] for i, C in enumerate((X, Y, Z))]
        rr[sel] = sum((v - a * np.round(v / a)) ** 2 for v in dd)
    inside = (ids >= 0) & (rr <= R**2)
    inside_p = {p: inside & (ids == p) for p in range(len(FCC_SITES))}
    dist = []
    for i, j, k in np.argwhere(cut & ~inside):
        p = ids[i, j, k]
        c = np.array([X[i, j, k], Y[i, j, k], Z[i, j, k]]) - sites[p]
        c -= a * np.round(c / a)
        r = math.sqrt(float((c**2).sum()))
        mp = (sites[p] + c * (2.0 * R - r) / r) / dx - 0.5
        nc = np.asarray(_nearest_cell(mp, inside_p[p]), np.float64)
        e = nc - mp
        e -= N * np.round(e / N)
        dist.append(math.sqrt(float((e**2).sum())))
    dist = np.asarray(dist)
    return {"N": int(N), "dx_lex": float(dx_lex), "n_particle_cells": int((ids >= 0).sum()), "n_cut": int(cut.sum()),
            "n_thin": int((cut & (frac < 1.0 / dx_lex)).sum()), "n_v2": int(len(dist)),
            "nearest_cells": {"median": float(np.median(dist)), "max": float(dist.max())}}


if __name__ == "__main__":
    folder = pathlib.Path(sys.argv[1]).resolve()
    c = json.loads((folder / "K5f12" / "config.json").read_text())
    g = np.load(folder / "K5f12" / "geometry.npz")
    res = {"K5": counts(g["ids"].astype(np.int64), g["frac"], c["N_used"], c["d"], c["dx"], c["dx_lex"])}
    line = [l for l in (HERE / "jobs" / "plan421.txt").read_text().splitlines() if l.startswith("K3_0_v16 ")][0]
    args = shlex.split(line.split(None, 1)[1])
    N, d, dx, _, _ = scan421.mesh(line.split(None, 1)[1])
    b = fcc_box(N, d, dx, sub=16)
    res["K3_q16"] = counts(b["ids"].astype(np.int64), b["frac"], N, d, dx, float(args[args.index("--dx_lex") + 1]))
    (HERE / "results" / "plan421_d0_counts.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
