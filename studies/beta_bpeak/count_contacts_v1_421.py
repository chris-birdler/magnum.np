#!/usr/bin/env python
"""
The contacts of K3 (jobs/plan421.txt, K3_0_v16: 104^3, q 16) and the V1 tie (thin_face, D0 of PLAN 4.21), written to
results/plan421_contacts_v1.json: the cell pairs of two particles that touch at a face, an edge or a corner (the run
geometry of geometry.fcc_box at q 16), split by the number of tied cells in the pair (0, 1, 2) and by the number of
outside cut cells (centre outside the sphere) in the pair; the masters of a different particle than their tied cell
(PLAN 4.21 "Contacts"; the code of tie_thin_face permits none).
    OMP_NUM_THREADS=1 python count_contacts_v1_421.py
"""
import json
import pathlib
import shlex

import numpy as np

import scan421
from field_terms_extra import tie_thin_face
from geometry import fcc_box, fcc_subcells

HERE = pathlib.Path(__file__).resolve().parent

if __name__ == "__main__":
    line = [l for l in (HERE / "jobs" / "plan421.txt").read_text().splitlines() if l.startswith("K3_0_v16 ")][0]
    N, d, dx, _, _ = scan421.mesh(line.split(None, 1)[1])
    args = shlex.split(line.split(None, 1)[1])
    dx_lex = float(args[args.index("--dx_lex") + 1])
    b = fcc_box(N, d, dx, sub=16)
    ids, frac = b["ids"], b["frac"]
    inside = fcc_box(N, d, dx, sub=0)["ids"] >= 0      # the cells with the centre inside (domains.geometry_from_config)
    tm = tie_thin_face(ids, frac, dx_lex, fcc_subcells(N, d, dx, 16)[2], 16)
    tied = np.zeros(N ** 3, bool)
    tied[tm.tied] = True
    tied = tied.reshape(N, N, N)
    pid = ids.reshape(-1)
    cross = int((pid[tm.tied] != pid[tm.src[:, 0]]).sum())
    res = {"run": "K3_0_v16", "N": int(N), "q": 16, "tie": tm.info, "masters_of_another_particle": cross,
           "pairs": {"face": 0, "edge": 0, "corner": 0}, "pairs_by_tied_cells": {"0": 0, "1": 0, "2": 0}}
    res["pairs_by_outside_cut_cells"] = {"0": 0, "1": 0, "2": 0}
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                k = abs(sx) + abs(sy) + abs(sz)
                if not k:
                    continue
                nb = np.roll(ids, (sx, sy, sz), axis=(0, 1, 2))
                sel = (ids >= 0) & (nb >= 0) & (ids != nb)
                res["pairs"][{1: "face", 2: "edge", 3: "corner"}[k]] += int(sel.sum())
                nt = tied[sel].astype(int) + np.roll(tied, (sx, sy, sz), axis=(0, 1, 2))[sel].astype(int)
                for v in (0, 1, 2):
                    res["pairs_by_tied_cells"][str(v)] += int((nt == v).sum())
                out = (ids >= 0) & ~inside
                no = out[sel].astype(int) + np.roll(out, (sx, sy, sz), axis=(0, 1, 2))[sel].astype(int)
                for v in (0, 1, 2):
                    res["pairs_by_outside_cut_cells"][str(v)] += int((no == v).sum())
    for key in ("pairs", "pairs_by_tied_cells", "pairs_by_outside_cut_cells"):
        res[key] = {k: v // 2 for k, v in res[key].items()}     # each pair is seen from both cells
    (HERE / "results" / "plan421_contacts_v1.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
