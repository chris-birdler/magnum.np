"""Location of the neighbour pairs with a large angle in a relaxed state (init.pt) or an end state (checkpoint.pt).

For each run directory: rebuild the geometry with fcc_box (main geometry.py, same sub as run_loops.surface_sub),
take the pairs of neighbour cells inside one particle (the pair masks of run_loops.diagnostics), and give for the
pairs with an angle above 60 and 90 degrees: the smaller f of the two cells, r/R of that cell (FCC site, minimum
image), and its number of magnetic neighbours in the same particle. Usage: pairs_loc.py <state> <run dir> ...
"""
import json
import math
import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from geometry import fcc_box, FCC_SITES  # noqa: E402


def sub_of(cfg):
    if cfg.get("surface", "voxel") != "fraction":
        return 0
    return cfg["staircase_q"] if cfg.get("staircase", "none") != "none" else 8


def analyse(run, state):
    cfg = json.load(open(run + "/config.json"))
    N, d, dx = cfg["N_used"], cfg["d"], cfg["dx"]
    g = fcc_box(N, d, dx, sub=sub_of(cfg))
    ids, f = g["ids"].astype(np.int64), g["frac"].astype(np.float64)
    m = torch.load(run + "/" + state, map_location="cpu", weights_only=False)["m"].double().numpy()
    a = N * dx
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rr = np.full((N, N, N), np.nan)
    for p, s in enumerate(FCC_SITES * a):
        dd = [c - sc for c, sc in zip((X, Y, Z), s)]
        dd = [v - a * np.round(v / a) for v in dd]
        r = np.sqrt(dd[0] ** 2 + dd[1] ** 2 + dd[2] ** 2) / (0.5 * d)
        rr[ids == p] = r[ids == p]
    nb = np.zeros((N, N, N), dtype=np.int64)
    for ax in range(3):
        for sh in (1, -1):
            nb += (ids >= 0) & (ids == np.roll(ids, sh, axis=ax))
    rows = []
    for ax in range(3):
        same = (ids >= 0) & (ids == np.roll(ids, 1, axis=ax))
        c = (m * np.roll(m, 1, axis=ax)).sum(-1)
        fo = np.roll(f, 1, axis=ax)
        ro = np.roll(rr, 1, axis=ax)
        no = np.roll(nb, 1, axis=ax)
        for i, j, k in zip(*np.nonzero(same & (c < 0.5))):
            if f[i, j, k] <= fo[i, j, k]:
                rows.append((math.degrees(math.acos(max(-1, min(1, c[i, j, k])))), f[i, j, k], fo[i, j, k],
                             rr[i, j, k], nb[i, j, k]))
            else:
                rows.append((math.degrees(math.acos(max(-1, min(1, c[i, j, k])))), fo[i, j, k], f[i, j, k],
                             ro[i, j, k], no[i, j, k]))
    rows = np.array(rows) if rows else np.zeros((0, 5))
    mag = ids >= 0
    surf = mag & (f < 1)
    out = {"run": run.split("/")[-1], "state": state, "n_cells": int(mag.sum()), "n_surface": int(surf.sum()),
           "f_surface_mean": float(f[surf].mean()) if surf.any() else None,
           "n_gt60": len(rows), "n_gt90": int((rows[:, 0] > 90).sum()) if len(rows) else 0,
           "max_angle": float(rows[:, 0].max()) if len(rows) else None}
    for lab, sel in (("gt60", rows[:, 0] > 60 if len(rows) else []), ("gt90", rows[:, 0] > 90 if len(rows) else [])):
        r = rows[sel] if len(rows) else rows
        if len(r):
            out[lab] = {"f_small_median": float(np.median(r[:, 1])), "f_small_max": float(r[:, 1].max()),
                        "share_f_small_lt_0.3": float((r[:, 1] < 0.3).mean()),
                        "share_both_f1": float(((r[:, 1] >= 1) & (r[:, 2] >= 1)).mean()),
                        "rR_min": float(r[:, 3].min()), "rR_median": float(np.median(r[:, 3])),
                        "nb_median": float(np.median(r[:, 4]))}
    return out


if __name__ == "__main__":
    state = sys.argv[1]
    for run in sys.argv[2:]:
        print(json.dumps(analyse(run, state)))
