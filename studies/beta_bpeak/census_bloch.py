#!/usr/bin/env python
"""
Census of Bloch points in the state files of run folders (Chris 2026-10-18, after the Bloch-point pair of K3_1_nx).
For each state file of a run (init.pt, ckpt_<i>_<stage>.pt, checkpoint.pt): the pairs of neighbour cells of one particle
above 60 / 90 / 120 / 150 degrees, their clusters, and the Bloch points: clusters with a topological charge |Q| > 0.5 on
the surface of a cube of half-edge 2 cells around the cluster centre that lies at least 95 % in one particle.
Restriction: a Bloch point closer than about 3 cells to the particle surface is not counted (no closed surface inside
the particle); the clusters at the surface (the corners of the voxel staircase) are counted separately.

    python census_bloch.py --out census.json runs/K1_0_nx runs/K1_1_nx …
"""
import argparse
import json
import pathlib
import re

import numpy as np
from scipy import ndimage

import analyze_defect as A

R_Q, IN_MIN = 2, 0.95


def states(d):
    order = lambda f: (0 if f.name == "init.pt" else 2 if f.name == "checkpoint.pt" else 1, f.name)
    return sorted([f for f in d.glob("*.pt") if f.name == "init.pt" or f.name == "checkpoint.pt"
                   or re.match(r"ckpt_\d\d_.+\.pt$", f.name)], key=order)


def census(m, ids):
    pa = A.pair_angles(m, ids)
    counts = {str(int(t)): int(sum(((ang > t) & same).sum() for ang, same in pa)) for t in (60, 90, 120, 150)}
    cells = np.zeros(ids.shape, bool)
    for a, (ang, same) in enumerate(pa):
        h = (ang > 60.0) & same
        cells |= h | np.roll(h, -1, axis=a)
    lab, n = ndimage.label(cells)
    bps, surface = [], 0
    for k in range(1, n + 1):
        pos = np.argwhere(lab == k)
        pid = int(np.bincount(ids[tuple(pos.T)][ids[tuple(pos.T)] >= 0]).argmax())
        cen = np.round(pos.mean(0)).astype(int)
        q, fin = A.charge(m, ids, cen, R_Q, pid)
        if fin < IN_MIN:
            surface += 1
            continue
        if abs(q) > 0.5:
            bps.append({"particle": pid, "centre": cen.tolist(), "Q": round(q, 2), "cells": len(pos)})
    return {"pairs_above": counts, "clusters": n, "surface_clusters": surface, "bloch_points": bps}


def main():
    import torch
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    res = {}
    for r in a.runs:
        d = pathlib.Path(r)
        cfg = json.loads((d / "config.json").read_text())
        ids, _, _ = A.geometry(cfg)
        res[d.name] = {}
        for f in states(d):
            m = torch.load(str(f), map_location="cpu")["m"].double().numpy()
            m /= np.maximum(np.linalg.norm(m, axis=-1, keepdims=True), 1e-30)
            c = census(m, ids)
            res[d.name][f.name] = c
            print("%s %s: pairs >60/90/120/150 %s, clusters %d (surface %d), Bloch points %d %s" % (
                d.name, f.name, "/".join(str(v) for v in c["pairs_above"].values()), c["clusters"],
                c["surface_clusters"], len(c["bloch_points"]), [(b["particle"], b["Q"]) for b in c["bloch_points"]]),
                flush=True)
    pathlib.Path(a.out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
