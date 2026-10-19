#!/usr/bin/env python
"""
Census of Bloch points in the state files of run folders (Chris 2026-10-18, after the Bloch-point pair of K3_1_nx).
For each state file of a run (init.pt, ckpt_<i>_<stage>.pt, checkpoint.pt): the pairs of neighbour cells of one particle
above 60 / 90 / 120 / 150 degrees, and the Bloch points: the elementary cubes of 2 x 2 x 2 cell centres, all in one
particle, with a topological charge |Q| > 0.5 (analyze_defect.bloch_points; exact integer charges, one cube per Bloch
point). Restriction: a Bloch point in a cube with a corner outside the particle (within about 1 cell of the surface)
is not counted.

    python census_bloch.py --out census.json runs/K1_0_nx runs/K1_1_nx …
"""
import argparse
import json
import pathlib
import re

import numpy as np

import analyze_defect as A



def states(d):
    order = lambda f: (0 if f.name == "init.pt" else 2 if f.name == "checkpoint.pt" else 1, f.name)
    return sorted([f for f in d.glob("*.pt") if f.name == "init.pt" or f.name == "checkpoint.pt"
                   or re.match(r"ckpt_\d\d_.+\.pt$", f.name)], key=order)


def census(m, ids):
    pa = A.pair_angles(m, ids)
    counts = {str(int(t)): int(sum(((ang > t) & same).sum() for ang, same in pa)) for t in (60, 90, 120, 150)}
    bps = [{"particle": p, "Q": q, "cube": list(c), "dist_to_surface_cells": round(d, 1)}
           for p, q, c, d in A.bloch_points(m, ids)]
    return {"pairs_above": counts, "bloch_points": bps}


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
            print("%s %s: pairs >60/90/120/150 %s, Bloch points %d %s" % (
                d.name, f.name, "/".join(str(v) for v in c["pairs_above"].values()), len(c["bloch_points"]),
                [(b["particle"], b["Q"], b["dist_to_surface_cells"]) for b in c["bloch_points"]]), flush=True)
    pathlib.Path(a.out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
