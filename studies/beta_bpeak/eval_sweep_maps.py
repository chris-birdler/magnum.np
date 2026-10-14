#!/usr/bin/env python
"""
Frequency sweep at 9 mT (PROTOKOLL 7.92), option 1 of Chris (2026-10-14, second step): does A 20 at 2 f lose in the same
zones as A 10 at f? Local, no GPU. Hypothesis (7.92): the eigen-frequencies of the loss scale with A, so the w(f) curve of
A 20 is that of A 10 shifted by the factor 2 in f. If the same motion carries the loss, the loss map of A 20 at 2 f must
look like the map of A 10 at f more than the map of A 20 at f does.

Data: the cut z = 0 (acc_slices.npz, stage amp00 = 9 mT) of the runs of eval_sweep_state.py (corner K1, seeds 0 … 2).
The A 10 box has 152 cells per edge (d/l_ex 304), the A 20 box 108 (216); the geometry in units of the box is the same
(same seed, stress axes and packing; the grain axes differ cell by cell). The A 10 map is averaged onto the A 20 grid
(block mean of the particle cells; a coarse cell counts if it is a particle cell of A 20 and holds A 10 particle cells).
Map: ln p, p = local loss density / its mean over the particle cells of the cut (p_rel). Second map (reported, not in
the reading): ln |m1|², the amplitude of the motion at the drive frequency.

Similarity r_s(k, f) = Pearson correlation over the common particle cells between the map of A 10 at f and the map of
A 20 at √2^s · f, s = 0, 1, 2, 3, for f = 15, 21, 30, 42 MHz (the frequencies where all 4 shifts exist). R_s(k) = mean
of r_s over these f.

Reading (fixed before the maps were compared, 2026-10-14):
  supports        R_2 is the largest of R_0 … R_3 in all 3 seeds;
  against         R_0 is the largest in at least 2 of 3 seeds (the loss pattern follows f, not f scaled with A);
  not determined  otherwise.

    python eval_sweep_maps.py > results/sweep9_maps.md      # also results/sweep9_maps.json
"""
import json
import math
import pathlib

import numpy as np

import eval_sweep_state as ES

POS = {15: 0, 21: 1, 30: 2, 42: 3, 60: 4, 85: 5, 120: 6}
INV = {v: f for f, v in POS.items()}
F_BASE = (15, 21, 30, 42)
SHIFTS = range(4)


def load(a, f, k, key):
    z = np.load(ES.folder(a, f, k) / "acc_slices.npz")
    v = z["amp00_" + key].astype(np.float64)
    return v, z["ids"]


def to_grid(v, ids, n_to):
    """Block mean of the particle cells of an n × n map onto an n_to × n_to grid (same box); count per coarse cell."""
    n = v.shape[0]
    idx = np.floor((np.arange(n) + 0.5) / n * n_to).astype(int)
    I, J = np.meshgrid(idx, idx, indexing="ij")
    m = ids >= 0
    acc, cnt = np.zeros((n_to, n_to)), np.zeros((n_to, n_to))
    np.add.at(acc, (I[m], J[m]), v[m])
    np.add.at(cnt, (I[m], J[m]), 1)
    return np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan), cnt


def corr(key, k, f, s):
    v10, i10 = load(10, f, k, key)
    v20, i20 = load(20, INV[POS[f] + s], k, key)
    c10, cnt = to_grid(v10, i10, v20.shape[0])
    ok = (i20 >= 0) & (cnt > 0)
    x, y = np.log(np.maximum(c10[ok], 1e-12)), np.log(np.maximum(v20[ok], 1e-12))
    return float(np.corrcoef(x, y)[0, 1]), int(ok.sum())


def main():
    res = {}
    L = ["# Loss maps in the cut: A 10 at f against A 20 at √2^s · f (9 mT, corner K1)", "",
         "r: Pearson correlation of the log maps over the common particle cells. R_s: mean over f = 15, 21, 30, 42 MHz. "
         "Reading fixed before the comparison (docstring).", ""]
    for key, lab in (("p_rel", "loss density p (reading)"), ("m1sq", "motion at f, |m1|² (reported only)")):
        L += ["## %s" % lab, "", "| seed | f (MHz) | " + " | ".join("s %d (A 20 at × %.2f)" % (s, math.sqrt(2) ** s)
                                                            for s in SHIFTS) + " |", "|---|---|" + "---|" * len(SHIFTS)]
        R = np.zeros((3, len(SHIFTS)))
        for k in ES.SEEDS:
            rows = []
            for f in F_BASE:
                r = [corr(key, k, f, s) for s in SHIFTS]
                rows.append([x[0] for x in r])
                L.append("| %d | %d | %s |" % (k, f, " | ".join("%.3f" % x[0] for x in r)))
            R[k] = np.mean(rows, 0)
            L.append("| %d | **R_s** | %s |" % (k, " | ".join("**%.3f**" % x for x in R[k])))
        best = [int(np.argmax(R[k])) for k in ES.SEEDS]
        if key == "p_rel":
            rd = ("supports" if all(b == 2 for b in best) else
                  "against" if sum(b == 0 for b in best) >= 2 else "not determined")
        else:
            rd = None
        res[key] = {"R": R.tolist(), "best_shift": best, "reading": rd}
        L += ["", "Largest R_s per seed: shift %s." % " / ".join(str(b) for b in best) +
              (" **Reading: %s**" % rd if rd else ""), ""]
    print("\n".join(L))
    (pathlib.Path(__file__).resolve().parent / "results" / "sweep9_maps.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
