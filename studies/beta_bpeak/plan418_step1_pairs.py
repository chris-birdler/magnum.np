"""PLAN 4.18 step 1, after the verdict: where are the large-angle pairs, and how large is the residual torque at the stop.

Reads diag_fields.npz of the step-1 runs (m_relaxed, ids, frac, h_relaxed_<term>). Pairs as eval_418.pair_angles:
neighbour cells of one particle, periodic. Pair classes: 0 = both cells full (f = 1), 1 = one cut cell, 2 = both cut
cells (0 < f < 1); CutPairExchange scales the classes 1 and 2. Tables:
1. per arm and seed: pairs of each class, number above 45, 60 and 90 deg, largest angle;
2. the pairs above 90 deg of T16 and B8: cells, f, the angle of the same pair and the tilt (eval_418.tilt_deg) of both
   cells in each arm of the same seed;
3. full-full pairs above 60 deg: distance (cells, Euclidean, periodic) of the nearer cell to the nearest cut cell;
4. residual torque |m x h| (A/m, sum of all terms) at the stop: largest value and f of its cell, 99th percentile and
   largest value in the cut cells and in the full cells.
Usage: plan418_step1_pairs.py <runs dir> > results/plan418_step1_pairs.md
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_418 import SEEDS, tilt_deg  # noqa: E402

ARMS = ["B8", "C8t", "C8h", "T4", "T16"]


def load(runs, arm, s):
    z = np.load(os.path.join(runs, "%s_s%d" % (arm, s), "diag_fields.npz"))
    return z


def pair_arrays(m, ids, f, ax):
    same = (ids >= 0) & (ids == np.roll(ids, 1, axis=ax))
    ang = np.degrees(np.arccos(np.clip((m * np.roll(m, 1, axis=ax)).sum(-1), -1.0, 1.0)))
    cut = (ids >= 0) & (f < 1.0)
    ncut = cut.astype(int) + np.roll(cut, 1, axis=ax).astype(int)
    return same, ang, ncut


def dist_to_cut(ids, f):
    N = ids.shape[0]
    cut = (ids >= 0) & (f < 1.0)
    d = distance_transform_edt(np.tile(~cut, (3, 3, 3)))     # periodic: tile and take the middle block
    return d[N:2 * N, N:2 * N, N:2 * N]


def main(runs):
    Z = {(a, s): load(runs, a, s) for a in ARMS for s in SEEDS}
    M = {k: z["m_relaxed"].astype(np.float64) for k, z in Z.items()}

    print("## 1. Pair angles by class (0 both full, 1 one cut cell, 2 both cut cells)\n")
    print("| arm | seed | class | pairs | > 45° | > 60° | > 90° | largest |")
    print("|---|---|---|---|---|---|---|---|")
    for s in SEEDS:
        for a in ARMS:
            ids, f = Z[a, s]["ids"], Z[a, s]["frac"].astype(np.float64)
            ang, nc = [], []
            for ax in range(3):
                same, an, ncut = pair_arrays(M[a, s], ids, f, ax)
                ang.append(an[same])
                nc.append(ncut[same])
            ang, nc = np.concatenate(ang), np.concatenate(nc)
            for k in (0, 1, 2):
                x = ang[nc == k]
                print("| %s | %d | %d | %d | %d | %d | %d | %.1f° |" % (
                    a, s, k, x.size, (x > 45).sum(), (x > 60).sum(), (x > 90).sum(), x.max()))

    print("\n## 2. Pairs above 90° of T16 and B8, the same pair in each arm\n")
    print("| run | cells | f | pair angle B8 / C8t / T4 / T16 | tilt cell 1 B8 / C8t / T4 / T16 | "
          "tilt cell 2 B8 / C8t / T4 / T16 |")
    print("|---|---|---|---|---|---|")
    arms4 = ["B8", "C8t", "T4", "T16"]
    for a0 in ("T16", "B8"):
        for s in SEEDS:
            ids, f = Z[a0, s]["ids"], Z[a0, s]["frac"].astype(np.float64)
            th = {a: tilt_deg(M[a, s], ids, f)[0] for a in arms4}
            for ax in range(3):
                same, ang, _ = pair_arrays(M[a0, s], ids, f, ax)
                for p in zip(*np.nonzero(same & (ang > 90.0))):
                    q = list(p)
                    q[ax] = (q[ax] - 1) % ids.shape[ax]
                    p, q = tuple(int(i) for i in p), tuple(q)
                    pa = [np.degrees(np.arccos(np.clip((M[a, s][p] * M[a, s][q]).sum(), -1, 1))) for a in arms4]
                    print("| %s_s%d | %s–%s | %.3f / %.3f | %s | %s | %s |" % (
                        a0, s, p, tuple(int(i) for i in q), f[p], f[q], " / ".join("%.1f" % x for x in pa),
                        " / ".join("%.1f" % th[a][p] for a in arms4), " / ".join("%.1f" % th[a][q] for a in arms4)))

    print("\n## 3. Full-full pairs above 60°: distance of the nearer cell to the nearest cut cell (cells)\n")
    print("| arm | pairs (4 seeds) | median | share at most 1.5 | largest | pairs above 90° |")
    print("|---|---|---|---|---|---|")
    for a in ARMS:
        ds, n90 = [], 0
        for s in SEEDS:
            ids, f = Z[a, s]["ids"], Z[a, s]["frac"].astype(np.float64)
            d = dist_to_cut(ids, f)
            full = (ids >= 0) & (f >= 1.0)
            for ax in range(3):
                same, ang, _ = pair_arrays(M[a, s], ids, f, ax)
                same &= full & np.roll(full, 1, axis=ax)
                ds.append(np.minimum(d, np.roll(d, 1, axis=ax))[same & (ang > 60.0)])
                n90 += int((same & (ang > 90.0)).sum())
        x = np.concatenate(ds)
        print("| %s | %d | %.1f | %.2f | %.1f | %d |" % (a, x.size, np.median(x), (x <= 1.5).mean(), x.max(), n90))

    print("\n## 4. Residual torque |m × h| at the stop (A/m)\n")
    print("| run | relax iterations | largest | f of its cell | cut cells: p99 | cut cells: largest | "
          "full cells: p99 | full cells: largest |")
    print("|---|---|---|---|---|---|---|---|")
    for a in ARMS:
        for s in SEEDS:
            z = Z[a, s]
            ids, f = z["ids"], z["frac"].astype(np.float64)
            h = sum(z["h_relaxed_" + str(t)].astype(np.float64) for t in z["term_names"])
            tq = np.linalg.norm(np.cross(M[a, s], h), axis=-1)
            mag = ids >= 0
            tq[~mag] = 0.0
            cut, full = mag & (f < 1.0), mag & (f >= 1.0)
            i = np.unravel_index(np.argmax(tq), tq.shape)
            it = json.load(open(os.path.join(runs, "%s_s%d" % (a, s), "summary.json")))["init"]["relax_iterations"]
            print("| %s_s%d | %d | %.2f | %.3f | %.2f | %.2f | %.2f | %.2f |" % (
                a, s, it, tq.max(), f[i], np.percentile(tq[cut], 99), tq[cut].max(),
                np.percentile(tq[full], 99), tq[full].max()))


if __name__ == "__main__":
    main(sys.argv[1])
