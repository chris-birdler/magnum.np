"""PLAN 4.18 step 3, after the verdict (go of Chris 2026-11-04): which cells form the large-angle pairs in the relaxed
states of C, D and Dk (4.16 geometry, 100 mT arms, seeds 921, 923, 925)?

Reads init.pt (the relaxed state of each run) and the geometry of the run (eval_418.geometry). The states during the
cycles are not available (no snapshots or checkpoints were kept), so the per-cycle counts of reading 3c cannot be
classified. Pairs as eval_418.pair_angles: neighbour cells of one particle, periodic. Pair types: full-full (both
f = 1), cut-full, cut-cut (0 < f < 1); "with a thin cell" if a cell of the pair has f below l_ex/dx (1/3). For the
full-full pairs: the distance (cells, Euclidean, periodic) of the nearer cell to the nearest cut cell.
Usage: plan418_step3_pairs.py <runs dir> > results/plan418_step3_pairs.md
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_418 import S3_SEEDS, S3_ARMS, THIN, geometry, load_run  # noqa: E402
from plan418_step1_pairs import dist_to_cut, pair_arrays  # noqa: E402


def pairs(m, ids, f):
    """Per pair above 60 deg: (angle, number of cut cells, thin cell in the pair, distance of the nearer cell to the
    nearest cut cell)."""
    dc = dist_to_cut(ids, f)
    thin = (ids >= 0) & (f < THIN)
    out = []
    for ax in range(3):
        same, ang, ncut = pair_arrays(m, ids, f, ax)
        sel = same & (ang > 60.0)
        th = thin | np.roll(thin, 1, axis=ax)
        d = np.minimum(dc, np.roll(dc, 1, axis=ax))
        for a, n, t, di in zip(ang[sel], ncut[sel], th[sel], d[sel]):
            out.append((float(a), int(n), bool(t), float(di)))
    return out


def main(runs):
    rows = {}
    geo = {}
    for s in S3_SEEDS:
        for a in S3_ARMS:
            cfg, m, rec = load_run(os.path.join(runs, "S3%s_s%d" % (a, s)))
            if s not in geo:
                geo[s] = geometry(cfg)
            ids, f = geo[s]
            rows[a, s] = pairs(m, ids, f)
            print("done", a, s, len(rows[a, s]), file=sys.stderr, flush=True)

    def kind(n, t):
        return {0: "full-full", 1: "cut-full", 2: "cut-cut"}[n] + (" (thin)" if t else "")

    kinds = ["full-full", "cut-full", "cut-full (thin)", "cut-cut", "cut-cut (thin)"]
    for lim in (90.0, 60.0):
        print("## Pairs above %.0f° in the relaxed states, by type\n" % lim)
        print("| arm | seed | all | " + " | ".join(kinds) + " | largest angle | type of the largest |")
        print("|---|---|---|" + "---|" * len(kinds) + "---|---|")
        for a in S3_ARMS:
            for s in S3_SEEDS:
                p = [x for x in rows[a, s] if x[0] > lim]
                cnt = {k: sum(1 for x in p if kind(x[1], x[2]) == k) for k in kinds}
                big = max(p, default=None)
                print("| %s | %d | %d | %s | %s | %s |" % (
                    a, s, len(p), " | ".join(str(cnt[k]) for k in kinds), "%.1f°" % big[0] if big else "–",
                    kind(big[1], big[2]) if big else "–"))
        print()
    print("## Full-full pairs above 60°: distance of the nearer cell to the nearest cut cell (cells)\n")
    print("| arm | seed | pairs | distance 1 | distance up to 2 | larger | largest angle |")
    print("|---|---|---|---|---|---|---|")
    for a in S3_ARMS:
        for s in S3_SEEDS:
            p = [x for x in rows[a, s] if x[1] == 0]
            print("| %s | %d | %d | %d | %d | %d | %s |" % (
                a, s, len(p), sum(1 for x in p if x[3] <= 1.0), sum(1 for x in p if x[3] <= 2.0),
                sum(1 for x in p if x[3] > 2.0), "%.1f°" % max(x[0] for x in p) if p else "–"))


if __name__ == "__main__":
    main(sys.argv[1])
