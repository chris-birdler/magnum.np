#!/usr/bin/env python
"""
PROTOKOLL 3.45: the T2 bias b of the stage values (s1_states.py --against: b = sum (q - q_ref) . q_ref / sum q_ref .
q_ref for each quantity and state). Reads T2 files of the B states (16 states, 3 repeats each); trees with bitwise
equal b in all states and quantities count once. For each quantity (stage3, stage4, stage5, dx), b in 1e-6:
  (a) the repeat noise: the difference of b between two repeats of one tree in one state, for the trees whose repeats
      are not bitwise equal (--noise; the CSR product with atomic additions); RMS and largest over all pairs, the
      range of the RMS over the states, and the floor for a difference of two medians of 3 repeats (normal noise:
      the median of 3 has the variance 0.449 sigma^2);
  (b) the change of b in each state between two trees (--pair OLD NEW); RMS over the states and largest;
  (c) for each tree: the RMS of b over the states (the T2 measure), the largest state and its share of the sum of
      squares, the mean of the signed b over the states +- SE (n 16), and e (RMS over the states);
  (d) group A (--group_a) against the other trees: the correlation of the signed b over the states within A,
      between A and the others, within the others; the RMS over the states of the group mean of b and their ratio
      A / others with a 95 % interval (bootstrap over the states, 5000 draws, seed 1); the number of states with one
      sign in all trees and in all trees of A;
  (e) for each pair: the ratio of the RMS of b of NEW and OLD, also without one state (--drop K0) in both trees.

    python t2_b_noise.py --tree NAME=FILE ... --noise NAME ... --pair OLD NEW ... --group_a NAME ... --drop 8
"""
import argparse
import math
import random
import statistics

import torch

QS = ("stage3", "stage4", "stage5", "dx")


def rms(v):
    return math.sqrt(sum(x * x for x in v) / len(v))


def corr(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    return (sum((a - mx) * (c - my) for a, c in zip(x, y))
            / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((c - my) ** 2 for c in y)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", nargs="+", required=True, help="NAME=FILE (T2 files of s1_states.py --against)")
    ap.add_argument("--noise", nargs="+", default=[], help="trees with repeats that are not bitwise equal")
    ap.add_argument("--pair", nargs=2, action="append", default=[], metavar=("OLD", "NEW"))
    ap.add_argument("--group_a", nargs="+", default=[], help="trees of group A for (d)")
    ap.add_argument("--drop", type=int, default=None, help="k0 of the state to drop in (e)")
    a = ap.parse_args()
    R, sig = {}, {}
    for t in a.tree:
        k, f = t.split("=", 1)
        r = torch.load(f, map_location="cpu", weights_only=False)
        s = repr(sorted((k0, q, v["b"]) for k0, st in r["states"].items() for q, v in st["q"].items() if "b" in v))
        if s in sig:
            print("tree %s: b bitwise equal to %s, not counted again" % (k, sig[s]))
            continue
        sig[s] = k
        R[k] = r
    ks = sorted(next(iter(R.values()))["states"])
    for r in R.values():
        if sorted(r["states"]) != ks:
            raise SystemExit("t2_b_noise: the trees have other states")
    print("%d distinct trees, %d states" % (len(R), len(ks)))
    B = lambda t, q: [R[t]["states"][k0]["q"][q]["b"] for k0 in ks]
    for q in QS:
        print("== %s, b in 1e-6" % q)
        if a.noise:
            d, per = [], {k0: [] for k0 in ks}
            for t in a.noise:
                for k0 in ks:
                    v = [rr[q]["b"] for rr in R[t]["states"][k0]["reps"]]
                    p = [v[i] - v[j] for i in range(len(v)) for j in range(i + 1, len(v))]
                    d += p
                    per[k0] += p
            pr = {k0: rms(v) for k0, v in per.items()}
            lo, hi = min(pr, key=pr.get), max(pr, key=pr.get)
            sigma = rms(d) / math.sqrt(2)
            print("(a) repeat noise, %d trees, %d pairs: RMS %.2f, largest %.2f; RMS for each state %.2f (k0 %d) ... "
                  "%.2f (k0 %d); floor for two medians of 3: %.2f"
                  % (len(a.noise), len(d), 1e6 * rms(d), 1e6 * max(abs(x) for x in d), 1e6 * pr[lo], lo, 1e6 * pr[hi],
                     hi, 1e6 * sigma * math.sqrt(2 * 0.449)))
        for o, n in a.pair:
            dd = [x - y for x, y in zip(B(n, q), B(o, q))]
            print("(b) change %s -> %s: RMS %.2f, largest %.2f" % (o, n, 1e6 * rms(dd), 1e6 * max(abs(x) for x in dd)))
        for t in R:
            v = B(t, q)
            i = max(range(len(v)), key=lambda j: abs(v[j]))
            m, se = statistics.mean(v), statistics.stdev(v) / math.sqrt(len(v))
            e = rms([R[t]["states"][k0]["q"][q]["e"] for k0 in ks])
            print("(c) %-6s RMS %6.2f; largest k0 %3d %+6.1f (%2.0f %% of the sum of squares); mean %+5.2f +- %.2f "
                  "(%+.2f SE); e %.3e" % (t, 1e6 * rms(v), ks[i], 1e6 * v[i], 100 * v[i] ** 2 / sum(x * x for x in v),
                                         1e6 * m, 1e6 * se, m / se, e))
        if a.group_a:
            A = [t for t in R if t in a.group_a]
            O = [t for t in R if t not in a.group_a]
            cA = [corr(B(x, q), B(y, q)) for i, x in enumerate(A) for y in A[i + 1:]]
            cX = [corr(B(x, q), B(y, q)) for x in A for y in O]
            cO = [corr(B(x, q), B(y, q)) for i, x in enumerate(O) for y in O[i + 1:]]
            mA = [statistics.mean(B(t, q)[i] for t in A) for i in range(len(ks))]
            mO = [statistics.mean(B(t, q)[i] for t in O) for i in range(len(ks))]
            rnd = random.Random(1)
            bs = []
            for _ in range(5000):
                idx = [rnd.randrange(len(ks)) for _ in ks]
                bs.append(rms([mA[i] for i in idx]) / rms([mO[i] for i in idx]))
            bs.sort()
            sgn = lambda ts: sum(len({x > 0 for x in (B(t, q)[i] for t in ts)}) == 1 for i in range(len(ks)))
            print("(d) %d trees in A, %d others; correlation within A %.2f ... %.2f (mean %.2f), A-others mean %.2f, "
                  "within the others mean %.2f; RMS of the group mean A %.2f, others %.2f, ratio %.2f (95 %% %.2f ... "
                  "%.2f); states with one sign: in all trees %d, in A %d"
                  % (len(A), len(O), min(cA), max(cA), statistics.mean(cA), statistics.mean(cX), statistics.mean(cO),
                     1e6 * rms(mA), 1e6 * rms(mO), rms(mA) / rms(mO), bs[124], bs[4874], sgn(list(R)), sgn(A)))
        if a.drop is not None:
            j = ks.index(a.drop)
            for o, n in a.pair:
                vo, vn = B(o, q), B(n, q)
                wo, wn = vo[:j] + vo[j + 1:], vn[:j] + vn[j + 1:]
                print("(e) RMS ratio %s / %s: %.2f; without k0 %d in both: %.2f" % (n, o, rms(vn) / rms(vo), a.drop,
                                                                                   rms(wn) / rms(wo)))


if __name__ == "__main__":
    main()
