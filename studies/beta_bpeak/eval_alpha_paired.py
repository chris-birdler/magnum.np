#!/usr/bin/env python
"""
Paired alpha comparison (PROTOKOLL 7.31): R300v01_* (alpha 0.01, 9 cycles) vs M300_dx3_* (alpha 0.1, 5 cycles),
the same relaxed state (init.pt), mesh dx 3, staircase surface, one job per amplitude (9 / 50 / 100 mT),
seeds 921 ... 924. Paired differences: mean ± SE, 95 % CI (Student t, df 3), two-sided p.

    python eval_alpha_paired.py > results/alpha_paired.md
"""
import math
import pathlib

import numpy as np
from scipy import stats

import analyze

RUNS = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = (921, 922, 923, 924)
AMPS = (9, 50, 100)


def one(name):
    _, rows = analyze.load_run(RUNS / name)
    r = rows[0]
    return r["w_loop"], r["b_peak"], r["balance"]


def summ(d):
    d = np.asarray(d, float)
    m, se = d.mean(), d.std(ddof=1) / math.sqrt(len(d))
    t = stats.t.ppf(0.975, len(d) - 1)
    p = 2 * stats.t.sf(abs(m / se), len(d) - 1) if se > 0 else float("nan")
    return "%+.3f ± %.3f [%+.2f, %+.2f], p %.3f" % (m, se, m - t * se, m + t * se, p)


def main():
    A = {(s, m): one("R300v01_s%d_b%d" % (s, m)) for s in SEEDS for m in AMPS}
    B = {(s, m): one("M300_dx3_s%d_b%d" % (s, m)) for s in SEEDS for m in AMPS}
    print("| seed | B | w α 0.1 | w α 0.01 | Δ ln w | balance α 0.1 / 0.01 |")
    print("|---|---|---|---|---|---|")
    for s in SEEDS:
        for m in AMPS:
            print("| %d | %d mT | %.3g | %.3g | %+.3f | %.2f / %.2f |" %
                  (s, m, B[(s, m)][0], A[(s, m)][0], math.log(A[(s, m)][0] / B[(s, m)][0]), B[(s, m)][2], A[(s, m)][2]))
    print("\n**Paired differences α 0.01 − α 0.1 (n = 4; mean ± SE [95 % CI], p)**\n")
    for m in AMPS:
        print("- Δ ln w at %d mT: %s" % (m, summ([math.log(A[(s, m)][0] / B[(s, m)][0]) for s in SEEDS])))

    def beta(D, s, lo, hi):
        return math.log(D[(s, hi)][0] / D[(s, lo)][0]) / math.log(D[(s, hi)][1] / D[(s, lo)][1])
    for lo, hi in ((9, 50), (50, 100), (9, 100)):
        b1 = [beta(B, s, lo, hi) for s in SEEDS]
        b2 = [beta(A, s, lo, hi) for s in SEEDS]
        print("- β %d→%d mT: α 0.1 %.2f ± %.2f, α 0.01 %.2f ± %.2f; Δβ %s" %
              (lo, hi, np.mean(b1), np.std(b1, ddof=1) / 2, np.mean(b2), np.std(b2, ddof=1) / 2,
               summ(np.array(b2) - np.array(b1))))
    bal = [A[k][2] for k in A]
    print("\nα 0.01: stages with balance outside 0.9 … 1.1: %d of %d" % (sum(abs(b - 1) > 0.1 for b in bal), len(bal)))


if __name__ == "__main__":
    main()
