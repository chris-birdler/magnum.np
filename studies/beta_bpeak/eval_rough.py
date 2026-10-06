#!/usr/bin/env python
"""
Surface roughness with the converged virgin state (PROTOKOLL 7.29). d/l_ex 300, base point, seeds 921 ... 924,
one job per amplitude (9 / 50 / 100 mT), paired by seed.
  alpha 0.1 : A = staircase 10 nm (M300_dx3_*), B = staircase 5 nm (M300_dx15_*), C = smoothed 10 nm (R300f_*)
  alpha 0.01: staircase 10 nm (R300v01_*) vs smoothed 10 nm (R300f01_*), 9 cycles
Δ ln w (paired, mean ± SE) and β of the segments 9→50, 50→100, 9→100 mT.

    python eval_rough.py > results/rough.md
"""
import math
import pathlib

import numpy as np

import analyze

RUNS = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = (921, 922, 923, 924)
AMPS = (9, 50, 100)
SETS = {0.1: {"A staircase 10 nm": "M300_dx3_s%d_b%d", "B staircase 5 nm": "M300_dx15_s%d_b%d",
              "C smoothed 10 nm": "R300f_s%d_b%d"},
        0.01: {"staircase 10 nm": "R300v01_s%d_b%d", "smoothed 10 nm": "R300f01_s%d_b%d"}}


def one(name):
    _, rows = analyze.load_run(RUNS / name)
    return rows[0]["w_loop"], rows[0]["b_peak"], rows[0]["balance"]


def ms(v):
    v = np.asarray(v, float)
    return v.mean(), v.std(ddof=1) / math.sqrt(len(v))


def beta(W, s, lo, hi):
    return math.log(W[(s, hi)][0] / W[(s, lo)][0]) / math.log(W[(s, hi)][1] / W[(s, lo)][1])


def main():
    for alpha, sets in SETS.items():
        if not all((RUNS / (p % (s, m)) / "DONE").exists() for p in sets.values() for s in SEEDS for m in AMPS):
            print("**α %g: not complete yet**\n" % alpha)
            continue
        D = {k: {(s, m): one(p % (s, m)) for s in SEEDS for m in AMPS} for k, p in sets.items()}
        keys = list(sets)
        ref = keys[0]
        print("**α %g** (reference: %s; paired over seeds %s, mean ± SE)\n" % (alpha, ref, ", ".join(map(str, SEEDS))))
        print("| comparison | Δ ln w 9 mT | Δ ln w 50 mT | Δ ln w 100 mT | Δβ 9→50 | Δβ 50→100 | Δβ 9→100 |")
        print("|---|---|---|---|---|---|---|")
        for k in keys[1:]:
            cells = []
            for m in AMPS:
                cells.append("%+.3f ± %.3f" % ms([math.log(D[k][(s, m)][0] / D[ref][(s, m)][0]) for s in SEEDS]))
            for lo, hi in ((9, 50), (50, 100), (9, 100)):
                cells.append("%+.2f ± %.2f" % ms([beta(D[k], s, lo, hi) - beta(D[ref], s, lo, hi) for s in SEEDS]))
            print("| %s − %s | %s |" % (k, ref, " | ".join(cells)))
        print("\n| set | β 9→50 | β 50→100 | β 9→100 | balance range |")
        print("|---|---|---|---|---|")
        for k in keys:
            b = [ms([beta(D[k], s, lo, hi) for s in SEEDS]) for lo, hi in ((9, 50), (50, 100), (9, 100))]
            bal = [D[k][(s, m)][2] for s in SEEDS for m in AMPS]
            print("| %s | %.2f ± %.2f | %.2f ± %.2f | %.2f ± %.2f | %.2f … %.2f |" %
                  (k, *b[0], *b[1], *b[2], min(bal), max(bal)))
        print()
    print("d/l_ex 300 (1 µm), base point, converged virgin state; α 0.1: 5 cycles, α 0.01: 9 cycles; cycles kept "
          "from 2 on. The α 0.01 runs start from the relaxed states of the α 0.1 runs.")


if __name__ == "__main__":
    main()
