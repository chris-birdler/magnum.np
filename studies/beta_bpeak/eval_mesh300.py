#!/usr/bin/env python
"""
mesh300 (PLAN 3e): dx 1.5 vs dx 3 at d/l_ex 300 with the converged virgin state, seeds 921 ... 924,
one job per amplitude (9 / 50 / 100 mT; 50 and 100 mT start from the init.pt of the 9 mT job of the
same mesh), 5 cycles, kept from cycle 2 (PROTOKOLL 7.25).
  A: ln w per seed, mesh error Δ ln w = ln w(dx 3) − ln w(dx 1.5), paired over the seeds (mean ± SE)
  B: β of the segments 9→50 and 50→100 mT per seed and mesh, paired Δβ
  C: rule of PLAN step 1b at 50 mT: compare with the clean d/l_ex 212 value (7.10: −0.06 ± 0.14)
  D: protocol check: mesh300 dx 3 (one job per amplitude from the virgin state) vs baseline B300 (the
     same seeds, amplitudes in sequence) at 9 / 50 / 100 mT

    python eval_mesh300.py > results/mesh300.md
"""
import math
import pathlib

import numpy as np

import analyze

RUNS = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = (921, 922, 923, 924)
AMPS = (9, 50, 100)


def one(name):
    _, rows = analyze.load_run(RUNS / name)
    r = rows[0]
    return r["w_loop"], r["b_peak"], r["balance"], r["lnw_err"]


def ms(v):
    v = np.asarray(v, float)
    return v.mean(), v.std(ddof=1) / math.sqrt(len(v))


def main():
    W = {}
    print("**A. Loss per seed and mesh** (w = W/K_d per cycle; balance w_dis/w_loop)\n")
    print("| seed | B | w dx 1.5 | w dx 3 | Δ ln w (dx 3 − dx 1.5) | balance dx 1.5 / dx 3 |")
    print("|---|---|---|---|---|---|")
    for s in SEEDS:
        for mT in AMPS:
            a = one("M300_dx15_s%d_b%d" % (s, mT))
            b = one("M300_dx3_s%d_b%d" % (s, mT))
            W[(s, mT, 1.5)], W[(s, mT, 3.0)] = a, b
            print("| %d | %d mT | %.4g | %.4g | %+.3f | %.2f / %.2f |" %
                  (s, mT, a[0], b[0], math.log(b[0] / a[0]), a[2], b[2]))
    print("\n| B | Δ ln w (paired, n = 4) | ratio w(dx 3)/w(dx 1.5) |")
    print("|---|---|---|")
    D50 = None
    for mT in AMPS:
        d = [math.log(W[(s, mT, 3.0)][0] / W[(s, mT, 1.5)][0]) for s in SEEDS]
        m, se = ms(d)
        if mT == 50:
            D50 = (m, se)
        print("| %d mT | %+.3f ± %.3f | %.2f |" % (mT, m, se, math.exp(m)))

    print("\n**B. β of the segments, per mesh (mean ± SE over the seeds) and paired Δβ**\n")
    print("| segment | β dx 1.5 | β dx 3 | Δβ (dx 3 − dx 1.5), paired |")
    print("|---|---|---|---|")
    for lo, hi in ((9, 50), (50, 100), (9, 100)):
        bb = {}
        for dx in (1.5, 3.0):
            bb[dx] = [math.log(W[(s, hi, dx)][0] / W[(s, lo, dx)][0]) / math.log(W[(s, hi, dx)][1] / W[(s, lo, dx)][1])
                      for s in SEEDS]
        d = np.array(bb[3.0]) - np.array(bb[1.5])
        print("| %d→%d mT | %.2f ± %.2f | %.2f ± %.2f | %+.2f ± %.2f |" %
              (lo, hi, *ms(bb[1.5]), *ms(bb[3.0]), *ms(d)))

    print("\n**C. Rule of PLAN step 1b at 50 mT** (same mesh error at both ends of d/l_ex 212 … 300)\n")
    ref = (-0.06, 0.14)
    diff = D50[0] - ref[0]
    se = math.hypot(D50[1], ref[1])
    print("d/l_ex 300: %+.3f ± %.3f; d/l_ex 212 (7.10, clean): %+.2f ± %.2f; difference %+.3f ± %.3f; "
          "|difference| + 1.645 SE = %.3f (rule: < 0.10)" % (*D50, *ref, diff, se, abs(diff) + 1.645 * se))

    print("\n**D. Protocol check, dx 3:** one job per amplitude from the virgin state (mesh300) vs amplitudes in "
          "sequence (baseline B300, same seeds)\n")
    print("| B | Δ ln w (mesh300 − baseline), paired n = 4 |")
    print("|---|---|")
    for mT in AMPS:
        d = []
        for s in SEEDS:
            _, rows = analyze.load_run(RUNS / ("B300_s%d" % s))
            r = min(rows, key=lambda r: abs(r["b_peak"] * 1500 - mT))
            d.append(math.log(W[(s, mT, 3.0)][0] / r["w_loop"]))
        print("| %d mT | %+.3f ± %.3f |" % (mT, *ms(d)))
    print("\nd/l_ex 300 (1 µm), base point, α 0.1, f = 30 MHz, converged virgin state, cycles from 2 on (5 cycles per "
          "job). Errors: SE over the 4 seeds (paired differences).")


if __name__ == "__main__":
    main()
