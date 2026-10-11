#!/usr/bin/env python
"""
Same-state mesh test at K3 (PROTOKOLL 7.63): dx 3 started from the coarse-grained relaxed dx 1.5 state ("dx 3, fine
state") against dx 1.5 (its own state) and dx 3 with its own relaxed state. Seeds 2031 … 2036, 9 / 50 / 100 mT,
steady-tail rule with ≥ 2 steady loops, finished reruns X_mc2 replace X.
Decision rule (fixed before the result, 7.63): mean ln[P(dx 3, fine state) / P(dx 1.5)] inside ±0.10 at 9, 50 and
100 mT → the dx 3 dynamics is correct and the mesh error comes from the relaxed state.

    python eval_cg.py > results/k3_same_state.md
"""
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import eval_k3fill as KF

R = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = range(2031, 2037)
AMPS = (9, 50, 100)
SETS = (("dx 1.5, own state", "MK3_dx15_s%d_b%d"), ("dx 3, fine state", "MK3cg_dx3_s%d_b%d"),
        ("dx 3, own state", "MK3_dx3_s%d_b%d"))


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    to_P = None
    L = {}
    for lab, fmt in SETS:
        for s in SEEDS:
            for b in AMPS:
                d = AD.rerun_dir(R / (fmt % (s, b))) or R / (fmt % (s, b))
                L[(lab, s, b)] = KF.loss(d)
                if to_P is None and (d / "config.json").exists():
                    import json
                    to_P = json.loads((d / "config.json").read_text())["units"]["Kd"] * 30e6 / 1e6
    print("# Same-state mesh test at K3 (A 20 pJ/m, d/l_ex 212)\n")
    print("Loss density P (W/cm³, geometric mean over the repeats), steady loops (≥ 2). n = repeats with a value.\n")
    print("| B (mT) | " + " | ".join(lab for lab, _ in SETS) + " |")
    print("|---|---|---|---|")
    for b in AMPS:
        cells = []
        for lab, _ in SETS:
            v = [math.log(L[(lab, s, b)][0] * to_P) for s in SEEDS if L[(lab, s, b)]]
            cells.append("%.1f (n %d)" % (math.exp(np.mean(v)), len(v)) if v else "–")
        print("| %d | %s |" % (b, " | ".join(cells)))
    print("\nPaired ratios per repeat (mean of ln ratio ± SE, as a factor; Welch-free paired t-test p):\n")
    print("| B (mT) | dx 3 fine state / dx 1.5 | dx 3 own state / dx 1.5 | dx 3 fine state / dx 3 own state |")
    print("|---|---|---|---|")
    verdict = []
    for b in AMPS:
        cells = []
        for a, c in ((SETS[1][0], SETS[0][0]), (SETS[2][0], SETS[0][0]), (SETS[1][0], SETS[2][0])):
            d = np.array([math.log(L[(a, s, b)][0] / L[(c, s, b)][0]) for s in SEEDS if L[(a, s, b)] and L[(c, s, b)]])
            if len(d) > 1:
                se = d.std(ddof=1) / math.sqrt(len(d))
                cells.append("× %.2f (ln %+.2f ± %.2f, n %d, p %.3f)" % (math.exp(d.mean()), d.mean(), se, len(d),
                                                                         stats.ttest_1samp(d, 0).pvalue))
                if a == SETS[1][0] and c == SETS[0][0]:
                    verdict.append(abs(d.mean()) <= 0.10)
            else:
                cells.append("n %d" % len(d))
        print("| %d | %s |" % (b, " | ".join(cells)))
    print("\nβ over 9 / 50 / 100 mT (power-law fit per repeat, mean ± SE):\n")
    for lab, _ in SETS:
        bt = []
        for s in SEEDS:
            pts = [L[(lab, s, b)] for b in AMPS]
            if all(pts):
                bt.append(np.polyfit(np.log([p[1] for p in pts]), np.log([p[0] for p in pts]), 1)[0])
        bt = np.array(bt)
        if len(bt) > 1:
            print("- %s: β %.2f ± %.2f (n %d)" % (lab, bt.mean(), bt.std(ddof=1) / math.sqrt(len(bt)), len(bt)))
    print("\nDecision rule (7.63): |mean ln ratio| ≤ 0.10 at all 3 amplitudes for dx 3 (fine state) / dx 1.5: %s"
          % ("met" if verdict and all(verdict) else "not met (%s)" % verdict))


if __name__ == "__main__":
    main()
