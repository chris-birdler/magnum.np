#!/usr/bin/env python
"""
Checks before the Sobol design (PROTOKOLL 7.24 / 7.26): d/l_ex 300, base point, alpha 0.1, converged
virgin state (init.pt of B300_*).
  A  frequency test at 9 and 20 mT: f/2 (F9_half_*), f (B300_*, first two stages), 2 f (F9_dbl_*),
     seeds 921 ... 923; exponent n of w per cycle ~ f^n (n = 0: rate-independent loss per cycle,
     n = 1: loss per cycle ~ f, viscous below resonance); beta(9 -> 20) at each f
  B  phase-resolved local response (phasemap.json of PH_s921, PH_s922; 16 snapshots per cycle, 3 cycles)
  C  ring-down (ringdown_*.csv of PH_*): H held at H max; decay of p_dis

    python eval_checks.py > results/checks.md
"""
import json
import math
import pathlib

import numpy as np

import analyze

RUNS = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = (921, 922, 923)
KD = 895e3


def stage_w(name, mT):
    _, rows = analyze.load_run(RUNS / name)
    r = min(rows, key=lambda r: abs(r["b_peak"] * 1500 - mT))
    return r["w_loop"], r["b_peak"], r["balance"]


def ms(v):
    v = np.asarray(v, float)
    return v.mean(), (v.std(ddof=1) / math.sqrt(len(v)) if len(v) > 1 else float("nan"))


def part_a():
    print("**A. Frequency test at small B** (w per cycle; n from a fit of ln w over ln f per seed; mean ± SE, n = 3)\n")
    print("| B | w f/2 | w f | w 2f | exponent n (w ~ f^n) | balance f/2 / f / 2f |")
    print("|---|---|---|---|---|---|")
    beta = {}
    for mT in (9, 20):
        ns, W = [], {0.5: [], 1.0: [], 2.0: []}
        bal = {0.5: [], 1.0: [], 2.0: []}
        for s in SEEDS:
            pts = {0.5: stage_w("F9_half_s%d" % s, mT), 1.0: stage_w("B300_s%d" % s, mT), 2.0: stage_w("F9_dbl_s%d" % s, mT)}
            x = np.log(list(pts))
            y = np.log([pts[k][0] for k in pts])
            ns.append(np.polyfit(x, y, 1)[0])
            for k in pts:
                W[k].append(pts[k][0])
                bal[k].append(pts[k][2])
                beta.setdefault((k, s), {})[mT] = pts[k][:2]
        print("| %d mT | %.3g | %.3g | %.3g | %.2f ± %.2f | %s |" %
              (mT, np.mean(W[0.5]), np.mean(W[1.0]), np.mean(W[2.0]), *ms(ns),
               " / ".join("%.2f…%.2f" % (min(bal[k]), max(bal[k])) for k in (0.5, 1.0, 2.0))))
    print("\n| f | β(9 → 20 mT), mean ± SE |")
    print("|---|---|")
    for k, lab in ((0.5, "f/2 (15 MHz)"), (1.0, "f (30 MHz)"), (2.0, "2f (60 MHz)")):
        b = [math.log(beta[(k, s)][20][0] / beta[(k, s)][9][0]) / math.log(beta[(k, s)][20][1] / beta[(k, s)][9][1])
             for s in SEEDS]
        print("| %s | %.2f ± %.2f |" % (lab, *ms(b)))


def part_b():
    print("\n**B. Phase-resolved local response** (16 snapshots per cycle, last 3 cycles; shares of the measured "
          "dissipation p_dis)\n")
    print("| run | B | R1 (local fundamental) | R_np (non-periodic) | local phase lag mean / spread [°] |")
    print("|---|---|---|---|---|")
    for s in (921, 922):
        f = RUNS / ("PH_s%d" % s) / "phasemap.json"
        if not f.exists():
            print("| PH_s%d | – | phasemap.json missing | | | | |" % s)
            continue
        for st in json.loads(f.read_text())["stages"]:
            print("| PH_s%d | %.0f mT | %.2f | %.3f | %.1f / %.1f |" %
                  (s, st["B_peak_mT"], st["R1"], st["R_np"],
                   math.degrees(st["lag_mean_rad"]), math.degrees(st["lag_spread_rad"])))
    print("\nThe harmonics 2 … 7 (R_harm in phasemap.json) are not usable: motion faster than T/16 aliases into "
          "them and the weight (k ω)² amplifies it (R_harm = 3.6 … 44 > 1 at 9 … 20 mT).")


def part_c():
    print("\n**C. Ring-down** (H held at H max after the last cycle; p_dis relative to the mean dissipation of "
          "the cycles before; E_rd = dissipated energy in the 2 periods of the hold / loss per cycle)\n")
    print("| run | B | p(0) | p(0.01 T) | p(0.05 T) | p(0.1 T) | p(0.5 T) | p(2 T) | E_rd / w_cycle |")
    print("|---|---|---|---|---|---|---|---|---|")
    for s in (921, 922):
        d = RUNS / ("PH_s%d" % s)
        cfg = json.loads((d / "config.json").read_text())
        T = 1.0 / cfg["freq"]
        for f in sorted(d.glob("ringdown_*.csv")):
            stage = f.stem.split("_", 2)[2]
            a = np.loadtxt(f, delimiter=",", comments="#")
            t, p = a[:, 0], a[:, 3]
            c = np.loadtxt(next(d.glob("samples_*_%s.csv" % stage)), delimiter=",", comments="#")
            last = c[c[:, 0] == c[:, 0].max()]
            pm = last[:-1, 8].mean()
            B = 1000.0 * np.abs(last[:, 7]).max()
            at = lambda x: p[np.argmin(np.abs(t - x * T))] / pm
            E = np.trapezoid(p, t) if hasattr(np, "trapezoid") else np.trapz(p, t)
            print("| PH_s%d | %.0f mT | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |" %
                  (s, B, at(0), at(0.01), at(0.05), at(0.1), at(0.5), at(t.max() / T), E / (pm * T)))


def main():
    part_a()
    part_b()
    part_c()
    print("\nd/l_ex 300, dx 3, base point, α 0.1, converged virgin state. A: B300 gives f; all runs start from the "
          "same init.pt per seed and run 9 then 20 mT (7 cycles, kept from cycle 2).")


if __name__ == "__main__":
    main()
