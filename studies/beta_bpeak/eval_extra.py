#!/usr/bin/env python
"""
Further local evaluations (PROTOKOLL 7.27):
  A  ring-down (PH_*): fit p(t) = p_inf + A exp(-t / tau) over the hold; drift of M_par during the hold
  B  Debye fit of the frequency test: w per cycle = c * x / (1 + x^2), x = omega tau, from f/2, f, 2 f
     per seed and amplitude (2 parameters, 3 points: 1 degree of freedom; tau is only roughly determined)
  C  events per mesh: kept cycles with w_dis/w_loop outside 0.9 ... 1.1, and drops of the number of
     large-angle pairs by >= 20 % between cycles, for mesh300 (dx 1.5 vs dx 3) and the baseline (dx 3)

    python eval_extra.py > results/extra.md
"""
import json
import math
import pathlib

import numpy as np
from scipy.optimize import curve_fit

RUNS = pathlib.Path(__file__).resolve().parent / "runs"


def part_a():
    print("**A. Ring-down fit** p(t) = p_inf + A exp(-t/τ) over the 2 periods of the hold (T = 33.3 ns)\n")
    print("| run | B | p(0)/p_cycle | τ [ns] | τ/T | p_inf/p_cycle | ΔM_par during the hold / M amplitude |")
    print("|---|---|---|---|---|---|---|")
    for s in (921, 922):
        d = RUNS / ("PH_s%d" % s)
        T = 1.0 / json.loads((d / "config.json").read_text())["freq"]
        for f in sorted(d.glob("ringdown_*.csv")):
            stage = f.stem.split("_", 2)[2]
            a = np.loadtxt(f, delimiter=",", comments="#")
            t, M, p = a[:, 0], a[:, 2], a[:, 3]
            c = np.loadtxt(next(d.glob("samples_*_%s.csv" % stage)), delimiter=",", comments="#")
            last = c[c[:, 0] == c[:, 0].max()][:-1]
            pm = last[:, 8].mean()
            Mamp = 0.5 * (last[:, 3].max() - last[:, 3].min())
            B = 1000.0 * np.abs(last[:, 7]).max()
            try:
                (pinf, A, tau), _ = curve_fit(lambda t, a, b, c: a + b * np.exp(-t / c), t, p,
                                              p0=(p[-1], p[0] - p[-1], 0.2 * T), maxfev=20000)
            except RuntimeError:
                pinf = A = tau = float("nan")
            print("| PH_s%d | %.0f mT | %.2f | %.1f | %.2f | %.3f | %+.3f |" %
                  (s, B, p[0] / pm, tau * 1e9, tau / T, pinf / pm, (M[-1] - M[0]) / Mamp))


def w_of(name, mT):
    import analyze
    _, rows = analyze.load_run(RUNS / name)
    r = min(rows, key=lambda r: abs(r["b_peak"] * 1500 - mT))
    return r["w_loop"]


def part_b():
    print("\n**B. Debye fit of the frequency test** w = c x/(1 + x²), x = ωτ (f/2, f, 2f = 15, 30, 60 MHz)\n")
    print("| seed | B | τ [ns] | f_relax = 1/(2πτ) [MHz] | residual rms of ln w |")
    print("|---|---|---|---|---|")
    F = np.array([15e6, 30e6, 60e6])
    taus = {}
    for mT in (9, 20):
        for s in (921, 922, 923):
            w = np.array([w_of("F9_half_s%d" % s, mT), w_of("B300_s%d" % s, mT), w_of("F9_dbl_s%d" % s, mT)])
            om = 2 * math.pi * F

            def model(om, lc, ltau):
                x = om * math.exp(ltau)
                return lc + np.log(x / (1 + x * x))
            try:
                (lc, ltau), _ = curve_fit(model, om, np.log(w), p0=(math.log(w.max() * 2), math.log(5e-9)), maxfev=20000)
                tau = math.exp(ltau)
                res = float(np.sqrt(np.mean((np.log(w) - model(om, lc, ltau)) ** 2)))
            except RuntimeError:
                tau, res = float("nan"), float("nan")
            taus.setdefault(mT, []).append(tau)
            print("| %d | %d mT | %.2f | %.0f | %.3f |" % (s, mT, tau * 1e9, 1 / (2 * math.pi * tau) / 1e6, res))
    for mT, v in taus.items():
        v = np.array(v)
        print("\n%d mT: τ = %.2f ± %.2f ns (mean ± SE, n = 3) → f_relax ≈ %.0f MHz" %
              (mT, v.mean() * 1e9, v.std(ddof=1) / math.sqrt(3) * 1e9, 1 / (2 * math.pi * v.mean()) / 1e6))


def events(name):
    s = json.loads((RUNS / name / "summary.json").read_text())
    nb, nd, nc = 0, 0, 0
    for st in s["stages"]:
        cy = st["cycles"]
        for i, c in enumerate(cy):
            if i >= 2:
                nc += 1
                if abs(c["w_dis"] / c["w_loop"] - 1) > 0.1:
                    nb += 1
            if i >= 1 and c["n_pairs_gt60"] < 0.8 * cy[i - 1]["n_pairs_gt60"]:
                nd += 1
    return nb, nd, nc


def part_c():
    print("\n**C. Events per mesh** (kept cycles with balance outside 0.9 … 1.1; drops of the number of large-angle "
          "pairs by ≥ 20 % from one cycle to the next, all cycles)\n")
    print("| set | runs | kept cycles | cycles with balance outside 0.9 … 1.1 | n60 drops |")
    print("|---|---|---|---|---|")
    for lab, names in (("mesh300 dx 1.5", sorted(p.name for p in RUNS.glob("M300_dx15_*"))),
                       ("mesh300 dx 3", sorted(p.name for p in RUNS.glob("M300_dx3_*"))),
                       ("baseline dx 3", ["B300_s%d" % s for s in range(921, 929)])):
        tot = np.array([events(n) for n in names]).sum(0)
        print("| %s | %d | %d | %d (%.0f %%) | %d |" % (lab, len(names), tot[2], tot[0], 100 * tot[0] / tot[2], tot[1]))


def main():
    part_a()
    part_b()
    part_c()


if __name__ == "__main__":
    main()
