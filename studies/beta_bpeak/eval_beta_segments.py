#!/usr/bin/env python
"""
beta per segment and from the fits, for each alpha (valid runs, d/l_ex 300, dx 3, base point).
Segments on the common amplitudes 9 / 20 / 50 / 100 / 150 mT (all runs have them).
alpha 0.1 : B300_* (clean, 8) + F_base_f_* with 9 mT from X100_* (spliced, 3)
alpha 0.03: A030_* (9 mT not valid: only 20 ... 150 mT; no converged 9 mT run exists)
alpha 0.01: A010c_* (clean, 3) + A010_* with 9 mT from X010_* (spliced, 3)

    python eval_beta_segments.py > results/beta_segments.md
"""
import math

import numpy as np

import analyze
import eval_relaxfix as E

RUNS = E.RUNS
COMMON = (9, 20, 50, 100, 150)


def pick(rows):
    """rows at the common amplitudes (nearest B_peak)."""
    out = []
    for mT in COMMON:
        r = min(rows, key=lambda r: abs(r["b_peak"] * 1500 - mT))
        out.append(r if abs(r["b_peak"] * 1500 / mT - 1) < 0.1 else None)
    return out


def seg(rows):
    v = []
    for a, b in zip(rows[:-1], rows[1:]):
        v.append(math.log(b["w_loop"] / a["w_loop"]) / math.log(b["b_peak"] / a["b_peak"]) if a and b else np.nan)
    return v


def group(alpha):
    runs = []
    if alpha == 0.1:
        runs += [analyze.load_run(RUNS / ("B300_s%d" % s))[1] for s in range(921, 929)]
        runs += [E.spliced(0.1, s)[0] for s in E.SEEDS]
    elif alpha == 0.01:
        runs += [analyze.load_run(RUNS / ("A010c_s%d" % s))[1] for s in (904, 905, 906)]
        runs += [E.spliced(0.01, s)[0] for s in E.SEEDS]
    else:
        analyze.ALLOW_UNCONVERGED = True
        for s in (901, 902, 903):
            rows = analyze.load_run(RUNS / ("A030_s%d" % s))[1]
            runs.append([r for r in rows if r["b_peak"] * 1500 > 12])        # 9 mT not valid
        analyze.ALLOW_UNCONVERGED = False
    return runs


def ms(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if len(v) < 2:
        return "–"
    return "%.2f ± %.2f" % (v.mean(), v.std(ddof=1) / math.sqrt(len(v)))


def main():
    hdr = ["9→20", "20→50", "50→100", "100→150"]
    print("| α | n | " + " | ".join("β %s mT" % h for h in hdr) + " | β fit 9…150 mT | β fit 10 mT | β fit 50 mT | β fit 100 mT |")
    print("|---|---|" + "---|" * 8)
    for alpha in (0.1, 0.03, 0.01):
        runs = group(alpha)
        S = np.array([seg(pick(r)) for r in runs])
        F = [E.fits(r) for r in runs]
        fit = [ms([f[k] for f in F]) for k in ("power", "beta@10", "beta@50", "beta@100")]
        if alpha == 0.03:
            fit = [ms([f["power"] for f in F]) + " (20…150)", "–", ms([f["beta@50"] for f in F]),
                   ms([f["beta@100"] for f in F])]
        print("| %g | %d | %s | %s |" % (alpha, len(runs), " | ".join(ms(S[:, i]) for i in range(4)), " | ".join(fit)))
    print()
    print("d/l_ex 300 (1 µm), dx 3, L_eff 12 l_ex, r_p 0, φ 0.65, f = 30 MHz. Mean ± SE over the runs. "
          "Segment: β = Δ ln w / Δ ln B between two amplitudes of one run. Fit: power law over all amplitudes "
          "and quadratic fit at 10 / 50 / 100 mT per run (eval_relaxfix.fits). α 0.03: the 9 mT stage is not "
          "valid (not converged relaxation), the fit covers 20 … 150 mT.")


if __name__ == "__main__":
    main()
