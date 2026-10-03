#!/usr/bin/env python
"""
beta(B_peak) with a valid 9 mT value (PROTOKOLL 7.15).

The 9 mT value comes from the relaxation test (X*: virgin state relaxed to
convergence, 9 and 20 mT). The values 20 ... 150 mT come from the runs with the
same seeds and settings but the old relaxation (F_base_f_* for alpha 0.1,
A010_* for alpha 0.01). In those runs only the first stage (9 mT) is disturbed
(PROTOKOLL 7.13); they are read with ALLOW_UNCONVERGED and their 9 mT stage is
dropped. Check of the splice: ln w at 20 mT of X* vs the old run, same seed.

Per realisation: beta per segment (adjacent amplitudes), power law over
9 ... 150 mT, and a quadratic fit ln w = c + beta0 x + k x^2 / 2 (x = ln b/b_ref),
beta(b) = beta0 + k x. Mean +- SE over the realisations (n given).

    python eval_relaxfix.py > results/relaxfix_beta.txt
"""
import math
import pathlib

import numpy as np

import analyze

RUNS = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = (901, 902, 903)
SETS = {0.1: ("F_base_f_s%d", "X100_s%d"), 0.01: ("A010_s%d", "X010_s%d")}
B9 = 0.008                       # b below this is the 9 mT stage (b = 0.006)
JS_T, KD, FREQ = 1.5, 895e3, 30e6
REF_MT = (10.0, 50.0, 100.0)


def mean_se(v):
    v = np.asarray(v, dtype=float)
    return v.mean(), (v.std(ddof=1) / math.sqrt(len(v)) if len(v) > 1 else float("nan")), len(v)


def fmt(v, nd=2):
    m, s, n = mean_se(v)
    return "%.*f ± %.*f" % (nd, m, nd, s)


def spliced(alpha, seed):
    old, new = (RUNS / (p % seed) for p in SETS[alpha])
    analyze.ALLOW_UNCONVERGED = True
    _, r_old = analyze.load_run(old)
    analyze.ALLOW_UNCONVERGED = False
    _, r_new = analyze.load_run(new)                       # converged: passes the gate
    r9 = [r for r in r_new if r["b_peak"] < B9]
    r20n = [r for r in r_new if r["b_peak"] > B9]
    rows = r9 + [r for r in r_old if r["b_peak"] > B9]
    r20o = [r for r in r_old if B9 < r["b_peak"] < 0.02]
    check = math.log(r20n[0]["w_loop"] / r20o[0]["w_loop"]) if r20n and r20o else float("nan")
    return rows, check


def fits(rows):
    b = np.array([r["b_peak"] for r in rows])
    lw = np.log([r["w_loop"] for r in rows])
    sw = np.array([r["lnw_err"] for r in rows])
    lb = np.log(b)
    seg = [(lw[i + 1] - lw[i]) / (lb[i + 1] - lb[i]) for i in range(len(b) - 1)]
    W = 1.0 / sw ** 2
    X1 = np.c_[np.ones_like(lb), lb]
    c1 = np.linalg.solve(X1.T @ (X1 * W[:, None]), X1.T @ (W * lw))
    chi1 = float(np.sum(W * (lw - X1 @ c1) ** 2) / (len(b) - 2))
    out = {"seg": seg, "power": float(c1[1]), "chi2_power": chi1}
    for mT in REF_MT:
        x = lb - math.log(mT / 1000.0 / JS_T)
        X2 = np.c_[np.ones_like(x), x, 0.5 * x ** 2]
        c2 = np.linalg.solve(X2.T @ (X2 * W[:, None]), X2.T @ (W * lw))
        out["beta@%g" % mT] = float(c2[1])
        out["curv"] = float(c2[2])
        out["P@%g" % mT] = math.exp(c2[0]) * KD * FREQ / 1e6     # W/cm^3 at b_ref
    return out


def main():
    print("beta(B_peak) with a valid 9 mT value: 9 mT from X* (converged relaxation), 20 ... 150 mT from")
    print("F_base_f_* (alpha 0.1) / A010_* (alpha 0.01), same seeds 901 ... 903. Mean ± SE over the realisations.")
    print("d/l_ex 300 (1 um), dx 3, L_eff 12 l_ex, r_p 0, phi 0.65, f = 30 MHz. B = B_peak.\n")
    for alpha in SETS:
        res, chk, bs = [], [], None
        for s in SEEDS:
            rows, c = spliced(alpha, s)
            res.append(fits(rows))
            chk.append(c)
            bs = [r["b_peak"] * JS_T * 1000 for r in rows]
        n = len(res)
        print("alpha = %g   (n = %d realisations)" % (alpha, n))
        print("  splice check, ln w(20 mT) X* - old run: %s  (single values %s)"
              % (fmt(chk, 3), " ".join("%+.3f" % c for c in chk)))
        print("  beta per segment:")
        for i in range(len(bs) - 1):
            print("    %5.0f -> %5.0f mT   %s" % (bs[i], bs[i + 1], fmt([r["seg"][i] for r in res])))
        print("  power law 9 ... 150 mT:  beta = %s   (chi2/dof %.1f ... %.1f)"
              % (fmt([r["power"] for r in res]), min(r["chi2_power"] for r in res),
                 max(r["chi2_power"] for r in res)))
        print("  quadratic fit: curvature d beta / d ln B = %s" % fmt([r["curv"] for r in res]))
        for mT in REF_MT:
            lp = np.log([r["P@%g" % mT] for r in res])
            m, se, _ = mean_se(lp)
            print("    B = %3.0f mT: beta = %s   P = %.3g W/cm^3 (x/ %.2f, geometric mean, factor = exp(SE))"
                  % (mT, fmt([r["beta@%g" % mT] for r in res]), math.exp(m), math.exp(se)))
        print()
    d = {a: [fits(spliced(a, s)[0]) for s in SEEDS] for a in SETS}
    print("alpha 0.01 - 0.1 (Welch, n = 3 + 3):")
    for k in ["power"] + ["beta@%g" % mT for mT in REF_MT]:
        a, b = (np.array([r[k] for r in d[x]]) for x in (0.01, 0.1))
        diff = a.mean() - b.mean()
        se = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
        print("  %-10s %+.2f ± %.2f  (%.1f SE)" % (k, diff, se, abs(diff) / se))


if __name__ == "__main__":
    main()
