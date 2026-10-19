#!/usr/bin/env python
"""
Test of the hypothesis (2026-10-04): loss = jump part (irreversible, Barkhausen-like, energy per jump
fixed, number of jumps ~ B -> w ~ B, alpha-independent) + damping part (reversible rotation, linear
viscous response -> w ~ B^2, ~ alpha):
    w(b) = c1 b + c2 b^2,  b = B_peak / Js,
weighted least squares on w (weight 1 / (w lnw_err)^2), per run, all valid runs (as eval_beta_segments).
Predictions: (i) the fit describes the data within the cycle errors; (ii) c2 changes strongly with
alpha; (iii) c1 is about the same at alpha 0.01 and 0.1.
For comparison the single power law w = c b^beta (also 2 parameters) on the same points.

    python eval_decompose.py > results/decompose.md
"""
import math

import numpy as np

import eval_beta_segments as G

KD, FREQ, JS_MT = 895e3, 30e6, 1500.0


def fit_run(rows):
    b = np.array([r["b_peak"] for r in rows])
    w = np.array([r["w_loop"] for r in rows])
    sw = w * np.array([r["lnw_err"] for r in rows])
    W = 1.0 / sw ** 2
    X = np.c_[b, b ** 2]
    cov = np.linalg.inv(X.T @ (X * W[:, None]))
    c = cov @ (X.T @ (W * w))
    chi_d = float(np.sum(W * (w - X @ c) ** 2)) / (len(b) - 2)
    Xp = np.c_[np.ones_like(b), np.log(b)]                     # power law, same weights in ln w
    Wl = 1.0 / np.array([r["lnw_err"] for r in rows]) ** 2
    cp = np.linalg.solve(Xp.T @ (Xp * Wl[:, None]), Xp.T @ (Wl * np.log(w)))
    chi_p = float(np.sum(Wl * (np.log(w) - Xp @ cp) ** 2)) / (len(b) - 2)
    share = {mT: c[0] * (mT / JS_MT) / (c[0] * (mT / JS_MT) + c[1] * (mT / JS_MT) ** 2) for mT in (10, 50, 100)}
    return c, np.sqrt(np.diag(cov)), chi_d, chi_p, share, len(b)


def ms(v):
    v = np.asarray(v, float)
    return v.mean(), v.std(ddof=1) / math.sqrt(len(v))


def main():
    res = {}
    print("| α | run | n points | c1 (jump, ∝ b) | c2 (damping, ∝ b²) | χ²/dof c1 b + c2 b² | χ²/dof power law "
          "| jump share of w at 10 / 50 / 100 mT |")
    print("|---|---|---|---|---|---|---|---|")
    for alpha in (0.1, 0.01):
        runs = G.group(alpha)
        R = []
        for i, rows in enumerate(runs):
            c, e, chi_d, chi_p, sh, n = fit_run(rows)
            R.append((c, e, chi_d, chi_p, sh))
            print("| %g | %d | %d | %.3g ± %.1g | %.3g ± %.1g | %.1f | %.1f | %.2f / %.2f / %.2f |" %
                  (alpha, i + 1, n, c[0], e[0], c[1], e[1], chi_d, chi_p, sh[10], sh[50], sh[100]))
        res[alpha] = R
    print("\n| α | n | c1 mean ± SE | c2 mean ± SE | median χ²/dof decomposition | median χ²/dof power law "
          "| jump share 10 / 50 / 100 mT |")
    print("|---|---|---|---|---|---|---|")
    for alpha, R in res.items():
        c1 = ms([r[0][0] for r in R])
        c2 = ms([r[0][1] for r in R])
        print("| %g | %d | %.3g ± %.2g | %.3g ± %.2g | %.1f | %.1f | %s |" %
              (alpha, len(R), *c1, *c2, np.median([r[2] for r in R]), np.median([r[3] for r in R]),
               " / ".join("%.2f ± %.2f" % ms([r[4][mT] for r in R]) for mT in (10, 50, 100))))
    (m1a, s1a), (m1b, s1b) = ms([r[0][0] for r in res[0.01]]), ms([r[0][0] for r in res[0.1]])
    (m2a, s2a), (m2b, s2b) = ms([r[0][1] for r in res[0.01]]), ms([r[0][1] for r in res[0.1]])
    q1 = m1a / m1b
    q2 = m2a / m2b
    print("\nRatio α 0.01 / α 0.1: c1 = %.2f ± %.2f, c2 = %.2f ± %.2f (error propagation of the SE)." %
          (q1, q1 * math.hypot(s1a / m1a, s1b / m1b), q2, q2 * math.hypot(s2a / m2a, s2b / m2b)))
    print("\nw in units of K_d per cycle; P = w K_d f. c1 b and c2 b² are the two parts of w. χ²/dof uses the "
          "cycle errors of ln w (both models have 2 parameters). Valid runs as eval_beta_segments: α 0.1 = 8 clean "
          "+ 3 spliced, α 0.01 = 3 clean + 3 spliced.")


if __name__ == "__main__":
    main()
