#!/usr/bin/env python
"""
Local Steinmetz exponent: power law w = c B^beta fitted over 3 neighbouring amplitudes (sliding
window), per run, weighted by the cycle error of ln w. 3 points = 1 degree of freedom: the fit
also tests whether the power law holds locally (chi^2). The window centre is the geometric mean
of its 3 amplitudes.
  Table A: common amplitudes 9 / 20 / 50 / 100 / 150 mT, all alpha (as eval_beta_segments).
  Table B: alpha 0.1, all 7 amplitudes 9 / 20 / 35 / 50 / 70 / 100 / 150 mT (finer windows).

    python eval_beta_windows.py > results/beta_windows.md
"""
import math

import numpy as np

import eval_beta_segments as G

ALL7 = (9, 20, 35, 50, 70, 100, 150)


def pick(rows, amps):
    out = []
    for mT in amps:
        r = min(rows, key=lambda r: abs(r["b_peak"] * 1500 - mT))
        out.append(r if abs(r["b_peak"] * 1500 / mT - 1) < 0.1 else None)
    return out


def window_fits(rows):
    """[(beta, se_cycle, chi2)] for each window of 3 neighbouring rows."""
    out = []
    for i in range(len(rows) - 2):
        w3 = rows[i:i + 3]
        if any(r is None for r in w3):
            out.append((np.nan, np.nan, np.nan))
            continue
        x = np.log([r["b_peak"] for r in w3])
        y = np.log([r["w_loop"] for r in w3])
        W = 1.0 / np.array([r["lnw_err"] for r in w3]) ** 2
        X = np.c_[np.ones(3), x]
        cov = np.linalg.inv(X.T @ (X * W[:, None]))
        c = cov @ (X.T @ (W * y))
        chi2 = float(np.sum(W * (y - X @ c) ** 2))
        out.append((float(c[1]), math.sqrt(cov[1, 1]), chi2))
    return out


def table(groups, amps):
    n_w = len(amps) - 2
    centres = [(amps[i] * amps[i + 1] * amps[i + 2]) ** (1 / 3) for i in range(n_w)]
    print("| α | n | " + " | ".join("β %d–%d–%d mT (≈ %.0f mT)" % (amps[i], amps[i + 1], amps[i + 2], centres[i])
                                  for i in range(n_w)) + " |")
    print("|---|---|" + "---|" * n_w)
    for alpha, runs, note in groups:
        R = [window_fits(pick(r, amps)) for r in runs]
        cells = []
        for i in range(n_w):
            b = np.array([r[i][0] for r in R])
            ok = np.isfinite(b)
            if ok.sum() < 2:
                cells.append("–")
                continue
            se_c = np.median([r[i][1] for r in R if np.isfinite(r[i][0])])
            chi = np.median([r[i][2] for r in R if np.isfinite(r[i][0])])
            cells.append("%.2f ± %.2f (cycle %.2f, χ² %.1f)" % (b[ok].mean(), b[ok].std(ddof=1) / math.sqrt(ok.sum()),
                                                            se_c, chi))
        print("| %g%s | %d | %s |" % (alpha, note, len(runs), " | ".join(cells)))


def main():
    A = [(0.1, G.group(0.1), ""), (0.03, G.group(0.03), " (no 9 mT)"), (0.01, G.group(0.01), "")]
    print("**A. Common amplitudes, all α** (3-point power-law fit per run, mean ± SE over the runs)\n")
    table(A, G.COMMON)
    print("\n**B. α 0.1, all 7 amplitudes** (B300_* and F_base_f_* + X100_*)\n")
    table([(0.1, G.group(0.1), "")], ALL7)
    print()
    print("d/l_ex 300 (1 µm), dx 3, L_eff 12 l_ex, r_p 0, φ 0.65, f = 30 MHz. ± = SE over the runs (realisation "
          "scatter). 'cycle' = median SE of β within one run from the cycle scatter of ln w. χ² = median of the fit "
          "χ² (1 degree of freedom; above ≈ 4 the 3 points are not on one power law within the cycle errors).")


if __name__ == "__main__":
    main()
