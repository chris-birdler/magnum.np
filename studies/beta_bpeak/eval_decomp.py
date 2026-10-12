#!/usr/bin/env python
"""
Loss decomposition (option A of Chris, 2026-10-12; method fixed before the result):
    P(B) = a1 · x + a2 · x² + a3 · x³,   x = B / 50 mT,   a_j ≥ 0
per design run (116 runs, steady rule ≥ 2 calm cycles with reruns, P in W/cm³ at the 7 nominal amplitudes), fitted
with relative errors (non-negative least squares of (model − P) / P). Physical reading of the terms (hypotheses to test):
    x¹  switching: every irreversible event dissipates a fixed energy, their number grows about ∝ B
    x²  linear viscous: reversible motion with a constant phase lag (damping), β = 2
    x³  Rayleigh: progressive depinning of the magnetisation from the pinning by the Herzer anisotropy, β = 3
Outputs:
- fit quality against the single power law (rms of ln P residuals), how often each term is determined (a_j > 2 SE,
  SE from the residual covariance of the unconstrained weighted fit);
- share of each term in P at 9 / 50 / 150 mT, for all runs and for the low-A and high-A thirds;
- factor effects on the shares (OLS, HC3, range effects as in the design analysis, with the lnd × lnQ term);
- damping check (centre C0 … C7, α 0.1 against α 0.01, steady rule ≥ 2; amplitudes missing at α 0.01 are left out):
  if the x² term is viscous, a2 scales with α (× 0.1), a1 and a3 do not.

    python eval_decomp.py > results/loss_decomposition.md      # also results/loss_decomposition.png
"""
import json
import math
import pathlib

import numpy as np
from scipy.optimize import nnls

import analyze
import analyze_design as AD

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE / "runs"
F_HZ = 30e6
B0 = 50.0
AMPS = np.array(AD.AMPS_MT)
REPORT_B = (9.0, 50.0, 150.0)
TERMS = ("x¹ switching", "x² viscous", "x³ Rayleigh")


def fit(B, P):
    """NNLS of the relative residual; returns coefficients a (W/cm³ at 50 mT), their SE (unconstrained), rms of ln
    residuals."""
    ok = np.isfinite(P) & (P > 0)
    B, P = B[ok], P[ok]
    x = B / B0
    X = np.column_stack([x, x ** 2, x ** 3]) / P[:, None]
    a, _ = nnls(X, np.ones(len(P)))
    model = (np.column_stack([x, x ** 2, x ** 3]) @ a)
    r = np.log(model / P)
    dof = max(len(P) - 3, 1)
    s2 = float(((X @ a - 1) ** 2).sum()) / dof
    try:
        se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    except np.linalg.LinAlgError:
        se = np.full(3, np.nan)
    return a, se, float(np.sqrt(np.mean(r ** 2))), len(P)


def shares(a, B):
    x = np.asarray(B, float) / B0
    t = np.array([a[0] * x, a[1] * x ** 2, a[2] * x ** 3])
    return t / t.sum(0)


def powerlaw_rms(B, P):
    ok = np.isfinite(P) & (P > 0)
    c = np.polyfit(np.log(B[ok]), np.log(P[ok]), 1)
    return float(np.sqrt(np.mean((np.polyval(c, np.log(B[ok])) - np.log(P[ok])) ** 2)))


def ols(y, D):
    import statsmodels.api as sm
    return sm.OLS(y, D).fit(cov_type="HC3")


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    recs, _ = AD.load_runs(RUNS, reruns=True)
    kd = json.loads((RUNS / recs[0]["name"] / "config.json").read_text())["units"]["Kd"]
    to_P = kd * F_HZ / 1e6
    P = np.array([[math.exp(r["lnw@%gmT" % m]) * to_P for m in AMPS] for r in recs])
    res = [fit(AMPS, p) for p in P]
    A = np.array([x[0] for x in res])
    SE = np.array([x[1] for x in res])
    rms = np.array([x[2] for x in res])
    rms_pl = np.array([powerlaw_rms(AMPS, p) for p in P])
    S = np.array([shares(a, REPORT_B) for a in A])           # [run, term, B]
    det = A > 2 * SE
    Areal = np.array([10.0 * ((1000.0 / math.exp(r["lnd"])) / 3.342) ** 2 for r in recs])
    q1, q2 = np.quantile(Areal, [1 / 3, 2 / 3])
    X = AD.coded(recs)
    D = np.column_stack([np.ones(len(X)), X, X[:, 0] * X[:, 3]])       # 1, lnQ, r_p, phi, lnd, lnd × lnQ

    L = ["# Loss decomposition P = a1·x + a2·x² + a3·x³ (x = B / 50 mT), 116 design runs", "",
         "Steady loops (≥ 2 calm cycles) with reruns, coarse mesh, P in W/cm³ per core volume at 30 MHz. Non-negative "
         "fit of the relative residual. Term names are hypotheses: x¹ switching, x² linear viscous, x³ Rayleigh.", ""]
    L += ["## Fit quality", "",
          "- rms of the ln P residual: 3 terms median %.3f (quartiles %.3f … %.3f); single power law median %.3f "
          "(%.3f … %.3f)." % (np.median(rms), *np.percentile(rms, [25, 75]), np.median(rms_pl), *np.percentile(rms_pl, [25, 75])),
          "- Terms present (a_j > 0) / determined (a_j > 2 SE): %s." % ", ".join(
              "%s %d / %d" % (n, (A[:, j] > 0).sum(), det[:, j].sum()) for j, n in enumerate(TERMS)), ""]
    L += ["## Share of each term in P (median over runs; quartiles)", "",
          "| group | term | 9 mT | 50 mT | 150 mT |", "|---|---|---|---|---|"]
    groups = (("all", np.ones(len(recs), bool)), ("low A (%.0f … %.1f pJ/m)" % (Areal.min(), q1), Areal <= q1),
              ("high A (%.1f … %.0f pJ/m)" % (q2, Areal.max()), Areal > q2))
    for gname, sel in groups:
        for j, n in enumerate(TERMS):
            L.append("| %s | %s | %s |" % (gname, n, " | ".join("%.2f (%.2f … %.2f)" % (np.median(S[sel, j, k]),
                                                                                          *np.percentile(S[sel, j, k], [25, 75]))
                                                              for k in range(3))))
    L += ["", "## Factor effects on the shares (OLS, HC3; range effect = 2 × coefficient; A 10 → 20 pJ/m is − lnd)", "",
          "| share | B | K_eff (lnQ) | r_p | φ | A (− lnd) | lnd × lnQ | R² |", "|---|---|---|---|---|---|---|---|"]
    eff = {}
    for j, n in enumerate(TERMS):
        for k, b in enumerate(REPORT_B):
            m = ols(S[:, j, k], D)
            c = 2 * m.params
            e = 2 * m.bse
            eff[(j, b)] = (c, e)
            L.append("| %s | %g mT | %+.2f ± %.2f | %+.2f ± %.2f | %+.2f ± %.2f | %+.2f ± %.2f | %+.2f ± %.2f | %.2f |" % (
                n, b, c[1], e[1], c[2], e[2], c[3], e[3], -c[4], e[4], c[5], e[5], m.rsquared))
    # damping check
    L += ["", "## Damping check (centre C0 … C7, steady loops ≥ 2 at both α)", ""]
    pairs = []
    for i in range(8):
        a1d, a0d = RUNS / ("C%d_a001" % i), RUNS / ("C%d" % i)
        x1, x0 = AD.run_record(a1d), AD.run_record(AD.rerun_dir(a0d) or a0d)
        if x1 is None or x0 is None:
            continue
        p1 = np.array([math.exp(x1["lnw@%gmT" % m]) * to_P if np.isfinite(x1["lnw@%gmT" % m]) else np.nan for m in AMPS])
        p0 = np.array([math.exp(x0["lnw@%gmT" % m]) * to_P if np.isfinite(x0["lnw@%gmT" % m]) else np.nan for m in AMPS])
        ok = np.isfinite(p1) & np.isfinite(p0)
        if ok.sum() < 4:
            continue
        f0, f1 = fit(AMPS[ok], p0[ok]), fit(AMPS[ok], p1[ok])
        pairs.append((i, ok.sum(), f0[0], f1[0]))
    L.append("| seed | amplitudes | a1 α0.1 → α0.01 | a2 α0.1 → α0.01 | a3 α0.1 → α0.01 |")
    L.append("|---|---|---|---|---|")
    for i, n, c0, c1 in pairs:
        L.append("| C%d | %d | %.1f → %.1f | %.1f → %.1f | %.1f → %.1f |" % (i, n, c0[0], c1[0], c0[1], c1[1], c0[2], c1[2]))
    if pairs:
        s0 = np.array([p[2] for p in pairs]).sum(0)
        s1 = np.array([p[3] for p in pairs]).sum(0)
        L.append("")
        L.append("Sum over the seeds, ratio α 0.01 / α 0.1: a1 × %.2f, a2 × %.2f, a3 × %.2f (viscous expectation for a2: "
                 "× 0.1)." % tuple(s1 / np.where(s0 > 0, s0, np.nan)))
    print("\n".join(L))
    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Bg = np.geomspace(9, 150, 40)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
    for ax, (gname, sel) in zip(axs[:2], groups[1:]):
        sh = np.array([shares(a, Bg) for a in A[sel]])
        med = np.median(sh, 0)
        ax.stackplot(Bg, med / med.sum(0), labels=TERMS, colors=["#1F4E79", "#B0B0B0", "#9B111E"], alpha=0.85)
        ax.set_xscale("log")
        ax.set_xticks([9, 20, 50, 100, 150])
        ax.set_xticklabels(["9", "20", "50", "100", "150"])
        ax.set_xlabel("B (mT)")
        ax.set_ylabel("share of P (median of the runs, normalised)")
        ax.set_title(gname)
        ax.legend(loc="lower right", fontsize=9)
    ax = axs[2]
    ax.scatter(Areal, S[:, 2, 1], c=[r["lnQ"] for r in recs], cmap="viridis", s=18)
    ax.set_xlabel("A at d = 1 µm (pJ/m)")
    ax.set_ylabel("share of the x³ term at 50 mT")
    ax.set_title("colour: ln Q (K_eff)")
    fig.tight_layout()
    fig.savefig(HERE / "results" / "loss_decomposition.png", dpi=130)


if __name__ == "__main__":
    main()
