#!/usr/bin/env python
"""
Further evaluations of the corner test, old against new (Chris 2026-10-18: evaluate all that is finished; descriptive,
no reading). The methods of the old main study on the 12 corner runs (new: <run>_nx; old: the partner of
eval_corners_nx.old_dir):
  decomposition  P(B) = a1 x + a2 x² + a3 x³, x = B / 50 mT, NNLS of the relative residual (eval_decomp.fit); the
                 shares of the three terms at 9, 50 and 150 mT;
  channels       R1 = fundamental share of the dissipation per stage (acc.json); β_f = β + d ln R1 / d ln B
                 (eval_channels: β = β_f − d ln R1 / d ln B);
  location       share of the loss in the hottest 10 % of the particle cells and in the outer shell r ≥ 0.8 of the
                 particles (acc.json, "where");
  jumps          share of the dissipation in jumps (p_dis residual above 5 robust sigma over a Fourier series up to
                 harmonic 8 of the drive phase, eval_jumps.analyse_cycle), mean over the cycles from 2 on.
Per quantity: mean over the 3 seeds per corner, old / new; and the paired difference new − old over the 12 runs
(mean ± SE).

    python eval_corners_fields.py > runs/corners_nx_fields.md
"""
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import conv_jobs as CJ
import eval_corners_nx as EN
import eval_decomp as ED
import eval_jumps as EJ
import units as U

R = pathlib.Path(__file__).resolve().parent / "runs"
B = np.asarray(AD.AMPS_MT, float)
X = np.log(B)
RUNS = ["%s_%d" % (c, s) for c in CJ.CORNERS_NX for s in range(3)]


def acc_stages(d):
    f = d / "acc.json"
    if not f.exists():
        return None
    a = json.loads(f.read_text())
    return a["stages"] if isinstance(a, dict) else a


def jumps(d):
    """Jump energy share per stage (mean over the cycles from 2 on)."""
    out = []
    for f in sorted(d.glob("samples_*.csv")):
        a = np.loadtxt(f, delimiter=",", comments="#")
        cyc = a[:, 0].astype(int)
        v = [EJ.analyse_cycle(a[cyc == k, 1], a[cyc == k, 8], 8, 5.0)[1] for k in sorted(set(cyc))[2:]]
        out.append(float(np.mean(v)) if v else float("nan"))
    return np.array(out[:7] + [float("nan")] * (7 - len(out)))


def record(d, lnF):
    lnP = EN.lnw(d) + lnF
    P = np.exp(lnP)
    a, _, rms, _ = ED.fit(B, P)
    sh = ED.shares(a, [9.0, 50.0, 150.0])
    acc = acc_stages(d)
    R1 = np.array([s.get("R1", np.nan) for s in acc]) if acc and len(acc) == 7 else np.full(7, np.nan)
    top = np.array([s["where"]["p top 10%"]["p_share"] for s in acc]) if acc and len(acc) == 7 else np.full(7, np.nan)
    shell = np.array([s["where"]["r>=0.8"]["p_share"] for s in acc]) if acc and len(acc) == 7 else np.full(7, np.nan)
    beta = float(np.polyfit(X, lnP, 1)[0]) if np.all(np.isfinite(lnP)) else float("nan")
    beta_f = beta + float(np.polyfit(X, np.log(R1), 1)[0]) if np.all(np.isfinite(R1)) else float("nan")
    return {"a": a, "sh9": sh[:, 0], "sh50": sh[:, 1], "sh150": sh[:, 2], "rms": rms, "R1": R1, "top": top,
            "shell": shell, "beta": beta, "beta_f": beta_f, "jump": jumps(d)}


def fmt(v, f="%.2f"):
    return " ".join(f % x if np.isfinite(x) else "–" for x in np.atleast_1d(v))


def main(root=R, arm=CJ.NX_SUFFIX):
    """arm: suffix of the new runs (default _nx, the corner test; _ad: option C, PROTOKOLL 7.118)."""
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    u = U.units()
    lnF = math.log(u["Kd"] * U.F_REL_DEFAULT * u["f_M"] / 1e6)
    D = {n: (record(EN.old_dir(root, n), lnF), record(root / (n + arm), lnF)) for n in RUNS}
    L = ["# Corner test: further evaluations, old against new (eval_corners_fields.py, descriptive)", "",
         "Old: DemagFieldPBC, atol 1e-5 (eval_corners_nx.old_dir); new: setup of PLAN 4.13. Mean over the 3 seeds of a "
         "corner; amplitudes 9, 20, 35, 50, 70, 100, 150 mT.", ""]
    for key, title, f in (("R1", "Fundamental share R1 of the dissipation (acc.json)", "%.2f"),
                          ("top", "Share of the loss in the hottest 10 % of the cells", "%.2f"),
                          ("shell", "Share of the loss in the outer shell r ≥ 0.8 (49 % of the volume)", "%.2f"),
                          ("jump", "Share of the dissipation in jumps (KH 8, 5 sigma)", "%.3f")):
        L += ["## " + title, "", "| corner | old | new |", "|---|---|---|"]
        for c in CJ.CORNERS_NX:
            o = np.nanmean([D["%s_%d" % (c, s)][0][key] for s in range(3)], 0)
            w = np.nanmean([D["%s_%d" % (c, s)][1][key] for s in range(3)], 0)
            L.append("| %s | %s | %s |" % (c, fmt(o, f), fmt(w, f)))
        d = np.array([D[n][1][key] - D[n][0][key] for n in RUNS])
        m = np.nanmean(d, 0)
        se = np.nanstd(d, 0, ddof=1) / np.sqrt(np.sum(np.isfinite(d), 0))
        L += ["", "Paired new − old (mean ± SE, n 12): " + "; ".join(
            "%g mT %+.3f ± %.3f" % (b, a, e) for b, a, e in zip(B, m, se)) + ".", ""]
    L += ["## Decomposition P = a1 x + a2 x² + a3 x³ (x = B / 50 mT): shares of the three terms", "",
          "| corner | at 9 mT old / new | at 50 mT old / new | at 150 mT old / new | rms of ln residual old / new |",
          "|---|---|---|---|---|"]
    for c in CJ.CORNERS_NX:
        rr = [D["%s_%d" % (c, s)] for s in range(3)]
        cell = lambda k, i: "%s / %s" % (fmt(np.mean([r[0][k] for r in rr], 0)), fmt(np.mean([r[1][k] for r in rr], 0)))
        L.append("| %s | %s | %s | %s | %.2f / %.2f |" % (c, cell("sh9", 0), cell("sh50", 0), cell("sh150", 0),
                                                         np.mean([r[0]["rms"] for r in rr]),
                                                         np.mean([r[1]["rms"] for r in rr])))
    L += ["", "## Exponents: β and the exponent of the fundamental channel β_f = β + d ln R1 / d ln B", "",
          "| corner | β old / new | β_f old / new |", "|---|---|---|"]
    for c in CJ.CORNERS_NX:
        rr = [D["%s_%d" % (c, s)] for s in range(3)]
        L.append("| %s | %.2f / %.2f | %.2f / %.2f |" % (
            c, np.nanmean([r[0]["beta"] for r in rr]), np.nanmean([r[1]["beta"] for r in rr]),
            np.nanmean([r[0]["beta_f"] for r in rr]), np.nanmean([r[1]["beta_f"] for r in rr])))
    print("\n".join(L))


if __name__ == "__main__":
    main()
