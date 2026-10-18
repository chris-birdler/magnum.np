#!/usr/bin/env python
"""
Solver scan (option (a) of the solver test, decision Chris 2026-10-18; rules fixed before the runs). The solver test
showed at K1, 150 mT, from the same state: ln w(DemagFieldPBC) − ln w(Newell) = +0.573 ± 0.018 (SE, n 3), the whole
old − new difference of the corner test there. The scan measures this solver effect at all 7 amplitudes for one seed
per corner (K1_0_nx … K4_0_nx): each stage once more with DemagFieldPBC in the cycles, from the same state as the
corner run (conv_jobs.py --newell solverscan). The Newell value is the corner run itself (the Newell restart of the
solver test reproduced it exactly).
Values per corner: Δ_k = ln w_FD − ln w_Newell at stage k (calm tail, steady-tail rule); β_FD and β_Newell = slope of
ln w over ln B at the 7 nominal amplitudes; Δβ_solver = β_FD − β_Newell.
Reading (mean Δβ_solver over the 4 corners, 90 % CI with t(3)):
  "β shifted by the solver"   the CI excludes 0 and |mean| > 0.10;
  "β robust to the solver"    the CI lies inside ± 0.10;
  "not determined"            all other cases, or fewer than 4 corners with all 7 values.
Information: Δ_k per corner and the mean per amplitude (± SE, n 4); the corner differences β(Ki) − β(K1) with FD and
with Newell; for comparison the corner test of the same seeds (old run with DemagFieldPBC, atol 1e-5, its own virgin
state): β_old − β_new and ln w_old − ln w_new.

    python eval_solverscan.py > runs/solverscan_report.md
"""
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import conv_jobs as CJ
import eval_corners_nx as EN
import eval_solvertest as ES

R = pathlib.Path(__file__).resolve().parent / "runs"
LIMIT = 0.10
X = np.log(np.asarray(AD.AMPS_MT))


def slope(y):
    return float(np.polyfit(X, y, 1)[0]) if np.all(np.isfinite(y)) else float("nan")


def reading(db):
    v = np.array([x for x in db if np.isfinite(x)])
    if len(v) < len(CJ.SCAN_RUNS):
        return "not determined", None
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.95, len(v) - 1) * se
    s = (v.mean(), se, v.mean() - h, v.mean() + h, len(v))
    if (s[2] > 0 or s[3] < 0) and abs(s[0]) > LIMIT:
        return "β shifted by the solver", s
    if -LIMIT < s[2] and s[3] < LIMIT:
        return "β robust to the solver", s
    return "not determined", s


def main(root=R):
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    n_st = len(AD.AMPS_MT)
    L = ["# Solver scan: DemagFieldPBC against Newell at all amplitudes from the same states (eval_solverscan.py)", "",
         "Δ_k = ln w_FD − ln w_Newell per stage (same state at the start of the stage).", "",
         "| corner run | " + " | ".join("%g mT" % b for b in AD.AMPS_MT) + " | β Newell | β FD | Δβ_solver | "
         "β old − β new (corner test) |", "|---|" + "---|" * (n_st + 4)]
    D, B, db = {}, {}, []
    for n in CJ.SCAN_RUNS:
        wn = np.array([ES.stage_value(root / n, "amp%02d" % k)[0] for k in range(n_st)])
        wf = np.array([ES.stage_value(root / CJ.scan_name(n, k), "amp%02d" % k)[0] for k in range(n_st)])
        D[n] = wf - wn
        B[n] = (slope(wn), slope(wf))
        db.append(B[n][1] - B[n][0])
        o = EN.lnw(EN.old_dir(root, n.replace(CJ.NX_SUFFIX, "")))
        bo, bn = EN.beta_pair(o, EN.lnw(root / n))
        L.append("| %s | " % n + " | ".join("%+.3f" % x for x in D[n]) + " | %.3f | %.3f | %+.3f | %+.3f |" % (
            B[n][0], B[n][1], B[n][1] - B[n][0], bo - bn))
    d = np.array([D[n] for n in CJ.SCAN_RUNS])
    m, se = np.nanmean(d, 0), np.nanstd(d, 0, ddof=1) / math.sqrt(len(d))
    L += ["", "Mean Δ_k over the corners (± SE, n %d): " % len(d) + "; ".join(
        "%g mT %+.3f ± %.3f" % (b, a, e) for b, a, e in zip(AD.AMPS_MT, m, se)) + "."]
    k1 = CJ.SCAN_RUNS[0]
    L += ["", "Corner differences β(Ki) − β(K1), Newell / FD: " + "; ".join(
        "%s %+.3f / %+.3f" % (n[:2], B[n][0] - B[k1][0], B[n][1] - B[k1][1]) for n in CJ.SCAN_RUNS[1:]) + "."]
    key, s = reading(db)
    L += ["", "Mean Δβ_solver: %s. **Reading: %s.**" % (
        "%+.3f ± %.3f (SE, n %d; 90 %% CI %+.3f … %+.3f)" % (s[0], s[1], s[4], s[2], s[3]) if s else "–", key)]
    print("\n".join(L))


if __name__ == "__main__":
    main()
