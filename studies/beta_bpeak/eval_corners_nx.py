#!/usr/bin/env python
"""
Corner test (decision Chris 2026-10-17: first the corners, then decide on a new main study; rules fixed before the
runs). The 12 corner runs of the main study (K1 … K4, 3 seeds each, dx 3 l_ex, 9 … 150 mT) again with the final
stray-field setting (newell_ext, 4 shells; conv_jobs.py --newell corners; runs <run>_nx). Only the solver changes.
Value per run: β = slope of ln w over ln B at the nominal amplitudes (main-fit definition; analyze_design.run_record
with the steady-tail rule ≥ 2 calm cycles). Paired per run: both fits use the amplitudes with a value in the old and the
new run, at least 6 of 7 (3 old runs miss one amplitude by the steady-tail rule); Δβ = β(new) − β(old).
Reading:
  "β robust to the solver"   the 90 % CI of the mean Δβ over the 12 runs lies inside ± 0.10, and no corner has
                             |mean Δβ| > 0.20 (n 3 per corner);
  "β shifted"                the 90 % CI of the mean Δβ excludes 0 and |mean Δβ| > 0.10;
  "not determined"           all other cases, or fewer than 9 runs with β in both.
Information: per corner β old and new (mean ± SE, n 3) and the corner differences β(Ki) − β(K1) old and new;
per amplitude the mean Δ ln w (new − old) with SE (n 12).

    python eval_corners_nx.py > runs/corners_nx_report.md
"""
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import conv_jobs as CJ

R = pathlib.Path(__file__).resolve().parent / "runs"
LIMIT, CORNER_LIMIT, MIN_RUNS = 0.10, 0.20, 9


def lnw(d):
    r = AD.run_record(d)
    return np.array([r["lnw@%gmT" % b] for b in AD.AMPS_MT]) if r else np.full(len(AD.AMPS_MT), np.nan)


MIN_AMPS = 6


def beta_pair(o, n):
    """(β old, β new) over the amplitudes with a value in both runs (at least MIN_AMPS), else NaN."""
    m = np.isfinite(o) & np.isfinite(n)
    if m.sum() < MIN_AMPS:
        return float("nan"), float("nan")
    x = np.log(np.asarray(AD.AMPS_MT))[m]
    return float(np.polyfit(x, o[m], 1)[0]), float(np.polyfit(x, n[m], 1)[0])


def stat(v):
    v = np.array([x for x in v if np.isfinite(x)])
    if len(v) < 2:
        return None
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.95, len(v) - 1) * se
    return v.mean(), se, v.mean() - h, v.mean() + h, len(v)


def reading(db, by_corner):
    s = stat(db)
    if s is None or s[4] < MIN_RUNS:
        return "not determined", s
    if LIMIT > s[3] and s[2] > -LIMIT and all(abs(np.nanmean(v)) <= CORNER_LIMIT for v in by_corner.values()):
        return "β robust to the solver", s
    if (s[2] > 0 or s[3] < 0) and abs(s[0]) > LIMIT:
        return "β shifted", s
    return "not determined", s


def main(root=R):
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    runs = [n.replace(CJ.NX_SUFFIX, "") for n, _ in CJ.corner_newell_lines()]
    W = {n: (lnw(root / n), lnw(root / (n + CJ.NX_SUFFIX))) for n in runs}
    B = {n: beta_pair(o, nw) for n, (o, nw) in W.items()}
    db = [B[n][1] - B[n][0] for n in runs]
    by_c = {c: [B[n][1] - B[n][0] for n in runs if n.startswith(c + "_")] for c in CJ.CORNERS_NX}
    key, s = reading(db, by_c)
    L = ["# Corner test: the 12 corner runs with the final stray-field setting (eval_corners_nx.py)", "",
         "β = slope of ln w over ln B, 9 … 150 mT (main-fit definition). Old: DemagFieldPBC; new: newell_ext, 4 shells. "
         "Δβ = new − old, paired per run.", "", "| run | β old | β new | Δβ |", "|---|---|---|---|"]
    for n in runs:
        L.append("| %s | %.3f | %.3f | %+.3f |" % (n, B[n][0], B[n][1], B[n][1] - B[n][0]))
    L += ["", "## Per corner (mean ± SE, n 3)", "", "| corner | β old | β new | Δβ |", "|---|---|---|---|"]
    co = {}
    for c in CJ.CORNERS_NX:
        bo = stat([B[n][0] for n in runs if n.startswith(c + "_")])
        bn = stat([B[n][1] for n in runs if n.startswith(c + "_")])
        d = stat(by_c[c])
        co[c] = (bo, bn)
        f = lambda x: "%.3f ± %.3f" % (x[0], x[1]) if x else "–"
        L.append("| %s | %s | %s | %s |" % (c, f(bo), f(bn), "%+.3f ± %.3f" % (d[0], d[1]) if d else "–"))
    L += ["", "Corner differences β(Ki) − β(K1), old / new: " + "; ".join(
        "%s %+.3f / %+.3f" % (c, co[c][0][0] - co["K1"][0][0], co[c][1][0] - co["K1"][1][0])
        for c in CJ.CORNERS_NX[1:] if co[c][0] and co[c][1] and co["K1"][0] and co["K1"][1]) + "."]
    L += ["", "Mean Δ ln w per amplitude (new − old, SE, n): " + "; ".join(
        "%g mT %s" % (b, "%+.3f ± %.3f (%d)" % (x[0], x[1], x[4]) if (x := stat([W[n][1][k] - W[n][0][k] for n in runs]))
                       else "–") for k, b in enumerate(AD.AMPS_MT)) + "."]
    L += ["", "Mean Δβ over the runs: %s. **Reading: %s.**" % (
        "%+.3f ± %.3f (SE, n %d; 90 %% CI %+.3f … %+.3f)" % (s[0], s[1], s[4], s[2], s[3]) if s else "–", key)]
    print("\n".join(L))


if __name__ == "__main__":
    main()
