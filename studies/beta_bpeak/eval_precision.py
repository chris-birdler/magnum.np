#!/usr/bin/env python
"""
Precision test of PLAN 4.13 (part 5, step 4; approved by Chris 2026-10-17; rules fixed before the runs). From the
init.pt of K1_0_nx, K2_0_nx, K3_0_nx and K4_0_nx (one seed per corner, n 4), the 9 mT stage only, all other settings as
the corner runs (conv_jobs.py --newell precision): <run>_p64 runs its cycles in fp64 (--precision_cycles double),
<run>_p32 in fp32; both start from the same init.pt.
Value per pair: Δ_p = ln P(p64) − ln P(p32) = ln w(p64) − ln w(p32), w = mean w_loop of the calm tail (steady-tail rule,
≥ 2 calm cycles). A pair counts only if both runs have a value.
Reading (90 % CI with t(n − 1)):
  "fp32 sufficient"      n ≥ 3 and the 90 % CI of the mean Δ_p lies inside ± 0.05;
  "fp32 not sufficient"  n ≥ 3, the 90 % CI excludes 0 and |mean Δ_p| > 0.05;
  "not determined"       all other cases.
Information: ε of the calm tails (p64, p32); wall time per cycle; ln w of the 9 mT stage of the corner run itself
(start in the same process after the relaxation) against p32.

    python eval_precision.py > runs/precision_nx_report.md
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import conv_jobs as CJ

R = pathlib.Path(__file__).resolve().parent / "runs"
LIMIT, MIN_PAIRS = 0.05, 3


def stage_tail(d, k=0):
    """(ln w, mean ε, wall per cycle) of the calm tail of measured stage k of run folder d; NaN without DONE or a
    tail < 2 cycles."""
    nan = (float("nan"),) * 3
    if not (d / "DONE").exists():
        return nan
    st = [s for s in json.loads((d / "summary.json").read_text())["stages"] if s.get("measure")]
    if len(st) <= k:
        return nan
    tail = analyze.steady_tail(st[k]["cycles"][analyze.first_measured(st[k]):], AD.BAL_EVENT)
    if len(tail) < 2:
        return nan
    return (math.log(np.mean([c["w_loop"] for c in tail])), float(np.mean([c["w_dis"] / c["w_loop"] - 1 for c in tail])),
            float(np.mean([c["wall_s"] for c in st[k]["cycles"]])))


def reading(dp):
    """(key, (mean, SE, lo, hi, n) or None) for the list of Δ_p (NaN = no pair)."""
    v = np.array([x for x in dp if np.isfinite(x)])
    if len(v) < MIN_PAIRS:
        return "not determined", None
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.95, len(v) - 1) * se
    s = (v.mean(), se, v.mean() - h, v.mean() + h, len(v))
    if -LIMIT < s[2] and s[3] < LIMIT:
        return "fp32 sufficient", s
    if (s[2] > 0 or s[3] < 0) and abs(s[0]) > LIMIT:
        return "fp32 not sufficient", s
    return "not determined", s


def main(root=R):
    L = ["# Precision test of PLAN 4.13: fp64 against fp32 in the drive cycles (eval_precision.py)", "",
         "9 mT stage from the same init.pt; Δ_p = ln w(p64) − ln w(p32) on the calm tails.", "",
         "| run | ln w p64 | ln w p32 | Δ_p | ε p64 | ε p32 | wall per cycle p64 / p32 (s) | ln w corner run − p32 |",
         "|---|---|---|---|---|---|---|---|"]
    dp = []
    for n in CJ.PRECISION_RUNS:
        a, b, c = stage_tail(root / (n + "_p64")), stage_tail(root / (n + "_p32")), stage_tail(root / n)
        dp.append(a[0] - b[0])
        L.append("| %s | %.4f | %.4f | %+.4f | %+.4f | %+.4f | %.0f / %.0f | %+.4f |" % (
            n, a[0], b[0], a[0] - b[0], a[1], b[1], a[2], b[2], c[0] - b[0]))
    key, s = reading(dp)
    L += ["", "Mean Δ_p: %s. **Reading: %s.**" % (
        "%+.4f ± %.4f (SE, n %d; 90 %% CI %+.4f … %+.4f)" % (s[0], s[1], s[4], s[2], s[3]) if s else "–", key)]
    print("\n".join(L))


if __name__ == "__main__":
    main()
