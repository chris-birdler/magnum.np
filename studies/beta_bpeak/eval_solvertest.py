#!/usr/bin/env python
"""
Solver test at 150 mT (decision Chris 2026-10-18; rules fixed before the runs). In the corner test the new runs of K1
had at 150 mT a loss about 45 % lower than the old runs at all 3 seeds (ln w old − new +0.57), with an energy balance
better than 1 % in both. Question: does the stray-field solver cause this? From the state of K1_s_nx after the 100 mT
stage (ckpt_05_amp05.pt), the 150 mT stage again: <run>_150fd with DemagFieldPBC in the cycles, <run>_150nx with the
final setting (newell_ext, 4 shells); all other settings as the corner run (conv_jobs.py --newell solvertest).
Value per seed: Δ_s = ln w(150fd) − ln w(150nx), w = mean w_loop of the calm tail of the 150 mT stage (steady-tail rule,
≥ 2 calm cycles). n 3 seeds, mean ± SE, 90 % CI with t(n − 1).
Reading:
  "solver effect"     the 90 % CI of the mean Δ_s excludes 0;
  "no solver effect"  the 90 % CI lies inside ± 0.10;
  "not determined"    all other cases, or fewer than 3 seeds with both values.
Information: the mean Δ_s against the old − new difference of the corner runs at 150 mT (+0.57 at K1, n 3); the 150 mT
stage of the corner run itself against <run>_150nx (the start from the state file); ε of the calm tails.

    python eval_solvertest.py > runs/solvertest_report.md
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
LIMIT, MIN_SEEDS, STAGE = 0.10, 3, "amp06"


def stage_value(d, name=STAGE):
    """(ln w, mean ε) of the calm tail of stage `name` that the run computed itself (not copied); NaN without DONE or a
    tail < 2 cycles."""
    nan = (float("nan"), float("nan"))
    if not (d / "DONE").exists():
        return nan
    st = [s for s in json.loads((d / "summary.json").read_text())["stages"]
          if s.get("name") == name and s.get("measure") and "copied_from" not in s]
    if not st:
        return nan
    tail = analyze.steady_tail(st[0]["cycles"][analyze.first_measured(st[0]):], AD.BAL_EVENT)
    if len(tail) < 2:
        return nan
    return math.log(np.mean([c["w_loop"] for c in tail])), float(np.mean([c["w_dis"] / c["w_loop"] - 1 for c in tail]))


def reading(ds):
    v = np.array([x for x in ds if np.isfinite(x)])
    if len(v) < MIN_SEEDS:
        return "not determined", None
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.95, len(v) - 1) * se
    s = (v.mean(), se, v.mean() - h, v.mean() + h, len(v))
    if s[2] > 0 or s[3] < 0:
        return "solver effect", s
    if -LIMIT < s[2] and s[3] < LIMIT:
        return "no solver effect", s
    return "not determined", s


def main(root=R):
    L = ["# Solver test at 150 mT: DemagFieldPBC against the final setting from the same state (eval_solvertest.py)", "",
         "| run | ln w fd | ln w newell | Δ_s (fd − newell) | ε fd | ε newell | ln w corner run − newell restart |",
         "|---|---|---|---|---|---|---|"]
    ds = []
    for n in CJ.SOLVER_RUNS:
        a, b, c = stage_value(root / (n + "_150fd")), stage_value(root / (n + "_150nx")), stage_value(root / n)
        ds.append(a[0] - b[0])
        L.append("| %s | %.4f | %.4f | %+.4f | %+.4f | %+.4f | %+.4f |" % (n, a[0], b[0], a[0] - b[0], a[1], b[1],
                                                                       c[0] - b[0]))
    key, s = reading(ds)
    L += ["", "Mean Δ_s: %s. **Reading: %s.** For comparison: ln w old − new of the corner runs K1 at 150 mT +0.57 "
          "(n 3)." % ("%+.3f ± %.3f (SE, n %d; 90 %% CI %+.3f … %+.3f)" % (s[0], s[1], s[4], s[2], s[3]) if s else "–",
                      key)]
    print("\n".join(L))


if __name__ == "__main__":
    main()
