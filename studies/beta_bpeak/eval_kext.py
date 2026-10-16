#!/usr/bin/env python
"""
Effect of the shell extrapolation on the loss at dx 1.5 l_ex (decision Chris 2026-10-16, option 1' of PROTOKOLL 7.108;
rules fixed before the runs). The shell test (results/shells_N0b.md) found: the Newell image sum converges as
1 / (2k + 1)², and 4 shells miss up to 5.75 % of the stray-field energy of the N0b F states (2.9 % of their block
averages on dx 1.5); the missing energy is the same on all meshes. Test: the dx 1.5 runs of N0b (T15 and I15, seeds
941 … 944) run again from their own relaxed state with the extrapolated kernel in the drive cycles
(--strayfield_cycles newell_ext; conv_jobs.py --newell kext; runs <run>_kext). Value per run:
d = ln P(kext) − ln P(4 shells) over the steady tail (eval_conv.run_record, the rules of the reports). The seed is the
unit: d_s = mean of the T15 and I15 values of seed s (n 4), mean ± SE and 90 % CI with t(3).
Reading:
  "no effect at dx 1.5"   the 90 % CI of mean d_s lies inside ± 0.10;
  "effect"                the 90 % CI excludes 0 and |mean d_s| > 0.10;
  "not determined"        all other cases, or fewer than 3 seeds with a value in both arms.
Information: the same per arm (n 4), the status and b_peak / b_target of the new runs.

    python eval_kext.py > runs/kext15_report.md
"""
import math

import numpy as np
from scipy import stats

import conv_jobs as CJ
import eval_conv as EC

LIMIT = 0.10


def reading(d):
    """(key, statistics) of the per-seed values d (NaN = no value)."""
    v = np.array([x for x in d if np.isfinite(x)])
    if len(v) < EC.MIN_N:
        return "not determined", {"n": len(v)}
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.95, len(v) - 1) * se
    ci = (v.mean() - h, v.mean() + h)
    st = {"n": len(v), "mean": float(v.mean()), "se": float(se), "ci": ci}
    if -LIMIT < ci[0] and ci[1] < LIMIT:
        return "no effect at dx 1.5", st
    if (ci[0] > 0 or ci[1] < 0) and abs(v.mean()) > LIMIT:
        return "effect", st
    return "not determined", st


def main():
    rec = {}
    for arm in CJ.KEXT_ARMS:
        for s in CJ.SEEDS:
            n = CJ.name(arm, s, CJ.PREFIX_NEWELL)
            rec[(arm, s)] = (EC.run_record(EC.R / n), EC.run_record(EC.R / (n + CJ.KEXT_SUFFIX)))
    L = ["# Effect of the shell extrapolation on the loss at dx 1.5 (eval_kext.py, PROTOKOLL 7.108)", "",
         "d = ln P(extrapolated kernel) − ln P(4 shells), steady tail, P in W/cm³; the same relaxed state.", "",
         "| run | P, 4 shells | P, extrapolated | d | status (new) | b_peak / b_target (new) |", "|---|---|---|---|---|---|"]
    d = {}
    for (arm, s), (a, b) in rec.items():
        d[(arm, s)] = b["lnP"] - a["lnP"]
        L.append("| %s | %.3f | %.3f | %+.4f | %s | %s |" % (CJ.name(arm, s, CJ.PREFIX_NEWELL), math.exp(a["lnP"]),
                                                            math.exp(b["lnP"]), d[(arm, s)], b["status"],
                                                            b["b_ratio"]))
    ds = [np.nanmean([d[(arm, s)] for arm in CJ.KEXT_ARMS]) if all(np.isfinite(d[(arm, s)]) for arm in CJ.KEXT_ARMS)
          else float("nan") for s in CJ.SEEDS]
    key, st = reading(ds)
    L += ["", "Per seed (mean of the arms): %s." % " / ".join("%+.4f" % x for x in ds)]
    for arm in CJ.KEXT_ARMS:
        _, sa = reading([d[(arm, s)] for s in CJ.SEEDS])
        if "mean" in sa:
            L.append("%s: %+.4f ± %.4f (SE, n %d; 90 %% CI %+.3f … %+.3f)." % (arm, sa["mean"], sa["se"], sa["n"],
                                                                         sa["ci"][0], sa["ci"][1]))
    if "mean" in st:
        L.append("Seed mean: %+.4f ± %.4f (SE, n %d; 90 %% CI %+.3f … %+.3f)." % (st["mean"], st["se"], st["n"],
                                                                             st["ci"][0], st["ci"][1]))
    L += ["", "**Reading: %s.**" % key]
    print("\n".join(L))


if __name__ == "__main__":
    main()
