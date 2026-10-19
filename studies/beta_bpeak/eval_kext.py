#!/usr/bin/env python
"""
Effect of the shell extrapolation on the loss at dx 1.5 l_ex (decision Chris 2026-10-16, option 1' of PROTOKOLL 7.108;
rules fixed before the runs; per arm, decision Chris 2026-10-16 before any result). The shell test
(results/shells_N0b.md) found: the Newell image sum converges as 1 / (2k + 1)², and with 4 shells E(4) is too small by
up to 5.75 % for the N0b F states. Estimated from check 2 (0.446 (E(4) − E(2)) / E∞): T15 0.05 / 2.8 / 1.9 / 3.1 %,
I15 0.6 / 1.0 / 0.9 / −0.1 % (s941 … s944). Test: the dx 1.5 runs of N0b (T15 and I15, seeds 941 … 944) run again
from their own relaxed state with the extrapolated kernel in the drive cycles (--strayfield_cycles newell_ext;
conv_jobs.py --newell kext; runs <run>_kext). Value per run: d = ln P(kext) − ln P(4 shells) over the steady tail
(eval_conv.run_record, the rules of the reports). Per arm (n 4 seeds): mean ± SE and 90 % CI with t(3).
Reading per arm:
  "no effect at dx 1.5"   the 90 % CI of the mean d lies inside ± 0.10;
  "effect"                the 90 % CI excludes 0 and |mean d| > 0.10;
  "not determined"        all other cases, or fewer than 3 seeds with a value in both runs.
T15 decides (its 4-shell error is like that of F); I15 has its own reading (an error 3 … 30 × smaller).
Information: status and b_peak / b_target of the old and the new runs (without --b_tol a shift of b_peak enters d:
I15 s942 … s944).

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
         "| run | P, 4 shells | P, extrapolated | d | status (old / new) | b_peak / b_target (old / new) |",
         "|---|---|---|---|---|---|"]
    d = {}
    for (arm, s), (a, b) in rec.items():
        d[(arm, s)] = b["lnP"] - a["lnP"]
        L.append("| %s | %.3f | %.3f | %+.4f | %s / %s | %s / %s |" % (
            CJ.name(arm, s, CJ.PREFIX_NEWELL), math.exp(a["lnP"]), math.exp(b["lnP"]), d[(arm, s)], a["status"],
            b["status"], a["b_ratio"], b["b_ratio"]))
    L.append("")
    keys = {}
    for arm in CJ.KEXT_ARMS:
        key, st = reading([d[(arm, s)] for s in CJ.SEEDS])
        keys[arm] = key
        if "mean" in st:
            L.append("%s: d = %+.4f ± %.4f (SE, n %d; 90 %% CI %+.3f … %+.3f). **Reading %s: %s.**" % (
                arm, st["mean"], st["se"], st["n"], st["ci"][0], st["ci"][1], arm, key))
        else:
            L.append("%s: %d seeds with a value. **Reading %s: %s.**" % (arm, st["n"], arm, key))
    L += ["", "**Decisive reading (T15): %s.**" % keys["T15"]]
    print("\n".join(L))


if __name__ == "__main__":
    main()
