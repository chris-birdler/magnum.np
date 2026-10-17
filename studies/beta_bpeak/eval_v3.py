#!/usr/bin/env python
"""
PLAN 4.12 V3 (PROTOKOLL 7.109; the decision after V1 delegated by Chris 2026-10-17): mesh convergence at 9 mT with the
final stray-field setting (DemagField(true_pbc=True, extrapolate=True), 4 shells; prefix N1b; the tighter tolerance of
the drive cycles). Rules fixed before the runs.

Value of a run: ln P of the calm tail (eval_conv.run_record, the rules of the B0b / N0b reports). Paired by seed:
Δ_a = ln P(T15) − ln P(F) (dx 1.5 − dx 0.75) and Δ_b = ln P(T3) − ln P(T15) (dx 3 − dx 1.5); mean ± SE, 90 % CI with
t(n − 1), n ≥ 3 seeds with a value in both arms. Reading (PLAN 4.11 with the "not determined" case of PLAN 4.12):
  converged at dx 1.5   the 90 % CI of Δ_a lies inside ± 0.10;
  not determined        the 90 % CI of Δ_a includes 0 (and does not lie inside ± 0.10), or n < 3;
  converging            the CI of Δ_a excludes 0, Δ_a has the sign of mean Δ_b, and |mean Δ_a| < |mean Δ_b|;
  not converging        the CI of Δ_a excludes 0, and Δ_a changes sign against mean Δ_b or |mean Δ_a| ≥ |mean Δ_b|.
The same for the I arms (ln P(I15) − ln P(F), ln P(I3) − ln P(I15)) as information. Also the reading of PLAN 4.10
(eval_conv.reading, SD_T). Information: μ_r = b_peak / h per arm (geometric mean) and the paired mesh steps of ln μ_r
(mean ± SE); ε = w_dis / w_loop − 1 of the F tails per seed (expected ≤ 0.03, an estimate from n 2, PROTOKOLL 7.108).

    python eval_v3.py > runs/v3_n1b_report.md
"""
import json
import math

import numpy as np

import analyze
import analyze_design as AD
import conv_jobs as CJ
import eval_conv as EC

PREFIX = CJ.PREFIX_FINAL


def conv_reading(lnP, a, b, c, seeds=CJ.SEEDS):
    """(key, Δ_a statistics, Δ_b statistics) with Δ_a = ln P(a) − ln P(b), Δ_b = ln P(c) − ln P(a)."""
    pa, pb = EC.paired(lnP, a, b, seeds), EC.paired(lnP, c, a, seeds)
    if "mean" not in pa or "mean" not in pb:
        return "not determined", pa, pb
    lo, hi = pa["ci"]
    if -EC.LIMIT < lo and hi < EC.LIMIT:
        return "converged at dx 1.5", pa, pb
    if lo <= 0 <= hi:
        return "not determined", pa, pb
    if np.sign(pa["mean"]) == np.sign(pb["mean"]) and abs(pa["mean"]) < abs(pb["mean"]):
        return "converging", pa, pb
    return "not converging", pa, pb


def tail(d):
    """(μ_r, ε) of the steady tail of the measured stage of one run folder; NaN without DONE."""
    if not (d / "DONE").exists():
        return float("nan"), float("nan")
    st = [x for x in json.loads((d / "summary.json").read_text())["stages"] if x.get("measure")][-1]
    cyc = analyze.steady_tail(st["cycles"][analyze.first_measured(st):], AD.BAL_EVENT)
    if len(cyc) < 2:
        return float("nan"), float("nan")
    return (float(np.mean([c["b_peak"] for c in cyc])) / st["h_amp"],
            float(np.mean([c["w_dis"] / c["w_loop"] - 1 for c in cyc])))


def main(root=EC.R):
    recs = EC.collect(root, PREFIX)
    text, _ = EC.report(recs)
    text = text.replace("# B0: mesh", "# %s (final stray-field setting, PLAN 4.12 V3): mesh" % PREFIX, 1)
    lnP = {k: r["lnP"] for k, r in recs.items()}
    L = [text, "", "## Convergence reading (PLAN 4.11 / 4.12, pre-registered)", ""]
    for lab, (a, b, c) in (("T arm (decisive)", ("T15", "F", "T3")), ("I arm (information)", ("I15", "F", "I3"))):
        key, pa, pb = conv_reading(lnP, a, b, c)
        L.append("- %s: Δ(1.5 − 0.75) %s; Δ(3 − 1.5) %s. **Reading: %s.**" % (lab, EC.fmt_p(pa), EC.fmt_p(pb), key))
    tl = {(arm, s): tail(root / CJ.name(arm, s, PREFIX)) for arm, _, _ in CJ.ARMS for s in CJ.SEEDS}
    mu = {k: math.log(v[0]) if np.isfinite(v[0]) else float("nan") for k, v in tl.items()}
    L += ["", "## Information: μ_r = b_peak / h and the energy balance", "",
          "| arm | μ_r (geometric mean) | SD of ln μ_r |", "|---|---|---|"]
    for arm, _, _ in CJ.ARMS:
        v = np.array([mu[(arm, s)] for s in CJ.SEEDS if np.isfinite(mu[(arm, s)])])
        L.append("| %s | %s | %s |" % (arm, "%.3f" % math.exp(v.mean()) if len(v) else "–",
                                        "%.4f" % v.std(ddof=1) if len(v) > 1 else "–"))
    for a, b in (("T3", "T15"), ("T15", "F"), ("I3", "I15"), ("I15", "F")):
        p = EC.paired(mu, a, b, CJ.SEEDS)
        L.append("- ln μ_r(%s) − ln μ_r(%s): %s" % (a, b, "%+.4f ± %.4f (SE, n %d)" % (p["mean"], p["se"], p["n"])
                                                    if "mean" in p else "n %d" % p["n"]))
    L.append("- ε of the F tails (s941 … s944): %s." % " / ".join(
        "%+.3f" % tl[("F", s)][1] if np.isfinite(tl[("F", s)][1]) else "–" for s in CJ.SEEDS))
    print("\n".join(L))


if __name__ == "__main__":
    main()
