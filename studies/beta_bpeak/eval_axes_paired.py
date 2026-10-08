#!/usr/bin/env python
"""
Paired comparison old vs new stress axes (PROTOKOLL 7.36, 7.37, 7.42): the same Sobol points and seeds (same
Herzer cubes and the same random initial state), computed with the old non-isotropic axis set (archive) and with
the isotropic tetrahedral set (new runs). Difference new - old of beta (windows, all) and ln w (10 / 50 / 100 mT),
mean +- SE over the pairs, Student t; the slope of the difference over r_p (the axis set acts through K_p = r_p K_eff,
so the difference should grow with r_p); points with r_p < 0.1 as a control (axes nearly irrelevant there).

    python eval_axes_paired.py > results/axes_paired.md
"""
import math
import pathlib

import numpy as np
from scipy import stats

import analyze_design as AD

HERE = pathlib.Path(__file__).resolve().parent
OLD = HERE / "runs" / "_archive_sobol_v2_oldaxes"
NEW = HERE / "runs"
OUTS = ["beta_all", "beta@9-35mT", "beta@35-70mT", "beta@70-150mT", "lnw@10mT", "lnw@50mT", "lnw@100mT"]


def main():
    old = {r["name"]: r for r in AD.load_runs(OLD)[0]}
    new = {r["name"]: r for r in AD.load_runs(NEW)[0]}
    names = sorted(set(old) & set(new))
    print("# Old vs new stress axes, paired by Sobol point and seed\n")
    print("Pairs: %d (%s). Not yet recomputed: %s.\n" % (len(names), ", ".join(names),
                                                       ", ".join(sorted(set(old) - set(new))) or "none"))
    print("| run | r_p | β all old → new | β 9-35 old → new | β 70-150 old → new | ln w 50 mT old → new | flagged old / new |")
    print("|---|---|---|---|---|---|---|")
    for n in names:
        o, w = old[n], new[n]
        print("| %s | %.2f | %.2f → %.2f | %.2f → %.2f | %.2f → %.2f | %.2f → %.2f | %s / %s |" %
              (n, o["r_p"], o["beta_all"], w["beta_all"], o["beta@9-35mT"], w["beta@9-35mT"], o["beta@70-150mT"],
               w["beta@70-150mT"], o["lnw@50mT"], w["lnw@50mT"], "yes" if o["flagged"] else "no",
               "yes" if w["flagged"] else "no"))
    rp = np.array([old[n]["r_p"] for n in names])
    print("\n## Difference new − old (mean ± SE, n = %d pairs)\n" % len(names))
    print("| output | mean ± SE | p (t) | slope over r_p ± SE | mean ± SE, r_p ≥ 1 (n) | mean of \\|Δ\\|, r_p < 0.1 (n) |")
    print("|---|---|---|---|---|---|")
    for k in OUTS:
        d = np.array([new[n][k] - old[n][k] for n in names])
        se = d.std(ddof=1) / math.sqrt(len(d))
        p = stats.ttest_1samp(d, 0.0).pvalue
        sl = stats.linregress(rp, d)
        hi, lo = d[rp >= 1.0], d[rp < 0.1]
        print("| %s | %+.3f ± %.3f | %.2f | %+.3f ± %.3f | %+.3f ± %.3f (%d) | %.3f (%d) |" %
              (k, d.mean(), se, p, sl.slope, sl.stderr, hi.mean(), hi.std(ddof=1) / math.sqrt(len(hi)), len(hi),
               np.abs(lo).mean() if len(lo) else float("nan"), len(lo)))
    print("\nln w differences are natural logs (+0.10 ≈ +10 % loss). The relaxed virgin state depends on the axes "
          "when r_p > 0, so a pair differs in the state as well as in the axes.")


if __name__ == "__main__":
    main()
