#!/usr/bin/env python
"""
Where do the loss and the factor effects come from? Exploratory evaluation (hints, not pre-registered) of the
snapshot and accumulator results of the α 0.1 Sobol runs (PROTOKOLL 7.46):
  structure  from domains.json: texture A_w/V (wall measure, 7.14), F90 (volume fraction rotated > 90° on the scale
             of a wall; one vortex per particle gives ≈ 4.3 %, distinct 180° walls more), l_C, of the relaxed virgin
             state (init) and at 9 / 50 / 150 mT (mean over the 4 quarter phases of the last cycle); switched volume
             per half cycle V_sw/V and its exponent over B (9 → 150 mT)
  loss       from acc.json (all samples): R1 (share of the local fundamental), share of the local loss in the core
             r/R < 0.5, the shell r/R ≥ 0.8, the texture top 10 %, the p top 10 %
  (1) range of each measure over the design and the factor effects on it (main effects, OLS HC3, coded factors)
  (2) does a measure carry the factor effect on β all? β all on the 4 factors with and without the (standardised)
      measure: the d/l_ex coefficient with and without, and the measure coefficient (OLS HC3)

    python eval_mechanism.py > results/mechanism.md
"""
import json
import math
import pathlib
import warnings

import numpy as np

import analyze_design as AD

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE / "runs"


def at(stages, mT, key="B_peak_mT"):
    return min(stages, key=lambda s: abs(s[key] / mT - 1))


def measures(name):
    d = RUNS / name
    dom = json.loads((d / "domains.json").read_text())
    acc = json.loads((d / "acc.json").read_text())["stages"]
    out = {}
    ini = dom["init"]
    out["A_w/V init"] = 1e3 * ini["A_w"] / ini["V"]
    out["F90 init"] = 100 * ini["F90"]
    out["l_C init"] = ini["l_C"]
    vsw = []
    for mT in (9, 50, 150):
        st = at(dom["stages"], mT)
        S = st["states"]
        out["A_w/V %d" % mT] = 1e3 * np.mean([s["A_w"] for s in S]) / S[0]["V"]
        out["F90 %d" % mT] = 100 * np.mean([s["F90"] for s in S])
        v = 100 * np.mean(st["half_cycle"]["V_sw"]) / S[0]["V"]
        out["V_sw %d" % mT] = v
        vsw.append(v)
    out["V_sw exponent"] = float(np.polyfit(np.log([9, 50, 150]), np.log(vsw), 1)[0])
    for mT in (9, 150):
        a = at(acc, mT)
        out["R1 %d" % mT] = a["R1"]
        out["core %d" % mT] = a["where"]["r<0.5"]["p_share"]
        out["shell %d" % mT] = a["where"]["r>=0.8"]["p_share"]
        out["texture top %d" % mT] = a["where"]["texture top 10%"]["p_share"]
        out["p top %d" % mT] = a["where"]["p top 10%"]["p_share"]
    return out


def ols(y, X):
    import statsmodels.api as sm
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sm.OLS(y, X).fit(cov_type="HC3")


def main():
    recs, _ = AD.load_runs(RUNS)
    M = [measures(r["name"]) for r in recs]
    keys = list(M[0].keys())
    X = AD.coded(recs)
    D = np.column_stack([np.ones(len(recs)), X])
    print("# Where the loss and the factor effects come from (exploratory, n = %d α 0.1 runs)\n" % len(recs))
    print("Main effects = range effect over the design (2 × coefficient of the coded factor), OLS with HC3 SE. "
          "Hints only: %d measures × 4 factors.\n" % len(keys))
    print("## 1. Measures over the design and the factor effects on them\n")
    print("| measure | median (min … max) | d/l_ex | ln Q_eff | r_p | φ |")
    print("|---|---|---|---|---|---|")
    for k in keys:
        y = np.array([m[k] for m in M])
        r = ols(y, D)
        cells = ["%+.3g ± %.2g%s" % (2 * r.params[j + 1], 2 * r.bse[j + 1], " *" if r.pvalues[j + 1] < 0.001 else "")
                 for j in range(4)]
        print("| %s | %.3g (%.3g … %.3g) | %s |" % (k, np.median(y), y.min(), y.max(), " | ".join(
            [cells[3], cells[0], cells[1], cells[2]])))
    print("\n\\* p < 0.001. Units: A_w/V 10⁻³/l_ex; F90, V_sw in % of the particle volume; l_C in l_ex; shares of "
          "the local loss (volume shares: core 0.125, shell 0.49, top 10 % 0.10).\n")
    f90 = np.array([m["F90 init"] for m in M])
    print("Runs with F90 (virgin state) above one vortex per particle (4.3 %%): %d of %d; maximum %.1f %%.\n" %
          ((f90 > 4.3).sum(), len(f90), f90.max()))
    print("## 2. Does a measure carry the d/l_ex effect on β all?\n")
    b = np.array([r["beta_all"] for r in recs])
    r0 = ols(b, D)
    print("β all on the 4 factors (linear): d/l_ex range effect %+.3f ± %.3f.\n" % (2 * r0.params[4], 2 * r0.bse[4]))
    print("| measure added | d/l_ex range effect on β all | measure: Δβ per 1 sd (p) |")
    print("|---|---|---|")
    for k in keys:
        y = np.array([m[k] for m in M])
        z = (y - y.mean()) / y.std()
        r = ols(b, np.column_stack([D, z]))
        print("| %s | %+.3f ± %.3f | %+.3f (%.2g) |" % (k, 2 * r.params[4], 2 * r.bse[4], r.params[5], r.pvalues[5]))
    print("\nA measure that carries the effect makes the d/l_ex coefficient shrink towards 0 when it is added. "
          "A measure can also be a consequence of the loss (for example the loss shares), not a cause.")


if __name__ == "__main__":
    main()
