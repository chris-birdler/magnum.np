#!/usr/bin/env python
"""
Exploratory evaluation of the stopped Sobol design v2 (PROTOKOLL 7.36): the finished α 0.1 runs only.
Too few runs for the 15-term model of PLAN 4.0: main effects only (5 terms), no decision rules, all results are
hints. The r_p effect contains an unknown part of the non-isotropic stress-axis set (eigenvalues of <u u^T>
0.253 / 0.285 / 0.461 instead of 1/3).
  table 1  per run: factors (real units), beta windows and beta_all, P at 50 mT (W/cm^3, core volume), flags
  table 2  per run from the accumulators: R1 and the shares of the local loss in core / shell / texture top 10 %
  table 3  main effects (range effect over the design range): Gamma GLM (mean, bootstrap) and OLS (HC3)
  table 4  corner replicates (2 seeds of K2 and K3)

    python eval_sobol_partial.py > results/sobol_partial.md
"""
import json
import math
import pathlib
import warnings

import numpy as np

import analyze_design as AD

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE / "runs" / "_archive_sobol_v2_oldaxes"     # default: the stopped design (old stress axes), archived 7.38
F_HZ = 30e6


def main():
    import argparse
    import statsmodels.api as sm
    global RUNS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=str(RUNS), help="run folder (default: archive of the stopped design)")
    ap.add_argument("--title", default="Sobol design v2, stopped: exploratory evaluation of the finished runs")
    ap.add_argument("--axes_note", default="**The r_p effect contains an unknown part of the non-isotropic "
                    "stress-axis set** (PROTOKOLL 7.36).")
    a = ap.parse_args()
    RUNS = pathlib.Path(a.runs)
    recs, skipped = AD.load_runs(RUNS)
    recs.sort(key=lambda r: r["name"])
    n = len(recs)
    kd = {r["name"]: json.loads((RUNS / r["name"] / "config.json").read_text())["units"]["Kd"] for r in recs}
    print("# %s\n" % a.title)
    print("Runs: %d finished at α 0.1 (%s); not usable: %s. Main effects only (n too small for the 15-term "
          "model); hints, no decision rules. %s\n" % (n, ", ".join(r["name"] for r in recs), "; ".join(skipped) or "none",
                                                     a.axes_note))
    print("## 1. Runs\n")
    print("| run | d/l_ex | L_eff/l_ex | r_p | φ | d/L_eff | β 9-35 | β 35-70 | β 70-150 | β all | P(50 mT) W/cm³ | flagged |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in recs:
        P50 = math.exp(r["lnw@50mT"]) * kd[r["name"]] * F_HZ / 1e6
        print("| %s | %.0f | %.0f | %.2f | %.3f | %.1f | %.2f | %.2f | %.2f | %.2f | %.0f | %s |" %
              (r["name"], math.exp(r["lnd"]), 1 / math.sqrt(math.exp(r["lnQ"])), r["r_p"], r["phi"], r["d_over_L"],
               r["beta@9-35mT"], r["beta@35-70mT"], r["beta@70-150mT"], r["beta_all"], P50, "yes" if r["flagged"] else ""))
    b = np.array([r["beta_all"] for r in recs])
    print("\nβ all: mean %.2f, sd %.2f, range %.2f … %.2f (n = %d). Base point of the earlier studies (d/l_ex 300, "
          "L_eff 12, r_p 0, φ 0.65): β_pow 1.76 ± 0.05 (n = 8, 7.17).\n" % (b.mean(), b.std(ddof=1), b.min(), b.max(), n))

    print("## 2. Where the loss sits (accumulators, all samples)\n")
    print("| run | R1 9 mT | R1 50 mT | R1 150 mT | core r/R<0.5 9 / 150 mT | shell r/R≥0.8 9 / 150 mT | texture top 10 % 9 / 150 mT | p top 10 % 9 / 150 mT |")
    print("|---|---|---|---|---|---|---|---|")
    for r in recs:
        try:
            st = json.loads((RUNS / r["name"] / "acc.json").read_text())["stages"]
        except OSError:
            continue
        def at(mT):
            return min(st, key=lambda s: abs(s["B_peak_mT"] / mT - 1))
        s9, s50, s150 = at(9), at(50), at(150)
        sh = lambda s, k: s["where"][k]["p_share"]
        print("| %s | %.2f | %.2f | %.2f | %.2f / %.2f | %.2f / %.2f | %.2f / %.2f | %.2f / %.2f |" %
              (r["name"], s9["R1"], s50["R1"], s150["R1"], sh(s9, "r<0.5"), sh(s150, "r<0.5"), sh(s9, "r>=0.8"),
               sh(s150, "r>=0.8"), sh(s9, "texture top 10%"), sh(s150, "texture top 10%"), sh(s9, "p top 10%"),
               sh(s150, "p top 10%")))
    print("\nVolume shares of the classes: core ≈ 0.125, shell ≈ 0.49, texture top 10 % and p top 10 % = 0.10.\n")

    print("## 3. Main effects (range effect over the design range, 5 terms: 1 + 4 linear factors)\n")
    X = AD.coded(recs)
    T = np.column_stack([np.ones(n), X])
    W = np.exp(np.array([[r["lnw@%gmT" % m] for m in AD.AMPS_MT] for r in recs]))
    fam = sm.families.Gamma(sm.families.links.Log())

    def glm_eff(idx):
        lnmu = []
        for k in range(W.shape[1]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                c = sm.GLM(W[idx, k], T[idx], family=fam).fit().params
            lnmu.append(c)
        C = np.array(lnmu)                          # [7, 5]: linear in the coded factors -> range effect = 2 c_j
        out = {}
        for o in AD.OUTPUTS:
            w = AD._OUT_W[o]
            for j, f in enumerate(AD.FACTORS):
                out[(o, f)] = float(2.0 * (w @ C[:, j + 1]))
        return out

    est = glm_eff(np.arange(n))
    rng = np.random.default_rng(20261007)
    boots = []
    for _ in range(1000):
        i = rng.integers(0, n, n)
        try:
            boots.append(glm_eff(i))
        except Exception:                           # noqa: BLE001
            pass
    print("Main fit: Gamma GLM (mean loss), 95 %% CI by bootstrap over the runs (%d resamples). Fit A: OLS of the "
          "per-run values with HC3. φ: the ln E[w] effect contains the dilution %.2f.\n" % (len(boots), math.log(0.69 / 0.55)))
    print("| output | factor | main fit (95 % CI) | fit A ± SE |")
    print("|---|---|---|---|")
    for o in ("beta_all", "beta@9-35mT", "beta@35-70mT", "beta@70-150mT", "lnw@10mT", "lnw@50mT", "lnw@100mT"):
        y = np.array([r[o] for r in recs])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = sm.OLS(y, T).fit(cov_type="HC3")
        for j, f in enumerate(AD.FACTORS):
            bb = np.array([x[(o, f)] for x in boots])
            lo, hi = np.percentile(bb, [2.5, 97.5])
            print("| %s | %s | %+.2f (%+.2f … %+.2f) | %+.2f ± %.2f |" % (o, f, est[(o, f)], lo, hi, 2 * m.params[j + 1],
                                                                       2 * m.bse[j + 1]))
    print("\nFactor ranges: ln Q_eff ↔ L_eff 30 … 12 l_ex; r_p 0 … 3; φ 0.55 … 0.69; ln d/l_ex 212 … 300 "
          "(at d = 1 µm: A 20 … 10 pJ/m).\n")

    print("## 4. Corner replicates\n")
    byname = {r["name"]: r for r in recs}
    for k in ("K1", "K2", "K3", "K4"):
        reps = [byname[x] for x in sorted(byname) if x.startswith(k + "_")]
        if len(reps) >= 2:
            print("- %s (%s): β all %s; ln w 50 mT %s" % (k, ", ".join(r["name"] for r in reps),
                                                     " / ".join("%.2f" % r["beta_all"] for r in reps),
                                                     " / ".join("%.2f" % r["lnw@50mT"] for r in reps)))
        elif reps:
            print("- %s: one seed only (%s)" % (k, reps[0]["name"]))


if __name__ == "__main__":
    main()
