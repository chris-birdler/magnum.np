#!/usr/bin/env python
"""
Loss channels, option O1 of Chris (2026-10-12, PROTOKOLL 7.71): local checks of the channel split of 7.70 with the
data at hand (no new runs).
  1. Identity: β = β_f − d ln R1 / d ln B (P_f = R1 · P). Per run (all 7 amplitudes with a calm accumulator window):
     β, β_f (power-law fit of R1 · P) and the shift −d ln R1 / d ln B (fit of ln R1); means ± SE for the A thirds.
     The A effect on β splits into the A effect on β_f and on the shift.
  2. Where the rest 1 − R1 sits: share of the non-f loss in each zone of acc.json (centre, rim, the 10 % cells with
     the highest loss, the 10 % with the largest texture), (p_share − R1 · p_fund_share) / (1 − R1), calm windows.
  3. Damping: centre runs C0 … C7 at α 0.1 and α 0.01 (same microstructures): R1, the jump energy A = RF / 2 and the
     loop-to-loop energy C = RN / 2 of the last 3 cycles (steady or not: α 0.01 has few calm windows; flagged).
  4. Factor effects on β_f and on the shift (OLS, HC3, range effect = 2 × coefficient, coded factors as the design
     analysis, with the lnd × lnQ term).
    python eval_channels.py > results/loss_channels_o1.md
"""
import json
import math

import numpy as np

import analyze
import analyze_design as AD
import deck_data as DD
import eval_jumps_time as EJ

AMPS = np.array(AD.AMPS_MT)
LNB = np.log(AMPS)
ZONES = ("r<0.5", "r>=0.8", "p top 10%", "texture top 10%")


def mse(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return "%.2f ± %.2f (n %d)" % (x.mean(), x.std(ddof=1) / math.sqrt(len(x)), len(x)) if len(x) > 1 else "–"


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    recs, _ = AD.load_runs(DD.RUNS, reruns=True)
    kd = json.loads((DD.RUNS / recs[0]["name"] / "config.json").read_text())["units"]["Kd"]
    to_P = kd * DD.F_HZ / 1e6
    P = np.array([[math.exp(r["lnw@%gmT" % m]) * to_P for m in AMPS] for r in recs])
    st = [DD.acc_stages(DD.RUNS / r["name"]) for r in recs]
    R1 = np.array([[s["R1"] if s["calm"] else np.nan for s in x] for x in st])
    A = np.array([DD.real(r)["A_pJm"] for r in recs])
    q1, q2 = np.quantile(A, [1 / 3, 2 / 3])
    groups = (("all", A > 0), ("low A", A <= q1), ("high A", A > q2))
    full = np.isfinite(R1).all(1)
    beta = np.array([np.polyfit(LNB, np.log(p), 1)[0] for p in P])
    bf = np.array([np.polyfit(LNB, np.log(R1[i] * P[i]), 1)[0] if full[i] else np.nan for i in range(len(recs))])
    shift = np.array([-np.polyfit(LNB, np.log(R1[i]), 1)[0] if full[i] else np.nan for i in range(len(recs))])

    L = ["# Loss channels: local checks (option O1, PROTOKOLL 7.71)", "",
         "Coarse mesh, α 0.1, 116 design runs (reruns used), calm accumulator windows only (7.70).", ""]
    L += ["## 1. β = β_f + shift (shift = −d ln R1 / d ln B), runs with all 7 amplitudes calm", "",
          "| group | β | β_f (loss at f) | shift | check β − β_f − shift |", "|---|---|---|---|---|"]
    for lab, sel in groups:
        s = sel & full
        L.append("| %s | %s | %s | %s | %.1e |" % (lab, mse(beta[s]), mse(bf[s]), mse(shift[s]),
                                                    np.nanmax(np.abs(beta[s] - bf[s] - shift[s]))))
    lo, hi = groups[1][1] & full, groups[2][1] & full
    d = lambda v: (np.nanmean(v[hi]) - np.nanmean(v[lo]),
                   math.sqrt(np.nanvar(v[hi], ddof=1) / hi.sum() + np.nanvar(v[lo], ddof=1) / lo.sum()))
    L += ["", "Difference high A − low A (Welch SE): β %+.2f ± %.2f, β_f %+.2f ± %.2f, shift %+.2f ± %.2f." % (
        *d(beta), *d(bf), *d(shift)), ""]

    L += ["## 2. Where the rest (1 − R1) sits: share of the non-f loss per zone (median, calm windows)", "",
          "| zone | volume | " + " | ".join("%g mT" % b for b in AMPS) + " |", "|---|---|" + "---|" * len(AMPS)]
    for z in ZONES:
        v = np.array([[(s["where"][z]["p_share"] - s["R1"] * s["where"][z]["p_fund_share"]) / (1 - s["R1"])
                       if s["calm"] else np.nan for s in x] for x in st])
        f = np.array([[s["where"][z]["p_fund_share"] if s["calm"] else np.nan for s in x] for x in st])
        vol = np.nanmedian([[s["where"][z]["vol_share"] for s in x] for x in st])
        L.append("| %s, rest | %.0f %% | %s |" % (z, 100 * vol, " | ".join("%.0f %%" % (100 * np.nanmedian(v[:, j]))
                                                                        for j in range(len(AMPS)))))
        L.append("| %s, at f | %.0f %% | %s |" % (z, 100 * vol, " | ".join("%.0f %%" % (100 * np.nanmedian(f[:, j]))
                                                                        for j in range(len(AMPS)))))

    L += ["", "## 3. Damping: centre C0 … C7 at α 0.1 and α 0.01 (last 3 cycles; * = not steady)", "",
          "| α | B (mT) | R1 median (calm n) | R1 median, all windows | A = RF/2 median | C = RN/2 median | steady |",
          "|---|---|---|---|---|---|---|"]
    for a in ("", "_a001"):
        rows = {}
        for i in range(8):
            d0 = DD.RUNS / ("C%d%s" % (i, a))
            dd = (AD.rerun_dir(d0) or d0) if a == "" else d0
            if not (dd / "acc.json").exists():
                continue
            sts = DD.acc_stages(dd)
            files = sorted(dd.glob("samples_*_amp*.csv"))
            for k, (s, f) in enumerate(zip(sts, files)):
                ser = EJ.series(f)
                m = EJ.measures(*ser) if ser is not None else None
                rows.setdefault(k, []).append((s["B_peak_mT"], s["R1"], s["calm"], m, EJ.steady(dd, k)))
        for k, v in sorted(rows.items()):
            calm = [x[1] for x in v if x[2]]
            L.append("| %s | %.0f | %s | %.2f | %.3f | %.3f | %d / %d |" % (
                "0.01" if a else "0.1", np.median([x[0] for x in v]),
                ("%.2f (%d)" % (np.median(calm), len(calm))) if calm else "– (0)", np.median([x[1] for x in v]),
                np.median([x[3]["RF"] / 2 for x in v if x[3]]), np.median([x[3]["RN"] / 2 for x in v if x[3]]),
                sum(x[4] for x in v), len(v)))

    L += ["", "## 4. Factor effects on β_f and on the shift (OLS, HC3, range effect; runs with all 7 amplitudes calm)", "",
          "| output | K_eff (lnQ) | r_p | φ | A (− lnd) | lnd × lnQ | R² | n |", "|---|---|---|---|---|---|---|---|"]
    import statsmodels.api as sm
    X = AD.coded(recs)
    D = np.column_stack([np.ones(len(X)), X, X[:, 0] * X[:, 3]])
    for name, y in (("β", beta), ("β_f", bf), ("shift", shift)):
        s = full
        m = sm.OLS(y[s], D[s]).fit(cov_type="HC3")
        c, e = 2 * m.params, 2 * m.bse
        L.append("| %s | %+.2f ± %.2f | %+.2f ± %.2f | %+.2f ± %.2f | %+.2f ± %.2f | %+.2f ± %.2f | %.2f | %d |" % (
            name, c[1], e[1], c[2], e[2], c[3], e[3], -c[4], e[4], c[5], e[5], m.rsquared, s.sum()))
    L += ["", "Note: the identity makes β, β_f and the shift dependent; the shift is not an independent cause.", ""]
    print("\n".join(L))


if __name__ == "__main__":
    main()
