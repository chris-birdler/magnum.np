#!/usr/bin/env python
"""
Design model against the two direct coarse-mesh tests of the A effect (review of the deck 2026-10-13, finding 1).

Direct tests (A 20 → 10 pJ/m at d = 1 µm, i.e. d/l_ex ≈ 212 → 300, coarse mesh dx 3):
  K1 pair    corner K1 (L_eff 30, K_σ 1.5 K_eff, φ 0.62), paired by seed, n 3 (eval_pair.py); design protocol (7
             amplitudes in one run).
  base point base microstructure (L_eff 12, r_p 0, φ 0.65), M212 against M300, 4 + 4 different seeds (Welch,
             eval_size.py); one job per amplitude, 9 / 50 / 100 mT, each from the relaxed virgin state.
Design model: Gamma GLM of analyze_design.py on all α 0.1 design runs, evaluated at the coded factors of the runs of
each test (φ = φ_vox, ln d/l_ex = ln d_lex_used), 1000 bootstrap resamples (90 % interval). Output: A effect on β =
β(A 20) − β(A 10) from the fit over 9 / 50 / 100 mT and over 9 … 150 mT, and ln P(A 20) − ln P(A 10) at 9 / 50 / 100 mT.

    python eval_model_points.py > results/model_points.md      # also results/model_points.json
"""
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import eval_pair as EP
import eval_size as ES

R = pathlib.Path(__file__).resolve().parent / "runs"
AMPS = list(AD.AMPS_MT)
N_BOOT = 1000


def coded_from_cfg(name):
    """Coded factors of one run folder from its config.json (the same definitions as analyze_design.run_record)."""
    cfg = json.loads((R / name / "config.json").read_text())
    L = cfg["Leff_lex"]
    raw = {"lnQ": math.log(cfg.get("Q_eff") or 1.0 / L ** 2), "r_p": cfg["r_p"], "phi": cfg.get("phi_vox", cfg["phi"]),
           "lnd": math.log(cfg.get("d_lex_used", cfg["d_lex"]))}
    lo_hi = AD.RANGES
    return np.array([2.0 * (raw[f] - lo_hi[f][0]) / (lo_hi[f][1] - lo_hi[f][0]) - 1.0 for f in AD.FACTORS])


def outputs(c, xa, xb):
    """A effect (A 20 = xa, A 10 = xb) from GLM coefficients c: Δβ over 9/50/100 and 9 … 150, Δ ln P at 9/50/100."""
    la, lb = AD.terms(xa[None, :]) @ c.T, AD.terms(xb[None, :]) @ c.T          # ln E[w] at the 7 amplitudes
    lnB = np.log(AMPS)
    out = {}
    for lab, amps in (("beta 9/50/100", (9, 50, 100)), ("beta 9-150", tuple(AMPS))):
        idx = [AMPS.index(b) for b in amps]
        out[lab] = float(np.polyfit(lnB[idx], la[0, idx], 1)[0] - np.polyfit(lnB[idx], lb[0, idx], 1)[0])
    for b in (9, 50, 100):
        out["lnP %d" % b] = float(la[0, AMPS.index(b)] - lb[0, AMPS.index(b)])
    return out


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    recs, _ = AD.load_runs(R, reruns=True)
    X = AD.coded(recs)
    T = AD.terms(X)
    W = np.exp(np.array([[r["lnw@%gmT" % m] for m in AMPS] for r in recs]))
    pts = {"K1 pair": (np.mean([coded_from_cfg("K1_%d" % k) for k in range(3)], 0),
                       np.mean([coded_from_cfg("K1A10_%d" % k) for k in range(3)], 0)),
           "base point": (np.mean([coded_from_cfg("M212_dx3_s%d_b9" % s) for s in range(931, 935)], 0),
                          np.mean([coded_from_cfg("M300_dx3_s%d_b9" % s) for s in range(921, 925)], 0))}
    c0 = AD.glm_coefs(T, W)
    est = {p: outputs(c0, *x) for p, x in pts.items()}
    rng = np.random.default_rng(20261013)
    boots = {p: [] for p in pts}
    for _ in range(N_BOOT):
        i = rng.integers(0, len(recs), len(recs))
        try:
            c = AD.glm_coefs(T[i], W[i])
        except Exception:                               # noqa: BLE001 - a non-converging resample is skipped
            continue
        for p, x in pts.items():
            boots[p].append(outputs(c, *x))
    # measured
    pair = json.loads((pathlib.Path(__file__).resolve().parent / "results" / "pair_K1A10.json").read_text())
    mesh = json.loads((pathlib.Path(__file__).resolve().parent / "results" / "mesh212.json").read_text())
    rec = EP.records()
    lnw = {a: {k: {b: rec[a][k]["lnw@%gmT" % b] for b in AMPS} for k in EP.SEEDS} for a in (20, 10)}
    d3 = [np.polyfit(np.log([9, 50, 100]), [lnw[20][k][b] for b in (9, 50, 100)], 1)[0] -
          np.polyfit(np.log([9, 50, 100]), [lnw[10][k][b] for b in (9, 50, 100)], 1)[0] for k in EP.SEEDS]
    s3 = EP.tci(d3)
    meas = {"K1 pair": {"beta 9/50/100": (s3[0], s3[2], s3[3], "n 3, paired"),
                        "beta 9-150": (pair["lever"]["9 … 150"]["mean"], pair["lever"]["9 … 150"]["lo"],
                                       pair["lever"]["9 … 150"]["hi"], "n 3, paired")},
            "base point": {"beta 9/50/100": (-mesh["lever"]["beta"]["dx3"]["mean"], -mesh["lever"]["beta"]["dx3"]["hi"],
                                             -mesh["lever"]["beta"]["dx3"]["lo"], "n 4 + 4, Welch")}}
    for b in (9, 50, 100):
        r = pair["ratio"]["%g" % b]
        meas["K1 pair"]["lnP %d" % b] = (-math.log(r["gm"]), -math.log(r["hi"]), -math.log(r["lo"]), "n 3, paired")
        lv = mesh["lever"]["lnP %d" % b]["dx3"]
        meas["base point"]["lnP %d" % b] = (-lv["mean"], -lv["hi"], -lv["lo"], "n 4 + 4, Welch")
    L = ["# Design model against the direct coarse-mesh tests of the A effect", "",
         "A effect = value at A 20 pJ/m minus value at A 10 pJ/m (d = 1 µm), coarse mesh. Model: Gamma GLM, %d design "
         "runs, %d bootstrap resamples, 90 %% interval. Measured: 90 %% interval (t; Welch for the base point)."
         % (len(recs), len(boots["K1 pair"])), "",
         "| point | output | measured | design model at the same factors |", "|---|---|---|---|"]
    res = {}
    for p in pts:
        res[p] = {}
        for k in ("beta 9/50/100", "beta 9-150", "lnP 9", "lnP 50", "lnP 100"):
            b = np.array([x[k] for x in boots[p]])
            lo, hi = np.percentile(b, [5, 95])
            m = meas[p].get(k)
            res[p][k] = {"model": [est[p][k], float(lo), float(hi)], "measured": list(m[:3]) if m else None}
            L.append("| %s | %s | %s | %+.2f (%+.2f … %+.2f) |" % (
                p, k, ("%+.2f (%+.2f … %+.2f, %s)" % m) if m else "not measured", est[p][k], lo, hi))
    res["coded"] = {p: {"A20": x[0].tolist(), "A10": x[1].tolist()} for p, x in pts.items()}
    L += ["", "Coded factors (lnQ, r_p, φ, ln d/l_ex): " + "; ".join(
        "%s A 20 %s, A 10 %s" % (p, np.round(x[0], 2).tolist(), np.round(x[1], 2).tolist()) for p, x in pts.items())]
    print("\n".join(L))
    (pathlib.Path(__file__).resolve().parent / "results" / "model_points.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
