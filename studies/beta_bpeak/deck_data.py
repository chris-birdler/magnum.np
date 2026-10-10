#!/usr/bin/env python
"""
Data for the management deck of the Sobol design (slide_decks/beta_powder_sobol/make_deck.py reads
results/deck_data.json). All numbers come from the run folders and results/sobol_analysis.json (PLAN 4.0 analysis).
  runs        per run: factors (reduced and real units at d = 1 µm), β windows and β all, P (W/cm³) at the 7 nominal
              amplitudes, flag, kind; the mechanism measures (switched volume, F90, loss shares)
  base        the base point (B300_s921 … s928): β all mean ± SE, P curve
  partial     main fit (Gamma GLM, 15 terms) predictions of β all and ln E[w] along each factor, the others averaged
              over the design points, 95 % bands by bootstrap over the runs
  per_amp     range effect of each factor on ln E[w] at each amplitude with 95 % CI
  thirds      runs in thirds by β all: mean P per amplitude with 95 % CI
  scatter     realisation scatter (sd of β all) at the centre (C0 … C7) and the corners
  alpha       α 0.01 / α 0.1 pairs at the centre (if present)

    python deck_data.py            # writes results/deck_data.json
"""
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import eval_mechanism as EM

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE / "runs"
F_HZ = 30e6
N_BOOT = 500
LEX_A10 = 3.342                       # nm, l_ex at A = 10 pJ/m and Js = 1.5 T


def real(r):
    """Real units at d = 1 µm (reading 1 of PLAN 4.0): A from d/l_ex, cube edge L_eff in nm."""
    dl = math.exp(r["lnd"])
    lex = 1000.0 / dl                                         # nm
    return {"d_lex": dl, "A_pJm": 10.0 * (lex / LEX_A10) ** 2, "Leff_lex": math.sqrt(math.exp(-r["lnQ"])),
            "Leff_nm": math.sqrt(math.exp(-r["lnQ"])) * lex, "r_p": r["r_p"], "phi": r["phi"]}


def main():
    recs, skipped = AD.load_runs(RUNS)
    kd = json.loads((RUNS / recs[0]["name"] / "config.json").read_text())["units"]["Kd"]
    to_P = kd * F_HZ / 1e6
    out = {"n_runs": len(recs), "skipped": skipped, "amps_mT": list(AD.AMPS_MT), "Kd": kd, "f_Hz": F_HZ}
    rows = []
    for r in recs:
        q = {"name": r["name"], "kind": "corner" if r["name"].startswith("K") else ("centre" if r["name"].startswith("C") else "sobol"),
             "flagged": bool(r["flagged"]), "beta_all": r["beta_all"],
             "beta_w": [r["beta@9-35mT"], r["beta@35-70mT"], r["beta@70-150mT"]],
             "P": [math.exp(r["lnw@%gmT" % m]) * to_P for m in AD.AMPS_MT]}
        q.update(real(r))
        q.update(EM.measures(r["name"]))
        rows.append(q)
    out["runs"] = rows

    # base point
    base = []
    for d in sorted(RUNS.glob("B300_s92?")):
        x = AD.run_record(d)
        if x is not None:
            base.append(x)
    bb = np.array([x["beta_all"] for x in base])
    out["base"] = {"n": len(base), "beta_mean": float(bb.mean()), "beta_se": float(bb.std(ddof=1) / math.sqrt(len(bb))),
                   "P": [float(np.mean([math.exp(x["lnw@%gmT" % m]) for x in base]) * to_P) for m in AD.AMPS_MT]}

    # main fit along each factor
    X = AD.coded(recs)
    T = AD.terms(X)
    W = np.exp(np.array([[r["lnw@%gmT" % m] for m in AD.AMPS_MT] for r in recs]))
    grid = np.linspace(-1, 1, 21)

    def curves(coefs):
        res = {}
        for j, f in enumerate(AD.FACTORS):
            ba, l10, l50, l100 = [], [], [], []
            for g in grid:
                Xg = X.copy()
                Xg[:, j] = g
                lnmu = AD.terms(Xg) @ coefs.T
                o = AD.outputs_from_lnmu(lnmu)
                ba.append(o["beta_all"].mean())
                l10.append(o["lnw@10mT"].mean())
                l50.append(o["lnw@50mT"].mean())
                l100.append(o["lnw@100mT"].mean())
            res[f] = np.array([ba, l10, l50, l100])
        return res

    def per_amp(coefs):
        res = {}
        for j, f in enumerate(AD.FACTORS):
            Th, Tl = AD.contrast(X, j)
            res[f] = (Th @ coefs.T - Tl @ coefs.T).mean(0)
        return res

    c0 = AD.glm_coefs(T, W)
    cur0, pa0 = curves(c0), per_amp(c0)
    rng = np.random.default_rng(20261010)
    cb, pb = [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(recs), len(recs))
        try:
            c = AD.glm_coefs(T[i], W[i])
        except Exception:                                      # noqa: BLE001
            continue
        cb.append(curves(c))
        pb.append(per_amp(c))
    lo_hi = {f: (RANGES_REAL(f, grid)) for f in AD.FACTORS}
    out["partial"] = {f: {"x_coded": grid.tolist(), "x_real": lo_hi[f],
                          "est": cur0[f].tolist(),
                          "lo": np.percentile([b[f] for b in cb], 2.5, axis=0).tolist(),
                          "hi": np.percentile([b[f] for b in cb], 97.5, axis=0).tolist()} for f in AD.FACTORS}
    out["per_amp"] = {f: {"est": pa0[f].tolist(),
                          "lo": np.percentile([b[f] for b in pb], 2.5, axis=0).tolist(),
                          "hi": np.percentile([b[f] for b in pb], 97.5, axis=0).tolist()} for f in AD.FACTORS}
    out["n_boot"] = len(cb)

    # thirds by β all
    b = np.array([r["beta_all"] for r in recs])
    P = W * to_P
    q1, q2 = np.quantile(b, [1 / 3, 2 / 3])
    groups = {"low": b <= q1, "mid": (b > q1) & (b <= q2), "high": b > q2}
    th = {}
    for k, g in groups.items():
        idx = np.where(g)[0]
        bs = [P[rng.choice(idx, len(idx))].mean(0) for _ in range(1000)]
        th[k] = {"n": int(g.sum()), "beta_mean": float(b[g].mean()), "P": P[g].mean(0).tolist(),
                 "lo": np.percentile(bs, 2.5, axis=0).tolist(), "hi": np.percentile(bs, 97.5, axis=0).tolist()}
    rb = []
    for _ in range(1000):
        il = rng.choice(np.where(groups["low"])[0], groups["low"].sum())
        ih = rng.choice(np.where(groups["high"])[0], groups["high"].sum())
        rb.append(P[ih].mean(0) / P[il].mean(0))
    th["ratio"] = {"est": (P[groups["high"]].mean(0) / P[groups["low"]].mean(0)).tolist(),
                   "lo": np.percentile(rb, 2.5, axis=0).tolist(), "hi": np.percentile(rb, 97.5, axis=0).tolist()}
    th["limits"] = [float(q1), float(q2)]
    out["thirds"] = th

    # realisation scatter
    sc = {}
    for tag in ("C", "K1", "K2", "K3", "K4"):
        v = [r["beta_all"] for r in recs if (r["name"].startswith(tag + "_") if tag != "C" else
                                             (r["name"].startswith("C") and "_" not in r["name"]))]
        if len(v) >= 2:
            sc[tag] = {"n": len(v), "sd": float(np.std(v, ddof=1)), "mean": float(np.mean(v))}
    out["scatter"] = sc

    # α 0.01 pairs (if present)
    out["alpha"] = AD.alpha_test(RUNS)
    r1 = {}
    for tag, suf in (("0.1", ""), ("0.01", "_a001")):
        v = []
        for i in range(8):
            f = RUNS / ("C%d%s" % (i, suf)) / "acc.json"
            if f.exists():
                st = json.loads(f.read_text())["stages"]
                v.append([min(st, key=lambda s: abs(s["B_peak_mT"] / m - 1))["R1"] for m in (9, 50, 150)])
        if v:
            v = np.array(v)
            r1[tag] = {"mT": [9, 50, 150], "mean": v.mean(0).tolist(), "se": (v.std(0, ddof=1) / math.sqrt(len(v))).tolist(),
                       "n": len(v)}
    out["alpha"]["R1"] = r1

    # checks
    acc_ok = all(all(s["check_ok"] for s in json.loads((RUNS / r["name"] / "acc.json").read_text())["stages"]) for r in recs)
    it = []
    for r in recs:
        for line in (RUNS / r["name"] / "stdout.log").read_text().splitlines():
            if "converged=True after" in line:
                it.append(int(line.split("after ")[1].split()[0]))
    out["checks"] = {"acc_sum_ok": acc_ok, "relax_iter_max": max(it), "relax_limit": 80000,
                     "n_flagged": int(sum(r["flagged"] for r in recs))}
    (HERE / "results" / "deck_data.json").write_text(json.dumps(out, indent=1))
    print("written results/deck_data.json: %d runs, %d bootstrap resamples, α pairs %d" %
          (len(recs), len(cb), out["alpha"]["n"]))


def RANGES_REAL(f, grid):
    """Real values of the coded grid: d/l_ex → A (pJ/m) at d = 1 µm; ln Q_eff → L_eff (l_ex); r_p; φ."""
    lo, hi = AD.RANGES[f]
    v = lo + (grid + 1) / 2 * (hi - lo)
    if f == "lnd":
        return [10.0 * ((1000.0 / math.exp(x)) / LEX_A10) ** 2 for x in v]
    if f == "lnQ":
        return [math.sqrt(math.exp(-x)) for x in v]
    return v.tolist()


if __name__ == "__main__":
    main()
