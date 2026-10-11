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
Since 2026-10-11 (PROTOKOLL 7.54, decision Chris): steady-tail rule with at least 2 calm cycles, finished reruns X_mc2
replace X; the main-fit numbers come from results/sobol_analysis_steady2(_mc2).json. New blocks:
  mesh        K3 and base point on both meshes: loss density per job, β power-law fits (eval_k3beta.compute)
  jumps       per amplitude: RF and BF medians of the Sobol runs (dx 3); K3 per mesh: RF, BF, R1; R1 of the Sobol runs
  events      event rate against the cycle number, energy sign of the events, one example run per cycle (S050, 9 mT)

    python deck_data.py            # writes results/deck_data.json
"""
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import eval_events as EV
import eval_jumps_time as EJ
import eval_k3beta as KB
import eval_k3fill as KF
import eval_mechanism as EM

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE / "runs"
F_HZ = 30e6
N_BOOT = 1000
LEX_A10 = 3.342                       # nm, l_ex at A = 10 pJ/m and Js = 1.5 T


def real(r):
    """Real units at d = 1 µm (reading 1 of PLAN 4.0): A from d/l_ex, cube edge L_eff in nm."""
    dl = math.exp(r["lnd"])
    lex = 1000.0 / dl                                         # nm
    return {"d_lex": dl, "A_pJm": 10.0 * (lex / LEX_A10) ** 2, "Leff_lex": math.sqrt(math.exp(-r["lnQ"])),
            "Leff_nm": math.sqrt(math.exp(-r["lnQ"])) * lex, "r_p": r["r_p"], "phi": r["phi"]}


def analysis_file():
    f = HERE / "results" / "sobol_analysis_steady2_mc2.json"
    return f if f.exists() else HERE / "results" / "sobol_analysis_steady2.json"


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2          # steady-tail rule, ≥ 2 calm cycles (7.54)
    recs, skipped = AD.load_runs(RUNS, reruns=True)
    kd = json.loads((RUNS / recs[0]["name"] / "config.json").read_text())["units"]["Kd"]
    to_P = kd * F_HZ / 1e6
    out = {"n_runs": len(recs), "skipped": skipped, "amps_mT": list(AD.AMPS_MT), "Kd": kd, "f_Hz": F_HZ}
    rows = []
    for r in recs:
        q = {"name": r["orig"], "folder": r["name"],
             "kind": "corner" if r["orig"].startswith("K") else ("centre" if r["orig"].startswith("C") else "sobol"),
             "flagged": bool(r["flagged"]), "beta_all": r["beta_all"],
             "beta_w": [r["beta@9-35mT"], r["beta@35-70mT"], r["beta@70-150mT"]],
             "P": [math.exp(r["lnw@%gmT" % m]) * to_P for m in AD.AMPS_MT]}
        q.update(real(r))
        q.update(EM.measures(r["orig"]))                 # snapshot measures from the original run (reruns have none)
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
                   "P": [float(np.nanmean([math.exp(x["lnw@%gmT" % m]) for x in base]) * to_P) for m in AD.AMPS_MT]}

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

    # the same weights as the pre-registered main fit (scatter model, PLAN 4.0; sobol_analysis.json)
    SA = json.loads(analysis_file().read_text())
    out["analysis_file"] = analysis_file().name
    sw = SA["scatter"]["lnw@50mT"]
    wts = np.array(sw["weights"]) if sw["any"] and len(sw["weights"]) == len(recs) else None
    out["weights_used"] = wts is not None
    c0 = AD.glm_coefs(T, W, wts)
    cur0, pa0 = curves(c0), per_amp(c0)
    rng = np.random.default_rng(20261010)
    cb, pb = [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(recs), len(recs))
        try:
            c = AD.glm_coefs(T[i], W[i], None if wts is None else wts[i])
        except Exception:                                      # noqa: BLE001
            continue
        cb.append(curves(c))
        pb.append(per_amp(c))
    lo_hi = {f: (RANGES_REAL(f, grid)) for f in AD.FACTORS}
    out["partial"] = {f: {"x_coded": grid.tolist(), "x_real": lo_hi[f],
                          "est": cur0[f].tolist(),
                          "lo": np.percentile([b[f] for b in cb], 2.5, axis=0).tolist(),
                          "hi": np.percentile([b[f] for b in cb], 97.5, axis=0).tolist()} for f in AD.FACTORS}
    # loss relative to the low end of A (coded d/l_ex = +1, the last grid point), interval per bootstrap sample
    rel = lambda c: c["lnd"][1:] - c["lnd"][1:, -1:]
    out["partial"]["lnd"]["rel_est"] = rel(cur0).tolist()
    out["partial"]["lnd"]["rel_lo"] = np.percentile([rel(b) for b in cb], 2.5, axis=0).tolist()
    out["partial"]["lnd"]["rel_hi"] = np.percentile([rel(b) for b in cb], 97.5, axis=0).tolist()
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
        v = [r["beta_all"] for r in recs if (r["orig"].startswith(tag + "_") if tag != "C" else
                                             (r["orig"].startswith("C") and "_" not in r["orig"]))]
        if len(v) >= 2:
            sc[tag] = {"n": len(v), "sd": float(np.std(v, ddof=1)), "mean": float(np.mean(v))}
    num = sum((v["n"] - 1) * v["sd"] ** 2 for v in sc.values())
    dof = sum(v["n"] - 1 for v in sc.values())
    sc["pooled"] = {"sd": math.sqrt(num / dof), "dof": dof, "n": sum(v["n"] for k, v in sc.items() if k != "pooled")}
    out["scatter"] = sc
    bw = np.array([r["beta_w"] for r in rows])
    out["beta_windows"] = {"mean": bw.mean(0).tolist(), "se": (bw.std(0, ddof=1) / math.sqrt(len(bw))).tolist()}

    # α 0.01 pairs (if present)
    # the α test as evaluated in PROTOKOLL 7.48 (all measured cycles): with the steady-tail rule some amplitudes of the
    # 8 pairs have no 2 calm cycles, and the pair means would be undefined
    analyze.STEADY_TAIL = None
    out["alpha"] = AD.alpha_test(RUNS)
    out["alpha"]["rule"] = "all measured cycles (PROTOKOLL 7.48)"
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
    db = []
    for i in range(8):
        a1, a0 = RUNS / ("C%d_a001" % i), RUNS / ("C%d" % i)
        if (a1 / "DONE").exists() and (a0 / "DONE").exists():
            x1, x0 = AD.run_record(a1), AD.run_record(a0)
            if x1 is not None and x0 is not None:
                db.append((x1["beta_all"] - x0["beta_all"], x1["beta_all"]))
    analyze.STEADY_TAIL = AD.BAL_EVENT
    if db:
        d_ = np.array([x[0] for x in db])
        out["alpha"]["dbeta"] = {"mean": float(d_.mean()), "se": float(d_.std(ddof=1) / math.sqrt(len(d_))), "n": len(d_),
                                 "max_pair": float(d_.max()), "beta_max_a001": float(max(x[1] for x in db))}

    # checks
    acc_ok = all(all(s["check_ok"] for s in json.loads((RUNS / r["orig"] / "acc.json").read_text())["stages"]) for r in recs)
    it = []
    for r in recs:
        for line in (RUNS / r["name"] / "stdout.log").read_text().splitlines():
            if "converged=True after" in line:
                it.append(int(line.split("after ")[1].split()[0]))
    out["checks"] = {"acc_sum_ok": acc_ok, "relax_iter_max": max(it), "relax_limit": 80000,
                     "n_flagged": int(sum(r["flagged"] for r in recs))}
    out["mesh"] = mesh_block(to_P)
    out["jumps"] = jumps_block(recs)
    out["events"] = events_block(recs, to_P)
    (HERE / "results" / "deck_data.json").write_text(json.dumps(out, indent=1))
    print("written results/deck_data.json: %d runs, %d bootstrap resamples, α pairs %d" %
          (len(recs), len(cb), out["alpha"]["n"]))


def mesh_block(to_P):
    """K3 and base point, both meshes: P (W/cm³) and realised B per job and seed; β power-law fits."""
    out = {}
    for lab, (fmt, seeds, amps) in KB.POINTS.items():
        key = "K3" if lab.startswith("K3") else "base"
        blk = {"amps": list(amps), "seeds": list(seeds), "P": {}, "B": {}}
        for dx in ("3", "15"):
            blk["P"][dx], blk["B"][dx] = {}, {}
            for s in seeds:
                v = [KF.loss(AD.rerun_dir(RUNS / (fmt % (dx, s, b))) or RUNS / (fmt % (dx, s, b))) for b in amps]
                blk["P"][dx][str(s)] = [x[0] * to_P if x else None for x in v]
                blk["B"][dx][str(s)] = [x[1] if x else None for x in v]
        out[key] = blk
    out["fits"] = KB.compute()
    for r in out["fits"]:
        for dx in ("3", "15"):
            r[dx]["per_s"] = {str(k): v for k, v in r[dx]["per_s"].items()}
    return out


def jumps_block(recs):
    """RF / BF medians of the Sobol runs per amplitude (last 3 calm cycles), R1 medians (acc.json); K3 per mesh."""
    out = {"amps": list(AD.AMPS_MT), "sobol": {"RF": [], "BF": [], "R1": []}, "K3": {}}
    M = {r["name"]: EJ.run_measures(RUNS / r["name"], len(AD.AMPS_MT)) for r in recs}
    for k in range(len(AD.AMPS_MT)):
        for key in ("RF", "BF"):
            v = [M[r["name"]][k][key] for r in recs if M[r["name"]][k] is not None and M[r["name"]][k]["steady"]]
            out["sobol"][key].append([float(np.median(v)), float(np.percentile(v, 25)), float(np.percentile(v, 75)), len(v)])
    r1 = []
    for r in recs:
        f = RUNS / r["orig"] / "acc.json"
        st = sorted(json.loads(f.read_text())["stages"], key=lambda x: x["B_peak_mT"])
        r1.append([x["R1"] for x in st])
    r1 = np.array(r1)
    out["sobol"]["R1"] = [[float(np.median(r1[:, k])), float(np.percentile(r1[:, k], 25)), float(np.percentile(r1[:, k], 75)),
                           len(r1)] for k in range(r1.shape[1])]
    for b in KF.AMPS:
        e = {}
        for dx in ("3", "15"):
            rf, bf, r1v = [], [], []
            for s in KF.SEEDS:
                d = KF.folder(dx, s, b)
                if (d / "DONE").exists():
                    m = EJ.run_measures(d, 1)[0]
                    if m is not None and m["steady"]:
                        rf.append(m["RF"])
                        bf.append(m["BF"])
                if (d / "acc.json").exists():
                    r1v.append(json.loads((d / "acc.json").read_text())["stages"][0]["R1"])
            e[dx] = {"RF": rf, "BF": bf, "R1": r1v}
        out["K3"][str(b)] = e
    return out


def events_block(recs, to_P):
    """Event rate per cycle number (3 groups of amplitudes), energy sign of the events in cycles 2 … 6, example S050."""
    st = {r["orig"]: json.loads((RUNS / r["orig"] / "summary.json").read_text())["stages"] for r in recs}
    out = {"groups": {}, "release": 0, "uptake": 0}
    for lab, ks in EV.GROUPS:
        out["groups"][lab] = [float(np.mean([EV.is_event(st[n][k]["cycles"][c]) for n in st for k in ks]))
                              for c in range(1, 7)]
        for n in st:
            for k in ks:
                for c in st[n][k]["cycles"][2:]:
                    if EV.is_event(c):
                        if c["w_loop"] <= 0 or c["w_dis"] > c["w_loop"]:
                            out["release"] += 1
                        else:
                            out["uptake"] += 1
    ex = json.loads((RUNS / "S050" / "summary.json").read_text())["stages"][0]["cycles"]
    out["example"] = {"run": "S050", "B_mT": 9, "P_loop": [c["w_loop"] * to_P for c in ex],
                      "P_dis": [c["w_dis"] * to_P for c in ex], "event": [bool(EV.is_event(c)) for c in ex]}
    return out


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
