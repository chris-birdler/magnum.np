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
import eval_cond as EC
import eval_size as ES

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE / "runs"
F_HZ = 30e6
N_BOOT = 1000
LEX_A10 = 3.342                       # nm, l_ex at A = 10 pJ/m and Js = 1.5 T


def acc_stages(d):
    """acc.json stages of run folder d, sorted by B, each with "calm": True if no cycle of the accumulator window
    (cycles_from … cycles_from + n_samples / 256 − 1) is an event cycle (eval_events.is_event, PROTOKOLL 7.53).
    acc_reduce.py sums all cycles of the window, so a window with an event cycle mixes the event into R1 and into the
    loss shares; the steady rule for P uses calm cycles only. Thus only calm windows are used (PROTOKOLL 7.70)."""
    acc = json.loads((d / "acc.json").read_text())["stages"]
    sm = {s["name"]: s for s in json.loads((d / "summary.json").read_text())["stages"]}
    out = []
    for a in sorted(acc, key=lambda s: s["B_peak_mT"]):
        c0, n = a["cycles_from"], a["n_samples"] // EJ.N_SAMPLES
        cyc = [c for c in sm[a["stage"]]["cycles"] if c0 <= c["cycle"] < c0 + n]
        assert len(cyc) == n, (d.name, a["stage"], c0, n)
        out.append(dict(a, calm=not any(EV.is_event(c) for c in cyc)))
    return out


def acc_calm(recs, value):
    """[run, amplitude] array of value(stage) from the calm accumulator windows, NaN for the others."""
    return np.array([[value(s) if s["calm"] else np.nan for s in acc_stages(RUNS / r["name"])] for r in recs])


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
    cur0, pa0, cc0, lc0 = curves(c0), per_amp(c0), cond_curves(c0, X, grid), loss_curves(c0, X)
    rng = np.random.default_rng(20261010)
    cb, pb, ccb, lcb = [], [], [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(recs), len(recs))
        try:
            c = AD.glm_coefs(T[i], W[i], None if wts is None else wts[i])
        except Exception:                                      # noqa: BLE001
            continue
        cb.append(curves(c))
        pb.append(per_amp(c))
        ccb.append(cond_curves(c, X, grid))
        lcb.append(loss_curves(c, X))
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
    out["cond"] = {"grid": grid.tolist(), "K_levels_coded": list(K_LEVELS), "A_levels_coded": list(A_LEVELS),
                   "A_real": RANGES_REAL("lnd", grid), "Leff_real": RANGES_REAL("lnQ", grid)}
    for k in cc0:
        arr = np.array([b[k] for b in ccb])
        out["cond"][k] = {"est": np.asarray(cc0[k]).tolist(), "lo": np.percentile(arr, 2.5, axis=0).tolist(),
                          "hi": np.percentile(arr, 97.5, axis=0).tolist()}
    out["map"] = map_grid(c0, X)
    out["loss_curves"] = {"amps": list(AD.AMPS_MT), "to_P": to_P}
    for k in lc0:                                            # geometric mean of E[w] over the design points
        arr = np.array([b[k] for b in lcb])
        out["loss_curves"][k] = {"est": lc0[k].tolist(), "lo": np.percentile(arr, 2.5, axis=0).tolist(),
                                 "hi": np.percentile(arr, 97.5, axis=0).tolist()}

    # thirds by β all
    b = np.array([r["beta_all"] for r in recs])
    P = W * to_P
    q1, q2 = np.quantile(b, [1 / 3, 2 / 3])
    groups = {"low": b <= q1, "mid": (b > q1) & (b <= q2), "high": b > q2}
    th = {}
    gm = lambda x: np.exp(np.log(x).mean(0))                # geometric mean, as the other loss curves of the deck
    for k, g in groups.items():
        idx = np.where(g)[0]
        bs = [gm(P[rng.choice(idx, len(idx))]) for _ in range(1000)]
        th[k] = {"n": int(g.sum()), "beta_mean": float(b[g].mean()), "P": gm(P[g]).tolist(),
                 "lo": np.percentile(bs, 2.5, axis=0).tolist(), "hi": np.percentile(bs, 97.5, axis=0).tolist()}
    rb = []
    for _ in range(1000):
        il = rng.choice(np.where(groups["low"])[0], groups["low"].sum())
        ih = rng.choice(np.where(groups["high"])[0], groups["high"].sum())
        rb.append(gm(P[ih]) / gm(P[il]))
    th["ratio"] = {"est": (gm(P[groups["high"]]) / gm(P[groups["low"]])).tolist(),
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
    sdl = []                                                  # loop-to-loop scatter of ln w in the steady tail (≥ 3)
    for r in recs:
        for st in json.loads((RUNS / r["name"] / "summary.json").read_text())["stages"]:
            skip = 2 if "h_amp_initial" in st else 1
            tail = analyze.steady_tail(st["cycles"][skip:], AD.BAL_EVENT)
            if len(tail) >= 3:
                sdl.append(float(np.std(np.log([c["w_loop"] for c in tail]), ddof=1)))
    out["loop_scatter"] = {"median": float(np.median(sdl)), "n": len(sdl)}
    out["mesh"] = mesh_block(to_P)
    out["jumps"] = jumps_block(recs)
    out["events"] = events_block(recs, to_P)
    out["beta_fits"] = beta_fits_block(recs)
    out["physics"] = physics_block(recs)
    out["channels"] = channels_block(recs, to_P)
    out["alpha_steady"] = alpha_block(to_P)
    out["conditioning"] = cond_block()
    (HERE / "results" / "deck_data.json").write_text(json.dumps(out, indent=1))
    print("written results/deck_data.json: %d runs, %d bootstrap resamples, α pairs %d" %
          (len(recs), len(cb), out["alpha"]["n"]))


K_LEVELS = (-1.0, 0.0, 1.0)            # coded ln Q_eff: K_eff 1.0 / 2.5 / 6.2 kJ/m³ (L_eff 30 / 19 / 12 l_ex)
A_LEVELS = (1.0, -1.0)                 # coded ln d/l_ex: A 10 / 20 pJ/m at d = 1 µm
IL = {f: j for j, f in enumerate(AD.FACTORS)}


def _outputs_at(coefs, X, **coded):
    """Main-fit outputs (β all, β windows, ln E[w] at 10 / 50 / 100 mT) with some factors set (coded), the other
    factors at the design points; mean over the design points (definition of the range effect)."""
    Xg = X.copy()
    for f, v in coded.items():
        Xg[:, IL[f]] = v
    o = AD.outputs_from_lnmu(AD.terms(Xg) @ coefs.T)
    return {k: float(np.mean(v)) for k, v in o.items()}


def cond_curves(coefs, X, grid):
    """Conditional curves (A × K_eff interaction, PROTOKOLL 7.66): along ln d/l_ex at 3 K_eff levels and along
    ln Q_eff at 2 A levels; the real-alloy path K_eff ∝ A⁻³ from (A 10 pJ/m, K_eff 6.2 kJ/m³) to K_eff 1.0 kJ/m³
    (A ≈ 18.4 pJ/m), inside the design."""
    out = {}
    keys = ("beta_all", "beta@9-35mT", "beta@35-70mT", "beta@70-150mT", "lnw@10mT", "lnw@50mT", "lnw@100mT")
    for kq in K_LEVELS:
        rows = [_outputs_at(coefs, X, lnd=g, lnQ=kq) for g in grid]
        for k in keys:
            out["A|%+.0f|%s" % (kq, k)] = [r[k] for r in rows]
        for k in ("lnw@10mT", "lnw@50mT", "lnw@100mT"):   # relative to d/l_ex 300 (A 10 pJ/m, d 1 µm): last point
            out["Arel|%+.0f|%s" % (kq, k)] = [r[k] - rows[-1][k] for r in rows]
    for ka in A_LEVELS:
        rows = [_outputs_at(coefs, X, lnQ=g, lnd=ka) for g in grid]
        for k in keys:
            out["K|%+.0f|%s" % (ka, k)] = [r[k] for r in rows]
    lo_d, hi_d = AD.RANGES["lnd"]
    lo_q, hi_q = AD.RANGES["lnQ"]
    path = []
    for s in np.linspace(0, 1, 11):                    # ln Q from its max down to its min, A from 10 pJ/m up
        lnq = hi_q - s * (hi_q - lo_q)                  # K_eff ∝ A⁻³ and A ∝ (d/l_ex)⁻² → ln Q = 6 ln(d/l_ex) + c
        lnd = hi_d - (hi_q - lnq) / 6.0
        path.append(_outputs_at(coefs, X, lnd=2 * (lnd - lo_d) / (hi_d - lo_d) - 1, lnQ=2 * (lnq - lo_q) / (hi_q - lo_q) - 1))
    for k in ("beta_all", "lnw@10mT", "lnw@50mT", "lnw@100mT"):
        out["alloy|%s" % k] = [r[k] for r in path]
    # conditional range effects on β all (physical direction): A 10 → 20 pJ/m = coded ln d/l_ex +1 → −1;
    # K_eff low → high = coded ln Q −1 → +1; the real-alloy path end − start
    for kq in K_LEVELS:
        b_ = out["A|%+.0f|beta_all" % kq]
        out["dA|%+.0f" % kq] = b_[0] - b_[-1]
    for ka in A_LEVELS:
        b_ = out["K|%+.0f|beta_all" % ka]
        out["dK|%+.0f" % ka] = b_[-1] - b_[0]
    out["dalloy"] = path[-1]["beta_all"] - path[0]["beta_all"]
    return out


def _mean_curve(coefs, X, **coded):
    """Mean over the design points of ln E[w] at the 7 amplitudes (geometric mean, the same averaging as the range
    effects), some factors set (coded)."""
    Xg = X.copy()
    for f, v in coded.items():
        Xg[:, IL[f]] = v
    return (AD.terms(Xg) @ coefs.T).mean(0)


def loss_curves(coefs, X):
    """Mean loss curves (W/cm³ after × to_P): each factor at its low and high end ('ends|f|-1/+1'), and the 4 corners of
    A × K_eff ('corner|lnd|lnQ'), the other factors at the design points."""
    out = {}
    for f in AD.FACTORS:
        for v in (-1.0, 1.0):
            out["ends|%s|%+.0f" % (f, v)] = _mean_curve(coefs, X, **{f: v})
    for xd in (-1.0, 1.0):
        for xq in (-1.0, 1.0):
            out["corner|%+.0f|%+.0f" % (xd, xq)] = _mean_curve(coefs, X, lnd=xd, lnQ=xq)
    return out


def map_grid(coefs, X, n=31):
    """β all and ln E[w] at 10 / 100 mT on a grid of A (ln d/l_ex) × K_eff (ln Q), r_p and φ at the design points."""
    g = np.linspace(-1, 1, n)
    res = {k: np.zeros((n, n)) for k in ("beta_all", "lnw@10mT", "lnw@100mT")}
    for i, xd in enumerate(g):
        for j, xq in enumerate(g):
            o = _outputs_at(coefs, X, lnd=xd, lnQ=xq)
            for k in res:
                res[k][j, i] = o[k]
    lnw7 = np.zeros((len(AD.AMPS_MT), n, n))
    for i, xd in enumerate(g):
        for j, xq in enumerate(g):
            lnw7[:, j, i] = _mean_curve(coefs, X, lnd=xd, lnQ=xq)
    return {"A_real": RANGES_REAL("lnd", g), "Leff_real": RANGES_REAL("lnQ", g), "lnw7_meancurve": lnw7.tolist(),
            **{k: v.tolist() for k, v in res.items()}}


WINDOWS5 = ((9, 20, 35), (20, 35, 50), (35, 50, 70), (50, 70, 100), (70, 100, 150))


def beta_fits_block(recs):
    """Power-law fits per run (ln w over ln B at the nominal amplitudes, as the design analysis): β 9 … 150 mT and
    in 5 windows of 3 amplitudes; mean ± SE, sd, share above 2; paired difference of the last and the first window."""
    lnB = np.log(AD.AMPS_MT)
    L = np.array([[r["lnw@%gmT" % m] for m in AD.AMPS_MT] for r in recs])
    out = {"windows": [], "n": len(recs)}
    for win in ((9, 150),) + WINDOWS5:
        idx = [i for i, m in enumerate(AD.AMPS_MT) if win[0] <= m <= win[-1]]
        b = np.array([np.polyfit(lnB[idx], L[k, idx], 1)[0] for k in range(len(recs))])
        out["windows"].append({"range": [win[0], win[-1]], "mean": float(b.mean()), "se": float(b.std(ddof=1) / math.sqrt(len(b))),
                               "sd": float(b.std(ddof=1)), "share_gt2": float((b > 2).mean())})
        if win == WINDOWS5[0]:
            first = b
        if win == WINDOWS5[-1]:
            d = b - first
            out["last_minus_first"] = {"mean": float(d.mean()), "se": float(d.std(ddof=1) / math.sqrt(len(d)))}
    return out


def physics_block(recs):
    """Loop shape and loss location (audit 2026-10-12, audit_20261012/timecurves.py functions): per run and amplitude
    μ_r and sin δ of the fundamental (kept cycles of the main rule); medians per A third; loss shares from acc.json;
    per-particle loss ratio by stress-axis angle (37° against 79°) for runs with r_p > 0.5."""
    import sys
    sys.path.insert(0, str(HERE / "audit_20261012"))
    import timecurves as TC                                     # noqa: E402
    R = [TC.analyse_run(r["orig"]) for r in recs]
    mu = np.array([[s["mu_r"] for s in x["stages"]] for x in R])
    sd = np.array([[s["sin_d"] for s in x["stages"]] for x in R])
    ident = np.array([[s["W_fund"] / s["W_loop"] for s in x["stages"]] for x in R])
    A = np.array([real(r)["A_pJm"] for r in recs])
    q1, q2 = np.quantile(A, [1 / 3, 2 / 3])
    out = {"amps": list(AD.AMPS_MT), "A_thirds": [float(q1), float(q2)], "identity_minmax": [float(ident.min()), float(ident.max())],
           "mu_median": np.median(mu, 0).tolist(), "sin_median": np.median(sd, 0).tolist(), "groups": {}}
    for lab, sel in (("low A", A <= q1), ("high A", A > q2)):
        out["groups"][lab] = {"n": int(sel.sum()), "A": [float(A[sel].min()), float(A[sel].max())],
                              "mu": np.median(mu[sel], 0).tolist(), "sin": np.median(sd[sel], 0).tolist(),
                              "sin_q": [np.percentile(sd[sel], 25, 0).tolist(), np.percentile(sd[sel], 75, 0).tolist()]}
    keys = ("r<0.5", "r>=0.8", "p top 10%", "texture top 10%")
    sh = {k: acc_calm(recs, lambda s, k=k: s["where"][k]["p_share"]) for k in keys}       # calm windows only
    vs = {k: acc_calm(recs, lambda s, k=k: s["where"][k]["vol_share"]) for k in keys}
    out["where_n"] = np.sum(np.isfinite(sh["r<0.5"]), 0).tolist()
    out["where"] = {k: np.nanmedian(v, 0).tolist() for k, v in sh.items()}
    out["where_q"] = {k: [np.nanpercentile(v, 25, 0).tolist(), np.nanpercentile(v, 75, 0).tolist()] for k, v in sh.items()}
    out["where_min"] = {k: np.nanmin(v, 0).tolist() for k, v in sh.items()}
    # share of the loss of the motion at f (fundamental) in each zone: what the snapshot maps show
    out["where_fund"] = {k: np.nanmedian(acc_calm(recs, lambda s, k=k: s["where"][k]["p_fund_share"]), 0).tolist()
                         for k in keys}
    out["where_volume"] = {"r<0.5": 0.125, "r>=0.8": 0.488}   # volume shares of a sphere: 0.5³ and 1 − 0.8³
    out["where_volume_cells"] = {k: float(np.nanmedian(v)) for k, v in vs.items()}   # the same on the cell grid
    # fraction of the runs in which the centre has a higher loss density than the rim, per amplitude
    ok = np.isfinite(sh["r<0.5"])
    out["centre_gt_rim"] = [float(np.mean((sh["r<0.5"][:, j] / 0.125 > sh["r>=0.8"][:, j] / 0.488)[ok[:, j]]))
                            for j in range(ok.shape[1])]
    # loss density in a region, per volume of particle material (W/cm³): share × P(core volume) / (φ × volume share)
    Pcore = np.array([[s["W_loop"] * F_HZ / 1e6 for s in x["stages"]] for x in R])
    phi = np.array([x["cfg"]["phi_vox"] for x in R])
    out["P_particle_median"] = np.median(Pcore / phi[:, None], 0).tolist()
    out["region_density"] = {k: np.nanmedian(sh[k] * Pcore / (phi[:, None] * out["where_volume"][k]), 0).tolist()
                             for k in ("r<0.5", "r>=0.8")}
    q = []
    for x in R:
        if x["cfg"]["r_p"] < 0.5:
            continue
        cp = x["cos_p"]
        hi, lo = np.abs(cp) > 0.5, np.abs(cp) < 0.5
        q.append([math.log(s["w_part"][hi].sum() / s["w_part"][lo].sum()) if np.all(s["w_part"] > 0) else np.nan
                  for s in x["stages"]])
    q = np.array(q)
    out["angle_ratio"] = {"n": int(len(q)), "mean": np.nanmean(q, 0).tolist(),
                          "se": (np.nanstd(q, 0, ddof=1) / np.sqrt(np.sum(np.isfinite(q), 0))).tolist()}
    p37, p79 = [], []                                           # per particle volume, W/cm³ (w_part in J/m³ per loop)
    for x in R:
        if x["cfg"]["r_p"] < 0.5:
            continue
        cp = x["cos_p"]
        hi, lo = np.abs(cp) > 0.5, np.abs(cp) < 0.5
        p37.append([s["w_part"][hi].mean() * F_HZ / 1e6 for s in x["stages"]])
        p79.append([s["w_part"][lo].mean() * F_HZ / 1e6 for s in x["stages"]])
    out["angle_P"] = {"37": np.median(p37, 0).tolist(), "79": np.median(p79, 0).tolist(),
                      "37_q": np.percentile(p37, [25, 75], axis=0).tolist(), "79_q": np.percentile(p79, [25, 75], axis=0).tolist()}
    return out


def alpha_block(to_P):
    """Damping test with the steady rule (≥ 2 calm cycles, decision Chris 2026-10-12): ln P(α 0.01) − ln P(α 0.1) per
    seed and amplitude (C0 … C7); amplitudes with fewer than 3 seeds have no value."""
    out = {"amps": list(AD.AMPS_MT), "mean": [], "se": [], "n": []}
    rows, abs_rows = [], []
    for i in range(8):
        a1, a0 = RUNS / ("C%d_a001" % i), RUNS / ("C%d" % i)
        if not ((a1 / "DONE").exists() and (a0 / "DONE").exists()):
            continue
        x1, x0 = AD.run_record(a1), AD.run_record(AD.rerun_dir(a0) or a0)
        if x1 is None or x0 is None:
            continue
        rows.append([x1["lnw@%gmT" % m] - x0["lnw@%gmT" % m] for m in AD.AMPS_MT])
        abs_rows.append(([x0["lnw@%gmT" % m] for m in AD.AMPS_MT], [x1["lnw@%gmT" % m] for m in AD.AMPS_MT]))
    v = np.array(rows, float)
    out["P_a01"], out["P_a001"] = [], []
    for k in range(len(AD.AMPS_MT)):
        ok = [r for r in abs_rows if np.isfinite(r[0][k]) and np.isfinite(r[1][k])]
        out["P_a01"].append(float(np.exp(np.mean([r[0][k] for r in ok]))) * to_P if len(ok) >= 3 else None)
        out["P_a001"].append(float(np.exp(np.mean([r[1][k] for r in ok]))) * to_P if len(ok) >= 3 else None)
    for k in range(len(AD.AMPS_MT)):
        c = v[:, k][np.isfinite(v[:, k])] if len(v) else np.array([])
        out["n"].append(int(len(c)))
        out["mean"].append(float(c.mean()) if len(c) >= 3 else None)
        out["se"].append(float(c.std(ddof=1) / math.sqrt(len(c))) if len(c) >= 3 else None)
    return out


def cond_block():
    """Conditioning test (PLAN 4.5 D, eval_cond.py): drift per seed and amplitude, decision, β design and late."""
    out = {"amps": list(EC.AMPS), "drift": {}, "events_cycles": []}
    V = {(s, b): EC.values(EC.tested_stage("cond_C%d_b%d" % (s, b))) for s in EC.SEEDS for b in EC.AMPS
         if (EC.R / ("cond_C%d_b%d" % (s, b)) / "DONE").exists()}
    for b in EC.AMPS:
        d = [math.log(V[(s, b)][1] / V[(s, b)][0]) for s in EC.SEEDS if (s, b) in V and None not in V[(s, b)][:2]]
        n, m, se, c90, c95, cls = EC.decide(d)
        out["drift"][str(b)] = {"per_seed": d, "n": n, "mean": m, "se": se, "ci90": c90, "decision": cls}
    for (s, b), v in V.items():
        out["events_cycles"] += v[3]
    kd = json.loads((EC.R / "cond_C0_b9" / "config.json").read_text())["units"]["Kd"]
    out["P_loops"] = {str(b): {str(s): [c["w_loop"] * kd * F_HZ / 1e6 for c in EC.tested_stage("cond_C%d_b%d" % (s, b))["cycles"]]
                               for s in EC.SEEDS if (s, b) in V} for b in EC.AMPS}
    return out


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
    r1 = acc_calm(recs, lambda s: s["R1"])                    # folder of the steady analysis, calm windows (7.70)
    out["sobol"]["R1"] = [[float(np.nanmedian(r1[:, k])), float(np.nanpercentile(r1[:, k], 25)),
                           float(np.nanpercentile(r1[:, k], 75)), int(np.isfinite(r1[:, k]).sum())]
                          for k in range(r1.shape[1])]
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
                    s0 = acc_stages(d)[0]
                    if s0["calm"]:
                        r1v.append(s0["R1"])
            e[dx] = {"RF": rf, "BF": bf, "R1": r1v}
        out["K3"][str(b)] = e
    return out


def channels_block(recs, to_P):
    """Loss channels per run and amplitude (PROTOKOLL 7.70; Chris 2026-10-12, options A + B + C):
      f      R1 = share of the cell loss from the local motion at the drive frequency (acc.json, calm windows only)
      jumps  A = energy above the smooth periodic curve of the box dissipation p_dis(t) = RF / 2 (exact: P − Ps sums
             to zero over the kept cycles, Ps = mean cycle with Fourier k ≤ 6); bursts shorter than ≈ T / 13
      other  1 − R1 − A: smooth local motion at other frequencies (harmonics, slow non-periodic motion)
      C      non-repeating energy (cycle to cycle) = RN / 2, part of the jumps and of the slow non-periodic motion
      B      energy share Eev and time share ev of the sample intervals in which n60 (neighbour pairs > 60°: vortex
             cores, Bloch points, steep walls) changes
    A, B, C from the last 3 calm cycles (eval_jumps_time.py, steady-tail rule), R1 from the calm accumulator window.
    Means over the runs (the three channels add to 1 in each run), SE, quartiles; per A third (as physics_block).
    Exponents: power-law fit 9 … 150 mT of the loss density in each channel, runs with all 7 amplitudes valid."""
    from scipy import stats
    n = len(AD.AMPS_MT)
    M = {r["name"]: EJ.run_measures(RUNS / r["name"], n) for r in recs}
    get = lambda key: np.array([[m[key] if m is not None and m["steady"] else np.nan for m in M[r["name"]]] for r in recs])
    R1 = acc_calm(recs, lambda s: s["R1"])
    J = get("RF") / 2.0
    C = get("RN") / 2.0
    ch = {"f": R1, "jumps": J, "other": 1.0 - R1 - J, "C": C, "Eev": get("Eev"), "ev": get("ev")}
    A = np.array([real(r)["A_pJm"] for r in recs])
    q1, q2 = np.quantile(A, [1 / 3, 2 / 3])
    P = np.array([[math.exp(r["lnw@%gmT" % m]) * to_P for m in AD.AMPS_MT] for r in recs])
    lnB = np.log(AD.AMPS_MT)
    beta = np.array([np.polyfit(lnB, np.log(p), 1)[0] for p in P])
    out = {"amps": list(AD.AMPS_MT), "groups": {}}
    for lab, sel in (("all", A > 0), ("low A", A <= q1), ("high A", A > q2)):
        g = {}
        for k, v in ch.items():
            x = v[sel]
            ok = np.isfinite(x)
            g[k] = {"mean": np.nanmean(x, 0).tolist(), "median": np.nanmedian(x, 0).tolist(),
                    "se": (np.nanstd(x, 0, ddof=1) / np.sqrt(ok.sum(0))).tolist(),
                    "q": [np.nanpercentile(x, 25, 0).tolist(), np.nanpercentile(x, 75, 0).tolist()],
                    "n": ok.sum(0).tolist()}
        # mean channel shares on the runs with all three channels valid at that amplitude (they add to 1)
        both = np.isfinite(R1[sel]) & np.isfinite(J[sel])
        g["stack"] = {k: [float(ch[k][sel][both[:, j], j].mean()) for j in range(n)] for k in ("f", "other", "jumps")}
        g["stack_n"] = both.sum(0).tolist()
        full = both.all(1)
        ex = {}
        for k in ("f", "other", "jumps"):
            e = np.array([np.polyfit(lnB, np.log(np.clip(ch[k][sel][i], 1e-6, None) * P[sel][i]), 1)[0]
                          for i in np.where(full)[0]])
            ex[k] = {"mean": float(e.mean()), "se": float(e.std(ddof=1) / math.sqrt(len(e))), "n": int(len(e))}
        ex["total"] = {"mean": float(beta[sel][full].mean()), "se": float(beta[sel][full].std(ddof=1) / math.sqrt(full.sum())),
                       "n": int(full.sum())}
        g["exponent"] = ex
        g["n_runs"] = int(sel.sum())
        out["groups"][lab] = g
    ok = np.isfinite(R1[:, 0]) & np.isfinite(R1[:, -1])
    rho, pv = stats.spearmanr(beta[ok], np.log(R1[ok, -1] / R1[ok, 0]))
    out["beta_vs_R1drop"] = {"rho": float(rho), "p": float(pv), "n": int(ok.sum())}   # partly by construction (7.71)
    # identity β = β_f + shift, β_f = exponent of R1 · P, shift = −d ln R1 / d ln B (PROTOKOLL 7.71)
    full = np.isfinite(R1).all(1)
    bf = np.array([np.polyfit(lnB, np.log(R1[i] * P[i]), 1)[0] if full[i] else np.nan for i in range(len(recs))])
    sh = np.array([-np.polyfit(lnB, np.log(R1[i]), 1)[0] if full[i] else np.nan for i in range(len(recs))])
    ms = lambda v: {"mean": float(np.mean(v)), "se": float(np.std(v, ddof=1) / math.sqrt(len(v))), "n": int(len(v))}
    dec = {}
    for lab, sel in (("all", A > 0), ("low A", A <= q1), ("high A", A > q2)):
        s = sel & full
        dec[lab] = {"beta": ms(beta[s]), "beta_f": ms(bf[s]), "shift": ms(sh[s])}
    lo, hi = (A <= q1) & full, (A > q2) & full
    dec["diff"] = {k: [float(np.mean(v[hi]) - np.mean(v[lo])),
                       float(math.sqrt(np.var(v[hi], ddof=1) / hi.sum() + np.var(v[lo], ddof=1) / lo.sum()))]
                   for k, v in (("beta", beta), ("beta_f", bf), ("shift", sh))}
    out["decomp"] = dec
    jr = (J[np.isfinite(J)])
    out["jumps_single_max"] = float(jr.max())
    # damping: centre C0 … C7 at α 0.01, R1 of all windows (few calm windows), median per amplitude
    r1a = []
    for i in range(8):
        d = RUNS / ("C%d_a001" % i)
        if (d / "acc.json").exists():
            r1a.append([s["R1"] for s in acc_stages(d)])
    r1a = np.array(r1a)
    out["alpha001_R1"] = {"median": np.median(r1a, 0).tolist(), "n": int(len(r1a))}
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
