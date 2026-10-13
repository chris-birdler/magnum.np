#!/usr/bin/env python
"""
Pair K1 / K1A10 (PLAN 4.8, request of Chris 2026-10-13). PLAN 4.8 fixed no analysis for the pair; the rules below
were written before the K1A10 values were read.

The pair: corner K1 (L_eff 30 l_ex, K_σ 1.5 K_eff, φ 0.62, α 0.1, seeds 2011 … 2013) at d/l_ex 212 (A 20 pJ/m at
d = 1 µm; the design runs K1_k, the rerun K1_k_mc2 where it exists, as in the design analysis) and at d/l_ex 300
(A 10 pJ/m; K1A10_k: ≤ 12 cycles, --min_calm 2, 512 samples, harmonics). Same seed, stress axes, packing and history;
the Herzer-cell grid scales with the box, so the cell-by-cell axes differ.
Values: ln w at the 7 nominal amplitudes (analyze_design.run_record: steady rule ≥ 2 calm cycles, interpolation in
ln B to the nominal B), P = w · K_d · f in W/cm³.
1. P per B for each A (geometric mean over the seeds), the ratio P(A 10) / P(A 20) paired by seed (geometric mean,
   90 % CI with t(n − 1)), the number of settled stages per B.
2. β (slope of ln P over ln B): 9 … 150 mT, 20 … 150 mT (without 9 mT), 3-point windows. Per A: per-seed fits (mean,
   SD, n; only seeds with all points). Δβ = β(A 10) − β(A 20) per seed, mean ± SE, 90 % CI with t(n − 1) (the seed is
   the unit; corrected 2026-10-13 after the review: before, the model with one offset per seed and one common slope on
   all values was the main number, its interval treats the amplitudes as independent and is too narrow; it is kept as
   a second column). Fewer than 3 seeds: no statement.
3. Prediction of the design model (Gamma GLM of analyze_design, all α 0.1 design runs) for the same change
   (coded factors of K1_k → K1A10_k, the other factors fixed): Δβ and Δ ln P, 1000 bootstrap resamples of the design
   runs (90 % percentile interval). Measured against predicted.
4. Loss channels per B and A (mean over the seeds, calm accumulator windows only): R1 = share of the cell loss at f,
   share of the loss in the 10 % of the cells with the highest loss, jump share RF / 2 of the last 3 calm cycles.
5. Harmonics (O2, PROTOKOLL 7.72; added 2026-10-13, item 1 of Chris): per B the shares of the cell loss at f, 2 f,
   3 f, 4 … 7 f and the rest (motion above 7 f or not repeating), mean ± SD over the seeds; the share of each part in
   the 10 % of the cells with the highest loss and in the rim r ≥ 0.8 R; the share of cells with Σ_k p_k > 1.001 p
   (Bessel check). A 10: K1A10_k; A 20: the C′ runs K1_k_f1 (same protocol) when they exist.

    python eval_pair.py > results/pair_K1A10.md      # also results/pair_K1A10.json (numbers for the deck)
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import deck_data as DD
import eval_alpha_fits as EAF
import eval_jumps_time as EJ

R = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = range(3)
AMPS = list(AD.AMPS_MT)
MIN_N = 3
FITS = (("9 … 150", tuple(AMPS)), ("20 … 150 (without 9 mT)", tuple(AMPS[1:])),
        ("9 / 20 / 35", (9.0, 20.0, 35.0)), ("20 / 35 / 50", (20.0, 35.0, 50.0)), ("35 / 50 / 70", (35.0, 50.0, 70.0)),
        ("50 / 70 / 100", (50.0, 70.0, 100.0)), ("70 / 100 / 150", (70.0, 100.0, 150.0)))
N_BOOT = 1000


def folder(a, k):
    d = R / ("K1_%d" % k if a == 20 else "K1A10_%d" % k)
    return (AD.rerun_dir(d) or d) if a == 20 else d


def records():
    """{A: {seed: run_record}}"""
    return {a: {k: AD.run_record(folder(a, k)) for k in SEEDS} for a in (20, 10)}


def lnP(rec, b, to_P):
    v = rec["lnw@%gmT" % b] if rec is not None else float("nan")
    return v + math.log(to_P) if np.isfinite(v) else float("nan")


def tci(v):
    v = np.asarray(v, float)
    if len(v) < MIN_N:
        return None
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.95, len(v) - 1) * se
    return v.mean(), se, v.mean() - h, v.mean() + h, len(v)


def slope_fe(rows):
    """eval_alpha_fits.fe (one offset per seed, one common slope) on rows (seed, B, y)."""
    return EAF.fe(rows)


def design_prediction(rec_pair, n_boot=N_BOOT, seed=20261013):
    """Δ outputs (A 10 − A 20) of the design GLM at the coded factors of the pair runs; (estimate, boot array)."""
    recs, _ = AD.load_runs(R, reruns=True)
    X = AD.coded(recs)
    T = AD.terms(X)
    W = np.exp(np.array([[r["lnw@%gmT" % m] for m in AMPS] for r in recs]))
    xa = AD.coded([rec_pair[20][k] for k in SEEDS])
    xb = AD.coded([rec_pair[10][k] for k in SEEDS])
    lnw_out = ["lnw_%g" % m for m in AMPS]

    def delta(c):
        oa, ob = AD.outputs_from_lnmu(AD.terms(xa) @ c.T), AD.outputs_from_lnmu(AD.terms(xb) @ c.T)
        d = {k: float(np.mean(ob[k] - oa[k])) for k in AD.OUTPUTS}
        la, lb = AD.terms(xa) @ c.T, AD.terms(xb) @ c.T
        d.update({lnw_out[i]: float(np.mean(lb[:, i] - la[:, i])) for i in range(len(AMPS))})
        sl = lambda L, idx: np.polyfit(np.log([AMPS[i] for i in idx]), L[:, idx].T, 1)[0]
        for lab, amps in FITS:
            idx = [AMPS.index(b) for b in amps]
            d["fit " + lab] = float(np.mean(sl(lb, idx) - sl(la, idx)))
        return d
    est = delta(AD.glm_coefs(T, W))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(recs), len(recs))
        try:
            boots.append(delta(AD.glm_coefs(T[i], W[i])))
        except Exception:                               # noqa: BLE001 - a non-converging resample is skipped
            pass
    return est, boots, len(recs), {"lnd_coded_A20": xa[:, AD.FACTORS.index("lnd")].tolist(),
                                   "lnd_coded_A10": xb[:, AD.FACTORS.index("lnd")].tolist()}


def nonuniformity(run, amp):
    """g = mean of |grad m̄| over the particle cells of the snapshot cut, in rad per particle diameter d: m̄ = m averaged
    over the 16 snapshots of the loop, normalised and smoothed inside each particle with sigma = L_eff / 2 (the texture
    measure of acc_reduce.py, there in 3D)."""
    import domains as DM
    z = np.load(R / run / "slices.npz")
    cfg = json.loads((R / run / "config.json").read_text())
    ids = z[amp + "_ids"].astype(np.int64)
    m0 = z[amp].astype(np.float32).mean(0)
    m0 /= np.maximum(np.linalg.norm(m0, axis=-1, keepdims=True), 1e-12)
    m0 = np.where((ids >= 0)[..., None], m0, 0.0)
    g = DM.grad_norm(DM.smooth_m(m0[:, :, None, :], ids[:, :, None], 0.5 * cfg["Leff_lex"] / cfg["dx_lex"]),
                     ids[:, :, None], cfg["dx_lex"])[:, :, 0]
    return float(g[ids >= 0].mean() * cfg.get("d_lex_used", cfg["d_lex"]))


def q_cut(run, amp):
    """Mean over the particle cells of the snapshot cut of q = Σ_j |m(t_{j+1}) − m(t_j)|² over the 16 snapshots."""
    z = np.load(R / run / "slices.npz")
    m = z[amp].astype(np.float64)
    q = ((np.roll(m, -1, axis=0) - m) ** 2).sum(axis=(0, 3))
    return float(q[z[amp + "_ids"] >= 0].mean())


def channels(a, k):
    """Per nominal B: (R1, top-10 % share, jump share) or NaN (calm windows / steady stages only)."""
    d = folder(a, k)
    out = {b: [np.nan, np.nan, np.nan] for b in AMPS}
    if (d / "acc.json").exists():
        for s in DD.acc_stages(d):
            b = min(AMPS, key=lambda x: abs(x / s["B_peak_mT"] - 1))
            if s["calm"]:
                out[b][0] = s["R1"]
                out[b][1] = s["where"]["p top 10%"]["p_share"]
    for i, m in enumerate(EJ.run_measures(d, len(AMPS))):
        if m is not None and m["steady"]:
            out[AMPS[i]][2] = m["RF"] / 2
    return out


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    rec = records()
    import json
    res = {"seeds": [folder(a, k).name for a in (20, 10) for k in SEEDS], "amps": AMPS, "ratio": {}, "beta": {},
           "lever": {}, "design": {}}
    u = json.loads((folder(20, 0) / "config.json").read_text())
    to_P = u["units"]["Kd"] * u["f_rel"] * u["units"]["f_M"] / 1e6
    L = ["# Pair K1 (A 20 pJ/m, d/l_ex 212) against K1A10 (A 10 pJ/m, d/l_ex 300)", "",
         "Corner K1: L_eff 30 l_ex (K_eff 1.0 kJ/m³ at A 20 pJ/m, d = 1 µm), K_σ 1.5 K_eff, φ 0.62, α 0.1, f 30 MHz, "
         "seeds 2011 … 2013 (K1_0 … K1_2). Steady rule ≥ 2 calm cycles. Paired by seed; 90 %% CI with t(n − 1); "
         "fewer than %d seeds: no statement. The cell-by-cell axes of the pair differ (the Herzer cells scale with "
         "the box)." % MIN_N, "",
         "Runs used: " + ", ".join(folder(a, k).name for a in (20, 10) for k in SEEDS), ""]
    # 1. P per B
    L += ["## Loss density per B", "",
          "| B (mT) | P at A 20 (W/cm³) | P at A 10 (W/cm³) | P(A 10) / P(A 20), paired (90 % CI) | settled A 20 / A 10 |",
          "|---|---|---|---|---|"]
    for b in AMPS:
        pa = [lnP(rec[20][k], b, to_P) for k in SEEDS]
        pb = [lnP(rec[10][k], b, to_P) for k in SEEDS]
        d = [y - x for x, y in zip(pa, pb) if np.isfinite(x) and np.isfinite(y)]
        s = tci(d)
        res["ratio"]["%g" % b] = ({"gm": math.exp(s[0]), "lo": math.exp(s[2]), "hi": math.exp(s[3]), "n": s[4]} if s
                                else {"n": len(d)})
        gm = lambda v: ("%.3g" % math.exp(np.nanmean(v))) if np.isfinite(v).any() else "–"
        L.append("| %g | %s | %s | %s | %d / %d |" % (
            b, gm(np.array(pa)), gm(np.array(pb)),
            ("× %.2f (%.2f … %.2f), n %d" % (math.exp(s[0]), math.exp(s[2]), math.exp(s[3]), s[4])) if s
            else "n %d: no statement" % len(d), int(np.isfinite(pa).sum()), int(np.isfinite(pb).sum())))
    # 2. β fits
    L += ["", "## β fits", "",
          "| fit (mT) | β at A 20, per seed (mean / SD, n) | β at A 10, per seed | Δβ = A 10 − A 20 per seed (90 % CI, t) | "
          "the same, model with seed offsets on all values |", "|---|---|---|---|---|"]
    lnw = {a: {k: {b: (rec[a][k]["lnw@%gmT" % b] if rec[a][k] is not None else np.nan) for b in AMPS} for k in SEEDS}
           for a in (20, 10)}
    meas = {}
    for lab, amps in FITS:
        cells = []
        for a in (20, 10):
            v = [np.polyfit(np.log(amps), [lnw[a][k][b] for b in amps], 1)[0] for k in SEEDS
                 if all(np.isfinite(lnw[a][k][b]) for b in amps)]
            cells.append("%.2f / %.2f (n %d)" % (np.mean(v), np.std(v, ddof=1), len(v)) if len(v) >= 2 else
                         "n %d" % len(v))
        rows = [(k, b, lnw[10][k][b] - lnw[20][k][b]) for k in SEEDS for b in amps
                if np.isfinite(lnw[10][k][b]) and np.isfinite(lnw[20][k][b])]
        r = slope_fe(rows) if len({x[0] for x in rows}) >= MIN_N else None
        meas[lab] = r
        res["beta"][lab] = {str(a): [float(x) for x in
                                     [np.polyfit(np.log(amps), [lnw[a][k][b] for b in amps], 1)[0] for k in SEEDS
                                      if all(np.isfinite(lnw[a][k][b]) for b in amps)]] for a in (20, 10)}
        dd = [np.polyfit(np.log(amps), [lnw[10][k][b] for b in amps], 1)[0] -
              np.polyfit(np.log(amps), [lnw[20][k][b] for b in amps], 1)[0] for k in SEEDS
              if all(np.isfinite(lnw[a][k][b]) for a in (20, 10) for b in amps)]
        sl = tci(dd)
        meas[lab] = sl
        if sl:
            res["lever"][lab] = {"mean": -sl[0], "se": sl[1], "lo": -sl[3], "hi": -sl[2], "n": sl[4],
                                 "per_seed": [-x for x in dd]}
        if r:
            h = stats.t.ppf(0.95, r[4]) * r[1]
            res.setdefault("lever_fe", {})[lab] = {"mean": -r[0], "se": r[1], "lo": -r[0] - h, "hi": -r[0] + h}
        L.append("| %s | %s | %s | %s | %s |" % (
            lab, cells[0], cells[1], ("%+.2f ± %.2f (%+.2f … %+.2f, n %d)" % sl) if sl else "no statement",
            EAF.fmt_fe(r) if r else "no statement"))
    # 3. design prediction
    est, boots, n_des, lnd = design_prediction(rec)
    L += ["", "## Measured against the prediction of the design model", "",
          "Design model: Gamma GLM (quadratic in the 4 coded factors) on %d α 0.1 design runs, evaluated at the coded "
          "factors of the pair runs (coded ln d/l_ex: A 20 %s, A 10 %s); %d bootstrap resamples (90 %% interval)." %
          (n_des, ", ".join("%.2f" % x for x in lnd["lnd_coded_A20"]), ", ".join("%.2f" % x for x in lnd["lnd_coded_A10"]),
           len(boots)), "",
          "| output | measured (paired, 90 % CI) | design model (90 % interval) |", "|---|---|---|"]
    for lab, _ in FITS:
        key = "fit " + lab
        b = np.array([x[key] for x in boots])
        m = meas[lab]
        res["design"][lab] = {"mean": -est[key], "lo": -float(np.percentile(b, 95)), "hi": -float(np.percentile(b, 5))}
        L.append("| Δβ %s mT | %s | %+.2f (%+.2f … %+.2f) |" % (
            lab, ("%+.2f (%+.2f … %+.2f)" % (m[0], m[2], m[3])) if m else "–",
            est[key], *np.percentile(b, [5, 95])))
    for b_ in AMPS:
        key = "lnw_%g" % b_
        bb = np.array([x[key] for x in boots])
        d = [lnw[10][k][b_] - lnw[20][k][b_] for k in SEEDS if np.isfinite(lnw[10][k][b_]) and np.isfinite(lnw[20][k][b_])]
        s = tci(d)
        res.setdefault("design_ratio", {})["%g" % b_] = {"est": math.exp(est[key]),
                                                         "lo": float(np.exp(np.percentile(bb, 5))),
                                                         "hi": float(np.exp(np.percentile(bb, 95)))}
        L.append("| P(A 10) / P(A 20) at %g mT | %s | × %.2f (%.2f … %.2f) |" % (
            b_, ("× %.2f (%.2f … %.2f)" % (math.exp(s[0]), math.exp(s[2]), math.exp(s[3]))) if s else "–",
            math.exp(est[key]), *np.exp(np.percentile(bb, [5, 95]))))
    res["P"] = {str(a): [[math.exp(lnP(rec[a][k], b, to_P)) for b in AMPS] for k in SEEDS] for a in (20, 10)}  # W/cm³
    # 4. channels
    ch = {a: {k: channels(a, k) for k in SEEDS} for a in (20, 10)}
    res["jump"] = {str(a): {"%g" % b: [ch[a][k][b][2] for k in SEEDS] for b in AMPS} for a in (20, 10)}
    L += ["", "## Loss channels and location (mean over the seeds with a value, n in brackets)", "",
          "| B (mT) | R1 A 20 | R1 A 10 | top-10 % share A 20 | top-10 % share A 10 | jump share A 20 | jump share A 10 |",
          "|---|---|---|---|---|---|---|"]
    for b in AMPS:
        cells = []
        for j in range(3):
            for a in (20, 10):
                v = [ch[a][k][b][j] for k in SEEDS if np.isfinite(ch[a][k][b][j])]
                cells.append("%.2f (%d)" % (np.mean(v), len(v)) if v else "–")
        L.append("| %g | %s |" % (b, " | ".join(cells)))
    L += ["", "R1: share of the cell loss at f (calm accumulator windows). Top-10 % share: share of the loss in the 10 % "
          "of the particle cells with the highest loss. Jump share: RF / 2 (energy above the smooth periodic curve, "
          "last 3 calm cycles)."]
    # 5. harmonics
    L += ["", "## Harmonics (O2): shares of the cell loss, mean ± SD over the seeds", ""]
    for a, pat in ((10, "K1A10_%d"), (20, "K1_%d_f1")):
        st = {b: [] for b in AMPS}
        for k in SEEDS:
            f = R / (pat % k) / "acc.json"
            if not (R / (pat % k) / "DONE").exists() or not f.exists():
                continue
            for s in json.loads(f.read_text())["stages"]:
                if "harm" not in s:
                    continue
                b = min(AMPS, key=lambda x: abs(x / s["B_peak_mT"] - 1))
                h, wt, rim = s["harm"], s["where"]["p top 10%"], s["where"]["r>=0.8"]
                st[b].append([h["R_k"][0], h["R_k"][1], h["R_k"][2], sum(h["R_k"][3:]), h["R_rest"],
                              wt["p_fund_share"], wt["p_harm_share"], wt["p_rest_share"],
                              rim["p_fund_share"], rim["p_harm_share"], rim["p_rest_share"], h["bessel_frac_over_1e-3"]])
        if not any(st.values()):
            L += ["A %d pJ/m: no runs with harmonics yet (%s)." % (a, pat % 0), ""]
            continue
        res.setdefault("harm", {})[str(a)] = {"%g" % b: np.array(st[b], float)[:, :5].tolist() for b in AMPS if st[b]}
        hk = {b: [] for b in AMPS}
        for k in SEEDS:
            f = R / (pat % k) / "acc.json"
            if (R / (pat % k) / "DONE").exists() and f.exists():
                for s in json.loads(f.read_text())["stages"]:
                    if "harm" in s:
                        hk[min(AMPS, key=lambda x: abs(x / s["B_peak_mT"] - 1))].append(s["harm"]["R_k"] + [s["harm"]["R_rest"]])
        res.setdefault("harm_k", {})[str(a)] = {"%g" % b: hk[b] for b in AMPS if hk[b]}
        L += ["A %d pJ/m (%s):" % (a, pat.replace("%d", "k")), "",
              "| B (mT) | n | f | 2 f | 3 f | 4 … 7 f | rest | top 10 %: f / harmonics / rest | rim r ≥ 0.8 R: f / "
              "harmonics / rest | cells with Bessel excess |", "|---|---|---|---|---|---|---|---|---|---|"]
        for b in AMPS:
            if not st[b]:
                continue
            x = np.array(st[b], float)
            m, sd = x.mean(0), (x.std(0, ddof=1) if len(x) > 1 else np.full(x.shape[1], np.nan))
            L.append("| %g | %d | %.2f ± %.2f | %.3f ± %.3f | %.3f ± %.3f | %.3f ± %.3f | %.2f ± %.2f | %.2f / %.2f / %.2f | "
                     "%.2f / %.2f / %.2f | %.0f %% |" % (b, len(x), m[0], sd[0], m[1], sd[1], m[2], sd[2], m[3], sd[3],
                                                     m[4], sd[4], m[5], m[6], m[7], m[8], m[9], m[10], 100 * m[11]))
        L.append("")
    L += ["Rest: motion above 7 f or motion that does not repeat from cycle to cycle. A zone share of the rest can exceed "
          "1 where the rest is small, because cells with a Bessel excess have a negative rest. The 4 … 7 f shares are "
          "flat in k (≈ 0.02 each): a step-like motion in each cycle and a slow drift within a cycle both give a flat "
          "spectrum; these data cannot separate them. The rim r ≥ 0.8 R holds ≈ 49 % of the volume."]
    # 6. motion in the snapshots and loss at f / outside f (added 2026-10-13 after the review of deck slide 16)
    L += ["", "## Motion in the snapshot cut and loss at f / outside f, ratio A 10 / A 20 per seed", "",
          "q = Σ over the 16 snapshots of one loop of |Δm|² per mesh cell, mean over the particle cells of the cut (the "
          "loss measure q of the deck; it resolves motion on the time scale T/16). Loss at f = R1 · P, outside f = "
          "(1 − R1) · P (calm accumulator windows; n < 3 where a window has an event).", "",
          "| B (mT) | q ratio per seed | q ratio of the means | P ratio (paired, geometric mean) | loss at f, ratio | loss outside f, ratio |",
          "|---|---|---|---|---|---|"]
    res["q"], res["f_split"] = {}, {}
    for b, amp in ((9, "amp00"), (50, "amp03"), (150, "amp06")):
        q = {a: [q_cut(folder(a, k).name if a == 10 else "K1_%d" % k, amp) for k in SEEDS] for a in (20, 10)}
        qr = [q[10][k] / q[20][k] for k in SEEDS]
        fs = []
        for k in SEEDS:
            r20, r10 = ch[20][k][b][0], ch[10][k][b][0]
            if np.isfinite(r20) and np.isfinite(r10):
                p20, p10 = math.exp(lnP(rec[20][k], b, to_P)), math.exp(lnP(rec[10][k], b, to_P))
                fs.append((r10 * p10 / (r20 * p20), (1 - r10) * p10 / ((1 - r20) * p20)))
        gmf = lambda v: math.exp(np.mean(np.log(v)))
        res["q"]["%g" % b] = {"per_seed": qr, "ratio_of_means": float(np.mean(q[10]) / np.mean(q[20]))}
        res["f_split"]["%g" % b] = {"f": gmf([x[0] for x in fs]) if fs else None,
                                    "out": gmf([x[1] for x in fs]) if fs else None, "n": len(fs)}
        L.append("| %d | %s | %.2f | %.2f | %s | %s |" % (
            b, " / ".join("%.2f" % x for x in qr), np.mean(q[10]) / np.mean(q[20]), res["ratio"]["%g" % b]["gm"],
            ("%.2f (n %d)" % (gmf([x[0] for x in fs]), len(fs))) if fs else "–",
            ("%.2f (n %d)" % (gmf([x[1] for x in fs]), len(fs))) if fs else "–"))
    print("\n".join(L))
    import eval_mechanism as EM
    res["g"] = {str(a): [[nonuniformity(folder(a, k).name if a == 10 else "K1_%d" % k, amp) for amp in
                          ("amp00", "amp03", "amp06")] for k in SEEDS] for a in (20, 10)}
    res["vsw"] = {str(a): [[EM.measures(folder(a, k).name if a == 10 else "K1_%d" % k)["V_sw %d" % b] for b in (9, 50, 150)]
                           for k in SEEDS] for a in (20, 10)}      # switched volume per half loop (% of V), 3D snapshots
    res["note"] = ("lever = beta(A 20) - beta(A 10), paired by seed; design = the same from the design GLM (90 % "
                   "bootstrap); vsw from the snapshot runs K1_k (the reruns K1_k_mc2 have no snapshots) and K1A10_k")
    (pathlib.Path(__file__).resolve().parent / "results" / "pair_K1A10.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
