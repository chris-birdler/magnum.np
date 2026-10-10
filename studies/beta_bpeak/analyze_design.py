#!/usr/bin/env python
"""
Analysis of the Sobol design v2 (PLAN 4.0, rules fixed before the start, PROTOKOLL 7.33).

Estimand: the mean loss E[w] over the realisations (a core averages ~1e9 particles); beta_core = d ln E[w] / d ln B.
  main fit   per nominal amplitude (9 ... 150 mT): Gamma GLM with log link of w on the 15 terms of the coded factors
             x = (ln Q_eff, r_p, phi, ln d/l_ex) (linear, quadratic, 2-factor interactions); ln E[w] at 10 / 50 / 100 mT
             and beta_core of the windows and over all amplitudes from the predicted means; SE and CI by a bootstrap
             over the runs
  fit A      OLS of the per-run ln w and beta with HC3 SE (estimates E[ln w], the typical run)
  fit B      Huber fit of the same (the bulk without rare states)
  no flags   main fit without the flagged runs: |w_dis/w_loop - 1| > 0.20 at any amplitude (state change in the
             measured cycles; PROTOKOLL 7.39)
  w_dis      main fit with w_dis (LLG dissipation) instead of w_loop
Range effect of a factor = mean over the design points of (output with the factor at its maximum) - (at its minimum).
Decision rules: primary outputs beta_all and lnw@50mT x 4 factors = 8 tests, Holm 5 %; class "lever" (significant after
Holm), "not a lever" (90 % CI inside +-0.15), else "not determined"; only hints for the other outputs.
Scatter model: ln(residual^2) of fit A on the 4 linear factors for the 2 primary outputs, Holm over the 4 slopes; if a
slope is significant, fit A and the main fit are repeated with the weights 1 / sigma^2(x).
alpha 0.01 test: pairs C#_a001 / C#: ln of the w ratio, and R = |m1|^2 ratio = 10 (R1'/R1) (w_dis'/w_dis) at 9 and 50 mT.
Mesh-check trigger: in the main fit |range effect| >= 0.15 on any beta or ln E[w] output with the 95 % CI excluding 0.

    python analyze_design.py [--runs runs] [--boot 1000] [--out results/sobol_analysis]
"""
import argparse
import json
import math
import pathlib
import re
import warnings

import numpy as np
from scipy import stats

import analyze

HERE = pathlib.Path(__file__).resolve().parent
JS_T = 1.5
AMPS_MT = (9.0, 20.0, 35.0, 50.0, 70.0, 100.0, 150.0)
WINDOWS = ((9.0, 20.0, 35.0), (35.0, 50.0, 70.0), (70.0, 100.0, 150.0))
LNW_MT = (10.0, 50.0, 100.0)
FACTORS = ("lnQ", "r_p", "phi", "lnd")
# coded range -1 ... +1 of each factor (design ranges, PLAN 4.0)
RANGES = {"lnQ": (math.log(1 / 30.0 ** 2), math.log(1 / 12.0 ** 2)), "r_p": (0.0, 3.0), "phi": (0.55, 0.69),
          "lnd": (math.log(212.0), math.log(300.0))}
OUTPUTS = ["beta@9-35mT", "beta@35-70mT", "beta@70-150mT", "beta_all", "lnw@10mT", "lnw@50mT", "lnw@100mT"]
PRIMARY = ("beta_all", "lnw@50mT")
MARGIN = 0.15
PATTERN = re.compile(r"^(S\d{3}|C\d|K\d_\d)$")
BAL_EVENT = 0.20      # a run is flagged if |w_dis/w_loop - 1| > BAL_EVENT at any amplitude (state change in the
#                       measured cycles; rule fixed 2026-10-07 from the earlier data, PROTOKOLL 7.39)


# --- data ---------------------------------------------------------------------------------------------------------

def run_record(d):
    """One run folder -> dict with the factors, ln w at the nominal amplitudes and the flags (None if not usable)."""
    d = pathlib.Path(d)
    if not (d / "DONE").exists():
        return None
    cfg, rows = analyze.load_run(d)
    if len(rows) < len(AMPS_MT):
        return None
    bm = np.array([1000.0 * JS_T * r["b_peak"] for r in rows])
    lw = np.array([math.log(r["w_loop"]) if r["w_loop"] > 0 else np.nan for r in rows])
    if not np.all(np.isfinite(lw)):
        return None
    # ln w at the nominal amplitudes: linear in ln B between the measured amplitudes (b is within ~4 % of the target)
    lnw = np.interp(np.log(AMPS_MT), np.log(bm), lw)
    lnw[0] = lw[0] + (lw[1] - lw[0]) / math.log(bm[1] / bm[0]) * math.log(AMPS_MT[0] / bm[0])
    lnw[-1] = lw[-1] + (lw[-1] - lw[-2]) / math.log(bm[-1] / bm[-2]) * math.log(AMPS_MT[-1] / bm[-1])
    Leff = cfg["Leff_lex"]
    rec = {"name": d.name, "alpha": cfg["alpha"], "seed": cfg["seed"],
           "lnQ": math.log(cfg.get("Q_eff") or 1.0 / Leff ** 2), "r_p": cfg["r_p"], "phi": cfg.get("phi_vox", cfg["phi"]),
           "lnd": math.log(cfg.get("d_lex_used", cfg["d_lex"])), "d_over_L": cfg.get("d_lex_used", cfg["d_lex"]) / Leff,
           "flagged": any(abs(r["balance"] - 1.0) > BAL_EVENT for r in rows),
           "flags": [analyze.flags(r) for r in rows],
           "w_dis": [r["w_dis"] for r in rows], "b_mT": bm.tolist()}
    for k, mT in enumerate(AMPS_MT):
        rec["lnw@%gmT" % mT] = float(lnw[k])
    lwd = np.log([r["w_dis"] for r in rows])                 # the same at the nominal amplitudes from w_dis
    lnwd = np.interp(np.log(AMPS_MT), np.log(bm), lwd)
    lnwd[0] = lwd[0] + (lwd[1] - lwd[0]) / math.log(bm[1] / bm[0]) * math.log(AMPS_MT[0] / bm[0])
    lnwd[-1] = lwd[-1] + (lwd[-1] - lwd[-2]) / math.log(bm[-1] / bm[-2]) * math.log(AMPS_MT[-1] / bm[-1])
    for k, mT in enumerate(AMPS_MT):
        rec["lnwdis@%gmT" % mT] = float(lnwd[k])
    rec.update(per_run_outputs(lnw))
    return rec


def per_run_outputs(lnw):
    """beta of the windows and over all amplitudes, ln w at 10 / 50 / 100 mT from ln w at the nominal amplitudes."""
    lb = np.log(AMPS_MT)
    out = {}
    for win in WINDOWS:
        idx = [AMPS_MT.index(m) for m in win]
        out["beta@%d-%dmT" % (win[0], win[-1])] = float(np.polyfit(lb[idx], lnw[idx], 1)[0])
    out["beta_all"] = float(np.polyfit(lb, lnw, 1)[0])
    for mT in LNW_MT:
        out["lnw@%gmT" % mT] = float(np.interp(math.log(mT), lb, lnw))
    return out


def load_runs(runs, alpha=0.1):
    recs, skipped = [], []
    for d in sorted(pathlib.Path(runs).iterdir()):
        if not d.is_dir() or not PATTERN.match(d.name):
            continue
        try:
            r = run_record(d)
        except Exception as e:                      # noqa: BLE001 - report, do not stop
            skipped.append("%s (%s)" % (d.name, e))
            continue
        if r is None:
            skipped.append("%s (%s)" % (d.name, "not DONE" if not (d / "DONE").exists() else
                                        "incomplete or w ≤ 0 at an amplitude (e.g. a negative mean loop area)"))
        elif abs(r["alpha"] - alpha) < 1e-12:
            recs.append(r)
    return recs, skipped


# --- design matrix ------------------------------------------------------------------------------------------------

def coded(recs):
    return np.array([[2.0 * (r[f] - RANGES[f][0]) / (RANGES[f][1] - RANGES[f][0]) - 1.0 for f in FACTORS] for r in recs])


def terms(X):
    """15 columns: 1, x_i, x_i^2, x_i x_j."""
    n, k = X.shape
    cols = [np.ones(n)] + [X[:, i] for i in range(k)] + [X[:, i] ** 2 for i in range(k)]
    cols += [X[:, i] * X[:, j] for i in range(k) for j in range(i + 1, k)]
    return np.column_stack(cols)


def contrast(X, j):
    """Range-effect contrast for factor j: mean over the design points of terms(x_j = +1) - terms(x_j = -1)."""
    hi, lo = X.copy(), X.copy()
    hi[:, j], lo[:, j] = 1.0, -1.0
    return terms(hi), terms(lo)


# --- main fit (Gamma GLM, log link) ----------------------------------------------------------------------------------

def glm_coefs(T, W, weights=None):
    """Gamma GLM with log link per amplitude: W [n, 7] (w > 0) -> coefficients [7, p]."""
    import statsmodels.api as sm
    fam = sm.families.Gamma(sm.families.links.Log())
    out = []
    for k in range(W.shape[1]):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = sm.GLM(W[:, k], T, family=fam, var_weights=weights).fit()
        out.append(res.params)
    return np.array(out)


def _slope_weights(idx):
    x = np.log(np.array(AMPS_MT))[list(idx)]
    return (x - x.mean()) / ((x - x.mean()) ** 2).sum()


_OUT_W = {}                                         # output -> weight vector over the 7 amplitudes (linear in ln w)
for _win in WINDOWS:
    _w = np.zeros(len(AMPS_MT))
    _w[[AMPS_MT.index(m) for m in _win]] = _slope_weights([AMPS_MT.index(m) for m in _win])
    _OUT_W["beta@%d-%dmT" % (_win[0], _win[-1])] = _w
_OUT_W["beta_all"] = _slope_weights(range(len(AMPS_MT)))
for _mT in LNW_MT:
    _w = np.zeros(len(AMPS_MT))
    _lb = np.log(AMPS_MT)
    _j = int(np.clip(np.searchsorted(_lb, math.log(_mT)), 1, len(AMPS_MT) - 1))
    _t = (math.log(_mT) - _lb[_j - 1]) / (_lb[_j] - _lb[_j - 1])
    _w[_j - 1], _w[_j] = 1.0 - _t, _t
    _OUT_W["lnw@%gmT" % _mT] = _w


def outputs_from_lnmu(lnmu):
    """lnmu [n, 7] (ln E[w] at the nominal amplitudes) -> dict output -> [n] (all outputs are linear in ln w:
    least-squares slopes and linear interpolation, same as per_run_outputs)."""
    return {k: lnmu @ _OUT_W[k] for k in OUTPUTS}


def range_effects_glm(coefs, X):
    eff = {}
    for j, f in enumerate(FACTORS):
        Th, Tl = contrast(X, j)
        oh, ol = outputs_from_lnmu(Th @ coefs.T), outputs_from_lnmu(Tl @ coefs.T)
        for o in OUTPUTS:
            eff[(o, f)] = float(np.mean(oh[o] - ol[o]))
    return eff


def main_fit(recs, n_boot=1000, seed=20261007, weights=None, key="lnw"):
    X = coded(recs)
    T = terms(X)
    W = np.exp(np.array([[r["%s@%gmT" % (key, m)] for m in AMPS_MT] for r in recs]))
    est = range_effects_glm(glm_coefs(T, W, weights), X)
    rng = np.random.default_rng(seed)
    boots, n_fail = [], 0
    n = len(recs)
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        try:
            c = glm_coefs(T[i], W[i], None if weights is None else weights[i])
            boots.append(range_effects_glm(c, X))
        except Exception:                           # noqa: BLE001 - a non-converging resample is counted
            n_fail += 1
    res = {}
    for key, v in est.items():
        b = np.array([bb[key] for bb in boots])
        se = float(b.std(ddof=1))
        res[key] = {"eff": v, "se": se, "ci95": [float(x) for x in np.percentile(b, [2.5, 97.5])],
                    "ci90": [float(x) for x in np.percentile(b, [5.0, 95.0])],
                    "p": float(2 * stats.norm.sf(abs(v) / se)) if se > 0 else float("nan")}
    return res, {"n_boot": len(boots), "n_fail": n_fail}


# --- fits A (OLS, HC3) and B (Huber) on the per-run outputs -------------------------------------------------------

def linear_fit(recs, kind="ols", weights=None):
    import statsmodels.api as sm
    X = coded(recs)
    T = terms(X)
    res = {}
    for o in OUTPUTS:
        y = np.array([r[o] for r in recs])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if kind == "ols":
                m = (sm.WLS(y, T, weights=weights) if weights is not None else sm.OLS(y, T)).fit(cov_type="HC3")
            else:
                m = sm.RLM(y, T, M=sm.robust.norms.HuberT()).fit()
        cov = np.asarray(m.cov_params())
        for j, f in enumerate(FACTORS):
            Th, Tl = contrast(X, j)
            c = (Th - Tl).mean(0)
            v = float(c @ m.params)
            se = float(math.sqrt(c @ cov @ c))
            res[(o, f)] = {"eff": v, "se": se, "ci95": [v - 1.96 * se, v + 1.96 * se],
                           "p": float(2 * stats.norm.sf(abs(v) / se)) if se > 0 else float("nan")}
        res[(o, "resid")] = (y - m.fittedvalues).tolist()
    return res


def scatter_model(recs, fitA):
    """ln(residual^2) of fit A on the 4 linear coded factors, for the primary outputs; Holm over the 4 slopes."""
    import statsmodels.api as sm
    X = coded(recs)
    D = np.column_stack([np.ones(len(recs)), X])
    out = {}
    for o in PRIMARY:
        r = np.array(fitA[(o, "resid")])
        z = np.log(r ** 2 + 1e-12)
        m = sm.OLS(z, D).fit()
        p = m.pvalues[1:]
        sig = holm(list(p))
        w = 1.0 / np.exp(m.fittedvalues)
        out[o] = {"slopes": dict(zip(FACTORS, m.params[1:].tolist())), "p": dict(zip(FACTORS, p.tolist())),
                  "significant": dict(zip(FACTORS, sig)), "any": bool(any(sig)), "weights": (w / w.mean()).tolist()}
    return out


# --- decision rules -----------------------------------------------------------------------------------------------

def holm(pvals, level=0.05):
    """Holm-Bonferroni step-down: list of booleans (rejected)."""
    order = np.argsort(pvals)
    m = len(pvals)
    rej = [False] * m
    for rank, i in enumerate(order):
        if pvals[i] <= level / (m - rank):
            rej[i] = True
        else:
            break
    return rej


def classify(main):
    keys = [(o, f) for o in PRIMARY for f in FACTORS]
    rej = holm([main[k]["p"] for k in keys])
    out = {}
    for k, r in zip(keys, rej):
        lo, hi = main[k]["ci90"]
        cls = "lever" if r else ("not a lever" if (-MARGIN < lo and hi < MARGIN) else "not determined")
        out[k] = cls
    return out


def mesh_trigger(main):
    hits = [(o, f) for (o, f), v in main.items() if abs(v["eff"]) >= MARGIN and (v["ci95"][0] > 0 or v["ci95"][1] < 0)]
    corners = sorted({"K2" if f in ("lnQ", "r_p") else "K3" for _, f in hits})
    return hits, corners


# --- alpha 0.01 test ----------------------------------------------------------------------------------------------

def alpha_test(runs):
    """Pairs C#_a001 / C#: ln w ratio per amplitude and ln R at 9 and 50 mT (R = 10 (R1'/R1)(w_dis'/w_dis))."""
    runs = pathlib.Path(runs)
    pairs = []
    for d in sorted(runs.glob("C?_a001")):
        p = runs / d.name[:2]
        if not ((d / "DONE").exists() and (p / "DONE").exists()):
            continue
        a, b = run_record(d), run_record(p)
        if a is None or b is None:
            continue
        rec = {"pair": d.name[:2], "lnw_ratio": [a["lnw@%gmT" % m] - b["lnw@%gmT" % m] for m in AMPS_MT]}
        try:
            A, B = (json.loads((x / "acc.json").read_text())["stages"] for x in (d, p))
            for mT in (9.0, 50.0):
                ia = int(np.argmin([abs(s["B_peak_mT"] / mT - 1) for s in A]))
                ib = int(np.argmin([abs(s["B_peak_mT"] / mT - 1) for s in B]))
                ka = int(np.argmin([abs(x / mT - 1) for x in a["b_mT"]]))
                kb = int(np.argmin([abs(x / mT - 1) for x in b["b_mT"]]))
                rec["lnR@%gmT" % mT] = math.log(10.0 * A[ia]["R1"] / B[ib]["R1"] * a["w_dis"][ka] / b["w_dis"][kb])
        except (OSError, KeyError, ValueError):
            pass
        pairs.append(rec)
    out = {"n": len(pairs), "pairs": pairs}
    if len(pairs) >= 3:
        L = np.array([p["lnw_ratio"] for p in pairs])
        out["lnw_ratio_mean"] = L.mean(0).tolist()
        out["lnw_ratio_se"] = (L.std(0, ddof=1) / math.sqrt(len(pairs))).tolist()
        for mT in (9.0, 50.0):
            v = np.array([p["lnR@%gmT" % mT] for p in pairs if "lnR@%gmT" % mT in p])
            if len(v) >= 3:
                m, se = float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v)))
                t = stats.t.ppf(0.975, len(v) - 1)
                ci = (m - t * se, m + t * se)
                hyp = [h for h, val in (("viscous R~1", 1.0), ("alpha-independent R~10", 10.0), ("resonance R~100", 100.0))
                       if ci[0] <= math.log(val) <= ci[1]]
                out["lnR@%gmT" % mT] = {"mean": m, "se": se, "ci95": list(ci), "n": len(v),
                                        "class": hyp[0] if len(hyp) == 1 else ("mixed" if not hyp else "ambiguous: " + ", ".join(hyp))}
    return out


# --- D3 check -----------------------------------------------------------------------------------------------------

def d3_check(recs, fitA):
    small = np.array([r["d_over_L"] < 10.0 for r in recs])
    out = {}
    for o in PRIMARY:
        r = np.array(fitA[(o, "resid")])
        if small.sum() >= 2 and (~small).sum() >= 2:
            t = stats.ttest_ind(r[small], r[~small], equal_var=False)
            out[o] = {"n_small": int(small.sum()), "diff": float(r[small].mean() - r[~small].mean()), "p": float(t.pvalue)}
    return out


# --- report -------------------------------------------------------------------------------------------------------

def analyse(recs, n_boot=1000, runs=None):
    fitA = linear_fit(recs, "ols")
    sc = scatter_model(recs, fitA)
    weights = np.array(sc["lnw@50mT"]["weights"]) if sc["lnw@50mT"]["any"] else None
    main, info = main_fit(recs, n_boot, weights=weights)
    if sc["beta_all"]["any"]:
        fitA = linear_fit(recs, "ols", weights=np.array(sc["beta_all"]["weights"]))
    fitB = linear_fit(recs, "huber")
    unflagged = [r for r in recs if not r["flagged"]]
    noflag = main_fit(unflagged, max(200, n_boot // 5))[0] if len(unflagged) >= 30 else None
    wdis = main_fit(recs, max(200, n_boot // 5), key="lnwdis")[0] if all("lnwdis@9mT" in r for r in recs) else None
    res = {"n_runs": len(recs), "n_flagged": len(recs) - len(unflagged), "boot": info, "main": main, "A": fitA,
           "B": fitB, "noflag": noflag, "wdis": wdis, "scatter": sc, "classes": classify(main), "mesh_trigger": mesh_trigger(main),
           "d3": d3_check(recs, fitA), "phi_dilution": math.log(RANGES["phi"][1] / RANGES["phi"][0])}
    if runs is not None:
        res["alpha"] = alpha_test(runs)
    return res


def report(res):
    L = ["# Sobol design v2: analysis (PLAN 4.0)", "",
         "Runs in the fits: %d (flagged %d); bootstrap %d resamples (%d failed). Range effect = output at factor max "
         "minus at min, other factors averaged over the design points. φ: the effect on ln E[w] contains the dilution "
         "%.3f (value without it in brackets)." % (res["n_runs"], res["n_flagged"], res["boot"]["n_boot"],
                                                   res["boot"]["n_fail"], res["phi_dilution"]), ""]
    L += ["## Primary outputs (8 tests, Holm 5 %)", "", "| output | factor | main fit (95 % CI) | 90 % CI | class | fit A | fit B |",
          "|---|---|---|---|---|---|---|"]
    for o in PRIMARY:
        for f in FACTORS:
            m, a, b = res["main"][(o, f)], res["A"][(o, f)], res["B"][(o, f)]
            dil = " [%+.3f]" % (m["eff"] - res["phi_dilution"]) if (f == "phi" and o.startswith("lnw")) else ""
            L.append("| %s | %s | %+.3f%s (%+.3f … %+.3f) | %+.3f … %+.3f | %s | %+.3f ± %.3f | %+.3f ± %.3f |" %
                     (o, f, m["eff"], dil, *m["ci95"], *m["ci90"], res["classes"][(o, f)], a["eff"], a["se"], b["eff"], b["se"]))
    L += ["", "## Secondary outputs (hints only)", "", "| output | factor | main fit (95 % CI) | fit A | fit B | without flagged runs | with w_dis |",
          "|---|---|---|---|---|---|---|"]
    for o in OUTPUTS:
        if o in PRIMARY:
            continue
        for f in FACTORS:
            m, a, b = res["main"][(o, f)], res["A"][(o, f)], res["B"][(o, f)]
            nf = res["noflag"][(o, f)]["eff"] if res["noflag"] else float("nan")
            wd = res["wdis"][(o, f)]["eff"] if res.get("wdis") else float("nan")
            L.append("| %s | %s | %+.3f (%+.3f … %+.3f) | %+.3f ± %.3f | %+.3f ± %.3f | %+.3f | %+.3f |" %
                     (o, f, m["eff"], *m["ci95"], a["eff"], a["se"], b["eff"], b["se"], nf, wd))
    hits, corners = res["mesh_trigger"]
    L += ["", "## Scatter model (fit A residuals, Holm over 4 slopes)", ""]
    for o, s in res["scatter"].items():
        L.append("- %s: slopes %s; p %s; significant: %s" % (
            o, ", ".join("%s %+.2f" % kv for kv in s["slopes"].items()), ", ".join("%s %.3f" % kv for kv in s["p"].items()),
            ", ".join(f for f, v in s["significant"].items() if v) or "none"))
    L += ["", "## Mesh-check trigger", "", "- hits: %s" % (", ".join("%s × %s" % h for h in hits) or "none"),
          "- corners to check: %s" % (", ".join(corners) or "none"), ""]
    L += ["## D3 check (fit A residuals, d/L_eff < 10 vs ≥ 10)", ""]
    for o, v in res["d3"].items():
        L.append("- %s: n(<10) = %d, Δ = %+.3f, p = %.3f" % (o, v["n_small"], v["diff"], v["p"]))
    if "alpha" in res:
        a = res["alpha"]
        L += ["", "## α 0.01 centre test (pairs: %d)" % a["n"], ""]
        if "lnw_ratio_mean" in a:
            L.append("- ln(w α0.01 / w α0.1) per amplitude: " + ", ".join(
                "%g mT %+.2f ± %.2f" % (m, x, s) for m, x, s in zip(AMPS_MT, a["lnw_ratio_mean"], a["lnw_ratio_se"])))
        for mT in (9.0, 50.0):
            k = "lnR@%gmT" % mT
            if k in a:
                v = a[k]
                L.append("- %s: R = %.1f (95 %% CI %.1f … %.1f, n = %d): %s" % (k, math.exp(v["mean"]), math.exp(v["ci95"][0]),
                                                                           math.exp(v["ci95"][1]), v["n"], v["class"]))
    return "\n".join(L) + "\n"


def to_json(res):
    def conv(x):
        if isinstance(x, dict):
            return {(" × ".join(k) if isinstance(k, tuple) else k): conv(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [conv(v) for v in x]
        if isinstance(x, (np.floating, np.integer)):
            return x.item()
        return x
    return json.dumps(conv(res), indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--out", default=str(HERE / "results" / "sobol_analysis"))
    a = ap.parse_args()
    recs, skipped = load_runs(a.runs)
    print("runs: %d usable, %d skipped%s" % (len(recs), len(skipped), (": " + "; ".join(skipped)) if skipped else ""))
    keys = ["name", "seed", "alpha"] + list(FACTORS) + ["d_over_L", "flagged"] + OUTPUTS
    with open(a.out + "_runs.csv", "w") as f:
        f.write(",".join(keys) + "\n")
        for r in recs:
            f.write(",".join(str(r[k]) for k in keys) + "\n")
    res = analyse(recs, a.boot, a.runs)
    pathlib.Path(a.out + ".md").write_text(report(res))
    pathlib.Path(a.out + ".json").write_text(to_json(res))
    print("written: %s.md / .json / _runs.csv" % a.out)


if __name__ == "__main__":
    main()
