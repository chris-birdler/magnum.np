#!/usr/bin/env python
"""
Evaluate runs of run_loops.py in reduced units:
w = W/K_d per cycle versus b_peak = B_peak/Js, and beta_eff = d ln w / d ln b.

    python analyze.py runs/*                          # tables, features, ranking, csv, plot
    python analyze.py --pair runs/T_L18 runs/T_L18_f2         # frequency separation
    python analyze.py --compare runs/M_L12_dx3 runs/M_L12_dx15  # two runs at the same b
    python analyze.py --skip 6 --compare runs/S_L12_dx3 runs/S_L12_dx15   # only the last cycles
    python analyze.py --cycles runs/S_L12_dx15                  # per-cycle table and drift

Per run and amplitude, w and b_peak are the MEAN over the cycles after the
first SKIP cycles (--skip, default 1: the first cycle is a transient).


Errors (target: 10 % scatter of w from cycle to cycle, S_TARGET):
  s_cyc     scatter of ln w from cycle to cycle, pooled over all amplitudes of
            one run (more degrees of freedom than one amplitude with 2 or 3
            cycles). Check: s_cyc <= S_TARGET.
  lnw_err   standard error of ln w at one amplitude = max(s_cyc, S_TARGET) / sqrt(n_measured)
            (conservative: with few amplitudes s_cyc has few degrees of freedom
            and can be much too small by chance)
  beta_err  from lnw_err of the two neighbour amplitudes (central difference);
            the error of b_peak is small and is ignored
  features  errors are interpolated like the features
  ranking   the seeds of the base point give the realisation scatter of each
            feature; see ranking()

Flags per amplitude:
  !   scatter above the target: also the last cycle differs by more than
      dW_tol (0.10) from the cycle before (information: with 10 % scatter
      about 1 stage in 4 has this flag)
  c   loop not closed: |M(end) - M(start)| > closure_tol (drift)
  *   energy balance |w_dis/w_loop - 1| > 2 % (information: w_dis comes from
      257 point samples per cycle and can miss short Barkhausen jumps)
  #   the largest number of large-angle neighbour pairs per cycle changes by
      more than 20 % (and by at least 3 pairs) between the measured cycles:
      vortex-core / Bloch-point events -> mesh dependent

Features per run (for the ranking of the parameters): beta at B_peak = 10,
50 and 100 mT from a local fit of ln w over ln b (see features()).
--w dis uses the LLG dissipation per cycle instead of the loop area. Both
agree on average; w_dis scatters less when the cycles are not periodic.

Comparison (--compare A B): ln w of B is interpolated (linear in ln b) to
the b values of A inside the common b range. A weighted fit
d(b) = ln w_A - ln w_B = c0 + d_beta * (ln b - <ln b>) gives the level
difference c0 (relative change of w) and the change of the slope d_beta, with
errors. Verdict against tolerances fixed before the runs (--tol_lnw, --tol_beta):
PASS if the 90 % interval of c0 and of d_beta lies inside the tolerance.

Frequency separation (--pair): W(f) = W_h + c f  ->  W_h = (f2 W1 - f1 W2)/(f2 - f1).
"""
import argparse
import csv
import json
import math
import pathlib
import sys

import warnings

import numpy as np

JS_T = 1.5            # Js of the unit system (units.py): b = B_peak / Js
B_REF_MT = (10.0, 50.0, 100.0)   # study range 10 ... 150 mT (Chris, 2026-09-30)
WKEY = "w_loop"       # per-cycle loss used for w: w_loop (oint H dB) or w_dis (LLG dissipation), --w
SKIP = 1            # cycles discarded at the start of each amplitude (--skip)
S_TARGET = 0.10
# parameters that define a run (the seed is the realisation); numerical ones are not ranked
PARAMS = ["d_lex", "phi", "Leff_lex", "Q_eff", "r_p", "dx_lex", "f_rel", "alpha", "direction"]
NUMERICAL = {"dx_lex", "f_rel", "alpha", "direction"}
# two-sided 95 % Student-t quantiles; the realisation scatter comes from few seeds
T975 = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57, 6: 2.45, 7: 2.36, 8: 2.31, 9: 2.26, 10: 2.23}


def run_key(cfg):
    """Parameters as given on the command line (Q_eff = None: derived from L_eff)."""
    return tuple(tuple(cfg[k]) if isinstance(cfg.get(k), list) else cfg.get(k) for k in PARAMS)


def param_value(cfg, k):
    v = cfg.get(k)
    if k == "Q_eff" and v is None:
        v = cfg.get("anisotropy", {}).get("Q_eff")
    return tuple(v) if isinstance(v, list) else v


def load_run(d):
    d = pathlib.Path(d)
    cfg = json.loads((d / "config.json").read_text())
    summ = json.loads((d / "summary.json").read_text())
    rows = []
    ss, dof = 0.0, 0
    for st in summ["stages"]:
        if not st.get("measure") or not st.get("complete") or not st["cycles"]:
            continue
        cyc = st["cycles"][SKIP:] if len(st["cycles"]) > SKIP else st["cycles"][-1:]
        w = np.array([c[WKEY] for c in cyc])
        if len(w) > 1 and np.all(w > 0):
            lw = np.log(w)
            ss += float(np.sum((lw - lw.mean()) ** 2))
            dof += len(w) - 1
        n60 = [c.get("n_pairs_gt60", 0) for c in cyc]
        rows.append({"stage": st["name"], "h_amp": st["h_amp"],
                     "b_peak": float(np.mean([c["b_peak"] for c in cyc])),
                     "w_loop": float(w.mean()),
                     "w_err": float(w.std(ddof=1) / math.sqrt(len(w))) if len(w) > 1 else float("nan"),
                     "w_dis": float(np.mean([c["w_dis"] for c in cyc])),
                     "w_loop_raw": float(np.mean([c["w_loop"] for c in cyc])),
                     "closure": cyc[-1]["closure"], "dW_rel": cyc[-1].get("dW_rel", float("nan")),
                     "steady": st.get("steady", False), "closed": st.get("closed", True),
                     "n_cycles": len(st["cycles"]), "n_meas": len(w),
                     "frac_gt60": float(max(c.get("frac_pairs_gt60", 0.0) for c in cyc)),
                     "n60_jump": bool(max(n60) > 1.2 * max(min(n60), 1) and max(n60) - min(n60) >= 3),
                     "mean_dt": float(np.mean([c["mean_dt"] for c in cyc]))})
    rows.sort(key=lambda r: r["b_peak"])
    cfg["done"] = (d / "DONE").exists()
    cfg["s_cyc"] = math.sqrt(ss / dof) if dof > 0 else float("nan")
    cfg["s_cyc_dof"] = dof
    for r in rows:
        s_use = cfg["s_cyc"] if cfg["s_cyc"] > S_TARGET else S_TARGET      # also for nan
        r["lnw_err"] = s_use / math.sqrt(r["n_meas"])
        r["balance"] = r["w_dis"] / r["w_loop_raw"] if r["w_loop_raw"] != 0 else float("nan")
    add_beta(rows, "w_loop", "beta", "lnw_err")
    return cfg, rows


def add_beta(rows, key, out_key, err_key=None):
    """Central difference of ln(key) over ln b; with err_key (error of ln key) also out_key_err."""
    n = len(rows)
    lb = [math.log(r["b_peak"]) for r in rows]
    lw = [math.log(r[key]) if r[key] > 0 else float("nan") for r in rows]
    for i in range(n):
        if n < 2:
            rows[i][out_key] = float("nan")
            if err_key:
                rows[i][out_key + "_err"] = float("nan")
            continue
        j0, j1 = max(i - 1, 0), min(i + 1, n - 1)
        rows[i][out_key] = (lw[j1] - lw[j0]) / (lb[j1] - lb[j0])
        if err_key:
            rows[i][out_key + "_err"] = math.hypot(rows[j0][err_key], rows[j1][err_key]) / abs(lb[j1] - lb[j0])


def features(rows):
    """beta at B_peak = 10 / 50 / 100 mT (the study range) and its error.
    Local weighted fit of ln w over ln b with all amplitudes within a factor
    2.2 of b_ref (at the ends of the range the fit is one-sided)."""
    f, e = {}, {}
    if len(rows) < 2:
        return f, e
    lb = np.log([r["b_peak"] for r in rows])
    lw = np.array([math.log(r["w_loop"]) if r["w_loop"] > 0 else np.nan for r in rows])
    sw = np.array([r["lnw_err"] for r in rows])
    for mT in B_REF_MT:
        k = "beta@%gmT" % mT
        lbr = math.log(mT / 1000.0 / JS_T)
        sel = np.isfinite(lw) & (np.abs(lb - lbr) <= math.log(2.2))
        if sel.sum() < 2 or not (lb.min() - 0.1 <= lbr <= lb.max() + 0.1):
            f[k] = e[k] = float("nan")
            continue
        x, y, wgt = lb[sel] - lbr, lw[sel], 1.0 / sw[sel] ** 2
        X = np.c_[np.ones_like(x), x]
        cov = np.linalg.inv(X.T @ (X * wgt[:, None]))
        coef = cov @ (X.T @ (wgt * y))
        f[k], e[k] = float(coef[1]), float(math.sqrt(cov[1, 1]))
    return f, e
    lb = np.log([r["b_peak"] for r in rows])
    be = np.array([r["beta"] for r in rows])
    bs = np.array([r["beta_err"] for r in rows])
    ok = np.isfinite(be)
    for b in B_REF:
        k = "beta@b=%g" % b
        lbr = math.log(b)
        if ok.sum() >= 2 and lb[ok][0] <= lbr <= lb[ok][-1]:
            f[k] = float(np.interp(lbr, lb[ok], be[ok]))
            e[k] = float(np.interp(lbr, lb[ok], bs[ok]))
        else:
            f[k] = e[k] = float("nan")
    if ok.any():
        i = int(np.nanargmax(be))
        f["beta_max"], e["beta_max"] = float(be[i]), float(bs[i])
        f["b_at_beta_max"], e["b_at_beta_max"] = float(rows[i]["b_peak"]), float("nan")
    return f, e


def label(cfg):
    return ("d/lex=%.1f phi=%.3f L_eff/lex=%g Q_eff=%.3g r_p=%g dx/lex=%g alpha=%g f/fM=%.3g dir=%s seed=%d" %
            (cfg.get("d_lex_used", cfg["d_lex"]), cfg["phi_vox"], cfg["Leff_lex"], param_value(cfg, "Q_eff"),
             cfg["r_p"], cfg["dx_lex"], cfg["alpha"], cfg["f_rel"], cfg["direction"], cfg["seed"]))


def flags(r):
    s = "" if (r["steady"] or not r["closed"]) else "!"
    if not r["closed"]:
        s += "c"
    if abs(r["balance"] - 1.0) > 0.02:
        s += "*"
    if r["n60_jump"]:
        s += "#"
    return s


def print_table(name, cfg, rows):
    print("\n== %s%s\n   %s" % (name, "" if cfg["done"] else "   (NOT FINISHED)", label(cfg)))
    print("   scatter per cycle s_cyc = %.3f (dof %d, target %.2f)%s" %
          (cfg["s_cyc"], cfg["s_cyc_dof"], S_TARGET, "  ABOVE TARGET" if cfg["s_cyc"] > S_TARGET else ""))
    print("   %-6s %9s %9s %6s %10s %8s %7s %6s %6s %8s %7s %3s %s" %
          ("stage", "h", "b_peak", "B[mT]", "w(%s)" % WKEY[2:], "w_err", "beta", "+-", "wd/w", "closure", "f>60",
           "n", "flags"))
    for r in rows:
        print("   %-6s %9.4g %9.4g %6.1f %10.4g %8.2g %7.3f %6.3f %6.3f %8.2g %7.2g %3d %s" %
              (r["stage"], r["h_amp"], r["b_peak"], 1000 * JS_T * r["b_peak"], r["w_loop"], r["w_err"],
               r["beta"], r["beta_err"],
               r["balance"], r["closure"], r["frac_gt60"], r["n_cycles"], flags(r)))
    f, e = features(rows)
    if f:
        print("   features: " + "  ".join("%s=%.3g+-%.2g" % (k, f[k], e[k]) if k != "b_at_beta_max"
                                          else "%s=%.3g" % (k, f[k]) for k in f))


def ranking(results, path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)      # nanmean of all-nan features
        return _ranking(results, path)


def _ranking(results, path):
    """
    Change of each feature against the base point, with errors.

    Base point = the largest group of runs that differ only in the seed.
      var_seed   variance of a feature over the base seeds (realisation + statistics)
      var_real   = max(var_seed - mean(stat. error^2 of the seeds), 0)
      error of one run      sigma_run^2  = stat. error^2 + var_real
      error of the base     sigma_base^2 = var_seed / n_seeds   (mean over the seeds)
      error of the change   sigma_delta  = sqrt(sigma_run^2 + sigma_base^2)
    Only runs that differ from the base in ONE physical parameter are ranked
    (controls change dx, f or alpha). sens = delta / delta ln(p) where defined.
    Sorted by the largest significance |delta| / sigma_delta over the features.
    """
    groups = {}
    for name, cfg, rows in results:
        if cfg.get("max_cycles") is not None or not cfg["done"] or len(rows) < 2:   # benchmark or unfinished
            continue
        groups.setdefault(run_key(cfg), []).append((name, cfg, rows))
    if not groups:
        print("\nranking: no finished runs")
        return []
    base_key, base = max(groups.items(), key=lambda kv: len(kv[1]))
    if len(base) < 2:
        print("\nranking: no base point with more than one seed -> realisation scatter unknown, "
              "the errors contain only the statistics of the cycles")
    bf = [features(rows) for _, _, rows in base]
    keys = [k for k in bf[0][0] if k != "b_at_beta_max"]
    t_crit = T975.get(len(base) - 1, 1.96) if len(base) > 1 else float("nan")
    base_mean, base_var, real_var = {}, {}, {}
    for k in keys:
        v = np.array([f.get(k, np.nan) for f, _ in bf])
        s = np.array([e.get(k, np.nan) for _, e in bf])
        base_mean[k] = float(np.nanmean(v))
        if np.sum(np.isfinite(v)) >= 2:
            var_seed = float(np.nanvar(v, ddof=1))
            real_var[k] = max(var_seed - float(np.nanmean(s ** 2)), 0.0)
            base_var[k] = var_seed / np.sum(np.isfinite(v))
        else:
            real_var[k] = 0.0
            base_var[k] = float(np.nanmean(s ** 2)) / max(np.sum(np.isfinite(v)), 1)

    print("\n== ranking against the base point (%s, %d seeds: %s)" %
          (", ".join("%s=%s" % kv for kv in zip(PARAMS, base_key)), len(base), " ".join(n for n, _, _ in base)))
    print("   realisation scatter (1 sigma): " +
          "  ".join("%s=%.3f" % (k, math.sqrt(real_var[k])) for k in keys))
    print("   significance: |d|/sig > t = %.2f (95 %%, Student-t with %d dof of the seed scatter); "
          "%d features x runs are tested, expect some false hits at this level" %
          (t_crit, len(base) - 1, len(keys)))
    out = []
    for key, runs in groups.items():
        if key == base_key:
            continue
        diff = [p for p, a, b in zip(PARAMS, key, base_key) if a != b]
        if len(diff) != 1 or diff[0] in NUMERICAL:
            continue
        p = diff[0]
        p_run, p_base = key[PARAMS.index(p)], base_key[PARAMS.index(p)]
        try:
            dlnp = math.log(p_run / p_base) if p_run > 0 and p_base > 0 else float("nan")
        except TypeError:
            dlnp = float("nan")
        p_run = p_run if isinstance(p_run, (int, float)) else float("nan")
        p_base = p_base if isinstance(p_base, (int, float)) else float("nan")
        for name, cfg, rows in runs:
            f, e = features(rows)
            rec = {"run": name, "param": p, "value": p_run, "base": p_base, "dlnp": dlnp, "z_max": 0.0}
            for k in keys:
                d = f.get(k, float("nan")) - base_mean[k]
                s = math.sqrt(e.get(k, float("nan")) ** 2 + real_var[k] + base_var[k])
                rec["d_" + k], rec["err_" + k] = d, s
                rec["sens_" + k] = d / dlnp if math.isfinite(dlnp) and dlnp != 0 else float("nan")
                if math.isfinite(d) and s > 0:
                    rec["z_max"] = max(rec["z_max"], abs(d) / s)
            out.append(rec)
    out.sort(key=lambda r: -r["z_max"])
    print("   %-14s %-8s %9s %9s " % ("run", "param", "value", "base") +
          " ".join("%-17s" % ("d " + k) for k in keys) + "  |d|/sig max")
    for r in out:
        print("   %-14s %-8s %9.4g %9.4g " % (r["run"], r["param"], r["value"], r["base"]) +
              " ".join("%+7.3f +- %-6.3f" % (r["d_" + k], r["err_" + k]) for k in keys) +
              "  %5.1f" % r["z_max"])
    if out:
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)
        print("   ranking csv: %s  (sens_* = change per unit of ln p)" % path)
    return out


def pair(d1, d2):
    c1, r1 = load_run(d1)
    c2, r2 = load_run(d2)
    f1, f2 = c1["f_rel"], c2["f_rel"]
    ignore = {"f_rel", "freq", "out", "max_cycles", "stages", "git_commit", "torch", "device", "gpu",
              "no_resume", "timer", "vti", "keep_stage_ckpt", "init_from", "s_cyc", "s_cyc_dof", "units"}
    diff = [k for k in c1 if k not in ignore and c1.get(k) != c2.get(k)]
    if diff:
        print("WARNING: runs differ in more than the frequency: %s" % diff)
    m2 = {r["stage"]: r for r in r2}
    out = []
    for a in r1:
        b = m2.get(a["stage"])
        if b is None:
            continue
        w_h = (f2 * a["w_loop"] - f1 * b["w_loop"]) / (f2 - f1)
        s1, s2 = a["w_loop"] * a["lnw_err"], b["w_loop"] * b["lnw_err"]
        s_h = math.hypot(f2 * s1, f1 * s2) / abs(f2 - f1)
        dyn = 1.0 - w_h / a["w_loop"]
        # dyn = f1 (w2/w1 - 1) / (f2 - f1); the errors of w1 and w2 are independent
        s_dyn = abs(f1 / (f2 - f1)) * math.hypot(s2 / a["w_loop"], b["w_loop"] * s1 / a["w_loop"] ** 2)
        out.append({"stage": a["stage"], "b_peak": 0.5 * (a["b_peak"] + b["b_peak"]),
                    "db_rel": abs(a["b_peak"] - b["b_peak"]) / a["b_peak"],
                    "w_f1": a["w_loop"], "w_f2": b["w_loop"], "w_h": w_h,
                    "lnw_h_err": s_h / w_h if w_h > 0 else float("nan"), "lnw_f1_err": a["lnw_err"],
                    "dyn_share_f1": dyn, "dyn_err": s_dyn})
    out.sort(key=lambda r: r["b_peak"])
    add_beta(out, "w_h", "beta_h", "lnw_h_err")
    add_beta(out, "w_f1", "beta_f1", "lnw_f1_err")
    print("\n== frequency separation  f1/fM=%.3g  f2/fM=%.3g" % (f1, f2))
    print("   %-6s %9s %7s %10s %10s %10s %14s %14s %14s" %
          ("stage", "b_peak", "db/b", "w(f1)", "w(f2)", "w_h", "dyn(f1)", "beta_h", "beta_f1"))
    for r in out:
        print("   %-6s %9.4g %7.2g %10.4g %10.4g %10.4g %6.2f +- %-4.2f %6.3f +- %-4.2f %6.3f +- %-4.2f" %
              (r["stage"], r["b_peak"], r["db_rel"], r["w_f1"], r["w_f2"], r["w_h"],
               r["dyn_share_f1"], r["dyn_err"], r["beta_h"], r["beta_h_err"], r["beta_f1"], r["beta_f1_err"]))
    return out


def compare(d1, d2, tol_lnw=0.15, tol_beta=0.2):
    c1, r1 = load_run(d1)
    c2, r2 = load_run(d2)
    lb2 = np.log([r["b_peak"] for r in r2])
    lw2 = np.log([r["w_loop"] for r in r2])
    e2 = np.array([r["lnw_err"] for r in r2])
    x, y, s = [], [], []
    print("\n== compare  A = %s\n             B = %s" % (d1, d2))
    print("   A: %s\n   B: %s" % (label(c1), label(c2)))
    print("   %-6s %9s %10s %10s %8s %6s" % ("stage", "b(A)", "ln w A", "ln w B(b)", "d", "+-"))
    for r in r1:
        lb = math.log(r["b_peak"])
        if not (lb2.min() <= lb <= lb2.max()) or r["w_loop"] <= 0:
            continue
        j = int(np.clip(np.searchsorted(lb2, lb), 1, len(lb2) - 1))
        t = (lb - lb2[j - 1]) / (lb2[j] - lb2[j - 1])
        yb = (1 - t) * lw2[j - 1] + t * lw2[j]
        sb = math.hypot((1 - t) * e2[j - 1], t * e2[j])
        d = math.log(r["w_loop"]) - yb
        sd = math.hypot(r["lnw_err"], sb)
        x.append(lb); y.append(d); s.append(sd)
        print("   %-6s %9.4g %10.4f %10.4f %+8.4f %6.3f" % (r["stage"], r["b_peak"], math.log(r["w_loop"]), yb, d, sd))
    if len(x) < 2:
        print("   fewer than 2 points in the common b range: no fit")
        return None
    x, y, s = np.array(x), np.array(y), np.array(s)
    wgt = 1.0 / s**2
    xm = float(np.sum(wgt * x) / np.sum(wgt))
    X = np.c_[np.ones_like(x), x - xm]
    cov = np.linalg.inv(X.T @ (X * wgt[:, None]))
    c0, db = cov @ (X.T @ (wgt * y))
    sc0, sdb = math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1])
    chi2 = float(np.sum(wgt * (y - X @ np.array([c0, db])) ** 2))
    z = 1.645
    ok_l = abs(c0) + z * sc0 < tol_lnw
    ok_b = abs(db) + z * sdb < tol_beta
    print("   level  c0     = %+.3f +- %.3f  (w_A/w_B = %.3f)  90%% within +-%.2f: %s" %
          (c0, sc0, math.exp(c0), tol_lnw, "PASS" if ok_l else "FAIL"))
    print("   slope  d_beta = %+.3f +- %.3f                     90%% within +-%.2f: %s" %
          (db, sdb, tol_beta, "PASS" if ok_b else "FAIL"))
    print("   chi2 = %.2f for %d dof (a large value: the difference is not linear in ln b)" % (chi2, len(x) - 2))
    return {"c0": c0, "c0_err": sc0, "d_beta": db, "d_beta_err": sdb, "chi2": chi2, "n": len(x),
            "pass": bool(ok_l and ok_b)}


def cycles(d):
    """Per-cycle table of every measured stage, and the drift of w over the kept
    cycles (relative change per cycle from a linear fit, with its error)."""
    d = pathlib.Path(d)
    summ = json.loads((d / "summary.json").read_text())
    print("\n== %s  (kept: cycles >= %d)" % (d, SKIP))
    for st in summ["stages"]:
        if not st.get("measure"):
            continue
        print("   %s  h = %.4g" % (st["name"], st["h_amp"]))
        print("     %3s %9s %11s %11s %7s %8s %5s" % ("c", "b_peak", "w_loop", "w_dis", "wd/w", "closure", "n60"))
        for c in st["cycles"]:
            print("     %3d %9.4g %11.4g %11.4g %7.3f %8.4f %5d%s" %
                  (c["cycle"], c["b_peak"], c["w_loop"], c["w_dis"],
                   c["w_dis"] / c["w_loop"] if c["w_loop"] else float("nan"), c["closure"],
                   c.get("n_pairs_gt60", 0), "" if c["cycle"] >= SKIP else "   (skipped)"))
        kept = st["cycles"][SKIP:]
        if len(kept) >= 3:
            for key in ("w_loop", "w_dis"):
                w = np.array([c[key] for c in kept])
                k = np.arange(len(w), dtype=float)
                A = np.c_[np.ones_like(k), k - k.mean()]
                coef, res, _, _ = np.linalg.lstsq(A, w, rcond=None)
                dof = len(w) - 2
                s2 = float(np.sum((w - A @ coef) ** 2)) / dof if dof > 0 else float("nan")
                se = math.sqrt(s2 / np.sum((k - k.mean()) ** 2)) if dof > 0 else float("nan")
                print("     %-6s mean %.4g  scatter %.1f %%  drift %+.2f +- %.2f %% per cycle" %
                      (key, w.mean(), 100 * w.std(ddof=1) / abs(w.mean()), 100 * coef[1] / w.mean(),
                       100 * se / abs(w.mean())))


def plot(results, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: no plot"); return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for name, cfg, rows in results:
        if not rows:
            continue
        b = [r["b_peak"] for r in rows]
        ax1.errorbar(b, [r["w_loop"] for r in rows], yerr=[r["w_loop"] * r["lnw_err"] for r in rows],
                     fmt="o-", ms=3, capsize=2, label=name)
        ax2.errorbar(b, [r["beta"] for r in rows], yerr=[r["beta_err"] for r in rows],
                     fmt="o-", ms=3, capsize=2, label=name)
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax1.set_xlabel("b = B_peak / Js"); ax1.set_ylabel("w = W / K_d per cycle")
    ax2.set_xlabel("b = B_peak / Js"); ax2.set_ylabel(r"$\beta_\mathrm{eff}$ = d ln w / d ln b")
    for v in (2.0, 3.0):
        ax2.axhline(v, color="0.7", lw=0.8, ls="--")
    ax1.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print("plot: %s" % path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--pair", nargs=2, metavar=("RUN_F1", "RUN_F2"))
    ap.add_argument("--compare", nargs=2, metavar=("RUN_A", "RUN_B"))
    ap.add_argument("--cycles", action="store_true", help="per-cycle tables and drift of the given runs")
    ap.add_argument("--skip", type=int, default=1, help="cycles discarded at the start of each amplitude")
    ap.add_argument("--w", choices=["loop", "dis"], default="loop", help="loss per cycle: loop area or LLG dissipation")
    ap.add_argument("--tol_lnw", type=float, default=0.15, help="tolerance of the level difference (ln w)")
    ap.add_argument("--tol_beta", type=float, default=0.2, help="tolerance of the slope difference")
    ap.add_argument("--csv", default="all_runs.csv")
    ap.add_argument("--ranking_csv", default="ranking.csv")
    ap.add_argument("--plot", default="beta_vs_b.png")
    a = ap.parse_args()
    global SKIP, WKEY
    SKIP = a.skip
    WKEY = "w_" + a.w

    if a.cycles:
        for d in a.runs:
            cycles(d)
        return
    if a.pair:
        pair(*a.pair)
        return
    if a.compare:
        compare(*a.compare, tol_lnw=a.tol_lnw, tol_beta=a.tol_beta)
        return

    results = []
    for d in a.runs:
        p = pathlib.Path(d)
        if not (p / "summary.json").exists():
            continue
        cfg, rows = load_run(p)
        results.append((p.name, cfg, rows))
        print_table(p.name, cfg, rows)
    if not results:
        sys.exit("no runs with summary.json found")

    keys = PARAMS + ["phi_vox", "d_lex_used", "seed", "s_cyc"]
    rk = next((list(rows[0].keys()) for _, _, rows in results if rows), [])
    with open(a.csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["run"] + keys + rk)
        for name, cfg, rows in results:
            for r in rows:
                w.writerow([name] + [param_value(cfg, k) for k in keys] + [r[k] for k in rk])
    print("\ncsv: %s" % a.csv)
    ranking(results, a.ranking_csv)
    plot(results, a.plot)


if __name__ == "__main__":
    main()
