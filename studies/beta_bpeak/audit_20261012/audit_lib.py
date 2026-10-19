"""Independent loader for the audit (2026-10-12). Reads summary.json / config.json directly.
No import of the study analysis code, so that the audit does not inherit its errors.

Rules implemented (owner decisions, PROTOKOLL 7.54 / PLAN 4.5):
- measured cycles of a stage: cycles from index 2 on if the drive was corrected after cycle 0
  ("h_amp_initial" in the stage), else from index 1 on (old protocol);
- event cycle: w_loop <= 0 or |w_dis / w_loop - 1| > 0.20;
- steady tail: the measured cycles after the last event cycle;
- a stage with fewer than MIN_CALM kept cycles is missing.
Loss density P = f * K_d * w (W/m^3), reported in W/cm^3 (1 W/cm^3 = 1e6 W/m^3).
"""
import json
import math
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
STUDY = HERE.parent
RUNS = STUDY / "runs"
TOL = 0.20
MIN_CALM = 2
AMPS = (9.0, 20.0, 35.0, 50.0, 70.0, 100.0, 150.0)


def read(run):
    run = pathlib.Path(run)
    if not run.is_absolute():
        run = RUNS / run
    cfg = json.loads((run / "config.json").read_text())
    summ = json.loads((run / "summary.json").read_text()) if (run / "summary.json").exists() else None
    return run, cfg, summ


def is_event(c, tol=TOL):
    return c["w_loop"] <= 0 or abs(c["w_dis"] / c["w_loop"] - 1.0) > tol


def measured(st):
    skip = 2 if "h_amp_initial" in st else 1
    return st["cycles"][skip:] if len(st["cycles"]) > skip else []


def tail(cyc, tol=TOL):
    ev = [i for i, c in enumerate(cyc) if is_event(c, tol)]
    return cyc[ev[-1] + 1:] if ev else list(cyc)


def p_factor(cfg):
    """W/cm^3 per unit of w (w = W / K_d)."""
    return cfg["freq"] * cfg["units"]["Kd"] / 1e6


def stages(run, rule="steady", min_calm=MIN_CALM, tol=TOL, key="w_loop", last=None):
    """List of dicts per measured stage: B_mT (mean of the kept cycles), w, P (W/cm^3), n_kept, missing, ...
    rule: 'steady' (tail after the last event) or 'all' (all measured cycles).
    last: if given, use only the last `last` kept cycles."""
    run, cfg, summ = read(run)
    out = []
    pf = p_factor(cfg)
    for si, st in enumerate(summ["stages"]):
        if not st.get("measure") or not st["cycles"]:
            continue
        cyc = measured(st)
        if not cyc:
            continue
        kept = tail(cyc, tol) if rule == "steady" else list(cyc)
        if last:
            kept = kept[-last:]
        miss = len(kept) < (min_calm if rule == "steady" else 1)
        use = kept if kept else cyc
        w = np.array([c[key] for c in use])
        wd = np.array([c["w_dis"] for c in use])
        b = np.array([c["b_peak"] for c in use])
        rec = {"si": si, "stage": st["name"], "B_mT": float(1000 * cfg["Js"] * b.mean()),
               "w": float("nan") if miss else float(w.mean()),
               "w_dis": float("nan") if miss else float(wd.mean()),
               "n_meas": len(cyc), "n_kept": len(kept) if not miss else len(kept),
               "missing": bool(miss), "n_cycles": len(st["cycles"]),
               "cyc_sd_ln": float(np.std(np.log(w), ddof=1)) if (len(w) > 1 and np.all(w > 0)) else float("nan"),
               "complete": st.get("complete", True), "b_target_mT": 1000 * cfg["Js"] * st.get("b_target", float("nan"))
               if st.get("b_target") is not None else float("nan")}
        rec["P"] = rec["w"] * pf
        out.append(rec)
    out.sort(key=lambda r: r["B_mT"])
    return cfg, summ, out


def nominal(B):
    return min(AMPS, key=lambda a: abs(math.log(B / a)))


def rerun_of(name):
    r = RUNS / (name + "_mc2")
    return r.name if (r / "DONE").exists() else None


def use_name(name, reruns=True):
    if reruns:
        r = rerun_of(name)
        if r:
            return r
    return name


def powerfit(lnB, lnP):
    """least squares ln P = a + beta ln B; returns beta, SE(beta) (from residuals, nan if dof 0), residuals."""
    lnB, lnP = np.asarray(lnB, float), np.asarray(lnP, float)
    X = np.c_[np.ones_like(lnB), lnB]
    coef, *_ = np.linalg.lstsq(X, lnP, rcond=None)
    res = lnP - X @ coef
    dof = len(lnB) - 2
    if dof > 0:
        s2 = float(res @ res) / dof
        se = math.sqrt(s2 / float(((lnB - lnB.mean()) ** 2).sum()))
    else:
        se = float("nan")
    return float(coef[1]), se, res, float(coef[0])


def tq(n, q=0.975):
    from scipy import stats
    return float(stats.t.ppf(q, n - 1)) if n > 1 else float("nan")


def mse(v):
    """mean, SE, n (finite values)."""
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    n = len(v)
    if n == 0:
        return float("nan"), float("nan"), 0
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan"), n


def fmt_ci(m, se, n, q=0.975):
    t = tq(n, q)
    return "%+.3f ± %.3f (%s %d%% CI %+.3f … %+.3f, n %d)" % (m, se, "t", int(round(200 * q - 100)), m - t * se, m + t * se, n)
