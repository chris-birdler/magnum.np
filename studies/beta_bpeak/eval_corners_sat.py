#!/usr/bin/env python
"""
Corner study of PLAN 4.15, phase 1 (readings fixed in PLAN 4.15 before the results): the runs K<c>_<k>_sat (static
saturation from the uniform start at 2 M_s, PROTOKOLL 7.132), n seeds per corner.
Per run (steady-tail rule with at least 2 calm cycles; ln w at the nominal amplitudes, analyze_design.run_record):
  β all     slope of ln w over ln B over the 7 amplitudes, at least 6 with a value (main-fit definition, 7.119 "full")
  windows   PLAN 4.0: 9-20-35, 35-50-70, 70-100-150 mT (analyze_design); PROTOKOLL 7.119: 9 … 50 mT (all 4),
            70 … 150 mT (all 3)
  ln P      at 9 / 50 / 100 mT (P in W/cm³ of the core)
A run has no value if its prepared state has a Bloch point after the last pass (PLAN 4.15 D7 (b)); the report lists it.
Per corner: mean ± SE, sample SD, 95 % CI (t distribution, n − 1 degrees of freedom). Pooled SD: square root of the
summed squared deviations from the corner means divided by the summed degrees of freedom. Levers β(Ki) − β(K1) with
the SE from the two corner SEs. Target of D1: half-width of the 95 % CI of β all at most 0.10 per corner; else the
smallest n with t(n − 1) · SD / √n ≤ 0.10 for the corner SD.
Flags per stage: B = a Bloch point at one sample or more (bloch_<i>_<stage>.csv); S9 ε and w (eval_corners_nx).
β all also without the stages with flag B, and without the stages with an S9 flag.
Bloch points of the preparation: the runs free after the first, second and third pass. Time per run: preparation
(init.wall_s plus the passes) and cycles. Information: Δβ (all, 9 … 50 mT) against the ad, nx and old runs of the
same microstructure (k = 0 … 2, the seeds of the old corners); the energy-terms check (S5).

    python eval_corners_sat.py [--runs DIR] [--json FILE] > results/corners_sat_report.md
"""
import argparse
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import conv_jobs as CJ
import eval_corners_ad as EA
import eval_corners_nx as EN
import units as U

HERE = pathlib.Path(__file__).resolve().parent
CORNERS = ("K1", "K2", "K3", "K4")
AMPS = list(AD.AMPS_MT)
OUTPUTS = ["beta_all", "beta@9-35mT", "beta@35-70mT", "beta@70-150mT", "beta_low", "beta_high",
           "lnP@9mT", "lnP@50mT", "lnP@100mT"]
LABEL = {"beta_all": "β 9 … 150 mT", "beta@9-35mT": "β 9-20-35 mT", "beta@35-70mT": "β 35-50-70 mT",
         "beta@70-150mT": "β 70-100-150 mT", "beta_low": "β 9 … 50 mT", "beta_high": "β 70 … 150 mT",
         "lnP@9mT": "ln P 9 mT", "lnP@50mT": "ln P 50 mT", "lnP@100mT": "ln P 100 mT"}
TARGET = 0.10                     # PLAN 4.15 D1: half-width of the 95 % CI of β all per corner


def sat_runs(root):
    """Names K<c>_<k>_sat with a DONE file, sorted by corner and k."""
    pat = AD.design_pattern(CJ.SAT_SUFFIX)
    ds = [d.name for d in root.iterdir() if d.is_dir() and pat.match(d.name) and d.name.startswith("K")
          and (d / "DONE").exists()]
    return sorted(ds, key=lambda n: (n[:2], int(n.split("_")[1])))


def beta_all(lnw):
    return EA.window_beta([lnw], "full")[0]


def outputs(lnw, rec, lnF):
    """The per-run outputs (OUTPUTS) from ln w at the nominal amplitudes and the run record."""
    o = {"beta_all": beta_all(lnw), "beta_low": EA.window_beta([lnw], "low")[0],
         "beta_high": EA.window_beta([lnw], "high")[0]}
    for k in ("beta@9-35mT", "beta@35-70mT", "beta@70-150mT"):
        o[k] = rec[k] if rec else float("nan")
    for mT in (9.0, 50.0, 100.0):
        o["lnP@%gmT" % mT] = float(lnw[AMPS.index(mT)] + lnF)
    return o


def stage_flags_b(d):
    """Flag B per stage index 0 … 6: a Bloch point at one sample or more (None if the count file is missing)."""
    out = []
    for i in range(len(AMPS)):
        f = sorted(d.glob("bloch_%02d_*.csv" % i))
        if not f:
            out.append(None)
            continue
        a = np.loadtxt(f[0], delimiter=",", comments="#", ndmin=2)
        out.append(bool(len(a) and (a[:, 2] + a[:, 3]).max() > 0))
    return out


def prep_record(d):
    """(passes, BP count after each pass, free at the end, preparation hours, cycle hours) from summary.json."""
    s = json.loads((d / "summary.json").read_text())
    i = s.get("init", {})
    passes = i.get("predemag_passes") or []
    nb = [len(p.get("bloch_points") or []) for p in passes]          # the Bloch points after each pass
    prep_h = (i.get("wall_s", 0.0) + sum(p.get("wall_s", 0.0) for p in passes)) / 3600.0
    cyc_h = sum(c.get("wall_s", 0.0) for st in s.get("stages", []) for c in st.get("cycles", [])) / 3600.0
    return len(passes), nb, bool(i.get("bloch_free", True)), prep_h, cyc_h


def cstat(v):
    """(mean, SE, SD, CI low, CI high, n) with the 95 % CI of the t distribution; None for fewer than 2 values."""
    v = np.array([x for x in v if np.isfinite(x)])
    if len(v) < 2:
        return None
    sd = v.std(ddof=1)
    se = sd / math.sqrt(len(v))
    h = stats.t.ppf(0.975, len(v) - 1) * se
    return float(v.mean()), float(se), float(sd), float(v.mean() - h), float(v.mean() + h), len(v)


def pooled_sd(groups):
    """Pooled SD over groups of values (NaN removed); (SD, degrees of freedom) or None."""
    ss, df = 0.0, 0
    for v in groups:
        v = np.array([x for x in v if np.isfinite(x)])
        if len(v) > 1:
            ss += float(((v - v.mean()) ** 2).sum())
            df += len(v) - 1
    return (math.sqrt(ss / df), df) if df else None


def n_needed(sd, target=TARGET, n_max=200):
    """Smallest n >= 2 with t(n − 1) · sd / √n <= target (None if above n_max)."""
    for n in range(2, n_max + 1):
        if stats.t.ppf(0.975, n - 1) * sd / math.sqrt(n) <= target:
            return n
    return None


def fmt(s, f="%.3f"):
    return ("%s ± %s" % (f % s[0], f % s[1])) if s else "–"


def main(root=None, json_out=None):
    root = pathlib.Path(root) if root else HERE / "runs"
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    u = U.units()
    lnF = math.log(u["Kd"] * U.F_REL_DEFAULT * u["f_M"] / 1e6)          # ln w -> ln P (W/cm³)
    runs = sat_runs(root)
    rows, novalue = {}, []
    for n in runs:
        d = root / n
        npass, nb, free, prep_h, cyc_h = prep_record(d)
        rec = AD.run_record(d)
        lnw = EN.lnw(d)
        fb = stage_flags_b(d)
        fs = [any(EN.stage_flags(s)) for s in EN.measured_stages(d)]
        lw_b = np.where([bool(x) for x in fb], np.nan, lnw) if len(fb) == len(lnw) else np.full(len(AMPS), np.nan)
        lw_s = np.where(fs, np.nan, lnw) if len(fs) == len(lnw) else np.full(len(AMPS), np.nan)
        r = {"corner": n[:2], "k": int(n.split("_")[1]), "passes": npass, "bp_after_pass": nb, "free": free,
             "prep_h": prep_h, "cycles_h": cyc_h, "flag_B": fb, "flag_S9": fs,
             "beta_all_noB": beta_all(lw_b), "beta_all_noS9": beta_all(lw_s), "lnP": (lnw + lnF).tolist()}
        r.update(outputs(lnw, rec, lnF))
        rows[n] = r
        if not free:
            novalue.append(n)
    valid = [n for n in runs if rows[n]["free"]]
    L = ["# Corner study of PLAN 4.15, phase 1 (eval_corners_sat.py)", "",
         "Runs K<c>_<k>_sat: static saturation from the uniform start at 2 M_s (PROTOKOLL 7.132). Mean ± SE over the "
         "seeds of a corner; SD = sample SD; CI = 95 % CI (t distribution). P in W/cm³.", "",
         "Runs with a value: %d of %d. Without a value (a Bloch point after the last pass, D7 (b)): %s." % (
             len(valid), len(runs), ", ".join(novalue) or "none"), ""]

    # --- per run --------------------------------------------------------------------------------------------------
    L += ["## Per run", "", "| run | passes | BPs after each pass | " + " | ".join(LABEL[k] for k in OUTPUTS) +
          " | stages with flag B | stages with an S9 flag | preparation / cycles (h) |", "|---|---|---|" +
          "---|" * len(OUTPUTS) + "---|---|---|"]
    for n in runs:
        r = rows[n]
        L.append("| %s | %d | %s | %s | %s | %s | %.2f / %.2f |" % (
            n, r["passes"], " / ".join(str(x) for x in r["bp_after_pass"]),
            " | ".join("%.3f" % r[k] if np.isfinite(r[k]) else "–" for k in OUTPUTS),
            ", ".join("%g" % AMPS[i] for i, x in enumerate(r["flag_B"]) if x) or "–",
            ", ".join("%g" % AMPS[i] for i, x in enumerate(r["flag_S9"]) if x) or "–", r["prep_h"], r["cycles_h"]))

    # --- per corner -----------------------------------------------------------------------------------------------
    S = {k: {c: cstat([rows[n][k] for n in valid if rows[n]["corner"] == c]) for c in CORNERS} for k in OUTPUTS}
    S["beta_all_noB"] = {c: cstat([rows[n]["beta_all_noB"] for n in valid if rows[n]["corner"] == c]) for c in CORNERS}
    S["beta_all_noS9"] = {c: cstat([rows[n]["beta_all_noS9"] for n in valid if rows[n]["corner"] == c])
                          for c in CORNERS}
    L += ["", "## Per corner (mean ± SE; SD; 95 % CI; n)", "", "| output | " + " | ".join(CORNERS) + " |",
          "|---|" + "---|" * len(CORNERS)]
    for k in OUTPUTS + ["beta_all_noB", "beta_all_noS9"]:
        lab = LABEL.get(k, {"beta_all_noB": "β 9 … 150 mT without flag B",
                            "beta_all_noS9": "β 9 … 150 mT without S9 flags"}.get(k))
        L.append("| %s | %s |" % (lab, " | ".join(
            "%s; SD %.3f; %.3f … %.3f; n %d" % (fmt(s), s[2], s[3], s[4], s[5]) if s else "–"
            for s in (S[k][c] for c in CORNERS))))
    P = {k: pooled_sd([[rows[n][k] for n in valid if rows[n]["corner"] == c] for c in CORNERS]) for k in OUTPUTS}
    L += ["", "Pooled SD over the corners: " + "; ".join(
        "%s %.3f (%d degrees of freedom)" % (LABEL[k], P[k][0], P[k][1]) for k in OUTPUTS if P[k]) + "."]

    # --- target and levers -----------------------------------------------------------------------------------------
    L += ["", "## Target of D1: half-width of the 95 %% CI of β 9 … 150 mT at most %.2f" % TARGET, "",
          "| corner | n | half-width | met | n that the corner SD needs |", "|---|---|---|---|---|"]
    tgt = {}
    for c in CORNERS:
        s = S["beta_all"][c]
        if s:
            hw = (s[4] - s[3]) / 2.0
            tgt[c] = {"n": s[5], "half_width": hw, "met": hw <= TARGET, "n_needed": n_needed(s[2])}
            L.append("| %s | %d | %.3f | %s | %s |" % (c, s[5], hw, "yes" if hw <= TARGET else "no",
                                                        tgt[c]["n_needed"] or "> 200"))
        else:
            L.append("| %s | – | – | – | – |" % c)
    if P["beta_all"]:
        L.append("")
        L.append("With the pooled SD %.3f: n %s per corner." % (P["beta_all"][0], n_needed(P["beta_all"][0])))
    L += ["", "## Levers β(Ki) − β(K1) (SE from the two corner SEs)", "",
          "| lever | " + " | ".join(LABEL[k] for k in OUTPUTS if k.startswith("beta")) + " |",
          "|---|" + "---|" * sum(k.startswith("beta") for k in OUTPUTS)]
    lev = {}
    for c in CORNERS[1:]:
        cells = []
        for k in (k for k in OUTPUTS if k.startswith("beta")):
            x, y = S[k][c], S[k]["K1"]
            if x and y:
                lev.setdefault(c, {})[k] = (x[0] - y[0], math.hypot(x[1], y[1]))
                cells.append("%+.3f ± %.3f" % lev[c][k])
            else:
                cells.append("–")
        L.append("| %s − K1 | %s |" % (c, " | ".join(cells)))

    # --- Bloch points and time --------------------------------------------------------------------------------------
    L += ["", "## Bloch points of the preparation and time", "", "| corner | runs | free after the first pass | "
          "second | third | stages with flag B | preparation h (median) | cycles h (median) |",
          "|---|---|---|---|---|---|---|---|"]
    for c in CORNERS:
        rc = [rows[n] for n in runs if rows[n]["corner"] == c]
        if not rc:
            continue
        free_after = [sum(1 for r in rc if len(r["bp_after_pass"]) > p and r["bp_after_pass"][p] == 0
                          or (len(r["bp_after_pass"]) <= p and r["free"])) for p in range(3)]
        nB = sum(bool(x) for r in rc for x in r["flag_B"])
        L.append("| %s | %d | %s | %s | %s | %d of %d | %.2f | %.2f |" % (
            c, len(rc), *("%d of %d" % (f, len(rc)) for f in free_after), nB, sum(len(r["flag_B"]) for r in rc),
            float(np.median([r["prep_h"] for r in rc])), float(np.median([r["cycles_h"] for r in rc]))))

    # --- information: paired with the old arms --------------------------------------------------------------------
    L += ["", "## Information: Δβ against the runs of the same microstructure (k = 0 … 2)", "",
          "| run | arm | Δβ 9 … 150 mT | Δβ 9 … 50 mT |", "|---|---|---|---|"]
    pairs = {}
    for n in valid:
        r = rows[n]
        if r["k"] > 2:
            continue
        base = "%s_%d" % (r["corner"], r["k"])
        lw = EN.lnw(root / n)
        for arm in ("ad", "nx", "old"):
            d = EA.arm_dir(root, base, arm)
            if not (d / "DONE").exists():
                continue
            lo = EN.lnw(d)
            db = [EA.window_beta([lw, lo], w) for w in ("full", "low")]
            pairs.setdefault(n, {})[arm] = [b[0] - b[1] for b in db]
            L.append("| %s | %s | %s | %s |" % (n, arm, *("%+.3f" % x if np.isfinite(x) else "–"
                                                         for x in pairs[n][arm])))

    # --- information: energy-terms check ---------------------------------------------------------------------------
    L += ["", "## Information: energy-terms check (S5) on the calm-tail cycles", "",
          "| amplitude | cycles | median \\|ε\\| | median \\|d\\| | max \\|d\\| |", "|---|---|---|---|---|"]
    bal = {k: [] for k in range(len(AMPS))}
    for n in valid:
        for k, s in enumerate(EN.measured_stages(root / n)):
            tail = {c["cycle"] for c in EN.calm_tail(s)}
            bal[k] += [x for x in EN.balance(s) if x[0] in tail]
    for k, b in enumerate(AMPS):
        a = np.array([x[1:] for x in bal[k]]) if bal[k] else None
        L.append("| %g mT | %d | %s |" % (b, len(bal[k]), " | ".join(
            ("%.4f" % np.median(np.abs(a[:, 0])), "%.4f" % np.median(np.abs(a[:, 2])), "%.4f" % np.abs(a[:, 2]).max())
            if a is not None else ("–",) * 3)))
    print("\n".join(L))
    if json_out:
        pathlib.Path(json_out).write_text(json.dumps({
            "runs": rows, "valid": valid, "no_value": novalue, "corner": S, "pooled_sd": P, "target": tgt,
            "levers": lev, "pairs": pairs}, indent=1, default=lambda x: None))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    main(a.runs, a.json)
