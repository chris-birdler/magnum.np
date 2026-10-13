#!/usr/bin/env python
"""
B0, step 1 of the mesh convergence series (PLAN 4.10). Rules fixed before the runs (2026-10-13).

Runs (conv_jobs.py): base microstructure at d/l_ex 96, 9 mT, seeds 941 … 944; per seed F (dx 0.75, virgin
relaxation), T15 / T3 (dx 1.5 / 3, start from the relaxed F state, then the production relaxation), I15 / I3 (dx 1.5 / 3,
independent virgin relaxation from the same seed).
Value of a run: ln P at 9 mT, P = loss density (W/cm³) from the mean w_loop of the calm tail (steady rule: ≥ 2 calm
cycles, calm = w_loop > 0 and |w_dis / w_loop − 1| ≤ 0.20). A run without a calm tail, without DONE or with a not
converged relaxation has no value and is listed.
Statistics (n = seeds with a value in both arms): paired difference per seed, mean ± SE, 90 % CI with t(n − 1); the SD
of the paired difference with its 90 % CI (χ², n − 1). Unpaired (Welch) differences of the arm means for the I arms.
Cycle noise: pooled SD of ln w_loop over the calm cycles of all runs (the floor of a paired SD).

Pre-registered reading (go / no-go for variant B; SD_T = SD of ln P(T15) − ln P(F), SD_I = SD of ln P(I15) − ln P(F)):
  not determined   fewer than 3 seeds with a value in F and T15;
  stop             SD_T ≥ 0.10: the transferred state does not remove the scatter → stop and re-scope;
  go, transfer works     SD_T < 0.10 and SD_I ≥ 0.15: the transferred state removes the state scatter → variant B
                   (4 amplitudes, separate go);
  go, scatter small      SD_T < 0.10 and SD_I < 0.15: the scatter is small at this point with either start → variant B;
                   the benefit of the transfer is not shown;
  go, benefit not determined   SD_T < 0.10 and fewer than 3 seeds with a value in F and I15 → variant B; the benefit
                   of the transfer is not determined.
ln P is taken at the measured peak flux density b_peak of the calm tail, without a β correction to 9 mT; the report
lists b_peak / b_target per run.
Information only (9 mT, one amplitude, no β): mean Δ(1.5 − 0.75), Δ(3 − 1.5), Δ(3 − 0.75) of the T arm with 90 % CI,
and whether the CI of Δ(1.5 − 0.75) lies inside ± 0.10 (the 4.10 criterion for "converged at dx 1.5", here at 9 mT only).

    python eval_conv.py > results/conv_b0.md      (on the instance, chain_conv.sh writes runs/conv_b0_report.md)
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import conv_jobs as CJ

R = pathlib.Path(__file__).resolve().parent / "runs"
MIN_N = 3
SD_GO = 0.10
SD_OLD = 0.15                       # lower end of the paired SD of the earlier mesh checks (0.15 … 0.46, PROTOKOLL 7.78)
LIMIT = 0.10
ARM_LAB = {"F": "F dx 0.75", "T15": "T dx 1.5", "T3": "T dx 3", "I15": "I dx 1.5", "I3": "I dx 3"}
PAIRS = (("T15", "F", "dx 1.5 − dx 0.75, transferred"), ("T3", "T15", "dx 3 − dx 1.5, transferred"),
         ("T3", "F", "dx 3 − dx 0.75, transferred"), ("I15", "F", "dx 1.5 − dx 0.75, independent (same seed)"),
         ("I3", "I15", "dx 3 − dx 1.5, independent (same seed)"), ("T15", "I15", "dx 1.5: transferred − independent"),
         ("T3", "I3", "dx 3: transferred − independent"))


def run_record(d):
    """{lnP, status, s_cyc, dof, iters, transfer, n_cycles, b_ratio} of one run folder. The steady rule (≥ 2 calm
    cycles) is set here for the call of analyze.load_run and restored afterwards."""
    rec = {"lnP": float("nan"), "status": "not run", "s_cyc": float("nan"), "dof": 0, "iters": None,
           "transfer": None, "n_cycles": 0, "b_ratio": float("nan")}
    if not d.exists():
        return rec
    if not (d / "summary.json").exists():
        rec["status"] = "started, no summary (refused or crashed)"
        return rec
    init = json.loads((d / "summary.json").read_text()).get("init", {})
    rec["iters"], rec["transfer"] = init.get("relax_iterations"), init.get("transfer")
    if not (d / "DONE").exists():
        rec["status"] = ("relaxation not converged (gate)" if init and not init.get("relax_converged", True)
                         else "not done")
        return rec
    old = analyze.STEADY_TAIL, analyze.STEADY_MIN
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    try:
        cfg, rows = analyze.load_run(d)
    except analyze.UnconvergedRun:
        rec["status"] = "relaxation not converged"
        return rec
    finally:
        analyze.STEADY_TAIL, analyze.STEADY_MIN = old
    rec["s_cyc"], rec["dof"] = cfg["s_cyc"], cfg["s_cyc_dof"]
    rec["n_cycles"] = rows[0]["n_cycles"] if rows else 0
    if rows and cfg.get("b_list"):
        rec["b_ratio"] = rows[0]["b_peak"] / float(cfg["b_list"][0])
    if not rows or not np.isfinite(rows[0]["w_loop"]) or rows[0]["w_loop"] <= 0:
        rec["status"] = "no calm tail (< 2 calm cycles)"
        return rec
    u = cfg["units"]
    rec["lnP"] = math.log(rows[0]["w_loop"] * u["Kd"] * cfg["f_rel"] * u["f_M"] / 1e6)
    rec["status"] = "settled"
    return rec


def collect(root=R):
    return {(arm, s): run_record(root / CJ.name(arm, s)) for arm, _, _ in CJ.ARMS for s in CJ.SEEDS}


def paired(lnP, a, b, seeds):
    """Statistics of ln P(a) − ln P(b) over the seeds with a value in both arms; None for fewer than MIN_N."""
    v = np.array([lnP[(a, s)] - lnP[(b, s)] for s in seeds
                  if np.isfinite(lnP.get((a, s), np.nan)) and np.isfinite(lnP.get((b, s), np.nan))])
    n = len(v)
    if n < MIN_N:
        return {"n": n}
    sd = v.std(ddof=1)
    se = sd / math.sqrt(n)
    h = stats.t.ppf(0.95, n - 1) * se
    sd_lo = math.sqrt((n - 1) * sd ** 2 / stats.chi2.ppf(0.95, n - 1))
    sd_hi = math.sqrt((n - 1) * sd ** 2 / stats.chi2.ppf(0.05, n - 1))
    return {"n": n, "mean": v.mean(), "se": se, "ci": (v.mean() - h, v.mean() + h), "sd": sd, "sd_ci": (sd_lo, sd_hi)}


def welch(lnP, a, b, seeds):
    """Mean ln P(a) − mean ln P(b) over all values of each arm (unpaired), SE, 90 % CI (Welch dof); None if < MIN_N."""
    x = np.array([lnP[(a, s)] for s in seeds if np.isfinite(lnP.get((a, s), np.nan))])
    y = np.array([lnP[(b, s)] for s in seeds if np.isfinite(lnP.get((b, s), np.nan))])
    if len(x) < MIN_N or len(y) < MIN_N:
        return None
    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y)
    se = math.sqrt(vx + vy)
    dof = (vx + vy) ** 2 / (vx ** 2 / (len(x) - 1) + vy ** 2 / (len(y) - 1)) if se > 0 else 1.0
    h = stats.t.ppf(0.95, dof) * se
    return {"mean": x.mean() - y.mean(), "se": se, "ci": (x.mean() - y.mean() - h, x.mean() - y.mean() + h),
            "nx": len(x), "ny": len(y)}


def reading(lnP, seeds):
    """(key, text) of the pre-registered reading."""
    t, i = paired(lnP, "T15", "F", seeds), paired(lnP, "I15", "F", seeds)
    if t["n"] < MIN_N:
        return "undetermined", "not determined: %d seeds with a value in F and T15 (need %d)" % (t["n"], MIN_N)
    if t["sd"] >= SD_GO:
        return "stop", ("stop: SD_T %.3f ≥ %.2f, the transferred state does not remove the scatter → stop and "
                        "re-scope" % (t["sd"], SD_GO))
    if i["n"] < MIN_N:
        return "go_undetermined", ("go: SD_T %.3f < %.2f, but only %d seeds with a value in F and I15: the benefit of "
                                   "the transfer is not determined → variant B (separate go)" % (t["sd"], SD_GO, i["n"]))
    if i["sd"] >= SD_OLD:
        return "go_transfer", ("go: SD_T %.3f < %.2f and SD_I %.3f ≥ %.2f, the transferred state removes the state "
                               "scatter → variant B (separate go)" % (t["sd"], SD_GO, i["sd"], SD_OLD))
    return "go_small", ("go: SD_T %.3f < %.2f and SD_I %.3f < %.2f: the scatter is small at this point with either "
                        "start → variant B (separate go); the benefit of the transfer is not shown" %
                        (t["sd"], SD_GO, i["sd"], SD_OLD))


def cycle_noise(recs):
    ss = sum(r["s_cyc"] ** 2 * r["dof"] for r in recs.values() if r["dof"] > 0 and np.isfinite(r["s_cyc"]))
    dof = sum(r["dof"] for r in recs.values() if r["dof"] > 0 and np.isfinite(r["s_cyc"]))
    return (math.sqrt(ss / dof), dof) if dof > 0 else (float("nan"), 0)


def fmt_p(p):
    if "mean" not in p:
        return "n %d: no statement" % p["n"]
    inside = "inside" if -LIMIT < p["ci"][0] and p["ci"][1] < LIMIT else "not inside"
    return ("%+.3f ± %.3f (90 %% CI %+.3f … %+.3f, %s ± %.2f) | %.3f (90 %% CI %.3f … %.3f) | %d" %
            (p["mean"], p["se"], p["ci"][0], p["ci"][1], inside, LIMIT, p["sd"], p["sd_ci"][0], p["sd_ci"][1], p["n"]))


def report(recs, seeds=CJ.SEEDS):
    lnP = {k: r["lnP"] for k, r in recs.items()}
    L = ["# B0: mesh convergence at 9 mT with a transferred state (PLAN 4.10, rules fixed before the runs)", "",
         "Base microstructure (L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1) at d/l_ex 96, 9 mT, seeds %s. F: dx 0.75 l_ex. "
         "T: start from the relaxed F state (block average), then relaxed again. I: independent virgin relaxation from "
         "the same seed. Value: ln P (P in W/cm³) of the calm tail (≥ 2 calm cycles). Errors: SE over the seeds, "
         "90 %% CI with t(n − 1); SD of the paired difference with its 90 %% CI (χ²)." %
         ", ".join(map(str, seeds)), ""]
    L += ["## Runs", "", "| arm | settled / runs | P (W/cm³), geometric mean | SD of ln P over the seeds | relaxation "
          "iterations (mean, range) | b_peak / b_target (range) | missing |", "|---|---|---|---|---|---|---|"]
    for arm, _, _ in CJ.ARMS:
        r = [recs[(arm, s)] for s in seeds]
        v = np.array([x["lnP"] for x in r if np.isfinite(x["lnP"])])
        it = [x["iters"] for x in r if x["iters"] is not None]
        miss = ["s%d: %s" % (s, recs[(arm, s)]["status"]) for s in seeds if recs[(arm, s)]["status"] != "settled"]
        br = [x["b_ratio"] for x in r if np.isfinite(x.get("b_ratio", np.nan))]
        L.append("| %s | %d / %d | %s | %s | %s | %s | %s |" % (
            ARM_LAB[arm], len(v), len(r), "%.3g" % math.exp(v.mean()) if len(v) else "–",
            "%.3f" % v.std(ddof=1) if len(v) > 1 else "–",
            "%.0f (%d … %d)" % (np.mean(it), min(it), max(it)) if it else "–",
            "%.4f … %.4f" % (min(br), max(br)) if br else "–", "; ".join(miss) or "none"))
    tr = [recs[(a, s)]["transfer"] for a in ("T15", "T3") for s in seeds if recs[(a, s)]["transfer"]]
    if tr:
        L += ["", "Transfer: %d runs, coarse particle cells without a fine sub-cell of the same particle: %d (max), zero "
              "vectors: %d (max)." % (len(tr), max(t["no_same_index_subcell"] for t in tr),
                                      max(t["zero_vector"] for t in tr))]
    sc, dof = cycle_noise(recs)
    L += ["", "Cycle noise (pooled SD of ln w_loop over the calm cycles): %.3f (dof %d)." % (sc, dof), ""]
    L += ["## Paired differences of ln P (same seed)", "",
          "| difference | mean ± SE (90 % CI) | SD of the difference (90 % CI) | n |", "|---|---|---|---|"]
    for a, b, lab in PAIRS:
        L.append("| %s | %s |" % (lab, fmt_p(paired(lnP, a, b, seeds))))
    L += ["", "## Unpaired differences of the arm means (Welch)", "", "| difference | mean ± SE (90 % CI) | n |",
          "|---|---|---|"]
    for a, b in (("I15", "F"), ("I3", "I15"), ("I3", "F")):
        w = welch(lnP, a, b, seeds)
        L.append("| %s − %s | %s |" % (ARM_LAB[a], ARM_LAB[b], "%+.3f ± %.3f (%+.3f … %+.3f) | %d / %d" % (
            w["mean"], w["se"], w["ci"][0], w["ci"][1], w["nx"], w["ny"]) if w else "too few values | –"))
    key, text = reading(lnP, seeds)
    L += ["", "## Reading (pre-registered)", "", "**%s**" % text, "",
          "Information only (9 mT, one amplitude, no β): Δ(1.5 − 0.75) of the T arm inside ± %.2f means dx 1.5 agrees "
          "with dx 0.75 at 9 mT within 10 %%; it is not the 4.10 criterion for β." % LIMIT]
    return "\n".join(L), key


def main():
    text, _ = report(collect())
    print(text)


if __name__ == "__main__":
    main()
