#!/usr/bin/env python
"""
Corner study of PLAN 4.15, phase 1: comparison with the corner runs of the old main study (request of Chris
2026-10-23). Descriptive: no reading rule was fixed before this evaluation.

Arms per microstructure k = 0 … 2 of each corner (the same seeds, thus the same microstructure):
  old  the corner run of the old main study: DemagFieldPBC, atol 1e-5, relative stop threshold of the relaxation,
       relaxed random start (the _mc2 rerun where it exists, eval_corners_nx.old_dir)
  nx   setup of PLAN 4.13: new stray-field solver, atol 1e-7, absolute stop threshold; the same start protocol
  sat  setup of PLAN 4.15: as nx, plus the preparation by static saturation
Thus nx − old is the change of the code at the same start protocol, sat − nx the change of the preparation, and
sat − old the sum of both.

Per run, the methods of eval_corners_origin.run_record, restricted to the measures that the old accumulators have
(no harmonics, no wall measure for the _mc2 reruns):
  β, β_f, shift  per window of eval_corners_ad (β_f from ln (R1 P); shift = β − β_f)
  ln P           per amplitude (P in W/cm³)
  R1             share of the loss at the drive frequency f
  c              loss share / volume share of the core (r < 0.5), the rim (r >= 0.8), the 10 % of the cells with the
                 largest texture and the 10 % with the largest loss
  jump           share of the dissipation in jumps (eval_corners_fields.jumps)
  events         per stage, the event cycles from the first accumulated cycle on. The old accumulators include them,
                 the new ones (--acc_harm) do not (PLAN 4.13 S3). The value of a stage (calm tail) never includes them.
Paired differences per corner (n 3) and over all corners (n 12): mean ± SE, t test against 0. Unpaired: β of all sat
runs with a value against the 3 old runs per corner, and the levers Ki − K1 of both (Welch). The change of a lever,
(Ki − K1)_sat − (Ki − K1)_old, is a sum of 4 corner means: SE from the 4 SEs, degrees of freedom by Welch and
Satterthwaite.

    python eval_corners_compare.py [--runs DIR] [--json FILE] > report.md
"""
import argparse
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import eval_corners_ad as EA
import eval_corners_fields as EF
import eval_corners_nx as EN
import eval_corners_sat as ES
import eval_k1_origin as EO
import units as U

HERE = pathlib.Path(__file__).resolve().parent
CORNERS = list(ES.CORNERS)
AMPS = list(AD.AMPS_MT)
ARMS = ("old", "nx", "sat")
K_PAIRED = (0, 1, 2)
DIFFS = [("nx", "old", "code (nx − old)"), ("sat", "nx", "preparation (sat − nx)"), ("sat", "old", "total (sat − old)")]
EXPONENTS = [("beta", "β"), ("beta_f", "β_f"), ("shift", "shift")]
PER_AMP = [("R1", "R1, share at f", "%.2f"), ("c core", "c core (r < 0.5)", "%.2f"), ("c rim", "c rim (r >= 0.8)", "%.2f"),
           ("c texture top 10 %", "c texture top 10 %", "%.2f"), ("c loss top 10 %", "c loss top 10 %", "%.2f"),
           ("jump", "share of the dissipation in jumps", "%.3f")]


def arm_dir(root, c, k, arm):
    base = "%s_%d" % (c, k)
    return root / (base + "_sat") if arm == "sat" else EA.arm_dir(root, base, arm)


def event_cycles(d, acc):
    """Per stage: the event cycles (analyze.steady_tail criterion) from the first accumulated cycle on."""
    out = []
    for s, a in zip(EN.measured_stages(d), acc):
        cyc = [c for c in s["cycles"] if c["cycle"] >= a["cycles_from"]]
        out.append(sum(c["w_loop"] <= 0 or abs(c["w_dis"] / c["w_loop"] - 1.0) > AD.BAL_EVENT for c in cyc))
    return out + [None] * (len(AMPS) - len(out))


def pad(v):
    return list(v) + [float("nan")] * (len(AMPS) - len(v))


def record(d, lnF):
    lnw = EN.lnw(d)
    acc = json.loads((d / "acc.json").read_text())["stages"]
    rec = {"lnP": (lnw + lnF).tolist()}
    lnR1 = np.log(pad([s["R1"] for s in acc]))
    for win in EA.WINDOWS:
        b, bf = EA.window_beta([lnw, lnw + lnR1], win)
        rec["beta " + win], rec["beta_f " + win], rec["shift " + win] = b, bf, b - bf
    rec["R1"] = np.exp(lnR1).tolist()
    for key, name in EO.REGIONS:
        rec["c " + name] = pad([s["where"][key]["p_share"] / s["where"][key]["vol_share"] for s in acc])
    rec["jump"] = EF.jumps(d).tolist()
    rec["events"] = event_cycles(d, acc)
    return rec


def mse(v):
    """(mean, SE, n, p of the t test against 0) over the finite values; None for fewer than 2."""
    v = np.array([x for x in v if x is not None and np.isfinite(x)], float)
    if len(v) < 2:
        return None
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), len(v), float(stats.ttest_1samp(v, 0.0).pvalue)


def lincomb(groups, coefs):
    """Σ c_i mean_i over independent groups: (value, SE, Welch-Satterthwaite dof, two-sided p)."""
    g = [np.array([x for x in v if np.isfinite(x)], float) for v in groups]
    t = [c * c * x.var(ddof=1) / len(x) for c, x in zip(coefs, g)]
    val = float(sum(c * x.mean() for c, x in zip(coefs, g)))
    se = math.sqrt(sum(t))
    dof = sum(t) ** 2 / sum(ti * ti / (len(x) - 1) for ti, x in zip(t, g))
    return val, se, float(dof), float(2 * stats.t.sf(abs(val) / se, dof))


def cell(s, f="%.3f", sign=True):
    if not s:
        return "–"
    return "%s ± %s" % ((f.replace("%", "%+") if sign else f) % s[0], f % s[1])


def pcell(s):
    return "%.3f" % s[3] if s else "–"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    root = pathlib.Path(a.runs)
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    u = U.units()
    lnF = math.log(u["Kd"] * U.F_REL_DEFAULT * u["f_M"] / 1e6)

    P = {c: {k: {arm: record(arm_dir(root, c, k, arm), lnF) for arm in ARMS} for k in K_PAIRED} for c in CORNERS}
    sat = [n for n in ES.sat_runs(root) if ES.prep_record(root / n)[2]]
    S = {n: record(root / n, lnF) for n in sat}
    bp = {"%s_%d" % (c, k): EA.bp_stages(root, "%s_%d" % (c, k)) for c in CORNERS for k in K_PAIRED}
    bp_nx = {n: v.get("nx") for n, v in bp.items()}
    bp_old = {n: v.get("old") for n, v in bp.items()}
    bp_sat = {"%s_%d" % (c, k): [i for i, f in enumerate(ES.stage_flags_b(arm_dir(root, c, k, "sat"))) if f]
              for c in CORNERS for k in K_PAIRED}

    def diffs(key, i=None):
        """{(hi, lo): {corner: [per k], "all": [...]}} of the paired differences of one measure."""
        out = {}
        for hi, lo, _ in DIFFS:
            dd = {c: [(P[c][k][hi][key] if i is None else P[c][k][hi][key][i]) -
                      (P[c][k][lo][key] if i is None else P[c][k][lo][key][i]) for k in K_PAIRED] for c in CORNERS}
            dd["all"] = [x for c in CORNERS for x in dd[c]]
            out[hi + "-" + lo] = dd
        return out

    out = {"runs_paired": P, "runs_sat": S, "bp_old": bp_old, "bp_nx": bp_nx, "bp_sat": bp_sat, "amps_mT": AMPS,
           "paired": {},
           "unpaired": {}}
    L = ["# Corner study, phase 1: comparison with the old main study (eval_corners_compare.py)", "",
         "Descriptive; no reading rule was fixed before this evaluation. Paired arms per microstructure k = 0 … 2 (the "
         "same seeds): old (old main study), nx (new code, the same random start), sat (new code, static saturation). "
         "nx − old: change of the code; sat − nx: change of the preparation. Mean ± SE; p: t test against 0 "
         "(n 3 per corner: 2 degrees of freedom). P in W/cm³.", ""]

    # --- exponents ------------------------------------------------------------------------------------------------
    L += ["## Exponents (paired runs)", ""]
    for w in EA.WINDOWS:
        L += ["### Window %s" % EA.WIN_LABEL[w], "", "| corner | measure | old | nx | sat | " +
              " | ".join("%s | p" % lab for _, _, lab in DIFFS) + " |", "|---|---|---|---|---|" + "---|---|" * len(DIFFS)]
        for key, lab in EXPONENTS:
            name = "%s %s" % (key, w)
            D = diffs(name)
            out["paired"][name] = {d: {c: mse(v) for c, v in D[d].items()} for d in D}
            for c in CORNERS + ["all"]:
                arm = [mse([P[cc][k][x][name] for cc in ([c] if c != "all" else CORNERS) for k in K_PAIRED])
                       for x in ARMS]
                L.append("| %s | %s | %s | %s |" % (c, lab, " | ".join(cell(s, sign=False) for s in arm), " | ".join(
                    "%s | %s" % (cell(out["paired"][name][hi + "-" + lo][c]), pcell(out["paired"][name][hi + "-" + lo][c]))
                    for hi, lo, _ in DIFFS)))
        L.append("")

    # --- loss density -----------------------------------------------------------------------------------------------
    L += ["## Factor on P (exp of the mean paired Δ ln P; SE of Δ ln P in brackets)", "",
          "| corner | comparison | " + " | ".join("%g mT" % b for b in AMPS) + " |", "|---|---|" + "---|" * len(AMPS)]
    out["paired"]["lnP"] = []
    for i in range(len(AMPS)):
        D = diffs("lnP", i)
        out["paired"]["lnP"].append({d: {c: mse(v) for c, v in D[d].items()} for d in D})
    for c in CORNERS + ["all"]:
        for hi, lo, lab in DIFFS:
            s = [out["paired"]["lnP"][i][hi + "-" + lo][c] for i in range(len(AMPS))]
            L.append("| %s | %s | %s |" % (c, lab, " | ".join("× %.2f (%.2f)" % (math.exp(x[0]), x[1]) if x else "–"
                                                          for x in s)))
    L.append("")

    # --- channels, location, jumps ------------------------------------------------------------------------------------
    L += ["## Channel, location and jumps per amplitude (paired runs; mean over the runs, n 3 per corner)", ""]
    for key, lab, f in PER_AMP:
        L += ["### " + lab, "", "| corner | arm | " + " | ".join("%g mT" % b for b in AMPS) + " |",
              "|---|---|" + "---|" * len(AMPS)]
        for c in CORNERS:
            for x in ARMS:
                v = np.array([P[c][k][x][key] for k in K_PAIRED], float)
                L.append("| %s | %s | %s |" % (c, x, " | ".join(f % m if np.isfinite(m) else "–"
                                                              for m in np.nanmean(v, axis=0))))
        out["paired"][key] = []
        for i in range(len(AMPS)):
            D = diffs(key, i)
            out["paired"][key].append({d: {c: mse(v) for c, v in D[d].items()} for d in D})
        for hi, lo, labd in DIFFS:
            L.append("| all | %s | %s |" % (labd, " | ".join(cell(out["paired"][key][i][hi + "-" + lo]["all"], f)
                                                          for i in range(len(AMPS)))))
        L.append("")

    L += ["### Event cycles in the old accumulators (sum over the 3 old runs per corner)", "",
          "| corner | " + " | ".join("%g mT" % b for b in AMPS) + " |", "|---|" + "---|" * len(AMPS)]
    for c in CORNERS:
        L.append("| %s | %s |" % (c, " | ".join(str(sum(P[c][k]["old"]["events"][i] or 0 for k in K_PAIRED))
                                              for i in range(len(AMPS)))))
    L += ["", "Bloch points. old (count in the start state and the end state): %s. nx (stages 0 … 6 with a BP in the "
          "census at the start or the end of the stage): %s. sat (stages with flag B): %s." % (
              ", ".join("%s %s / %s" % (n, *v) for n, v in bp_old.items() if v and any(v)) or "none",
              ", ".join("%s %s" % (n, v) for n, v in bp_nx.items() if v) or "none",
              ", ".join("%s %s" % (n, v) for n, v in bp_sat.items() if v) or "none"), ""]

    # --- unpaired: all sat runs against the 3 old runs per corner ----------------------------------------------------
    L += ["## All sat runs against the 3 old runs per corner (unpaired)", "",
          "| corner | measure | old (n 3) | sat (all) | sat − old | dof | p |", "|---|---|---|---|---|---|---|"]
    G = {}
    for key, lab in EXPONENTS:
        name = key + " full"
        for c in CORNERS:
            G[(name, c, "old")] = [P[c][k]["old"][name] for k in K_PAIRED]
            G[(name, c, "sat")] = [S[n][name] for n in sat if n.startswith(c + "_")]
            r = lincomb([G[(name, c, "sat")], G[(name, c, "old")]], [1, -1])
            arm = [ES.cstat(G[(name, c, x)]) for x in ("old", "sat")]
            out["unpaired"]["%s %s" % (name, c)] = {"old": arm[0], "sat": arm[1], "diff": r}
            L.append("| %s | %s 9 … 150 mT | %s | %s (n %d) | %+.3f ± %.3f | %.1f | %.4f |" % (
                c, lab, cell(arm[0][:2], sign=False), cell(arm[1][:2], sign=False), arm[1][5], *r))
    L += ["", "### Levers Ki − K1, β 9 … 150 mT (Welch)", "",
          "| lever | old (n 3 + 3) | p | sat (all runs) | p | change sat − old | dof | p |", "|---|---|---|---|---|---|---|---|"]
    name = "beta full"
    for c in CORNERS[1:]:
        lo = lincomb([G[(name, c, "old")], G[(name, "K1", "old")]], [1, -1])
        ls = lincomb([G[(name, c, "sat")], G[(name, "K1", "sat")]], [1, -1])
        ch = lincomb([G[(name, c, "sat")], G[(name, "K1", "sat")], G[(name, c, "old")], G[(name, "K1", "old")]],
                     [1, -1, -1, 1])
        out["unpaired"]["lever %s" % c] = {"old": lo, "sat": ls, "change": ch}
        L.append("| %s − K1 | %+.3f ± %.3f | %.4f | %+.3f ± %.3f | %.4f | %+.3f ± %.3f | %.1f | %.4f |" % (
            c, lo[0], lo[1], lo[3], ls[0], ls[1], ls[3], *ch))
    L += ["", "### Change of the levers Ki − K1 from the paired runs, split (β 9 … 150 mT; n 3 per corner; Welch)", "",
          "| lever | " + " | ".join(lab for _, _, lab in DIFFS) + " |", "|---|" + "---|" * len(DIFFS)]
    D = diffs(name)
    for c in CORNERS[1:]:
        parts = {hi + "-" + lo: lincomb([D[hi + "-" + lo][c], D[hi + "-" + lo]["K1"]], [1, -1]) for hi, lo, _ in DIFFS}
        out["unpaired"]["lever parts %s" % c] = parts
        L.append("| %s − K1 | %s |" % (c, " | ".join("%+.3f ± %.3f (p %.3f)" % (parts[k][0], parts[k][1], parts[k][3])
                                                     for k in parts)))
    L.append("")

    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1, default=float))
    print("\n".join(L))


if __name__ == "__main__":
    main()
