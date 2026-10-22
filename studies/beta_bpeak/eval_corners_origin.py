#!/usr/bin/env python
"""
Corner study of PLAN 4.15, phase 1: where and in which form the loss occurs at the 4 corners, and which of these
measures changes in the same field window as β (Chris 2026-10-22, options 1 and 2; descriptive, it names no
mechanism). Methods of the main study, per run K<c>_<k>_sat with a value (D7 (b)):
  β, β_f, shift  β: slope of ln P over ln B in the windows 9 … 150, 9 … 50 and 70 … 150 mT
                 (eval_corners_ad.window_beta); β_f: the same slope of ln (R1 P), the loss at the drive frequency f;
                 shift = β − β_f = −d ln R1 / d ln B (identity of PROTOKOLL 7.71)
  channels       R1 = share at f, R_harm = share at 2f … 7f, R_rest = the remainder (acc.json, calm cycles; 7.70, 7.72)
  location       concentration c = loss share / volume share of the core (r < 0.5), the rim (r >= 0.8), the 10 % of
                 the cells with the largest texture and the 10 % with the largest loss (acc.json; 7.69)
  jumps          share of the dissipation in jumps (eval_corners_fields.jumps; 7.119)
  walls          at 9, 50 and 150 mT (the stages with snapshots): wall area per box volume a_w = A_w / V, the
                 wall paths s_0 and s_1 of the 2 half cycles, ln (X / V) with X = A_w s_mean² (domains.json;
                 eval_k1_origin.run_record; 7.125). Per box volume, because the boxes of the corners differ in size
                 and P is per box volume. The wall measure holds within ± 12 % only if the walls are at least
                 5 L_eff apart (7.14); this is not checked here
Per corner: mean ± SE over the runs, n. Differences Ki − K1 with the SE from the two corner SEs. All K1 runs count
(Chris 2026-10-22: K1_4 and K1_8 are valid realisations of K1).
The loss of each particle is in eval_k1_particles.py (only the runs with the key "particles" in acc.json).

    python eval_corners_origin.py [--runs DIR] [--json FILE] > report.md
"""
import argparse
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import eval_corners_ad as EA
import eval_corners_fields as EF
import eval_corners_nx as EN
import eval_corners_sat as ES
import eval_k1_origin as EO
import units as U

HERE = pathlib.Path(__file__).resolve().parent
CORNERS = ES.CORNERS
AMPS = list(AD.AMPS_MT)
WALL_STAGES = {0: "9 mT", 3: "50 mT", 6: "150 mT"}
PER_AMP = [("lnP", "ln P (P in W/cm³)"), ("R1", "R1, share at f"), ("R_harm", "R_harm, share at 2f … 7f"),
           ("R_rest", "R_rest, share outside f … 7f"), ("c core", "c core (r < 0.5)"), ("c rim", "c rim (r >= 0.8)"),
           ("c texture top 10 %", "c texture top 10 %"), ("c loss top 10 %", "c loss top 10 %"),
           ("jump", "share of the dissipation in jumps")]
WALLS = [("a_w", "a_w = A_w / V, 10⁻³ / l_ex", "%.3f"), ("s0", "s_0 (H max <-> H min), l_ex", "%.2f"),
         ("s1", "s_1 (H = 0 falling <-> rising), l_ex", "%.2f"), ("lnXv", "ln (X / V), X / V in l_ex", "%.2f")]
EXPONENTS = [("beta", "β"), ("beta_f", "β_f"), ("shift", "shift")]


def run_record(d, lnF):
    lnw = EN.lnw(d)
    acc = json.loads((d / "acc.json").read_text())["stages"]
    cfg = json.loads((d / "config.json").read_text())
    V = (cfg["N_used"] * cfg["dx_lex"]) ** 3                 # box volume, l_ex³ (P is per box volume)
    rec = {"lnP": (lnw + lnF).tolist()}
    lnR1 = np.log([s["R1"] for s in acc])
    for win in EA.WINDOWS:
        b, bf = EA.window_beta([lnw, lnw + lnR1], win)
        rec["beta " + win], rec["beta_f " + win], rec["shift " + win] = b, bf, b - bf
    for i in range(len(AMPS)):
        o = EO.run_record(d, i) if i in WALL_STAGES else None
        if o is None:
            s = acc[i]
            o = {"R1": s["harm"]["R_k"][0], "R_harm": float(sum(s["harm"]["R_k"][1:])), "R_rest": s["harm"]["R_rest"]}
            for key, name in EO.REGIONS:
                o["c " + name] = s["where"][key]["p_share"] / s["where"][key]["vol_share"]
        for k in ("R1", "R_harm", "R_rest") + tuple("c " + name for _, name in EO.REGIONS):
            rec.setdefault(k, []).append(o[k])
        if i in WALL_STAGES:
            o["a_w"], o["lnXv"] = 1e3 * o["A_w"] / V, o["lnX"] - math.log(V)
            for k, _, _ in WALLS:
                rec.setdefault(k, []).append(o[k])
    rec["jump"] = EF.jumps(d).tolist()
    return rec


def stat(v):
    v = np.array([x for x in v if np.isfinite(x)])
    return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), len(v)) if len(v) > 1 else None


def cell(s, f, sign=False):
    return ("%s ± %s" % ((f.replace("%", "%+") if sign else f) % s[0], f % s[1])) if s else "–"


def diff(a, b):
    return (a[0] - b[0], math.sqrt(a[1] ** 2 + b[1] ** 2), 0) if a and b else None


def table(rows, cols, title, f):
    """rows: {corner: [values of the runs per column]} -> markdown lines with mean ± SE and Ki − K1."""
    S = {c: [stat(v) for v in rows[c]] for c in CORNERS}
    L = ["### " + title, "", "| corner | n | " + " | ".join(cols) + " |", "|---|---|" + "---|" * len(cols)]
    for c in CORNERS:
        n = max((s[2] for s in S[c] if s), default=0)
        L.append("| %s | %d | %s |" % (c, n, " | ".join(cell(s, f) for s in S[c])))
    for c in CORNERS[1:]:
        L.append("| %s − K1 | | %s |" % (c, " | ".join(cell(diff(a, b), f, True) for a, b in zip(S[c], S["K1"]))))
    return L + [""], S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    root = pathlib.Path(a.runs)
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    u = U.units()
    lnF = math.log(u["Kd"] * U.F_REL_DEFAULT * u["f_M"] / 1e6)
    runs = [n for n in ES.sat_runs(root) if (root / n / "acc.json").exists() and (root / n / "domains.json").exists()]
    novalue = [n for n in runs if not ES.prep_record(root / n)[2]]
    runs = [n for n in runs if n not in novalue]
    R = {n: run_record(root / n, lnF) for n in runs}
    by = {c: [n for n in runs if n.startswith(c + "_")] for c in CORNERS}
    L = ["# Corner study, phase 1: where and in which form the loss occurs (eval_corners_origin.py)", "",
         "Runs: %s. Without a value (D7 (b)): %s. Mean ± SE over the runs of a corner; Ki − K1 with the SE from the two "
         "corner SEs. c = loss share / volume share of the region (1: the mean density). a_w = A_w / V in 1/l_ex "
         "(V: box volume), s in l_ex, X = A_w s_mean². The texture is smoothed with sigma L_eff / 2 (acc_reduce.py); "
         "L_eff differs between the corners." % (", ".join("%s n %d" % (c, len(by[c])) for c in CORNERS), ", ".join(novalue) or "none"),
         ""]
    out = {"runs": R, "corners": {}}
    cols = [EA.WIN_LABEL[w] for w in EA.WINDOWS]
    L += ["## Exponents per window: β = β_f + shift", ""]
    for k, lab in EXPONENTS:
        t, S = table({c: [[R[n]["%s %s" % (k, w)] for n in by[c]] for w in EA.WINDOWS] for c in CORNERS}, cols, lab,
                     "%.3f")
        L += t
        out["corners"][k] = S
    L += ["## Channels, location and jumps per amplitude", ""]
    amp_cols = ["%g mT" % b for b in AMPS]
    for k, lab in PER_AMP:
        t, S = table({c: [[R[n][k][i] for n in by[c]] for i in range(len(AMPS))] for c in CORNERS}, amp_cols, lab,
                     "%.3f")
        L += t
        out["corners"][k] = S
    L += ["## Walls at 9, 50 and 150 mT", ""]
    for k, lab, f in WALLS:
        t, S = table({c: [[R[n][k][j] for n in by[c]] for j in range(len(WALL_STAGES))] for c in CORNERS},
                     list(WALL_STAGES.values()), lab, f)
        L += t
        out["corners"][k] = S
    print("\n".join(L))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
