#!/usr/bin/env python
"""
Check of the measuring method of B0b (go of Chris 2026-10-15, PROTOKOLL 7.102; exploratory, after the pre-registered
reading of eval_conv.py). B0b: base microstructure, d/l_ex 96, 9 mT, seeds 941 … 944, atol 1e-6; arms F (dx 0.75),
I15 / I3 (dx 1.5 / 3, independent relaxation), T15 / T3 (dx 1.5 / 3, start from the relaxed F state).

Per run: number of cycles, event cycles (7.53), steady tail, P per tail cycle (W/cm³, f K_d w_loop), its cycle SD,
w_dis / w_loop and closure in the tail, RN and RF of the last 3 cycles (eval_jumps_time.measures: share of the
dissipation that changes from cycle to cycle / that is not smooth), relaxed energy E / K_d, relaxation iterations,
neighbour pairs above 60° in the relaxed state. Per arm: geometric mean of P, SD of ln P over the seeds, median cycle
SD. Reproducibility floor: F s941 and s942 repeated their stage after the restart of the instance (2026-10-14 18:47 UTC)
from init.pt; stdout.log holds the cycle lines of both attempts.

    python eval_b0b_check.py > results/b0b_check.md      # also results/b0b_check.json
"""
import json
import math
import pathlib
import re

import numpy as np

import analyze
import analyze_design as AD
import eval_jumps_time as EJ

HERE = pathlib.Path(__file__).resolve().parent
R = HERE / "runs"
ARMS = ("F", "I15", "I3", "T15", "T3")
SEEDS = (941, 942, 943, 944)
K_D = 1.5 ** 2 / (2 * 4e-7 * math.pi) / 1e6          # J/cm³
F_HZ = 30e6
LINE = re.compile(r"^\[amp00 c(\d+)\].*w_loop=([0-9.eE+-]+)\s+w_dis=([0-9.eE+-]+)\s+closure=([0-9.eE+-]+)")


def attempts(d):
    """Cycle lines {cycle: (w_loop, w_dis, closure)} per start of run_loops in stdout.log."""
    out = []
    for line in (d / "stdout.log").read_text().splitlines():
        if "magnum.np 2.2.0" in line:
            out.append({})
        m = LINE.match(line)
        if m and out:
            out[-1][int(m.group(1))] = tuple(float(x) for x in m.groups()[1:])
    return out


def run(d):
    s = json.loads((d / "summary.json").read_text())
    st = [x for x in s["stages"] if x.get("measure")][0]
    meas = st["cycles"][2:] if "h_amp_initial" in st else st["cycles"][1:]
    tail = analyze.steady_tail(meas, AD.BAL_EVENT)
    ev = [c["cycle"] for c in meas if c["w_loop"] <= 0 or abs(c["w_dis"] / c["w_loop"] - 1) > AD.BAL_EVENT]
    P = np.array([F_HZ * K_D * c["w_loop"] for c in tail])
    m = EJ.measures(*EJ.series(sorted(d.glob("samples_*_amp*.csv"))[0]))
    i = s.get("init", {})
    return {"cycles": len(st["cycles"]), "events": ev, "tail": len(tail), "P": float(np.exp(np.log(P).mean())),
            "P_cycles": P.tolist(), "cycle_sd": float(P.std(ddof=1) / P.mean()),
            "wd_w": [c["w_dis"] / c["w_loop"] for c in tail], "closure_max": max(c["closure"] for c in tail),
            "RN": float(m["RN"]), "RF": float(m["RF"]), "E": i.get("E_over_Kd"), "relax_it": i.get("relax_iterations"),
            "pairs60": i.get("n_pairs_gt60")}


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    res = {"runs": {}, "arms": {}, "repro": {}}
    L = ["# B0b: check of the measuring method (exploratory, after the pre-registered reading)", "",
         "| run | cycles | events | tail | P (W/cm³) per tail cycle | cycle SD | w_dis / w_loop (tail) | largest closure "
         "(tail) | RN | RF | E / K_d relaxed | relaxation iterations | pairs > 60° |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for a in ARMS:
        for s in SEEDS:
            d = R / ("B0b_%s_s%d" % (a, s))
            if not (d / "DONE").exists():
                L.append("| %s | missing |" % d.name)
                continue
            x = res["runs"][d.name] = run(d)
            L.append("| %s | %d | %s | %d | %s | %.1f %% | %s | %.0e | %.3f | %.3f | %.6f | %s | %s |" % (
                d.name, x["cycles"], x["events"] or "–", x["tail"], " / ".join("%.2f" % p for p in x["P_cycles"]),
                100 * x["cycle_sd"], " / ".join("%.3f" % v for v in x["wd_w"]), x["closure_max"], x["RN"], x["RF"],
                x["E"], x["relax_it"], x["pairs60"]))
    L += ["", "| arm | runs | P geometric mean (W/cm³) | SD of ln P over the seeds | median cycle SD | median RN | ln P per "
          "seed (941 / 942 / 943 / 944) |", "|---|---|---|---|---|---|---|"]
    for a in ARMS:
        xs = {s: res["runs"].get("B0b_%s_s%d" % (a, s)) for s in SEEDS}
        v = [x for x in xs.values() if x]
        lp = np.log([x["P"] for x in v])
        res["arms"][a] = {"n": len(v), "P_gm": float(np.exp(lp.mean())), "sd_lnP": float(lp.std(ddof=1)) if len(v) > 1
                          else float("nan"), "lnP": {s: (math.log(x["P"]) if x else None) for s, x in xs.items()}}
        L.append("| %s | %d | %.2f | %.3f | %.1f %% | %.3f | %s |" % (
            a, len(v), np.exp(lp.mean()), res["arms"][a]["sd_lnP"], 100 * np.median([x["cycle_sd"] for x in v]),
            np.median([x["RN"] for x in v]), " / ".join("%.3f" % math.log(x["P"]) if x else "–" for x in xs.values())))
    L += ["", "Reproducibility floor: F s941 and s942 after the restart (from init.pt) against the first attempt (the first "
          "attempt relaxed in the same process; 7.59: another start path gives differences at the solver tolerance). "
          "stdout.log prints w_loop and w_dis to 4 digits and the closure to 2 digits.", "",
          "| run | cycles of the first attempt | largest abs. w_loop ratio − 1 | largest abs. w_dis ratio − 1 | "
          "largest abs. closure ratio − 1 |", "|---|---|---|---|---|"]
    for s in (941, 942):
        at = attempts(R / ("B0b_F_s%d" % s))
        if len(at) < 2:
            continue
        a0, a1 = at[0], at[-1]
        com = sorted(set(a0) & set(a1))
        r = [max(abs(a1[c][j] / a0[c][j] - 1) for c in com) for j in range(3)]
        res["repro"][s] = {"cycles": com, "w_loop": r[0], "w_dis": r[1], "closure": r[2]}
        L.append("| B0b_F_s%d | %s | %.1e | %.1e | %.1e |" % (s, com, *r))
    print("\n".join(L))
    (HERE / "results" / "b0b_check.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
