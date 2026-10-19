#!/usr/bin/env python
"""
Regularity of the motion per frequency, particle size and B (question of Chris 2026-10-14: what happens at 60 MHz?).
Per run and stage, from the last 3 cycles (samples_*.csv, eval_jumps_time.series / measures): RN = non-repeating share
Σ|P − Pm| / ΣP (P the box dissipation, Pm its mean cycle), RF = non-smooth share Σ|P − Ps| / ΣP (Ps the smooth part of
Pm, Fourier terms ≤ 6), loop closure |ΔM| / M_peak (mean of the 3 cycles), ε span = largest minus smallest
w_dis / w_loop − 1 of the 3 cycles. Mean of the 3 seeds. 9 mT: all 8 frequencies; all B: 15 / 30 / 60 MHz (C′).

    python eval_sweep_regime.py > results/sweep9_regime.md
"""
import json
import pathlib

import numpy as np

import eval_jumps_time as EJ
import eval_sweep_state as ES
R = pathlib.Path(__file__).resolve().parent / "runs"
AMPS = (9, 20, 35, 50, 70, 100, 150)
def stage_files(d):
    return sorted(d.glob("samples_*_amp*.csv"))
def run_stats(d):
    summ = json.loads((d / "summary.json").read_text())["stages"]
    out = []
    for k, f in enumerate(stage_files(d)):
        s = EJ.series(f)
        if s is None:
            out.append(None); continue
        m = EJ.measures(*s)
        cyc = summ[k]["cycles"][-3:]
        eps = [c["w_dis"] / c["w_loop"] - 1 for c in cyc]
        out.append({"RN": m["RN"], "RF": m["RF"], "clo": float(np.mean([c["closure"] for c in cyc])),
                    "eps_span": float(max(eps) - min(eps))})
    return out
rows = {}
for a in (10, 20):
    for f in ES.F_MHZ:
        for k in ES.SEEDS:
            rows[(a, f, k)] = run_stats(ES.folder(a, f, k))
print("9 mT, mean of 3 seeds: RN (non-repeating share), RF (non-smooth share), closure, eps span of the last 3 cycles")
for a, lab in ((10, "large"), (20, "small")):
    for f in ES.F_MHZ:
        v = [rows[(a, f, k)][0] for k in ES.SEEDS if rows[(a, f, k)][0]]
        print("%-5s %3d MHz  RN %.3f  RF %.3f  closure %.4f  eps span %.3f   (n %d)" % (lab, f, np.mean([x["RN"] for x in v]),
              np.mean([x["RF"] for x in v]), np.mean([x["clo"] for x in v]), np.mean([x["eps_span"] for x in v]), len(v)))
print("\nC' runs, all B: RN per B (mean of 3 seeds)")
for a, lab in ((10, "large"), (20, "small")):
    for f in (15, 30, 60):
        rn = []
        for i in range(7):
            v = [rows[(a, f, k)][i]["RN"] for k in ES.SEEDS if len(rows[(a, f, k)]) > i and rows[(a, f, k)][i]]
            rn.append(np.mean(v) if v else np.nan)
        clo = []
        for i in range(7):
            v = [rows[(a, f, k)][i]["clo"] for k in ES.SEEDS if len(rows[(a, f, k)]) > i and rows[(a, f, k)][i]]
            clo.append(np.mean(v) if v else np.nan)
        print("%-5s %2d MHz RN: %s | closure: %s" % (lab, f, " ".join("%.3f" % x for x in rn), " ".join("%.4f" % x for x in clo)))
