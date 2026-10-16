#!/usr/bin/env python
"""
Energy balance on dx 0.75 with the reference solver (decision Chris 2026-10-16, PROTOKOLL 7.108; rules fixed before the
runs). In the stage 2 F runs of PLAN 4.11 the steady tail has ε = w_dis / w_loop − 1 of +0.170 (s941), +0.091 (s942),
+0.054 (s943), +0.032 (s944), although the loops close (closure 0.0000); on dx 1.5 ε is −0.006 … +0.028, on dx 3
+0.000 … +0.003. In a closed periodic loop ε must be 0. Both energies are trapezoid sums over the samples of a cycle. Test: the F runs of s941 and s942 again from their relaxed state with 1024 samples per cycle
instead of 256, all else equal (conv_jobs.py --newell eps).
Per seed, over the measured steady tail (analyze rules, ≥ 2 calm cycles): ε_256, ε_1024 (mean ε) and
Δ = ln w_loop(1024) − ln w_loop(256) (means of the tail).
Reading per seed: "sampling" if |ε_1024| ≤ 0.03 and |ε_1024| ≤ |ε_256| / 3 (the sampled energy integrals carry the
error); "not sampling" if |ε_1024| ≥ 2/3 |ε_256| (the error is elsewhere: time integration or the model); else "partly".
Overall: the class of both seeds if they agree, else "mixed". Loss value: "w_loop unaffected" if |Δ| ≤ 0.02 in both
seeds, else the Δ values are the size of the sampling error of the stage 2 F values.

    python eval_eps075.py > runs/eps075_report.md
"""
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD

R = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = (941, 942)


def tail(d):
    """(mean ε, mean w_loop, n cycles) over the steady tail of the measured stage."""
    st = [x for x in json.loads((d / "summary.json").read_text())["stages"] if x.get("measure")][-1]
    cyc = analyze.steady_tail(st["cycles"][analyze.first_measured(st):], AD.BAL_EVENT)
    if len(cyc) < 2:
        return float("nan"), float("nan"), len(cyc)
    return (float(np.mean([c["w_dis"] / c["w_loop"] - 1 for c in cyc])), float(np.mean([c["w_loop"] for c in cyc])),
            len(cyc))


def klass(e256, e1024):
    if not (np.isfinite(e256) and np.isfinite(e1024)):
        return "not determined"
    if abs(e1024) <= 0.03 and abs(e1024) <= abs(e256) / 3:
        return "sampling"
    if abs(e1024) >= 2 * abs(e256) / 3:
        return "not sampling"
    return "partly"


def main():
    L = ["# Energy balance on dx 0.75 with 1024 samples per cycle (eval_eps075.py, PROTOKOLL 7.108)", "",
         "| seed | ε, 256 samples | ε, 1024 samples | ln w_loop(1024) − ln w_loop(256) | tail cycles | class |",
         "|---|---|---|---|---|---|"]
    cl, dl = [], []
    for s in SEEDS:
        a, b = tail(R / ("N0b_F_s%d" % s)), tail(R / ("N0b_F_s%d_smp1024" % s))
        d = math.log(b[1] / a[1]) if np.isfinite(a[1]) and np.isfinite(b[1]) else float("nan")
        c = klass(a[0], b[0])
        cl.append(c)
        dl.append(d)
        L.append("| %d | %+.3f | %+.3f | %+.4f | %d / %d | %s |" % (s, a[0], b[0], d, a[2], b[2], c))
    overall = cl[0] if len(set(cl)) == 1 else "mixed"
    loss = ("w_loop unaffected" if all(np.isfinite(x) and abs(x) <= 0.02 for x in dl)
            else "w_loop changes by up to %.3f in ln" % max(abs(x) for x in dl if np.isfinite(x)) if any(np.isfinite(dl)) else
            "not determined")
    L += ["", "**Reading: %s; loss value: %s.**" % (overall, loss)]
    print("\n".join(L))


if __name__ == "__main__":
    main()
