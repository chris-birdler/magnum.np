#!/usr/bin/env python
"""
Energy balance against the number of samples per cycle, from the existing sample files (decision Chris 2026-10-16,
PROTOKOLL 7.108; no new simulation). Each measured cycle has 257 rows (start + 256 samples). The script computes
W_loop = trapz(H, B) and W_dis = trapz(p_dis, t) with all rows (256 samples), every 2nd row (128) and every 4th row
(64), over the steady tail of the run (analyze rules). Check of the method: with 256 samples the script must give the
values of summary.json (the CSV has 10 digits).
Values per run: ε_N = mean(W_dis / W_loop) − 1 over the tail, Δ_N = ln mean W_loop(N) − ln mean W_loop(256), and the
surplus S_N = f (mean W_dis − mean W_loop) in W/cm³ (the same unit as P). No reading: information for check 1.

    python eval_subsample.py <runs folder> > results/eps_subsample.md
"""
import json
import math
import pathlib
import sys

import numpy as np

import analyze
import analyze_design as AD

ARMS = ("F", "T15", "I15", "T3", "I3")
SEEDS = (941, 942, 943, 944)
STRIDES = (1, 2, 4)


def trapz(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1])))


def run(d):
    """{stride: (ε, ln mean W_loop, S in W/cm³)}, max rel. deviation from summary.json, tail cycle numbers."""
    cfg = json.loads((d / "config.json").read_text())
    summ = json.loads((d / "summary.json").read_text())
    si, st = [(i, x) for i, x in enumerate(summ["stages"]) if x.get("measure")][-1]
    tail = analyze.steady_tail(st["cycles"][analyze.first_measured(st):], AD.BAL_EVENT)
    a = np.loadtxt(d / ("samples_%02d_%s.csv" % (si, st["name"])), delimiter=",")
    out, dev = {}, 0.0
    for k in STRIDES:
        wl, wd = [], []
        for c in tail:
            r = a[a[:, 0] == c["cycle"]]
            assert len(r) == cfg["samples"] + 1, (d, c["cycle"], len(r))
            r = r[::k]
            W_loop, W_dis = trapz(r[:, 2], r[:, 7]), trapz(r[:, 8], r[:, 1])    # oint H dB, int p_dis dt
            if k == 1:
                dev = max(dev, abs(W_loop / c["W_loop"] - 1), abs(W_dis / c["W_dis"] - 1))
            wl.append(W_loop)
            wd.append(W_dis)
        wl, wd = np.array(wl), np.array(wd)
        out[k] = (float(np.mean(wd / wl) - 1), math.log(wl.mean()), cfg["freq"] * (wd.mean() - wl.mean()) / 1e6)
    return out, dev, [c["cycle"] for c in tail]


def main():
    root = pathlib.Path(sys.argv[1])
    L = ["# Energy balance against the samples per cycle (eval_subsample.py, PROTOKOLL 7.108)", "",
         "Existing sample files, steady tail. ε_N = W_dis / W_loop − 1 (mean over the tail cycles); "
         "Δ_N = ln W_loop(N) − ln W_loop(256); S_N = f (W_dis − W_loop) in W/cm³. "
         "N = 256 (all samples), 128 (every 2nd), 64 (every 4th).", "",
         "| run | tail cycles | ε_256 | ε_128 | ε_64 | Δ_128 | Δ_64 | S_256 | S_128 | S_64 | check |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    worst = 0.0
    for pre in ("N0b", "B0b"):
        for arm in ARMS:
            for s in SEEDS:
                d = root / ("%s_%s_s%d" % (pre, arm, s))
                if not (d / "summary.json").exists():
                    continue
                v, dev, cyc = run(d)
                worst = max(worst, dev)
                L.append("| %s | %s | %+.3f | %+.3f | %+.3f | %+.4f | %+.4f | %+.3f | %+.3f | %+.3f | %.1e |" % (
                    d.name, ",".join(map(str, cyc)), v[1][0], v[2][0], v[4][0], v[2][1] - v[1][1], v[4][1] - v[1][1],
                    v[1][2], v[2][2], v[4][2], dev))
    L += ["", "Check of the method: largest relative deviation of W_loop and W_dis (256 samples) from summary.json: "
          "%.1e." % worst]
    print("\n".join(L))


if __name__ == "__main__":
    main()
