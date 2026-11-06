#!/usr/bin/env python
"""
Numbers of gate 1 of PLAN 4.21 (PROTOKOLL 7.146 part 6), written to results/plan421_gate1.json:
  half_cycles  per K5 run (stdout logs of the folder): the preparation half cycles that ended, h of the last one, and
               the half cycles that ended with a Bloch point
  stall        the samples of diag_stall421.py (K5f12d_diag.log) before and in the stall (the stall: state time at or
               after the time of the first sample in which |s_j| falls below 1e-3 and stays below): medians and
               ranges of the last step, the smallest |s_j|, the largest master field, the share of rejected steps
  weights      the V2 tied cells whose largest trilinear weight is at most 0.5 (then s_j can be 0), at K5 (the run
               geometry of K5f12) and at K3 at q 16 (geometry.fcc_box of jobs/plan421.txt; a few minutes)
    python eval_gate421.py runs/plan421_gate1
"""
import json
import pathlib
import re
import statistics as st
import sys

import numpy as np

import scan421
from field_terms_extra import tie_ghost_mirror
from geometry import fcc_box

HERE = pathlib.Path(__file__).resolve().parent
HC = re.compile(r"\[predemag\] half cycle (\d+): h ([0-9.e+-]+) .* end \+(\d+) -(\d+)")
NUM = r" +([-0-9.e+naninf]+)"


def half_cycles(log):
    rows = [HC.search(l) for l in log.read_text().splitlines()]
    rows = [(int(m[1]), float(m[2]), int(m[3]) + int(m[4])) for m in rows if m]
    return {"n_ended": len(rows), "h_last": rows[-1][1] if rows else None,
            "n_ended_with_bloch": sum(r[2] > 0 for r in rows)}


def stall(diag):
    rows = []
    for l in diag.read_text().splitlines():
        if " wall " not in l:
            continue
        g = lambda k: float(re.search(k + NUM, l).group(1))
        rows.append({"t": g(r" t"), "dt": g("dt last"), "dt_min": g(" min"), "ns": g(r"min\|s_j\|"),
                     "hm": g("masters"), "ha": g(" all"), "tries": g("tries"), "rej": g("rejected"),
                     "cell": int(re.search(r"tied cell (-?\d+)", l).group(1))})
    rows = [r for r in rows if r["tries"] > 10]
    k = next(i for i in range(len(rows)) if all(r["ns"] < 1e-3 for r in rows[i:]))
    out = {}
    for name, part in (("before", rows[:k]), ("stall", rows[k:])):
        q = lambda x: {"median": st.median(r[x] for r in part), "min": min(r[x] for r in part),
                       "max": max(r[x] for r in part)}
        out[name] = {"n_samples": len(part), "dt_last_s": q("dt"), "dt_smallest_s": min(r["dt_min"] for r in part),
                     "min_s_j": q("ns"), "max_h_masters_A_per_m": q("hm"), "max_h_all_A_per_m": q("ha"),
                     "rejected_share": sum(r["rej"] for r in part) / sum(r["tries"] for r in part),
                     "cells_of_min_s_j": sorted({r["cell"] for r in part if r["ns"] < 1e-3})}
    out["before"]["n_samples_s_j_below_0.1"] = sum(r["ns"] < 0.1 for r in rows[:k])
    out["stall"]["state_time_s"] = [rows[k]["t"], rows[-1]["t"]]
    return out


def weights(folder):
    c = json.loads((folder / "K5f12" / "config.json").read_text())
    g = np.load(folder / "K5f12" / "geometry.npz")
    cases = [("K5", c["N_used"], c["d"], c["dx"], g["ids"].astype(np.int64), g["frac"])]
    line = [l for l in (HERE / "jobs" / "plan421.txt").read_text().splitlines() if l.startswith("K3_0_v16 ")][0]
    N, d, dx, _, _ = scan421.mesh(line.split(None, 1)[1])
    b = fcc_box(N, d, dx, sub=16)
    cases.append(("K3_q16", N, d, dx, b["ids"].astype(np.int64), b["frac"]))
    res = {}
    for name, N, d, dx, ids, frac in cases:
        wmax = np.asarray(tie_ghost_mirror(ids, frac, N, d, dx).w).max(1)
        res[name] = {"N": int(N), "n_tied": int(len(wmax)), "n_wmax_le_0.5": int((wmax <= 0.5).sum()),
                     "share": float((wmax <= 0.5).mean()), "wmax_median": float(np.median(wmax))}
    return res


if __name__ == "__main__":
    folder = pathlib.Path(sys.argv[1]).resolve()
    res = {"half_cycles": {n: half_cycles(folder / (n + ".log")) for n in ("K5v", "K5f", "K5f12", "K5f12d")},
           "stall": stall(folder / "K5f12d_diag.log"), "weights": weights(folder)}
    (HERE / "results" / "plan421_gate1.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
