#!/usr/bin/env python
"""
Stage 4 of PLAN 11: the bundle of speed patches on the prepared states of the main study. For each state the arms
old, null (perturbed start, old code) and new each ran 4 drive cycles from the same start (--max_cycles 4); each arm
folder has summary.json and m.pt of step_dump.py (m after the 1024 sample intervals).
- ln P: ln w_loop of cycle 3 (P = f K_d w_loop). D = ln P(arm) − ln P(old).
- Event cycle (β(B_peak) record 7.54): w_loop at most 0 or |w_dis / w_loop − 1| above 0.20, in cycle 3.
- State change, two readings: (a) PLAN 11 as approved: magnetic cells with m · m0 below 0.9 between the start and the
  end; (b) proposal of 2026-10-22: cells with m_arm · m_old below 0.9 at the end (arm against old).
- Deciding criterion: the states without an event cycle and without a state change in all arms (at least 6); the 90 %
  t-interval of the mean of D_new inside ± δ (δ 0.01). Fail if |mean D_new| is above δ; not determined otherwise.

    python s4_eval.py --root DIR --states K1_0_sat ... --init_dir DIR --json OUT
"""
import argparse
import json
import math

import numpy as np
import torch
from scipy import stats

DELTA = 0.01


def load(d, init, n_end):
    s = json.load(open(d + "/summary.json"))
    cyc = s["stages"][0]["cycles"]
    c = cyc[3]
    m = torch.load(d + "/m.pt", map_location="cpu")["m"][n_end].double()
    m0 = torch.load(init, map_location="cpu", weights_only=False)["m"].double()
    mag = m0.norm(dim=-1) > 0
    eps = c["w_dis"] / c["w_loop"] - 1 if c["w_loop"] != 0 else float("inf")
    return {"lnP": math.log(c["w_loop"]) if c["w_loop"] > 0 else float("nan"), "w_loop": c["w_loop"], "eps": eps,
            "event_cycle": bool(c["w_loop"] <= 0 or abs(eps) > 0.20), "n_cycles": len(cyc),
            "w_loop_cycles": [x["w_loop"] for x in cyc], "steps": [x["steps"] for x in cyc],
            "start_end_cells": int(((m * m0).sum(-1)[mag] < 0.9).sum()), "m": m, "mag": mag}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--states", nargs="+", required=True)
    ap.add_argument("--init_dir", required=True)
    ap.add_argument("--n_end", type=int, default=1024)
    ap.add_argument("--json", required=True)
    a = ap.parse_args()
    rows = []
    for st in a.states:
        init = "%s/%s/init.pt" % (a.init_dir, st)
        arms = {k: load("%s/%s/%s" % (a.root, st, k), init if k != "null" else "%s/%s/null/init.pt" % (a.root, st),
                        a.n_end) for k in ("old", "null", "new")}
        old = arms["old"]
        r = {"state": st}
        for k in ("null", "new"):
            x = arms[k]
            dot = (x["m"] * old["m"]).sum(-1)[old["mag"]]
            r["div_cells_" + k] = int((dot < 0.9).sum())
            r["div_rms_" + k] = float((x["m"] - old["m"]).norm(dim=-1)[old["mag"]].pow(2).mean().sqrt())
            r["D_" + k] = x["lnP"] - old["lnP"]
        for k, x in arms.items():
            r.update({k + "_" + q: x[q] for q in ("lnP", "w_loop", "eps", "event_cycle", "start_end_cells",
                                                   "w_loop_cycles", "steps")})
        r["event_any"] = any(x["event_cycle"] for x in arms.values())
        r["change_a"] = any(x["start_end_cells"] > 0 for x in arms.values())
        r["change_b"] = r["div_cells_null"] > 0 or r["div_cells_new"] > 0
        rows.append(r)
    res = {"rows": rows, "delta": DELTA}
    for reading in ("a", "b"):
        use = [r for r in rows if not r["event_any"] and not r["change_" + reading]]
        out = {"n": len(use), "states": [r["state"] for r in use]}
        if len(use) >= 2:
            D = np.array([r["D_new"] for r in use])
            Dn = np.array([r["D_null"] for r in use])
            mean, se = float(D.mean()), float(D.std(ddof=1) / math.sqrt(len(D)))
            t90 = float(stats.t.ppf(0.95, len(D) - 1))
            lo, hi = mean - t90 * se, mean + t90 * se
            out.update({"mean_D_new": mean, "SE": se, "t90": t90, "interval90": [lo, hi],
                        "p_t_test_zero": float(stats.ttest_1samp(D, 0.0).pvalue),
                        "rms_ratio_new_null": float(math.sqrt((D ** 2).mean() / max((Dn ** 2).mean(), 1e-30)))})
            if len(use) < 6:
                out["result"] = "not determined: fewer than 6 states"
            elif -DELTA <= lo and hi <= DELTA:
                out["result"] = "pass"
            elif abs(mean) > DELTA:
                out["result"] = "FAIL"
            else:
                out["result"] = "not determined"
        else:
            out["result"] = "not determined: fewer than 6 states"
        res["reading_" + reading] = out
    for r in rows:
        print("%-9s D_new %+.2e  D_null %+.2e  events old/null/new %s/%s/%s  start-end cells %d  divergence cells null %d new %d"
              "  RMS |dm| null %.1e new %.1e" % (r["state"], r["D_new"], r["D_null"], r["old_event_cycle"],
                                                  r["null_event_cycle"], r["new_event_cycle"], r["old_start_end_cells"],
                                                  r["div_cells_null"], r["div_cells_new"], r["div_rms_null"],
                                                  r["div_rms_new"]))
    for reading, name in (("a", "PLAN 11 as approved (start against end)"), ("b", "proposal (arm against old)")):
        o = res["reading_" + reading]
        print("reading %s, %s: n %d, %s" % (reading, name, o["n"], o["result"]) +
              ("; mean D_new %+.2e ± %.1e (SE), 90 %% interval [%+.2e, %+.2e], t-test p %.2f, RMS ratio new/null %.2f"
               % (o["mean_D_new"], o["SE"], *o["interval90"], o["p_t_test_zero"], o["rms_ratio_new_null"])
               if "mean_D_new" in o else ""))
    json.dump(res, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
