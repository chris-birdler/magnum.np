#!/usr/bin/env python
"""
S3 of PLAN 11: one drive cycle of one production state, arms old, null runs and new. Each arm folder has the
summary.json of run_loops (--max_cycles 1) and m.pt of step_dump.py (m after the 256 sample intervals of the cycle).
ln P = ln w_loop of the cycle (P = f K_d w_loop; f and K_d are the same in all arms). D = ln P(arm) − ln P(old).
State change: magnetic cells (|m0| above 0) with m · m0 below 0.9 between the start state (init.pt of the arm) and the
end of the cycle (criterion of check_events.py, PROTOKOLL 3.26).
Criteria: |D_new| at most 1e-4 (deciding) and at most 10 × max(RMS(D_null), 1e-6); calibration: RMS(D_null) at most
1e-5. Writes a JSON with all numbers.

    python s3_eval.py --old DIR --null DIR DIR DIR --new DIR --init_old INIT --init_null INIT INIT INIT --json OUT
"""
import argparse
import json
import math

import torch


def arm(d, init):
    s = json.load(open(d + "/summary.json"))
    c = s["stages"][0]["cycles"][0]
    m = torch.load(d + "/m.pt", map_location="cpu")["m"][256].double()
    m0 = torch.load(init, map_location="cpu", weights_only=False)["m"].double()
    mag = m0.norm(dim=-1) > 0
    dot = (m * m0).sum(-1)[mag]
    return {"dir": d, "w_loop": c["w_loop"], "w_dis": c["w_dis"], "steps": c["steps"], "bloch": c.get("bloch"),
            "lnP": math.log(c["w_loop"]), "state_change_cells": int((dot < 0.9).sum()), "min_m_dot_m0": float(dot.min()),
            "cells_below_0.97": int((dot < 0.97).sum()), "cells_below_0.99": int((dot < 0.99).sum())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--null", nargs="+", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--init_old", required=True)
    ap.add_argument("--init_null", nargs="+", required=True)
    ap.add_argument("--json", required=True)
    a = ap.parse_args()
    old = arm(a.old, a.init_old)
    nulls = [arm(d, i) for d, i in zip(a.null, a.init_null)]
    new = arm(a.new, a.init_old)
    # divergence of each arm from the old arm at the end of the cycle (proposal of 2026-10-22: the start-to-end
    # criterion marks the first cycle of a prepared state in every arm)
    m_old = torch.load(a.old + "/m.pt", map_location="cpu")["m"][256].double()
    mag = torch.load(a.init_old, map_location="cpu", weights_only=False)["m"].double().norm(dim=-1) > 0
    for x in nulls + [new]:
        mx = torch.load(x["dir"] + "/m.pt", map_location="cpu")["m"][256].double()
        dot = (mx * m_old).sum(-1)[mag]
        x["div_cells_0.9"] = int((dot < 0.9).sum())
        x["div_cells_0.9999"] = int((dot < 0.9999).sum())
        x["div_rms_dm"] = float((mx - m_old).norm(dim=-1)[mag].pow(2).mean().sqrt())
    D_null = [x["lnP"] - old["lnP"] for x in nulls]
    D_new = new["lnP"] - old["lnP"]
    rms = math.sqrt(sum(d * d for d in D_null) / len(D_null))
    change = any(x["state_change_cells"] > 0 for x in [old, new] + nulls)
    div = any(x["div_cells_0.9"] > 0 for x in nulls + [new])
    res = {"old": old, "null": nulls, "new": new, "D_null": D_null, "D_new": D_new, "RMS_D_null": rms,
           "calibration_ok": rms <= 1e-5, "state_change": change, "divergence": div,
           "pass": (not change) and abs(D_new) <= 1e-4 and abs(D_new) <= 10 * max(rms, 1e-6),
           "pass_proposal": (not div) and abs(D_new) <= 1e-3 and abs(D_new) <= 10 * max(rms, 1e-6)}
    for name, x in [("old", old)] + [("null %d" % k, x) for k, x in enumerate(nulls, 1)] + [("new", new)]:
        print("%-7s w_loop %.10e  steps %5d  state-change cells %d  min m.m0 %.4f" % (
            name, x["w_loop"], x["steps"], x["state_change_cells"], x["min_m_dot_m0"]))
    print("D_null %s  RMS %.2e (calibration %s)" % (" ".join("%+.2e" % d for d in D_null), rms,
                                                     "ok" if res["calibration_ok"] else "NOT ok: Chris sets the bound"))
    for name, x in [("null %d" % k, x) for k, x in enumerate(nulls, 1)] + [("new", new)]:
        print("%-7s against old at the end: cells below 0.9 %d, below 0.9999 %d, RMS |dm| %.2e" % (
            name, x["div_cells_0.9"], x["div_cells_0.9999"], x["div_rms_dm"]))
    print("proposal (arm against old at the end, |D_new| at most 1e-3 and 10 RMS(D_null)): %s"
          % ("not usable (divergence)" if div else ("pass" if res["pass_proposal"] else "FAIL")))
    print("D_new %+.2e  state change %s  RESULT S3 %s" % (D_new, change,
                                                          "not usable (state change)" if change else
                                                          ("pass" if res["pass"] else "FAIL")))
    json.dump(res, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
