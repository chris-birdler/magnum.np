#!/usr/bin/env python
"""
Test of protocol S (Chris 2026-10-18, Bloch points; rules fixed before the runs). The 4 corner runs with the most Bloch
points (K2_0, K2_1, K4_1, K3_1; conv_jobs.satstart_lines, names <run>_nx_ps): the start uniform along the drive direction
(no random start, no relaxation before), AC demagnetization from saturation and the relaxation, AC from twice the drive
of the 150 mT stage of <run>_nx (300 mT) and the relaxation; if a Bloch point remains, one further pass from 300 mT and
the relaxation; then the 9 mT stage.
Reading:
  "clean"      in all 4 runs: no Bloch point after the last relaxation, the remanence |M| / M_s (box average) below
               0.01, and no Bloch point at any sample of the 9 mT stage;
  "not clean"  all other cases (the report names the failed criterion).
Information per run: Bloch points after each pass (after its relaxation) and their maximum during the pass; time and
number of half cycles per pass; relaxation iterations per pass; ln w of the 9 mT stage against the corner run <run>_nx
(random start, relaxation) and against the AC-300 test run <run>_nx_pd300 (option C, attempt 0).

    python eval_protocol_s.py > runs/protocol_s_report.md
"""
import json
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import conv_jobs as CJ
import eval_precision as EP

R = pathlib.Path(__file__).resolve().parent / "runs"
REM_MAX = 0.01


def record(d):
    """Facts of one test run; None without DONE."""
    if not (d / "DONE").exists():
        return None
    s = json.loads((d / "summary.json").read_text())
    i = s["init"]
    passes = i.get("predemag_passes", [])
    st = [x for x in s["stages"] if x.get("measure")][0]
    return {"passes": [{"h0": p["h0"], "n_half": p["n_half"], "bp_after": len(p["bloch_points"]),
                        "bp_max": max([h["n_bloch_max"] for h in p["halves"]] or [0]),
                        "relax_it": p["relax_iterations"], "wall_h": p["wall_s"] / 3600} for p in passes],
            "bloch_after": len(i["predemag"]["bloch_points"]),
            "remanence": float(np.linalg.norm(i["predemag"]["M_over_Ms"])),
            "bloch_9mT_max": max(c.get("bloch", {}).get("max", 0) for c in st["cycles"]),
            "converged": bool(i.get("relax_converged")), "lnw": EP.stage_tail(d)[0]}


def main(root=R):
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    L = ["# Protocol S: start saturated, AC from saturation, AC from 300 mT (eval_protocol_s.py)", "",
         "Per pass: h0 (H/M_s), half cycles, Bloch points max during the pass / after its relaxation, relaxation "
         "iterations, time (h). Pass 0: from saturation; pass 1: from 300 mT; pass 2: the further pass from 300 mT.", "",
         "| run | pass 0 | pass 1 | pass 2 | BP after the last relaxation | remanence M/M_s | BP max at 9 mT "
         "| ln w 9 mT − corner run | ln w 9 mT − AC-300 test |", "|---|" + "---|" * 8]
    bad = []
    for n in CJ.PD_RUNS:
        r = record(root / (n + "_ps"))
        if r is None:
            bad.append("%s not done" % n)
            L.append("| %s | – | – | – | – | – | – | – | – |" % n)
            continue
        cells = ["%.4g, %d, %d / %d, %d, %.2f" % (p["h0"], p["n_half"], p["bp_max"], p["bp_after"], p["relax_it"],
                                                  p["wall_h"]) for p in r["passes"]]
        cells += ["–"] * (3 - len(cells))
        ref_c = EP.stage_tail(root / n)[0]
        pd = root / (n + "_pd300")
        ref_p = EP.stage_tail(pd)[0] if (pd / "DONE").exists() else float("nan")
        L.append("| %s | %s | %d | %.4f | %d | %+.3f | %+.3f |" % (
            n, " | ".join(cells), r["bloch_after"], r["remanence"], r["bloch_9mT_max"], r["lnw"] - ref_c,
            r["lnw"] - ref_p))
        if not r["converged"]:
            bad.append("%s: a relaxation did not converge" % n)
        if r["bloch_after"]:
            bad.append("%s: %d Bloch points after the last relaxation" % (n, r["bloch_after"]))
        if not r["remanence"] < REM_MAX:
            bad.append("%s: remanence %.4f" % (n, r["remanence"]))
        if r["bloch_9mT_max"]:
            bad.append("%s: Bloch points at 9 mT" % n)
    L += ["", "**Reading: %s.**%s" % ("clean" if not bad else "not clean", (" " + "; ".join(bad)) if bad else "")]
    print("\n".join(L))


if __name__ == "__main__":
    main()
