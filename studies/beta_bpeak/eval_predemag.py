#!/usr/bin/env python
"""
Test of the AC demagnetization after the virgin relaxation (Chris 2026-10-18, Bloch points; rules fixed before the
runs). The 4 corner runs with the most Bloch points (K2_0_nx, K2_1_nx, K4_1_nx, K3_1_nx) from their relaxed virgin
state: AC demagnetization from (a) the drive of their 150 mT stage (pd150), (b) twice that (pd300), (c) saturation
(pdsat), then the relaxation again and the 9 mT stage (conv_jobs.predemag_lines).
Reading per variant:
  "clean"      in all 4 runs: no Bloch point after the second relaxation, the remanence |M| / M_s (box average) below
               0.01, and no Bloch point at any sample of the 9 mT stage;
  "not clean"  all other cases (the report names the failed criterion).
Information: Bloch points of the virgin state, maximum during the demagnetization, after the second relaxation;
iterations of the second relaxation; time of the demagnetization; ln w of the 9 mT stage against the corner run itself
(its virgin state with the Bloch points).

    python eval_predemag.py > runs/predemag_report.md
"""
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import conv_jobs as CJ
import eval_precision as EP

R = pathlib.Path(__file__).resolve().parent / "runs"
VARIANTS = (("pd150", "150 mT"), ("pd300", "2 × 150 mT"), ("pdsat", "saturation"))
REM_MAX = 0.01


def record(d):
    """Facts of one test run; None without DONE."""
    if not (d / "DONE").exists():
        return None
    s = json.loads((d / "summary.json").read_text())
    i, pd = s["init"], s["init"].get("predemag", {})
    st = [x for x in s["stages"] if x.get("measure")][0]
    return {"bloch_virgin": sum(i.get("bloch_virgin", [0, 0])),
            "bloch_demag_max": max([h["n_bloch_max"] for h in pd.get("halves", [])] or [0]),
            "bloch_after": len(pd.get("bloch_points", [])),
            "remanence": float(np.linalg.norm(pd.get("M_over_Ms", [np.nan] * 3))),
            "bloch_9mT_max": max(c.get("bloch", {}).get("max", 0) for c in st["cycles"]),
            "relax2_it": pd.get("relax_iterations"), "demag_s": pd.get("wall_s", float("nan")),
            "n_half": pd.get("n_half"), "h0": pd.get("h0"), "lnw": EP.stage_tail(d)[0]}


def main(root=R):
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    L = ["# AC demagnetization after the virgin relaxation: Bloch points and remanence (eval_predemag.py)", "",
         "| run | variant | h0 | half cycles | BP virgin | BP max in the demag. | BP after relax 2 | remanence M/M_s | BP max at 9 mT "
         "| relax 2 iterations | demag. time (h) | ln w 9 mT − corner run |", "|---|" + "---|" * 11]
    for v, lab in VARIANTS:
        bad = []
        for n in CJ.PD_RUNS:
            r = record(root / ("%s_%s" % (n, v)))
            if r is None:
                bad.append("%s not done" % n)
                L.append("| %s | %s | – | – | – | – | – | – | – | – | – | – |" % (n, lab))
                continue
            ref = EP.stage_tail(root / n)[0]
            L.append("| %s | %s | %.4g | %s | %d | %d | %d | %.4f | %d | %s | %.2f | %+.3f |" % (
                n, lab, r["h0"], r["n_half"], r["bloch_virgin"], r["bloch_demag_max"], r["bloch_after"],
                r["remanence"], r["bloch_9mT_max"], r["relax2_it"], r["demag_s"] / 3600, r["lnw"] - ref))
            if r["bloch_after"]:
                bad.append("%s: %d Bloch points after relax 2" % (n, r["bloch_after"]))
            if not r["remanence"] < REM_MAX:
                bad.append("%s: remanence %.4f" % (n, r["remanence"]))
            if r["bloch_9mT_max"]:
                bad.append("%s: Bloch points at 9 mT" % n)
        L.append("")
        L.append("**Variant %s (%s): %s.**%s" % (v, lab, "clean" if not bad else "not clean",
                                                 (" " + "; ".join(bad)) if bad else ""))
        L.append("")
    print("\n".join(L))


if __name__ == "__main__":
    main()
