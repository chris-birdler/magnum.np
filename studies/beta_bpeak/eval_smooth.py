#!/usr/bin/env python
"""
Test of the Bloch-point cut (Chris 2026-10-19; rules fixed before the runs). The relaxed virgin state of each corner
run with Bloch points (conv_jobs.SM_RUNS, init.pt of <run>) is loaded; run_loops --bloch_smooth 3 cuts each Bloch point
by a tube to the particle surface and relaxes again (at most 3 steps, tube radius 2, 3, 4 cells); then the 9 mT stage
runs (<run>_sm, conv_jobs.smooth_lines).
Per run "clean": no Bloch point after the last relaxation and no Bloch point at any sample of the 9 mT stage.
Reading (agreed with Chris: the Bloch point is gone after the relaxation and does not come back in the 9 mT stage, in at
least 3 cases):
  "works"          at least 3 runs clean (the report gives the number of clean runs of all runs);
  "does not work"  fewer than 3 runs clean.
Information: Bloch points of the loaded state, after each cut and after each relaxation; tube lengths; relaxation
iterations; energy change by the cut and the relaxation; ln w of the 9 mT stage against the corner run.

    python eval_smooth.py > runs/smooth_report.md
"""
import json
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import conv_jobs as CJ
import eval_precision as EP

R = pathlib.Path(__file__).resolve().parent / "runs"
MIN_CLEAN = 3


def record(d):
    """Facts of one test run; None without DONE."""
    if not (d / "DONE").exists():
        return None
    s = json.loads((d / "summary.json").read_text())
    i = s["init"]
    steps = i.get("bloch_smooth", [])
    st = [x for x in s["stages"] if x.get("measure")][0]
    a = np.loadtxt(d / "bloch_00_amp00.csv", delimiter=",", comments="#", ndmin=2)
    n9 = a[:, 2] + a[:, 3]
    final = len(steps[-1]["bloch_points"]) if steps else sum(i.get("bloch_loaded", [0, 0]))
    return {"loaded": sum(i.get("bloch_loaded", [0, 0])), "steps": steps, "final": final,
            "max9": int(n9.max()), "share9": float((n9 > 0).mean()), "converged": bool(i.get("relax_converged")),
            "E_loaded": i.get("E_over_Kd"), "E_end": steps[-1]["E_over_Kd"] if steps else None,
            "lnw": EP.stage_tail(d)[0], "cycles": len(st["cycles"])}


def main(root=R):
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    L = ["# Bloch-point cut: tube to the particle surface, then relaxation (eval_smooth.py)", "",
         "Per step: Bloch points before the cut / after the cut / after the relaxation (relaxation iterations; tube "
         "lengths in cells).", "",
         "| run | loaded | step 0 | step 1 | step 2 | after the last relaxation | 9 mT: max (share of samples) "
         "| E/K_d loaded → end | ln w 9 mT − corner run |", "|---|" + "---|" * 8]
    clean, done = [], 0
    for n in CJ.SM_RUNS:
        r = record(root / (n + "_sm"))
        if r is None:
            L.append("| %s | – | – | – | – | not done | – | – | – |" % n)
            continue
        done += 1
        cells = ["%d / %d / %d (%d it; %s)" % (len(s["bloch_points_before"]), s["n_bloch_after_cut"],
                                               len(s["bloch_points"]), s["relax_iterations"],
                                               ", ".join("%.0f" % t["length_cells"] for t in s["tubes"]))
                 for s in r["steps"]]
        cells += ["–"] * (3 - len(cells))
        ok = r["final"] == 0 and r["max9"] == 0 and r["converged"]
        if ok:
            clean.append(n)
        e = ("%.6g → %.6g" % (r["E_loaded"], r["E_end"])) if r["E_end"] is not None else "–"
        L.append("| %s | %d | %s | %d | %d (%.0f %%) | %s | %+.3f |" % (
            n, r["loaded"], " | ".join(cells), r["final"], r["max9"], 100 * r["share9"], e,
            r["lnw"] - EP.stage_tail(root / n)[0]))
    L += ["", "Clean: %d of %d runs done (%s)." % (len(clean), done, ", ".join(clean) or "none"), "",
          "**Reading: %s.**" % ("works" if len(clean) >= MIN_CLEAN else
                                "does not work" if done == len(CJ.SM_RUNS) else "open (runs not done)")]
    print("\n".join(L))


if __name__ == "__main__":
    main()
