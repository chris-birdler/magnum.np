#!/usr/bin/env python
"""
Test of the start "each particle uniform along its own random direction" (Chris 2026-10-18, Bloch points, way 2; rules
fixed before the runs). The 4 corner runs with the most Bloch points (conv_jobs.PD_RUNS) relax from this start
(--init_mode particle) instead of random m on 6 l_ex blocks, then the 9 mT stage (conv_jobs.initmode_lines, <run>_ip).
Reading:
  "clean"      in all 4 runs: no Bloch point after the virgin relaxation and at no sample of the 9 mT stage;
  "not clean"  all other cases (the report names the run).
Information: Bloch points of the block start (the corner run) and of the particle start; relaxation iterations;
energy E / (K_d V) and |M| / M_s of the relaxed state; ln w of the 9 mT stage against the corner run.

    python eval_initmode.py > runs/initmode_report.md
"""
import json
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import conv_jobs as CJ
import eval_precision as EP

R = pathlib.Path(__file__).resolve().parent / "runs"


def init_facts(d):
    s = json.loads((d / "summary.json").read_text())
    i = s["init"]
    bp = i.get("bloch_virgin")
    st = [x for x in s["stages"] if x.get("measure")]
    return {"bp": sum(bp) if bp is not None else None, "it": i.get("relax_iterations"), "E": i.get("E_over_Kd"),
            "M": float(np.linalg.norm(i.get("M_over_Ms", [np.nan] * 3))),
            "bp9": max(c.get("bloch", {}).get("max", 0) for c in st[0]["cycles"]) if st else None}


def main(root=R):
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    census = json.loads((root / "census_bloch_nx_v2.json").read_text()) if (root / "census_bloch_nx_v2.json").exists() else {}
    L = ["# Start with each particle uniform: Bloch points (eval_initmode.py)", "",
         "| run | BP block start (corner run) | BP particle start | relaxation iterations block / particle | "
         "E/(K_d V) block / particle | M/M_s particle | BP max at 9 mT | ln w 9 mT particle − block |",
         "|---|---|---|---|---|---|---|---|"]
    bad = []
    for n in CJ.PD_RUNS:
        d = root / (n + "_ip")
        if not (d / "DONE").exists():
            bad.append("%s not done" % n)
            L.append("| %s | – | – | – | – | – | – | – |" % n)
            continue
        a, b = init_facts(root / n), init_facts(d)
        bpb = len(census.get(n, {}).get("init.pt", {}).get("bloch_points", [])) if census else a["bp"]
        L.append("| %s | %s | %d | %s / %s | %.5f / %.5f | %.4f | %d | %+.3f |" % (
            n, bpb, b["bp"], a["it"], b["it"], a["E"], b["E"], b["M"], b["bp9"],
            EP.stage_tail(d)[0] - EP.stage_tail(root / n)[0]))
        if b["bp"] or b["bp9"]:
            bad.append("%s: %d Bloch points after the relaxation, %d max at 9 mT" % (n, b["bp"], b["bp9"]))
    L += ["", "**Reading: %s.**%s" % ("clean" if not bad else "not clean", (" " + "; ".join(bad)) if bad else "")]
    print("\n".join(L))


if __name__ == "__main__":
    main()
