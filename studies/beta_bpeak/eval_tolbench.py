#!/usr/bin/env python
"""
Cost of the cycle tolerance and of the final stray-field setting on the mesh of the main study (dx 3; decision Chris
2026-10-17: measure, then decide; conv_jobs.py --newell tolbench). No reading: the numbers go to Chris.
Per corner (K1_0 d/l_ex 212, K2_0 d/l_ex 300) and variant (a5: atol 1e-5 as the main study; a6: atol 1e-6; tight:
atol 1e-7, rtol 1e-6; the same relaxed state): wall time and steps per cycle (cycles 1 and 2; cycle 0 is excluded, it
contains the transient and, for a6 and tight, the kernel), mean step, ε = w_dis / w_loop − 1 and P of cycles 1 and 2.
Reference for the solver cost: the 9 mT stage of the old corner run (DemagFieldPBC, the same corner and seed): wall
time per step.

    python eval_tolbench.py > runs/tolbench_report.md
"""
import json
import pathlib

import numpy as np

import conv_jobs as CJ

R = pathlib.Path(__file__).resolve().parent / "runs"


def cycles(d):
    if not (d / "DONE").exists() and not (d / "summary.json").exists():
        return None, None
    S = json.loads((d / "summary.json").read_text())
    cfg = json.loads((d / "config.json").read_text())
    st = [x for x in S["stages"] if x.get("measure")][0]
    return st["cycles"], cfg


def main(root=R):
    L = ["# Cost of the cycle tolerance and of the final stray-field setting on dx 3 (eval_tolbench.py)", "",
         "| run | wall per cycle (s) | steps per cycle | mean step (ps) | ε (cycles 1, 2) | P (W/cm³) | GPU |",
         "|---|---|---|---|---|---|---|"]
    for c in CJ.TB_CORNERS:
        base = None
        for v, _ in CJ.TB_VARIANTS:
            cyc, cfg = cycles(root / ("TB_%s_%s" % (c, v)))
            if not cyc or len(cyc) < 3:
                L.append("| TB_%s_%s | not complete | | | | | |" % (c, v))
                continue
            w = np.mean([x["wall_s"] for x in cyc[1:3]])
            s = np.mean([x["steps"] for x in cyc[1:3]])
            base = base or (w, s)
            L.append("| TB_%s_%s | %.0f (%+.0f %%) | %.0f (%+.0f %%) | %.3f | %s | %.3f | %s |" % (
                c, v, w, 100 * (w / base[0] - 1), s, 100 * (s / base[1] - 1),
                1e12 * np.mean([x["mean_dt"] for x in cyc[1:3]]),
                " / ".join("%+.3f" % (x["w_dis"] / x["w_loop"] - 1) for x in cyc[1:3]),
                cfg["freq"] * np.mean([x["W_loop"] for x in cyc[1:3]]) / 1e6, cfg.get("gpu", "–")))
        old, ocfg = cycles(root / c)
        new, ncfg = cycles(root / ("TB_%s_a5" % c))
        if old and new and len(new) >= 3:
            po = np.mean([x["wall_s"] / x["steps"] for x in old[1:3]])
            pn = np.mean([x["wall_s"] / x["steps"] for x in new[1:3]])
            L.append("| %s, old run (DemagFieldPBC, 9 mT stage) | – | – | – | – | – | %s |" % (c, ocfg.get("gpu", "–")))
            L.append("")
            L.append("%s: wall time per step %.2f ms (DemagFieldPBC, old run) against %.2f ms (final setting, a5): "
                     "%+.0f %%." % (c, 1e3 * po, 1e3 * pn, 100 * (pn / po - 1)))
            L.append("")
    print("\n".join(L))


if __name__ == "__main__":
    main()
