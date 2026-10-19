#!/usr/bin/env python
"""
Sensitivity check for the rule gap of PROTOKOLL 7.44: S029 has a negative mean w_loop at 9 mT (a state change in the
kept cycles released stored energy). Here S029 enters the analysis with the stationary cycles after the event at
9 mT (cycles 3 ... 6; all other amplitudes as usual), and the pre-registered analysis is repeated (bootstrap 1000).
Output: the 8 primary range effects and classes with and without S029.

    python eval_s029.py > results/s029_sensitivity.md
"""
import json
import pathlib

import numpy as np

import analyze
import analyze_design as AD

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE / "runs"
_orig = analyze.load_run


def load_run_s029(d):
    cfg, rows = _orig(d)
    if pathlib.Path(d).name == "S029":
        summ = json.loads((pathlib.Path(d) / "summary.json").read_text())
        st = [s for s in summ["stages"] if s.get("measure")][0]
        cyc = st["cycles"][3:]                                   # after the event in cycles 1 ... 2
        r = rows[0]
        r["w_loop"] = float(np.mean([c["w_loop"] for c in cyc]))
        r["w_loop_raw"] = r["w_loop"]
        r["w_dis"] = float(np.mean([c["w_dis"] for c in cyc]))
        r["balance"] = r["w_dis"] / r["w_loop_raw"]
        r["n_meas"] = len(cyc)
    return cfg, rows


def main():
    base, _ = AD.load_runs(RUNS)
    analyze.load_run = load_run_s029
    withs, skipped = AD.load_runs(RUNS)
    analyze.load_run = _orig
    s = [r for r in withs if r["name"] == "S029"][0]
    print("# S029 sensitivity (rule gap, PROTOKOLL 7.44)\n")
    print("S029 at 9 mT with cycles 3 … 6: w_loop %.3g (reduced units), flagged %s; runs %d → %d (skipped now: %s)\n" %
          (np.exp(s["lnw@9mT"]), s["flagged"], len(base), len(withs), skipped or "none"))
    r0, r1 = AD.analyse(base, 1000), AD.analyse(withs, 1000)
    print("| output | factor | without S029 (n = %d) | with S029 (n = %d) | class without / with |" % (len(base), len(withs)))
    print("|---|---|---|---|---|")
    for o in AD.PRIMARY:
        for f in AD.FACTORS:
            a, b = r0["main"][(o, f)], r1["main"][(o, f)]
            print("| %s | %s | %+.3f (%+.3f … %+.3f) | %+.3f (%+.3f … %+.3f) | %s / %s |" %
                  (o, f, a["eff"], *a["ci95"], b["eff"], *b["ci95"], r0["classes"][(o, f)], r1["classes"][(o, f)]))
    print("\nMesh-check corners: without %s, with %s." % (r0["mesh_trigger"][1], r1["mesh_trigger"][1]))


if __name__ == "__main__":
    main()
