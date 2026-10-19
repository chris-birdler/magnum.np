#!/usr/bin/env python
"""
Are the loops of the main data periodic? (option 1 of Chris 2026-10-14, PROTOKOLL 7.100). The stop rule of a stage (≥ 2
calm cycles at the end, |w_dis / w_loop − 1| ≤ 0.2, at most 12 cycles) checks the energy balance of each cycle, not
whether the cycles repeat. At 60 MHz it accepted cycles of a state that still drifts (7.99). This script checks every
stage of the runs of the main analysis (analyze_design.load_runs with reruns, α 0.1 design, 30 MHz) and of C′.

Per stage, from the last 3 complete cycles (the cycles of the steady tail): RN = share of the box dissipation that
changes from cycle to cycle (eval_jumps_time.measures), loop closure |ΔM| / M_peak (mean of the 3 cycles, summary.json),
relative cycle-to-cycle SD of w_loop and the drift (w of the last minus w of the first of the 3 cycles, relative).

Classes (fixed before the data were looked at, 2026-10-14):
  periodic       RN ≤ 0.05 and closure ≤ 0.001
  not periodic   RN > 0.10 or closure > 0.002
  marginal       otherwise

    python eval_periodic.py > results/periodic.md      # also results/periodic.json
"""
import json
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import eval_jumps_time as EJ

R = pathlib.Path(__file__).resolve().parent / "runs"
AMPS = list(AD.AMPS_MT)
CPRIME = ["%s_%d%s" % (b, k, s) for b, ss in (("K1", ("_fhalf", "_f1", "_f2")), ("K1A10", ("_fhalf", "", "_f2")))
          for k in range(3) for s in ss]


def classify(rn, clo):
    if rn <= 0.05 and clo <= 0.001:
        return "periodic"
    if rn > 0.10 or clo > 0.002:
        return "not periodic"
    return "marginal"


def stages(d):
    summ = json.loads((d / "summary.json").read_text())["stages"]
    files = sorted(d.glob("samples_*_amp*.csv"))
    out = []
    for k, f in enumerate(files):
        s = EJ.series(f)
        cyc = summ[k]["cycles"][-3:] if k < len(summ) else []
        if s is None or len(cyc) < 3:
            out.append(None)
            continue
        m = EJ.measures(*s)
        w = np.array([c["w_loop"] for c in cyc])
        clo = float(np.mean([c["closure"] for c in cyc]))
        out.append({"RN": float(m["RN"]), "clo": clo, "sd": float(w.std(ddof=1) / abs(w.mean())),
                    "drift": float((w[-1] - w[0]) / abs(w.mean())), "cls": classify(m["RN"], clo)})
    return out


def table(names, title):
    D = {n: stages(R / n) for n in names}
    L = ["## %s (%d runs)" % (title, len(names)), "",
         "| B (mT) | periodic | marginal | not periodic | median RN | median closure | median cycle SD of w | "
         "median abs. drift |", "|---|---|---|---|---|---|---|---|"]
    out = {}
    L2 = ["", "Added after the classes were looked at: is the loss per cycle stationary? Not stationary: |drift| over the "
          "3 cycles > 5 % or cycle SD of w > 5 %.", "", "| B (mT) | w not stationary | of them also not periodic |",
          "|---|---|---|"]
    for i, b in enumerate(AMPS):
        v = [D[n][i] for n in names if len(D[n]) > i and D[n][i]]
        c = {k: sum(x["cls"] == k for x in v) for k in ("periodic", "marginal", "not periodic")}
        ns = [x for x in v if abs(x["drift"]) > 0.05 or x["sd"] > 0.05]
        out["%g" % b] = {"n": len(v), **c, "w_not_stationary": len(ns)}
        L2.append("| %g | %d of %d | %d |" % (b, len(ns), len(v), sum(x["cls"] == "not periodic" for x in ns)))
        L.append("| %g | %d | %d | %d | %.3f | %.4f | %.1f %% | %.1f %% |" % (
            b, c["periodic"], c["marginal"], c["not periodic"], np.median([x["RN"] for x in v]),
            np.median([x["clo"] for x in v]), 100 * np.median([x["sd"] for x in v]),
            100 * np.median([abs(x["drift"]) for x in v])))
    return L + L2, out, D


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    recs, _ = AD.load_runs(R, reruns=True)
    names = [r["name"] for r in recs]
    L = ["# Periodicity of the evaluated loops (last 3 cycles of each stage)", "",
         "Classes fixed before the data were looked at: periodic RN ≤ 0.05 and closure ≤ 0.001; not periodic RN > 0.10 "
         "or closure > 0.002; marginal otherwise. RN: share of the dissipation that changes from cycle to cycle.", ""]
    l1, o1, D1 = table(names, "Main analysis, α 0.1 design at 30 MHz")
    l2, o2, D2 = table(CPRIME, "C′ (corner K1, 15 / 30 / 60 MHz, both particle sizes)")
    L += l1 + [""] + l2
    res = {"main": o1, "cprime": o2,
           "per_run": {n: [x["cls"] if x else None for x in v] for n, v in {**D1, **D2}.items()}}
    print("\n".join(L))
    (pathlib.Path(__file__).resolve().parent / "results" / "periodic.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
