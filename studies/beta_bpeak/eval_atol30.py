#!/usr/bin/env python
"""
Tolerance twins of the main run at 30 MHz (decision Chris 2026-10-14, option 1 with patch; reading fixed in PROTOKOLL
7.101 before the runs). Runs K1_k_L6 (k = 0 … 2): the production line of K1_k with atol 1e-6 and an own relaxation with
atol 1e-6 (atol is an INIT_KEY, 7.88). Reference: the main runs K1_k (atol 1e-5, the same line).

Question: is the loop non-closure |M(end) − M(start)| / M_peak of the main run physical, or an error of the time
integrator?
Per run: c = median of the closure over the measured cycles (cycles 2 … 6 of all 7 stages, no calm filter).
Per seed k: r_k = c(K1_k_L6) / c(K1_k).
Reading: "physical" if r_k ≥ 0.5 in at least 2 of 3 seeds; "numerical" if r_k < 0.5 in at least 2 of 3 seeds; else
"not determined" (a missing or incomplete run, or c(K1_k) = 0).
Controls (reported, not part of the reading): (1) the relaxed start state of K1_k_L6 against K1_k (E / K_d, M / M_s);
(2) w_loop of cycle 0 at 9 mT: "same start" if |w(L6) / w(K1_k) − 1| ≤ 1e-3, else the comparison is statistical only.
Also reported per B: the median closure, the median w_dis / w_loop of the measured cycles (atol 1e-5 makes w_dis too
large, 7.88 and 7.101) and P = f K_d w_loop (mean of the measured cycles, W/cm³).

    python eval_atol30.py > runs/atol30_report.md      # also results/atol30.json
"""
import json
import math
import pathlib

import numpy as np

import analyze_design as AD

R = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = range(3)
AMPS = list(AD.AMPS_MT)
K_D = 1.5 ** 2 / (2 * 4e-7 * math.pi) / 1e6          # J/cm³
F_HZ = 30e6
CTRL_TOL = 1e-3


def measured(st):
    """The cycles that analyze.py measures: from 2 on after a drive correction, else from 1 on."""
    return st["cycles"][2:] if "h_amp_initial" in st else st["cycles"][1:]


def load(d):
    """(measured stages, init record) of a finished run with all 7 stages complete, else None."""
    if not (d / "DONE").exists() or not (d / "summary.json").exists():
        return None
    s = json.loads((d / "summary.json").read_text())
    st = [x for x in s["stages"] if x.get("measure")]
    if len(st) != len(AMPS) or not all(x.get("complete") and x.get("cycles") for x in st):
        return None
    return st, s.get("init", {})


def closure_median(st):
    return float(np.median([c["closure"] for s in st for c in measured(s)]))


def reading(r):
    """r: per seed the ratio c(L6) / c(K1_k), None if not available."""
    ok = [x for x in r if x is not None and math.isfinite(x)]
    if sum(x >= 0.5 for x in ok) >= 2:
        return "physical"
    if sum(x < 0.5 for x in ok) >= 2:
        return "numerical"
    return "not determined"


def per_b(st):
    """Per stage: median closure, median w_dis / w_loop, P (W/cm³) of the measured cycles."""
    out = []
    for s in st:
        cyc = measured(s)
        out.append({"closure": float(np.median([c["closure"] for c in cyc])),
                    "wd_w": float(np.median([c["w_dis"] / c["w_loop"] for c in cyc if c["w_loop"] > 0])),
                    "P": float(F_HZ * K_D * np.mean([c["w_loop"] for c in cyc]))})
    return out


def main():
    res, ratios = {"seeds": {}}, []
    L = ["# Tolerance twins of the main run at 30 MHz (reading fixed in PROTOKOLL 7.101)", "",
         "c: median closure |M(end) − M(start)| / M_peak over the measured cycles (2 … 6) of all 7 stages.", "",
         "| seed | c atol 1e-5 (K1_k) | c atol 1e-6 (K1_k_L6) | r | start E/K_d rel. diff. | start max abs. diff. "
         "of M/M_s | w_loop cycle 0 at 9 mT, L6 / 1e-5 − 1 | control |", "|---|---|---|---|---|---|---|---|"]
    tabs = {}
    for k in SEEDS:
        a, b = load(R / ("K1_%d" % k)), load(R / ("K1_%d_L6" % k))
        if a is None or b is None:
            ratios.append(None)
            L.append("| %d | %s | %s | – | – | – | – | – |" % (k, "missing" if a is None else "%.6f" % closure_median(a[0]),
                                                             "missing" if b is None else "%.6f" % closure_median(b[0])))
            res["seeds"][k] = None
            continue
        c5, c6 = closure_median(a[0]), closure_median(b[0])
        r = c6 / c5 if c5 > 0 else None
        ratios.append(r)
        i5, i6 = a[1], b[1]
        dE = abs(i6["E_over_Kd"] / i5["E_over_Kd"] - 1) if i5.get("E_over_Kd") else float("nan")
        dM = float(np.max(np.abs(np.subtract(i6.get("M_over_Ms", [np.nan] * 3), i5.get("M_over_Ms", [np.nan] * 3)))))
        dw = a[0][0]["cycles"][0]["w_loop"]
        dw = b[0][0]["cycles"][0]["w_loop"] / dw - 1 if dw > 0 else float("nan")
        ctrl = "same start" if abs(dw) <= CTRL_TOL else "start differs (statistical comparison)"
        tabs[k] = (per_b(a[0]), per_b(b[0]))
        res["seeds"][k] = {"c5": c5, "c6": c6, "r": r, "start_dE_rel": dE, "start_dM": dM, "w0_rel": dw,
                           "control": ctrl, "per_b_1e5": tabs[k][0], "per_b_1e6": tabs[k][1]}
        L.append("| %d | %.6f | %.6f | %s | %.2e | %.2e | %+.2e | %s |" % (k, c5, c6, "%.2f" % r if r is not None else "–",
                                                                          dE, dM, dw, ctrl))
    rd = reading(ratios)
    res["reading"] = rd
    L += ["", "**Reading: %s** (r ≥ 0.5 physical, r < 0.5 numerical; at least 2 of 3 seeds)" % rd, ""]
    if tabs:
        L += ["Per B, mean over the %d complete seeds (reported, not part of the reading):" % len(tabs), "",
              "| B (mT) | median closure 1e-5 | median closure 1e-6 | median w_dis / w_loop 1e-5 | median w_dis / w_loop "
              "1e-6 | P 1e-5 (W/cm³) | P 1e-6 (W/cm³) | P ratio 1e-6 / 1e-5 |", "|---|---|---|---|---|---|---|---|"]
        for i, bm in enumerate(AMPS):
            g = lambda j, key: float(np.mean([tabs[k][j][i][key] for k in tabs]))
            pr = float(np.mean([tabs[k][1][i]["P"] / tabs[k][0][i]["P"] for k in tabs]))
            L.append("| %g | %.6f | %.6f | %.4f | %.4f | %.3g | %.3g | %.3f |" % (
                bm, g(0, "closure"), g(1, "closure"), g(0, "wd_w"), g(1, "wd_w"), g(0, "P"), g(1, "P"), pr))
    print("\n".join(L))
    (pathlib.Path(__file__).resolve().parent / "results" / "atol30.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
