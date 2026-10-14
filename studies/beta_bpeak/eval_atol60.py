#!/usr/bin/env python
"""
Long tolerance test at 60 MHz (option 2 of Chris 2026-10-14; reading fixed in PROTOKOLL 7.100 before the runs).
Runs: K1_k_sw060L6 (atol 1e-6) and K1_k_sw060L5 (atol 1e-5), k = 0 … 2: the small particle (A 20 runs) at 60 MHz,
9 mT only, exactly 40 cycles. Reference: the 9 mT stage of K1_k_f2 (C′, atol 1e-5, stop rule).

Per run, from cycle 2 on: for each window of 3 cycles ending at cycle c: RN (share of the box dissipation that changes
from cycle to cycle), mean closure, and w(c). c_settle = first cycle from which every later window has RN ≤ 0.05,
closure ≤ 0.001 and w within ±5 % of the mean of the last 10 cycles. Class: "periodic" if c_settle exists, "stays
irregular" if the last 10 cycles have median RN > 0.10 or median closure > 0.002, else "marginal".
Reading over the 3 seeds:
  tolerance artefact  with 1e-6 periodic with c_settle ≤ 12 in ≥ 2 of 3 seeds, with 1e-5 not periodic by cycle 12 in
                      ≥ 2 of 3 seeds;
  settles slowly      both tolerances periodic in ≥ 2 of 3 seeds, and c_settle > 12 in ≥ 2 of 3 seeds for one of them;
  stays irregular     both tolerances "stays irregular" in ≥ 2 of 3 seeds;
  otherwise           not determined.
Loss value: P = mean of the last 10 cycles (W/cm³); "the C′ value holds" if within ±10 % of the 9 mT value of K1_k_f2
for both tolerances in ≥ 2 of 3 seeds, "the C′ value is off" if outside ±10 % in ≥ 2 of 3 seeds, else not determined.

    python eval_atol60.py > runs/atol60_report.md      # also results/atol60.json
"""
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import eval_sweep as SW

R = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = range(3)
FIRST = 2
LAST = 10
K_D = 1.5 ** 2 / (2 * 4e-7 * math.pi) / 1e6          # J/cm³
F_HZ = 60e6


def windows(rn, clo, w, last=LAST):
    """(c_settle or None, class) from per-window RN and closure and per-cycle w (all indexed by the end cycle of the
    window, from FIRST on)."""
    rn, clo, w = np.asarray(rn, float), np.asarray(clo, float), np.asarray(w, float)
    ref = w[-last:].mean()
    ok = (rn <= 0.05) & (clo <= 0.001) & (np.abs(w / ref - 1.0) <= 0.05)
    settle = None
    for i in range(len(ok)):
        if ok[i:].all():
            settle = i
            break
    if settle is not None:
        return settle + FIRST, "periodic"
    if np.median(rn[-last:]) > 0.10 or np.median(clo[-last:]) > 0.002:
        return None, "stays irregular"
    return None, "marginal"


def reading(r6, r5):
    """r6, r5: per seed (c_settle, class) for atol 1e-6 and 1e-5."""
    per6 = [c == "periodic" and s is not None and s <= 12 for s, c in r6]
    not5 = [not (c == "periodic" and s is not None and s <= 12) for s, c in r5]
    if sum(per6) >= 2 and sum(not5) >= 2:
        return "tolerance artefact"
    both_per = sum(c == "periodic" for _, c in r6) >= 2 and sum(c == "periodic" for _, c in r5) >= 2
    slow = any(sum(c == "periodic" and s is not None and s > 12 for s, c in rr) >= 2 for rr in (r6, r5))
    if both_per and slow:
        return "settles slowly"
    if sum(c == "stays irregular" for _, c in r6) >= 2 and sum(c == "stays irregular" for _, c in r5) >= 2:
        return "stays irregular"
    return "not determined"


def value_reading(p6, p5, ref):
    hold = [abs(a / r - 1) <= 0.10 and abs(b / r - 1) <= 0.10 for a, b, r in zip(p6, p5, ref)]
    off = [abs(a / r - 1) > 0.10 or abs(b / r - 1) > 0.10 for a, b, r in zip(p6, p5, ref)]
    if sum(hold) >= 2:
        return "the C′ value holds"
    if sum(off) >= 2:
        return "the C′ value is off"
    return "not determined"


def run_series(d):
    """Per-window RN and closure and per-cycle w from cycle FIRST on."""
    st = json.loads((d / "summary.json").read_text())["stages"][0]
    cyc = {c["cycle"]: c for c in st["cycles"]}
    a = np.loadtxt(sorted(d.glob("samples_*_amp*.csv"))[0], delimiter=",", comments="#")      # the 9 mT stage
    P = {c: a[a[:, 0] == c][1:, 8] for c in cyc}
    n = int(np.bincount([len(v) for v in P.values()]).argmax())
    rn, clo, w = [], [], []
    for c in sorted(cyc):
        if c < FIRST:
            continue
        win = [c - 2, c - 1, c]
        if not all(x in P and len(P[x]) == n for x in win):
            continue
        X = np.array([P[x] for x in win])
        rn.append(float(np.abs(X - X.mean(0)).sum() / X.sum()))
        clo.append(float(np.mean([cyc[x]["closure"] for x in win])))
        w.append(cyc[c]["w_loop"])
    return rn, clo, w


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    res, L = {}, ["# Long tolerance test at 60 MHz, 9 mT, small particle (reading fixed in PROTOKOLL 7.100)", "",
                  "| seed | atol | class | c_settle | P last 10 cycles (W/cm³, mean ± SD) | C′ value (W/cm³) |",
                  "|---|---|---|---|---|---|"]
    r = {"L6": [], "L5": []}
    pv = {"L6": [], "L5": []}
    ref = []
    for k in SEEDS:
        ref.append(F_HZ * K_D * SW.w9(R / ("K1_%d_f2" % k)))
        for tag, tol in (("L6", "1e-6"), ("L5", "1e-5")):
            d = R / ("K1_%d_sw060%s" % (k, tag))
            rn, clo, w = run_series(d)
            s, c = windows(rn, clo, w)
            p = F_HZ * K_D * np.array(w[-LAST:])
            r[tag].append((s, c))
            pv[tag].append(float(p.mean()))
            res["%d|%s" % (k, tag)] = {"class": c, "c_settle": s, "P_mean": float(p.mean()), "P_sd": float(p.std(ddof=1)),
                                       "RN": rn, "closure": clo, "w": w}
            L.append("| %d | %s | %s | %s | %.1f ± %.1f | %.1f |" % (k, tol, c, s if s is not None else "–", p.mean(),
                                                                    p.std(ddof=1), ref[-1]))
    rd, vr = reading(r["L6"], r["L5"]), value_reading(pv["L6"], pv["L5"], ref)
    res["reading"], res["value_reading"] = rd, vr
    L += ["", "**Reading (motion): %s**" % rd, "", "**Reading (loss value): %s**" % vr]
    print("\n".join(L))
    (pathlib.Path(__file__).resolve().parent / "results" / "atol60.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
