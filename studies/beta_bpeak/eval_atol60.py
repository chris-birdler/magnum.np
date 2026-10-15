#!/usr/bin/env python
"""
Long tolerance test at 60 MHz (option 2 of Chris 2026-10-14; reading fixed in PROTOKOLL 7.100 before the runs).
Runs: K1_k_sw060L6 (atol 1e-6) and K1_k_sw060L5 (atol 1e-5), k = 0 … 2: the small particle (A 20 runs) at 60 MHz,
9 mT only, exactly 40 cycles. Reference: the 9 mT stage of K1_k_f2 (C′, atol 1e-5, stop rule). L5 starts from
K1_k/init.pt; L6 relaxes again with atol 1e-6 from the same seed (atol is an INIT_KEY, PROTOKOLL 7.101).

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


def start_diff(k, tag):
    """Start control (reported only, PROTOKOLL 7.101): L6 relaxes itself, L5 starts from K1_k/init.pt. Relative
    difference of the relaxed energy E/K_d and of w_loop of cycle 0 against L5."""
    if tag == "L5":
        return "reference"
    s6, s5 = [json.loads((R / ("K1_%d_sw060%s" % (k, t)) / "summary.json").read_text()) for t in ("L6", "L5")]
    dE = s6["init"]["E_over_Kd"] / s5["init"]["E_over_Kd"] - 1
    dw = s6["stages"][0]["cycles"][0]["w_loop"] / s5["stages"][0]["cycles"][0]["w_loop"] - 1
    return "%+.1e / %+.1e" % (dE, dw)


def exploratory(res):
    """Exploratory, added after the readings (2026-10-15, PROTOKOLL 7.102): (1) the cycle-by-cycle difference L6
    against L5 (does it grow with the cycle number?); (2) the cycle from which the running mean of P (from cycle 2)
    stays within ±5 % of the mean of the last 10 cycles; (3) the loss channels R_k, the rest and two zone shares of L6,
    L5 and C′ K1_k_f2 (C′ kept other cycles: a window difference, not a tolerance difference)."""
    L = ["", "## Exploratory (not pre-registered; added 2026-10-15 after the readings)", "",
         "| seed | largest abs. w_loop(L6) / w_loop(L5) − 1 in cycles 0–9 / 10–19 / 20–29 / 30–39 | largest abs. closure "
         "ratio − 1 in the same blocks | running mean of P within ±5 % from cycle (L6 / L5) |", "|---|---|---|---|"]
    ex = res["exploratory"] = {"trajectory": {}, "channels": {}}
    for k in SEEDS:
        s6, s5 = [json.loads((R / ("K1_%d_sw060%s" % (k, t)) / "summary.json").read_text())["stages"][0]["cycles"]
                  for t in ("L6", "L5")]
        w6, w5 = np.array([c["w_loop"] for c in s6]), np.array([c["w_loop"] for c in s5])
        c6, c5 = np.array([c["closure"] for c in s6]), np.array([c["closure"] for c in s5])
        blk = [slice(0, 10), slice(10, 20), slice(20, 30), slice(30, 40)]
        dw = [float(np.max(np.abs(w6[b] / w5[b] - 1))) for b in blk]
        dc = [float(np.max(np.abs(c6[b] / np.maximum(c5[b], 1e-30) - 1))) for b in blk]

        def settle(w):
            ref = w[-LAST:].mean()
            run = np.cumsum(w[FIRST:]) / np.arange(1, len(w) - FIRST + 1)
            bad = np.nonzero(np.abs(run / ref - 1) > 0.05)[0]
            return int(bad[-1] + 1 + FIRST) if len(bad) else FIRST
        ex["trajectory"][k] = {"dw_blocks": dw, "dclosure_blocks": dc, "settle_L6": settle(w6), "settle_L5": settle(w5)}
        L.append("| %d | %s | %s | %d / %d |" % (k, " / ".join("%.1e" % x for x in dw), " / ".join("%.1e" % x for x in dc),
                                                 settle(w6), settle(w5)))
    L += ["", "| seed | run | cycles kept | R1 | R2 | R3 | R4 … R7 | rest | loss share, texture top 10 % | loss share, "
          "r ≥ 0.8 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for k in SEEDS:
        for name, lab in (("K1_%d_sw060L6" % k, "L6"), ("K1_%d_sw060L5" % k, "L5"), ("K1_%d_f2" % k, "C′ (window)")):
            st = json.loads((R / name / "acc.json").read_text())["stages"][0]
            h, wh = st["harm"], st["where"]
            ex["channels"][name] = {"R_k": h["R_k"], "rest": h["R_rest"], "cycles": len(h["cycles_kept"]),
                                    "texture": wh["texture top 10%"]["p_share"], "rim": wh["r>=0.8"]["p_share"]}
            L.append("| %d | %s | %d | %.3f | %.3f | %.3f | %.3f | %+.3f | %.3f | %.3f |" % (
                k, lab, len(h["cycles_kept"]), h["R_k"][0], h["R_k"][1], h["R_k"][2], sum(h["R_k"][3:]), h["R_rest"],
                wh["texture top 10%"]["p_share"], wh["r>=0.8"]["p_share"]))
    L += ["", "rest = 1 − Σ R_k; a negative rest is the per-cell over-completeness of the harmonic sum (7.90)."]
    return L


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    res, L = {}, ["# Long tolerance test at 60 MHz, 9 mT, small particle (reading fixed in PROTOKOLL 7.100)", "",
                  "| seed | atol | class | c_settle | P last 10 cycles (W/cm³, mean ± SD) | C′ value (W/cm³) | start: E/K_d "
                  "and w_loop of cycle 0, relative to L5 |", "|---|---|---|---|---|---|---|"]
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
            L.append("| %d | %s | %s | %s | %.1f ± %.1f | %.1f | %s |" % (k, tol, c, s if s is not None else "–", p.mean(),
                                                                         p.std(ddof=1), ref[-1], start_diff(k, tag)))
    rd, vr = reading(r["L6"], r["L5"]), value_reading(pv["L6"], pv["L5"], ref)
    res["reading"], res["value_reading"] = rd, vr
    L += ["", "**Reading (motion): %s**" % rd, "", "**Reading (loss value): %s**" % vr]
    L += exploratory(res)
    print("\n".join(L))
    (pathlib.Path(__file__).resolve().parent / "results" / "atol60.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
