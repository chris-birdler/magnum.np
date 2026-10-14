#!/usr/bin/env python
"""
Frequency sweep at 9 mT (decision Chris 2026-10-14, option B; PROTOKOLL 7.89; rules fixed before the runs).

Question: is the loss at 9 mT the low-frequency side of a damped eigen-motion whose frequency is lower at A 10 than at
A 20 (hypothesis of slide "What sets β")? Data: corner K1, A 20 (K1_k) and A 10 (K1A10_k), seeds k = 0 … 2, 9 mT,
f = 15 / 21 / 30 / 42 / 60 / 85 / 120 / 240 MHz (15 / 30 / 60 from the first stage of the C′ runs, the others from the
sweep runs *_swNNN), the same start state. Value: w = loss per cycle (w_loop / K_d, steady rule ≥ 2 calm cycles).
Method checks: the control runs K1_0_sw030 / K1A10_0_sw030 must give the w of the C′ runs at 30 MHz (relative
difference < 1e-3); energy balance of the calm cycles.
Per A: geometric mean of w over the seeds at each f (and the per-seed values); local exponent n = 1 + d ln w / d ln f
between neighbours (P ∝ f^n); the frequency of the largest w: f_max = vertex of a parabola in ln f through the largest
mean and its two neighbours (an end point: "at the edge"). Paired loss ratio P(A 10) / P(A 20) at each f.
Reading (fixed now):
  supports   w of A 10 has an interior maximum (both neighbours lower) and f_max(A 20) > f_max(A 10) or the A 20
             maximum lies at or above 240 MHz while w(A 20) still rises: an eigen-frequency that is lower at A 10;
  against    both w rise without an interior maximum up to 240 MHz with n ≤ 2 above 60 MHz (no resonance in the range),
             or the maxima lie at the same f within one grid step (the eigen-frequency does not depend on A);
  otherwise  not determined.

    python eval_sweep.py > results/sweep9.md      # also results/sweep9.json
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD

R = pathlib.Path(__file__).resolve().parent / "runs"
F_MHZ = (15.0, 21.2132, 30.0, 42.4264, 60.0, 84.8528, 120.0, 240.0)
SEEDS = range(3)


def folder(a, f, k):
    base = "K1A10_%d" % k if a == 10 else "K1_%d" % k
    if abs(f - 15.0) < 1e-6:
        return R / (base + "_fhalf")
    if abs(f - 30.0) < 1e-6:
        return R / base if a == 10 else R / (base + "_f1")
    if abs(f - 60.0) < 1e-6:
        return R / (base + "_f2")
    return R / ("%s_sw%03d" % (base, round(f)))


def w9(d):
    """w at 9 mT (the first measured stage, steady rule) or NaN."""
    if not (d / "DONE").exists():
        return float("nan")
    _, rows = analyze.load_run(d)
    rows = [r for r in rows if 1500.0 * r["b_peak"] < 15.0]
    return float(rows[0]["w_loop"]) if rows and np.isfinite(rows[0]["w_loop"]) and rows[0]["w_loop"] > 0 else float("nan")


def f_max(f, w):
    """Vertex of a parabola in ln f through the largest w and its neighbours; (f_max in MHz, at the edge?)."""
    i = int(np.nanargmax(w))
    if i in (0, len(w) - 1) or not np.all(np.isfinite(w[i - 1:i + 2])):
        return float(f[i]), True
    x, y = np.log(f[i - 1:i + 2]), np.log(w[i - 1:i + 2])
    c = np.polyfit(x, y, 2)
    return (float(math.exp(np.clip(-c[1] / (2 * c[0]), x[0], x[2]))) if c[0] < 0 else float(f[i])), False


def reading(f, m10, m20):
    """Pre-registered reading from the mean curves w(f) of A 10 and A 20."""
    f10, e10 = f_max(f, m10)
    f20, e20 = f_max(f, m20)
    step = math.log(f[1] / f[0])
    rise20_end = m20[-1] > m20[-2]
    if not e10 and ((not e20 and f20 > f10 and math.log(f20 / f10) > step) or (e20 and f20 >= f[-1] and rise20_end)):
        return "supports", f10, f20
    n_hi = [1 + math.log(m[j + 1] / m[j]) / math.log(f[j + 1] / f[j]) for m in (m10, m20)
            for j in range(len(f) - 1) if f[j] >= 60.0]
    if (e10 and e20 and all(n <= 2.0 for n in n_hi)) or (not e10 and not e20 and abs(math.log(f20 / f10)) <= step):
        return "against", f10, f20
    return "not determined", f10, f20


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    f = np.array(F_MHZ)
    W = {a: np.array([[w9(folder(a, x, k)) for x in F_MHZ] for k in SEEDS]) for a in (20, 10)}
    L = ["# Frequency sweep at 9 mT: corner K1, A 20 and A 10 pJ/m (PROTOKOLL 7.89)", "",
         "w = loss per cycle (steady rule ≥ 2 calm cycles), 3 seeds, the same start state. P ∝ f · w.", ""]
    # control
    ctl = []
    for a, base in ((20, "K1_0"), (10, "K1A10_0")):
        c, ref = w9(R / (base + "_sw030")), w9(folder(a, 30.0, 0))
        ctl.append((a, c, ref))
        L.append("Control A %d: w(sweep, 30 MHz) = %.6g, w(C′) = %.6g, relative difference %.2g"
                 % (a, c, ref, abs(c / ref - 1) if np.isfinite(c) and np.isfinite(ref) else float("nan")))
    L += ["", "| f (MHz) | w A 20 (geometric mean; seeds) | w A 10 | P(A 10) / P(A 20), paired (90 % CI) |",
          "|---|---|---|---|"]
    gm = {a: np.exp(np.nanmean(np.log(W[a]), 0)) for a in W}
    ratio = []
    for j, x in enumerate(F_MHZ):
        d = np.log(W[10][:, j]) - np.log(W[20][:, j])
        d = d[np.isfinite(d)]
        if len(d) >= 3:
            h = stats.t.ppf(0.95, len(d) - 1) * d.std(ddof=1) / math.sqrt(len(d))
            ratio.append((math.exp(d.mean()), math.exp(d.mean() - h), math.exp(d.mean() + h)))
            rc = "× %.2f (%.2f … %.2f)" % ratio[-1]
        else:
            ratio.append(None)
            rc = "n %d" % len(d)
        L.append("| %.0f | %.3g (%s) | %.3g (%s) | %s |" % (
            x, gm[20][j], " / ".join("%.3g" % v for v in W[20][:, j]), gm[10][j],
            " / ".join("%.3g" % v for v in W[10][:, j]), rc))
    L += ["", "Local exponent n of P ∝ f^n between neighbours (mean curves):", ""]
    for a in (20, 10):
        n = [1 + math.log(gm[a][j + 1] / gm[a][j]) / math.log(f[j + 1] / f[j]) for j in range(len(f) - 1)]
        L.append("- A %d: %s" % (a, " / ".join("%.1f" % x for x in n)))
    key, f10, f20 = reading(f, gm[10], gm[20])
    L += ["", "Frequency of the largest w: A 10 %.0f MHz%s, A 20 %.0f MHz%s." % (
        f10, " (at the edge)" if f_max(f, gm[10])[1] else "", f20, " (at the edge)" if f_max(f, gm[20])[1] else ""),
        "", "**Reading (pre-registered): %s**" % key]
    print("\n".join(L))
    (pathlib.Path(__file__).resolve().parent / "results" / "sweep9.json").write_text(json.dumps({
        "f_MHz": list(F_MHZ), "w": {str(a): W[a].tolist() for a in W}, "gm": {str(a): gm[a].tolist() for a in gm},
        "ratio": ratio, "f_max": {"10": f10, "20": f20}, "reading": key, "control": ctl}, indent=1))


if __name__ == "__main__":
    main()
