#!/usr/bin/env python
"""
Are the events in the measured cycles part of the steady physics, or the end of a transient (PROTOKOLL 7.55)?
Event cycle: w_loop <= 0 or |w_dis / w_loop - 1| > 0.20 (7.54). Data: summary.json of the 116 α 0.1 runs, 7 amplitudes,
cycles 1 … 6 (cycle 0 holds the drive ramp and the drive correction).

Energy per cycle: the stored energy changes by ΔE = W_loop − W_dis (input work minus dissipation). ΔE < 0: the state
gives energy away (it goes down to a lower minimum). In a steady but not periodic state (intermittent jumps), the
stored energy goes down and up again: over many cycles ΣΔE ≈ 0 and uptake cycles (ΔE > 0) are as frequent as release
cycles. In a transient, ΔE < 0 dominates and the event rate decreases with the cycle number.

Tests
  T1  event rate against the cycle number (1 … 6), per amplitude group
  T2  sign of the energy change in the event cycles (release w_dis > w_loop, or uptake), and ΣΔE over cycles 2 … 6
      relative to the loop loss
  T3  after an event (cycles 2 … 5, one cycle before and after): change of the loop loss and of n60 (number of
      neighbour pairs > 60°: vortex cores, Bloch points). Lower and staying lower: the state went to a new, simpler
      minimum.
  T4  closure |M(end) − M(start)| / M_peak of the calm cycles: a periodic orbit has closure at the numerical floor
      (1.5e-7 in a small test run, PROTOKOLL 7.52); larger values mean that the state still moves from cycle to cycle
  T5  late events (cycle >= 4 after >= 2 calm cycles): closure of the 2 calm cycles before, against the calm cycles of
      run-amplitudes without any event in cycles 2 … 6

    python eval_events.py > results/events.md
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze_design as AD

HERE = pathlib.Path(__file__).resolve().parent
R = HERE / "runs"
TOL = AD.BAL_EVENT
GROUPS = (("9 mT", (0,)), ("20 … 35 mT", (1, 2)), ("50 … 150 mT", (3, 4, 5, 6)))


def is_event(c):
    return c["w_loop"] <= 0 or abs(c["w_dis"] / c["w_loop"] - 1.0) > TOL


def main():
    recs, _ = AD.load_runs(R)
    stages = {r["name"]: json.loads((R / r["name"] / "summary.json").read_text())["stages"] for r in recs}
    print("# Events in the measured cycles: physics or transient? (exploratory)\n")
    print("n = %d α 0.1 runs × 7 amplitudes; event cycle: w_loop ≤ 0 or |w_dis / w_loop − 1| > %.2f.\n" % (len(recs), TOL))

    print("## T1. Event rate against the cycle number (share of the run-amplitudes)\n")
    print("| amplitudes | n | " + " | ".join("cycle %d" % c for c in range(1, 7)) + " |")
    print("|---|---|" + "---|" * 6)
    for lab, ks in GROUPS:
        rate = [np.mean([is_event(stages[n][k]["cycles"][c]) for n in stages for k in ks]) for c in range(1, 7)]
        print("| %s | %d | %s |" % (lab, len(stages) * len(ks), " | ".join("%.1f %%" % (100 * x) for x in rate)))

    print("\n## T2. Energy sign of the event cycles (cycles 2 … 6) and ΣΔE\n")
    print("| amplitudes | event cycles | release (w_dis > w_loop) | uptake (w_dis < w_loop) | ΣΔE over cycles 2 … 6, "
          "median / loop loss (run-amplitudes with an event) | the same without an event |")
    print("|---|---|---|---|---|---|")
    for lab, ks in GROUPS:
        rel = upt = 0
        s_ev, s_calm = [], []
        for n in stages:
            for k in ks:
                cyc = stages[n][k]["cycles"][2:]
                wl = np.median([c["w_loop"] for c in cyc])
                ev = [c for c in cyc if is_event(c)]
                rel += sum(c["w_loop"] <= 0 or c["w_dis"] > c["w_loop"] for c in ev)
                upt += sum(c["w_loop"] > 0 and c["w_dis"] < c["w_loop"] for c in ev)
                dE = sum(c["w_loop"] - c["w_dis"] for c in cyc) / wl if wl > 0 else float("nan")
                (s_ev if ev else s_calm).append(dE)
        print("| %s | %d | %d | %d | %+.2f (n %d) | %+.3f (n %d) |" % (lab, rel + upt, rel, upt, np.nanmedian(s_ev),
                                                                       len(s_ev), np.nanmedian(s_calm), len(s_calm)))

    print("\n## T3. Before and after an event (event in cycle 2 … 5; one cycle before and after)\n")
    print("| amplitudes | events | loop loss after / before: median | share lower after | n60 after − before: median "
          "| share with fewer n60 after | share with more n60 after |")
    print("|---|---|---|---|---|---|---|")
    for lab, ks in GROUPS:
        rw, dn = [], []
        for n in stages:
            for k in ks:
                cyc = stages[n][k]["cycles"]
                for c in range(2, 6):
                    if is_event(cyc[c]) and cyc[c - 1]["w_loop"] > 0 and cyc[c + 1]["w_loop"] > 0:
                        rw.append(cyc[c + 1]["w_loop"] / cyc[c - 1]["w_loop"])
                        dn.append(cyc[c + 1]["n_pairs_gt60"] - cyc[c - 1]["n_pairs_gt60"])
        rw, dn = np.array(rw), np.array(dn)
        print("| %s | %d | %.2f | %.0f %% | %+.0f | %.0f %% | %.0f %% |" % (lab, len(rw), np.median(rw), 100 * np.mean(rw < 1),
                                                                          np.median(dn), 100 * np.mean(dn < 0),
                                                                          100 * np.mean(dn > 0)))
    calm_dn = []
    for n in stages:
        for k in range(7):
            cyc = stages[n][k]["cycles"]
            if not any(is_event(c) for c in cyc[1:]):
                calm_dn += [cyc[c + 1]["n_pairs_gt60"] - cyc[c - 1]["n_pairs_gt60"] for c in range(2, 6)]
    calm_dn = np.array(calm_dn)
    print("\nReference, run-amplitudes without an event in cycles 1 … 6 (same cycle pairs): n60 after − before median "
          "%+.0f, fewer %.0f %%, more %.0f %% (n %d).\n" % (np.median(calm_dn), 100 * np.mean(calm_dn < 0),
                                                           100 * np.mean(calm_dn > 0), len(calm_dn)))

    print("## T4. Closure of the calm cycles (cycles 2 … 6, no event)\n")
    print("| amplitudes | calm cycles | closure median | 90 % quantile | share > 1e-4 |")
    print("|---|---|---|---|---|")
    for lab, ks in GROUPS:
        cl = np.array([c["closure"] for n in stages for k in ks for c in stages[n][k]["cycles"][2:] if not is_event(c)])
        print("| %s | %d | %.1e | %.1e | %.0f %% |" % (lab, len(cl), np.median(cl), np.percentile(cl, 90),
                                                      100 * np.mean(cl > 1e-4)))

    print("\n## T5. Late events: closure of the 2 calm cycles before, against calm run-amplitudes\n")
    late, ref = [], []
    for n in stages:
        for k in range(7):
            cyc = stages[n][k]["cycles"]
            evs = [is_event(c) for c in cyc]
            if not any(evs[2:]):
                ref += [c["closure"] for c in cyc[4:]]
            for c in range(4, 7):
                if evs[c] and not evs[c - 1] and not evs[c - 2]:
                    late.append((n, AD.AMPS_MT[k], c, cyc[c - 2]["closure"], cyc[c - 1]["closure"],
                                 cyc[c]["w_dis"] / cyc[c]["w_loop"] if cyc[c]["w_loop"] > 0 else float("nan")))
    lc = np.array([x[3:5] for x in late]).ravel()
    ref = np.array(ref)
    print("Late events: %d. Closure of the 2 calm cycles before: median %.1e (90 %% quantile %.1e). Calm run-amplitudes, "
          "cycles 4 … 6: median %.1e (90 %% quantile %.1e), n %d.\n" % (len(late), np.median(lc), np.percentile(lc, 90),
                                                                       np.median(ref), np.percentile(ref, 90), len(ref)))
    mean_late = np.array([0.5 * (x[3] + x[4]) for x in late])
    print("Mann-Whitney (mean closure of the 2 cycles before each late event, n %d, against the calm cycles): p = %.1e.\n"
          % (len(mean_late), stats.mannwhitneyu(mean_late, ref, alternative="greater").pvalue))
    print("| run | B (mT) | cycle | closure 2 cycles before | 1 cycle before | balance of the event |")
    print("|---|---|---|---|---|---|")
    for x in sorted(late, key=lambda x: (x[1], x[0])):
        print("| %s | %g | %d | %.1e | %.1e | %.2f |" % x)

    print("\n## T6. Mesh pairs: event rate in the cycles 1 … 4 at dx 3 and dx 1.5 (K3 and base, 9 / 50 / 100 mT)\n")
    print("| mesh | jobs | " + " | ".join("cycle %d" % c for c in range(1, 5)) + " | events that release energy |")
    print("|---|---|" + "---|" * 4 + "---|")
    for dx in ("3", "15"):
        jobs = [json.loads((d / "summary.json").read_text())["stages"][0]["cycles"]
                for d in sorted(R.glob("M*_dx%s_s*_b*" % dx)) if (d / "DONE").exists() and d.name.split("_b")[1] in ("9", "50", "100")]
        rate = [np.mean([is_event(cy[c]) for cy in jobs]) for c in range(1, 5)]
        ev = [cy[c] for cy in jobs for c in range(1, 5) if is_event(cy[c])]
        rel = sum(c["w_loop"] <= 0 or c["w_dis"] > c["w_loop"] for c in ev)
        print("| dx %s | %d | %s | %d of %d |" % ("1.5" if dx == "15" else dx, len(jobs),
                                                " | ".join("%.0f %%" % (100 * x) for x in rate), rel, len(ev)))


if __name__ == "__main__":
    main()
