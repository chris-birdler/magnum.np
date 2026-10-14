#!/usr/bin/env python
"""
Frequency sweep at 9 mT (PROTOKOLL 7.91), option 1 of Chris (2026-10-14): is the steady state at 9 mT the same at all
frequencies? Local, no GPU. Question: do the dips of w(f) (A 10 at 42 MHz, A 20 at 85 MHz) and the seed scatter of A 20
at 42 MHz come with another steady state (state selection during the start-up), or with the same state?

Runs (corner K1, seeds k = 0 … 2, all from the same relaxed start state of the seed, 9 mT only): 15 / 30 / 60 MHz = the
first stage of C′ (K1_k_fhalf / K1_k_f1 / K1_k_f2, K1A10_k_fhalf / K1A10_k / K1A10_k_f2), 21 / 42 / 85 / 120 / 240 MHz =
the sweep (K1_k_swNNN, K1A10_k_swNNN). Kept cycles: acc.json "cycles_kept" (the calm cycles of the accumulators).

Markers per run (mean over the kept cycles; spread = largest minus smallest cycle mean of the kept cycles):
  Mp0 … Mp3  mean of M along the field in each of the 4 particles, / Ms (the domain state of each particle)
  Mx, My, Mz mean of the box magnetisation, / Ms
  n60        number of neighbour cell pairs with an angle above 60° (vortex cores, Bloch points, steep walls)
Start-up: number of event cycles (|w_dis / w_loop − 1| > 0.2) from cycle 1 (cycle 0 is the ramp from zero field) to
the first kept cycle, and the energy that the state gives away in them, Σ (W_dis − W_loop) / K_d (the unit of w).

Reading (fixed before the markers were looked at, 2026-10-14):
  Two runs of the same A and seed are in the same state if every marker differs by at most 3 × the larger of the two
  spreads + a floor (1e-4 for M / Ms, 1 for n60). Else they are in different states.
  state selection  at the dip frequency (A 10: 42 MHz, A 20: 85 MHz) the run differs from both neighbour frequencies
                   in at least 2 of 3 seeds, and at the other frequencies a run differs from both neighbours in at most
                   1 of 3 seeds;
  same state       at the dip frequency the run is in the same state as at least one neighbour in at least 2 of 3 seeds;
  otherwise        not determined.
The scatter of A 20 at 42 MHz is reported with the same markers (no reading).

    python eval_sweep_state.py > results/sweep9_state.md      # also results/sweep9_state.json
"""
import json
import math
import pathlib

import numpy as np

R = pathlib.Path(__file__).resolve().parent / "runs"
MU0 = 4e-7 * math.pi
F_MHZ = (15, 21, 30, 42, 60, 85, 120, 240)
CPRIME = {15: "fhalf", 30: "f", 60: "f2"}
SEEDS = range(3)
DIP = {10: 42, 20: 85}
EV = 0.2
FLOOR = {"M": 1e-4, "n60": 1.0}
MARKERS = ("Mp0", "Mp1", "Mp2", "Mp3", "Mx", "My", "Mz", "n60")


def folder(a, f, k):
    if f in CPRIME:
        fr = CPRIME[f]
        if a == 10:
            return R / ("K1A10_%d" % k if fr == "f" else "K1A10_%d_%s" % (k, fr))
        return R / ("K1_%d_%s" % (k, "f1" if fr == "f" else fr))
    return R / ("%s_%d_sw%03d" % ("K1A10" if a == 10 else "K1", k, f))


def run(a, f, k):
    d = folder(a, f, k)
    cfg = json.loads((d / "config.json").read_text())
    Ms = cfg["Js"] / MU0
    st = json.loads((d / "summary.json").read_text())["stages"][0]
    acc = json.loads((d / "acc.json").read_text())["stages"][0]
    kept = acc["harm"]["cycles_kept"]
    a_ = np.loadtxt(d / ("samples_00_%s.csv" % st["name"]), delimiter=",", comments="#")
    col = {"Mx": 4, "My": 5, "Mz": 6, "Mp0": 9, "Mp1": 10, "Mp2": 11, "Mp3": 12, "n60": 13}
    per = {m: [] for m in MARKERS}
    for c in kept:
        rows = a_[a_[:, 0] == c][1:]                                   # first row = end of the previous cycle
        for m, j in col.items():
            per[m].append(rows[:, j].mean() / (1.0 if m == "n60" else Ms))
    cyc = {c["cycle"]: c for c in st["cycles"]}
    first = min(kept)
    ev = [c for i, c in cyc.items() if 1 <= i < first and abs(c["w_dis"] / c["w_loop"] - 1) > EV]
    Kd = st["cycles"][0]["W_loop"] / st["cycles"][0]["w_loop"]
    return {"name": d.name, "kept": kept, "w": float(np.mean([cyc[c]["w_loop"] for c in kept])),
            "mean": {m: float(np.mean(v)) for m, v in per.items()},
            "spread": {m: float(np.max(v) - np.min(v)) for m, v in per.items()},
            "events": len(ev), "release": float(sum(c["W_dis"] - c["W_loop"] for c in ev) / Kd),
            "eps": float(np.mean([cyc[c]["w_dis"] / cyc[c]["w_loop"] - 1 for c in kept])),
            "closure": float(np.mean([cyc[c]["closure"] for c in kept]))}


def same(r1, r2):
    """(same state?, the marker with the largest difference relative to its limit, that ratio)."""
    worst, ratio = None, 0.0
    for m in MARKERS:
        lim = 3 * max(r1["spread"][m], r2["spread"][m]) + FLOOR["n60" if m == "n60" else "M"]
        q = abs(r1["mean"][m] - r2["mean"][m]) / lim
        if q > ratio:
            worst, ratio = m, q
    return ratio <= 1.0, worst, ratio


def main():
    D = {(a, f, k): run(a, f, k) for a in (20, 10) for f in F_MHZ for k in SEEDS}
    L = ["# Steady state at 9 mT across the frequency sweep (corner K1, 3 seeds)", "",
         "Markers: mean over the kept cycles. Same state: every marker within 3 × the larger spread + floor. Reading "
         "fixed before the markers were looked at (docstring).", ""]
    res = {"runs": {"%d|%d|%d" % key: v for key, v in D.items()}, "neighbours": {}}
    for a in (20, 10):
        L += ["## A %d pJ/m" % a, "", "| f (MHz) | seed | w (10⁻⁷) | start-up events (energy given away, 10⁻⁷) | n60 "
              "| Mp0 … Mp3 (10⁻³ Ms) | ⟨M⟩ (10⁻³ Ms) | same state as the lower / higher neighbour (largest ratio, "
              "marker) |", "|---|---|---|---|---|---|---|---|"]
        for i, f in enumerate(F_MHZ):
            for k in SEEDS:
                r = D[(a, f, k)]
                cmp = []
                for j in (i - 1, i + 1):
                    if 0 <= j < len(F_MHZ):
                        s, m, q = same(r, D[(a, F_MHZ[j], k)])
                        res["neighbours"]["%d|%d|%d|%d" % (a, f, F_MHZ[j], k)] = [s, m, q]
                        cmp.append("%s (%.1f, %s)" % ("same" if s else "differs", q, m))
                    else:
                        cmp.append("–")
                L.append("| %d | %d | %.2f | %d (%.2g) | %.1f | %s | %s | %s |" % (
                    f, k, r["w"] * 1e7, r["events"], r["release"] * 1e7 if r["events"] else 0.0, r["mean"]["n60"],
                    " / ".join("%.2f" % (1e3 * r["mean"]["Mp%d" % p]) for p in range(4)),
                    " / ".join("%.2f" % (1e3 * r["mean"][m]) for m in ("Mx", "My", "Mz")), " / ".join(cmp)))
        L.append("")
    # reading
    def differs_both(a, f, k):
        i = F_MHZ.index(f)
        return all(not res["neighbours"]["%d|%d|%d|%d" % (a, f, F_MHZ[j], k)][0] for j in (i - 1, i + 1)
                   if 0 <= j < len(F_MHZ))
    L += ["## Reading", ""]
    res["reading"] = {}
    for a in (20, 10):
        fd = DIP[a]
        n_dip = sum(differs_both(a, fd, k) for k in SEEDS)
        n_same = sum(any(res["neighbours"]["%d|%d|%d|%d" % (a, fd, F_MHZ[j], k)][0]
                         for j in (F_MHZ.index(fd) - 1, F_MHZ.index(fd) + 1)) for k in SEEDS)
        other = {f: sum(differs_both(a, f, k) for k in SEEDS) for f in F_MHZ if f != fd}
        if n_dip >= 2 and all(v <= 1 for v in other.values()):
            rd = "state selection"
        elif n_same >= 2:
            rd = "same state"
        else:
            rd = "not determined"
        res["reading"][str(a)] = {"dip_MHz": fd, "differs_both_at_dip": n_dip, "same_as_a_neighbour_at_dip": n_same,
                                  "differs_both_elsewhere": other, "reading": rd}
        L.append("- A %d (dip at %d MHz): differs from both neighbours in %d of 3 seeds; same as a neighbour in %d of 3. "
                 "Elsewhere (seeds that differ from both neighbours): %s. **%s**"
                 % (a, fd, n_dip, n_same, ", ".join("%d MHz %d" % kv for kv in other.items()), rd))
    # added after the reading (not part of it): each run against the 30 MHz state of the same A and seed
    L += ["", "## Added after the reading: each run against the state at 30 MHz (same A and seed)", "",
          "Entry: same / differs (largest ratio, marker). Not part of the fixed reading.", "",
          "| A | seed | " + " | ".join("%d MHz" % f for f in F_MHZ) + " |", "|---|---|" + "---|" * len(F_MHZ)]
    res["vs30"] = {}
    for a in (20, 10):
        for k in SEEDS:
            cells = []
            for f in F_MHZ:
                sm, m, q = same(D[(a, f, k)], D[(a, 30, k)])
                res["vs30"]["%d|%d|%d" % (a, f, k)] = [sm, m, q]
                cells.append("–" if f == 30 else "%s (%.1f, %s)" % ("same" if sm else "differs", q, m))
            L.append("| %d | %d | %s |" % (a, k, " | ".join(cells)))
    # added after the reading (not part of it): is w(f) of A 20 the curve of A 10 shifted in frequency?
    # Grid positions in steps of √2 (15 MHz = 0 … 120 MHz = 6, 240 MHz = 8). For a shift of s steps: SD over the common
    # positions of ln w20(p + s) − ln w10(p) (geometric means of the seeds; a constant factor is allowed), and the SD of
    # the same difference for each seed pair (the pair has the same seed, stress axes and packing).
    pos = {15: 0, 21: 1, 30: 2, 42: 3, 60: 4, 85: 5, 120: 6, 240: 8}
    lw = {(a, f): np.log([D[(a, f, k)]["w"] for k in SEEDS]) for a in (20, 10) for f in F_MHZ}
    L += ["", "## Added after the reading: w(f) of A 20 against the curve of A 10 shifted in frequency", "",
          "Shift s in steps of √2 (s 2: factor 2). SD of ln w20(√2^s f) − ln w10(f) over the common frequencies: mean "
          "curves; per seed pair (mean ± SD of 3). Small SD: the curves have the same shape after the shift. "
          "Factor: exp of the mean difference.", "",
          "| shift s (factor) | common frequencies | SD, mean curves | SD, per seed pair | factor w20 / w10 |",
          "|---|---|---|---|---|"]
    res["shift"] = {}
    inv = {v: f for f, v in pos.items()}
    for sh in range(0, 4):
        common = [f for f in F_MHZ if pos[f] + sh in inv]
        dm = np.array([lw[(20, inv[pos[f] + sh])].mean() - lw[(10, f)].mean() for f in common])
        dk = np.array([[lw[(20, inv[pos[f] + sh])][k] - lw[(10, f)][k] for f in common] for k in SEEDS])
        sdk = dk.std(axis=1, ddof=1)
        res["shift"][str(sh)] = {"n": len(common), "sd_mean": float(dm.std(ddof=1)), "sd_seed": sdk.tolist(),
                                 "factor": float(math.exp(dm.mean()))}
        L.append("| %d (× %.2f) | %d | %.2f | %.2f ± %.2f | × %.2f |" % (sh, math.sqrt(2) ** sh, len(common),
                                                                     dm.std(ddof=1), sdk.mean(), sdk.std(ddof=1),
                                                                     math.exp(dm.mean())))
    print("\n".join(L))
    (pathlib.Path(__file__).resolve().parent / "results" / "sweep9_state.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
