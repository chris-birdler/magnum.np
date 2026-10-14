#!/usr/bin/env python
"""
Frequency test C′ (PLAN 4.8, analysis fixed before the runs; PROTOKOLL 7.72, 7.75, 7.83).

Runs: corner K1 at A 20 pJ/m (K1_k_fhalf / K1_k_f1 / K1_k_f2) and at A 10 pJ/m (K1A10_k_fhalf / K1A10_k (the pair runs
at f) / K1A10_k_f2), seeds k = 0 … 2 (2011 … 2013), f/2 = 15 MHz, f = 30 MHz, 2 f = 60 MHz, the same sample interval
T(f)/512, harmonics k ≤ 7 in the accumulators, steady rule ≥ 2 calm cycles (≤ 12 cycles).
Values: ln w at the 7 nominal amplitudes (analyze_design.run_record), w = loss per cycle (w_loop / K_d).
1. Per A and B: w(f/2) / w(f) and w(2 f) / w(f), paired by seed (geometric mean, 90 % CI with t(n − 1), n).
   A loss per cycle that does not depend on f gives 1 (rate-independent); a loss ∝ f per cycle gives 0.5 / 2.
2. Rate-independent share w_h / w(f) with w = w_h + c · f: from f and 2 f: 2 − w(2f)/w(f); from f/2 and f:
   2 · w(f/2)/w(f) − 1 (mean over the seeds).
3. β (slope of ln w over ln B, 9 … 150 mT) per A and frequency (per seed, mean ± SD) and the A effect on β
   β(A 20) − β(A 10) at each frequency, paired by seed (90 % CI with t). Question: does the A effect on β come from
   the rate-dependent part? Also the A effect on β of w_h = 2 w(f) − w(2 f) where w_h > 0 at all B.
4. Harmonics (7.72 reading rule): share of the cell loss at f_drive (R1), at 2 … 7 × f_drive and the rest, per A, B and
   frequency (mean over the seeds, calm cycles). If the shares depend on B only, the response is a static
   nonlinearity; if they move with f, a rate effect.
5. Bessel check (7.83): share of the cells with Σ_k p_k > 1.001 p per frequency. Same sample interval at all f: if the
   excess comes from motion between the samples it is the same at all f; if it comes from a drift within a cycle it
   grows with the cycle time (largest at f/2).
8. Net over-completeness of the split (PROTOKOLL 7.90): e = max(Σ_{k ≤ 7} R_k − 1, 0) = max(−R_rest, 0) per stage, 3D
   (acc.json). It is a lower bound of the per-cell excess x = Σ max(Σ_k p_k − p, 0) / Σ p (cells with a positive rest
   compensate). x itself is known only in the cut plane (acc_slices.npz), and only where the mean rest of the plane is
   > 0 in files written before 2026-10-14: acc_reduce.relative_slice stored a plane with a mean rest ≤ 0 as zeros. For
   such a plane the stored signed mean gives a lower bound of x in the plane: −rest_mean / p_mean.
9. Energy balance and rest (PROTOKOLL 7.90): ε = w_dis / w_loop − 1 of the kept cycles. If the excess of the measured
   dissipation over the loop area sits in the rest, the rest relative to w_loop without it is R_rest (1 + ε) − ε. A
   periodic state (loop closes, ε the same in each cycle) with ε > 0 breaks the energy balance: then ε is an error of
   the computed trajectory or of the sampled integral, not a loss. Per run at 9 mT: ε per cycle, closure, excess power
   (W_dis − W_loop) f.

    python eval_freq.py > results/freq_Cprime.md      # also results/freq_Cprime.json
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import acc_reduce as AR
import analyze
import analyze_design as AD

R = pathlib.Path(__file__).resolve().parent / "runs"
AMPS = list(AD.AMPS_MT)
SEEDS = range(3)
FREQ = ("fhalf", "f", "f2")
FLAB = {"fhalf": "f/2", "f": "f", "f2": "2 f"}
MIN_N = 3
E_TOL = 0.02          # net over-completeness above this value: the split of the stage is not reliable (7.90)


def folder(a, fr, k):
    if a == 10:
        return R / ("K1A10_%d" % k if fr == "f" else "K1A10_%d_%s" % (k, fr))
    return R / ("K1_%d_%s" % (k, "f1" if fr == "f" else fr))


def tci(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if len(v) < MIN_N:
        return None
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.95, len(v) - 1) * se
    return v.mean(), se, v.mean() - h, v.mean() + h, len(v)


def lnw_all():
    out = {}
    for a in (20, 10):
        for fr in FREQ:
            for k in SEEDS:
                r = AD.run_record(folder(a, fr, k))
                out[(a, fr, k)] = np.array([r["lnw@%gmT" % b] for b in AMPS]) if r else np.full(len(AMPS), np.nan)
    return out


def beta(lnw, amps=AMPS):
    idx = [AMPS.index(b) for b in amps]
    y = lnw[idx]
    if not np.all(np.isfinite(y)):
        return float("nan")
    return float(np.polyfit(np.log(amps), y, 1)[0])


def harm(a, fr):
    """{B: [R1, R_2..7, rest, bessel frac]} mean over the seeds (calm cycles, acc.json 'harm')."""
    out = {b: [] for b in AMPS}
    for k in SEEDS:
        f = folder(a, fr, k) / "acc.json"
        if not f.exists():
            continue
        for s in json.loads(f.read_text())["stages"]:
            if "harm" not in s:
                continue
            b = min(AMPS, key=lambda x: abs(x / s["B_peak_mT"] - 1))
            h = s["harm"]
            out[b].append([h["R_k"][0], sum(h["R_k"][1:]), h["R_rest"], h["bessel_frac_over_1e-3"]])
    return {b: (np.mean(v, 0), len(v)) for b, v in out.items() if v}


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    W = lnw_all()
    res = {"ratio": {}, "wh": {}, "beta": {}, "A_effect": {}, "harm": {}, "bessel": {}}
    L = ["# Frequency test C′: corner K1 at A 20 and A 10 pJ/m, f/2 / f / 2 f (15 / 30 / 60 MHz)", "",
         "Steady rule ≥ 2 calm cycles. w = loss per cycle. Paired by seed (3 seeds), 90 %% CI with t(n − 1); fewer than %d "
         "seeds: no statement. A 20: K1_k_fhalf / K1_k_f1 / K1_k_f2; A 10: K1A10_k_fhalf / K1A10_k / K1A10_k_f2." % MIN_N,
         "", "Settled values (seeds with a value, of 3) per A and frequency at 9 / 20 / 35 / 50 / 70 / 100 / 150 mT:", ""]
    for a in (20, 10):
        for fr in FREQ:
            n = [int(sum(np.isfinite(W[(a, fr, k)][i]) for k in SEEDS)) for i in range(len(AMPS))]
            L.append("- A %d, %s: %s" % (a, FLAB[fr], " / ".join(map(str, n))))
    # 1. ratios
    L += ["", "## 1. Loss per cycle relative to f, paired by seed", "",
          "| A | B (mT) | w(f/2) / w(f) | w(2 f) / w(f) | rate-independent share, from f / 2 f | from f/2 / f |",
          "|---|---|---|---|---|---|"]
    for a in (20, 10):
        for i, b in enumerate(AMPS):
            cells, rr = [], {}
            for fr in ("fhalf", "f2"):
                d = [W[(a, fr, k)][i] - W[(a, "f", k)][i] for k in SEEDS]
                s = tci(d)
                rr[fr] = s
                cells.append("× %.2f (%.2f … %.2f, n %d)" % (math.exp(s[0]), math.exp(s[2]), math.exp(s[3]), s[4])
                             if s else "no statement")
            wh2 = 2 - math.exp(rr["f2"][0]) if rr["f2"] else float("nan")
            whh = 2 * math.exp(rr["fhalf"][0]) - 1 if rr["fhalf"] else float("nan")
            res["ratio"]["%d|%g" % (a, b)] = {fr: ([math.exp(x) for x in (rr[fr][0], rr[fr][2], rr[fr][3])] + [rr[fr][4]]
                                                    if rr[fr] else None) for fr in rr}
            res["wh"]["%d|%g" % (a, b)] = {"from_f_2f": wh2, "from_fhalf_f": whh}
            L.append("| %d | %g | %s | %s | %.2f | %.2f |" % (a, b, cells[0], cells[1], wh2, whh))
    # 3. β per frequency and the A effect
    L += ["", "## 3. β (9 … 150 mT) per frequency and the A effect β(A 20) − β(A 10), paired by seed", "",
          "| frequency | β at A 20 (mean / SD, n) | β at A 10 | A effect on β (90 % CI) |", "|---|---|---|---|"]
    for fr in FREQ:
        b20 = [beta(W[(20, fr, k)]) for k in SEEDS]
        b10 = [beta(W[(10, fr, k)]) for k in SEEDS]
        d = [x - y for x, y in zip(b20, b10)]
        s = tci(d)
        res["beta"][fr] = {"20": b20, "10": b10}
        res["A_effect"][fr] = list(s) if s else None
        f2 = lambda v: "%.2f / %.2f (n %d)" % (np.nanmean(v), np.nanstd(v, ddof=1), int(np.isfinite(v).sum())) \
            if np.isfinite(v).sum() >= 2 else "n %d" % int(np.isfinite(v).sum())
        L.append("| %s | %s | %s | %s |" % (FLAB[fr], f2(np.array(b20)), f2(np.array(b10)),
                                           ("%+.2f ± %.2f (%+.2f … %+.2f, n %d)" % s) if s else "no statement"))
    # A effect of the rate-independent part w_h = 2 w(f) − w(2 f) and the rate-dependent part w(f) − w_h
    L += ["", "A effect on β of the parts of w(f) (from f and 2 f; w_h = 2 w(f) − w(2 f), w_d = w(f) − w_h = w(2 f) − "
          "w(f)), per seed where the part is > 0 at all 7 B:", ""]
    for part in ("w_h", "w_d"):
        bb = {}
        for a in (20, 10):
            bb[a] = []
            for k in SEEDS:
                wf, w2 = np.exp(W[(a, "f", k)]), np.exp(W[(a, "f2", k)])
                v = 2 * wf - w2 if part == "w_h" else w2 - wf
                bb[a].append(beta(np.log(v)) if np.all(v > 0) else float("nan"))
        d = [x - y for x, y in zip(bb[20], bb[10])]
        s = tci(d)
        res["A_effect"][part] = {"beta20": bb[20], "beta10": bb[10], "effect": list(s) if s else None}
        L.append("- %s: β at A 20 %s, at A 10 %s; A effect %s" % (
            part, " / ".join("%.2f" % x for x in bb[20]), " / ".join("%.2f" % x for x in bb[10]),
            ("%+.2f (%+.2f … %+.2f, n %d)" % (s[0], s[2], s[3], s[4])) if s else "no statement (n < 3 or a part ≤ 0)"))
    # 4. harmonics
    L += ["", "## 4. Share of the cell loss at the drive frequency (R1), at 2 … 7 × the drive frequency and the rest "
          "(mean of the seeds, calm cycles)", "",
          "| A | B (mT) | R1 at f/2 / f / 2 f | 2 … 7 × at f/2 / f / 2 f | rest at f/2 / f / 2 f |", "|---|---|---|---|---|"]
    H = {(a, fr): harm(a, fr) for a in (20, 10) for fr in FREQ}
    for a in (20, 10):
        for b in AMPS:
            cells = []
            for j in range(3):
                cells.append(" / ".join(("%.2f" % H[(a, fr)][b][0][j]) if b in H[(a, fr)] else "–" for fr in FREQ))
            res["harm"]["%d|%g" % (a, b)] = {fr: (H[(a, fr)][b][0].tolist() if b in H[(a, fr)] else None) for fr in FREQ}
            L.append("| %d | %g | %s |" % (a, b, " | ".join(cells)))
    # 5. Bessel
    L += ["", "## 5. Cells with Σ_k p_k > 1.001 p (Bessel check), mean over seeds and B", "",
          "| A | f/2 | f | 2 f |", "|---|---|---|---|"]
    for a in (20, 10):
        vals = []
        for fr in FREQ:
            v = [H[(a, fr)][b][0][3] for b in H[(a, fr)]]
            vals.append(float(np.mean(v)) if v else float("nan"))
        res["bessel"][str(a)] = dict(zip(FREQ, vals))
        L.append("| %d | %s |" % (a, " | ".join("%.0f %%" % (100 * v) if np.isfinite(v) else "–" for v in vals)))
    # 5b. A effect on ln P at each frequency: P(A 10) / P(A 20), paired by seed (added 2026-10-14 after the review)
    L += ["", "## 5b. Loss ratio P(A 10) / P(A 20) at each frequency, paired by seed, 90 % CI", "",
          "| B (mT) | 15 MHz | 30 MHz | 60 MHz |", "|---|---|---|---|"]
    res["A_ratio"] = {}
    for i, b in enumerate(AMPS):
        cells, row = [], []
        for fr in FREQ:
            s = tci([W[(10, fr, k)][i] - W[(20, fr, k)][i] for k in SEEDS])
            cells.append(("× %.2f (%.2f … %.2f)" % (math.exp(s[0]), math.exp(s[2]), math.exp(s[3]))) if s else "–")
            row.append([math.exp(s[0]), math.exp(s[2]), math.exp(s[3])] if s else None)
        res["A_ratio"]["%g" % b] = row
        L.append("| %g | %s |" % (b, " | ".join(cells)))
    # 6. exponent of P in f: P ∝ f^n, n = 1 + d ln w / d ln f (w = loss per cycle), from f/2 → f and f → 2 f
    L += ["", "## 6. Exponent n of the loss density in f (P ∝ f^n), paired by seed, 90 % CI", "",
          "n = 1: the loss per cycle does not depend on f (hysteresis); n = 2: loss per cycle ∝ f (viscous).", "",
          "| A | B (mT) | n from 15 → 30 MHz | n from 30 → 60 MHz |", "|---|---|---|---|"]
    res["n_f"] = {}
    for a in (20, 10):
        for i, b in enumerate(AMPS):
            cells, row = [], []
            for lo, hi in (("fhalf", "f"), ("f", "f2")):
                s = tci([1 + (W[(a, hi, k)][i] - W[(a, lo, k)][i]) / math.log(2.0) for k in SEEDS])
                cells.append(("%.2f (%.2f … %.2f)" % (s[0], s[2], s[3])) if s else "no statement")
                row.append([s[0], s[2], s[3]] if s else None)
            res["n_f"]["%d|%g" % (a, b)] = row
            L.append("| %d | %g | %s | %s |" % (a, b, *cells))
    # 7. where the loss sits at each frequency: concentration of the total loss in the 10 % with the highest g and at
    # the rim (acc.json "where", calm cycles), mean over the seeds
    L += ["", "## 7. Concentration c of the total loss per zone and frequency (mean over the seeds)", "",
          "| A | B (mT) | 10 % with the highest g: f/2 / f / 2 f | rim r ≥ 0.8 R: f/2 / f / 2 f |", "|---|---|---|---|"]
    res["zones"] = {}
    for a in (20, 10):
        for b in AMPS:
            cells = []
            for zk in ("texture top 10%", "r>=0.8"):
                vals = []
                for fr in FREQ:
                    v = []
                    for k in SEEDS:
                        f = folder(a, fr, k) / "acc.json"
                        for s in json.loads(f.read_text())["stages"]:
                            if min(AMPS, key=lambda x: abs(x / s["B_peak_mT"] - 1)) == b:
                                w = s["where"][zk]
                                v.append(w["p_share"] / w["vol_share"])
                    vals.append(float(np.mean(v)) if v else float("nan"))
                res["zones"]["%d|%g|%s" % (a, b, zk)] = vals
                cells.append(" / ".join("%.1f" % x for x in vals))
            L.append("| %d | %g | %s | %s |" % (a, b, *cells))
    # 8. net over-completeness of the split, 3D (acc.json): e = max(Σ_k R_k − 1, 0) per stage (PROTOKOLL 7.90)
    L += ["", "## 8. Net over-completeness e = max(Σ_{k ≤ 7} R_k − 1, 0) of the frequency split, 3D, all seeds and B", "",
          "e > 0: the parts at f … 7f add up to more than the loss (negative rest). e is a lower bound of the per-cell "
          "excess. Largest e (seed, nominal B) and number of stages with e > %.2f." % E_TOL, "",
          "| A | f/2 | f | 2 f |", "|---|---|---|---|"]
    res["overcomplete"] = {}
    for a in (20, 10):
        cells, row = [], {}
        for fr in FREQ:
            e = []
            for k in SEEDS:
                f = folder(a, fr, k) / "acc.json"
                for s in json.loads(f.read_text())["stages"] if f.exists() else []:
                    if "harm" in s:
                        e.append((max(-s["harm"]["R_rest"], 0.0), k, min(AMPS, key=lambda x: abs(x / s["B_peak_mT"] - 1))))
            if not e:
                cells.append("–")
                continue
            emax = max(e)
            if emax[0] == 0.0:                          # no stage with a negative rest: no location to name
                emax = (0.0, None, None)
            row[fr] = {"max": emax[0], "seed": emax[1], "B_mT": emax[2], "n_over": sum(x[0] > E_TOL for x in e),
                       "n": len(e)}
            cells.append(("%.3f (seed %d, %.0f mT)" % emax if emax[1] is not None else "0") +
                         ", %d of %d" % (row[fr]["n_over"], len(e)))
        res["overcomplete"][str(a)] = row
        L.append("| %d | %s |" % (a, " | ".join(cells)))
    # per-cell excess x in the cut plane: planes with a mean rest > 0 (x known), planes with a mean rest ≤ 0 (lower bound)
    L += ["", "Per-cell excess x = Σ max(Σ_k p_k − p, 0) / Σ p in the cut plane z = 0 (2 of the 4 particles): largest x "
          "of the planes with a mean rest > 0 (number of planes); planes with a mean rest ≤ 0 (stored as zeros): number "
          "and largest lower bound −rest_mean / p_mean.", "", "| A | f/2 | f | 2 f |", "|---|---|---|---|"]
    res["cut_excess"] = {}
    for a in (20, 10):
        cells, row = [], {}
        for fr in FREQ:
            x, lb = [], []
            for k in SEEDS:
                z = np.load(folder(a, fr, k) / "acc_slices.npz")
                ids = z["ids"] >= 0
                for key in z.files:
                    if key.endswith("_rest_rel"):
                        st = key[:-len("_rest_rel")]
                        if AR.plane_known(z, st, "rest"):
                            rest, pp = AR.slice_values(z, st, "rest"), AR.slice_values(z, st, "p")
                            x.append(float(np.maximum(-rest[ids], 0).sum() / pp[ids].sum()))
                        else:
                            lb.append(-float(z[st + "_rest_mean"]) / float(z[st + "_p_mean"]))
            row[fr] = {"x_max": max(x), "n_x": len(x), "n_zeroed": len(lb), "lb_max": max(lb) if lb else None}
            cells.append("%.4f (%d); %d zeroed%s" % (max(x), len(x), len(lb), (", ≥ %.4f" % max(lb)) if lb else ""))
        res["cut_excess"][str(a)] = row
        L.append("| %d | %s |" % (a, " | ".join(cells)))
    # 9. energy balance and rest: ε = w_dis / w_loop − 1 of the kept cycles, rest without ε = R_rest (1 + ε) − ε
    L += ["", "## 9. Energy balance ε = w_dis / w_loop − 1 and the rest without ε, R_rest (1 + ε) − ε (mean over the seeds)",
          "", "| A | f | " + " | ".join("%g mT" % b for b in AMPS) + " |", "|---|---|" + "---|" * len(AMPS)]
    res["balance"] = {}
    det = []
    for a in (20, 10):
        for fr in FREQ:
            out = {b: [] for b in AMPS}
            for k in SEEDS:
                d = folder(a, fr, k)
                summ = {st["name"]: st for st in json.loads((d / "summary.json").read_text())["stages"]}
                fhz = json.loads((d / "config.json").read_text())["freq"]
                for s in json.loads((d / "acc.json").read_text())["stages"]:
                    h, cyc = s["harm"], {c["cycle"]: c for c in summ[s["stage"]]["cycles"]}
                    kept = [cyc[c] for c in h["cycles_kept"]]
                    eps = [c["w_dis"] / c["w_loop"] - 1 for c in kept]
                    e = float(np.mean(eps))
                    b = min(AMPS, key=lambda x: abs(x / s["B_peak_mT"] - 1))
                    out[b].append((h["R_rest"], h["R_rest"] * (1 + e) - e, e))
                    if b == AMPS[0]:
                        det.append((a, FLAB[fr], k, min(eps), max(eps), float(np.mean([c["closure"] for c in kept])),
                                    h["R_rest"], h["R_rest"] * (1 + e) - e,
                                    float(np.mean([c["W_dis"] - c["W_loop"] for c in kept])) * fhz * 1e-6))
            res["balance"]["%d|%s" % (a, fr)] = {"%g" % b: np.mean(out[b], 0).tolist() for b in AMPS}
            L.append("| %d | %s | %s |" % (a, FLAB[fr], " | ".join("%.3f → %.3f (ε %+.3f)" % tuple(np.mean(out[b], 0))
                                                                   for b in AMPS)))
    L += ["", "At %g mT per run: ε of the kept cycles (min … max), closure |ΔM| / M_peak (mean), rest, rest without ε, "
          "excess power (W_dis − W_loop) f in W/cm³:" % AMPS[0], "",
          "| A | f | seed | ε min … max | closure | rest | rest without ε | excess power (W/cm³) |",
          "|---|---|---|---|---|---|---|---|"]
    L += ["| %d | %s | %d | %+.3f … %+.3f | %.4f | %+.3f | %+.3f | %+.3f |" % x for x in det]
    res["balance_9mT"] = det
    # 10. the A effect on β from its two ends (question of Chris 2026-10-14): β(A 20) − β(A 10) ≈ [ln r(9) − ln r(150)] /
    # ln(150/9), r = P(A 10) / P(A 20) of section 5b; compared with the fitted A effect of section 3
    L += ["", "## 10. A effect on β from the two ends (9 and 150 mT), r = P(A 10) / P(A 20)", "",
          "| f | r at 9 mT | r at 150 mT | from the 9 mT end | from the 150 mT end | sum | fitted (90 % CI) |",
          "|---|---|---|---|---|---|---|"]
    res["A_effect_ends"] = {}
    lb = math.log(AMPS[-1] / AMPS[0])
    for j, fr in enumerate(FREQ):
        r9, r150 = res["A_ratio"]["%g" % AMPS[0]][j][0], res["A_ratio"]["%g" % AMPS[-1]][j][0]
        e9, e150 = math.log(r9) / lb, -math.log(r150) / lb
        fit = res["A_effect"][fr]
        res["A_effect_ends"][fr] = {"r9": r9, "r150": r150, "end9": e9, "end150": e150}
        L.append("| %s | %.2f | %.2f | %+.2f | %+.2f | %+.2f | %+.2f (%+.2f … %+.2f) |" % (
            FLAB[fr], r9, r150, e9, e150, e9 + e150, fit[0], fit[2], fit[3]))
    print("\n".join(L))
    (pathlib.Path(__file__).resolve().parent / "results" / "freq_Cprime.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
