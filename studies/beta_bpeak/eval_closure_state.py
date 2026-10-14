#!/usr/bin/env python
"""
Loop closure and settling of the main run, and the same measures for C′ (PROTOKOLL 7.101; questions of Chris,
2026-10-14). Recomputes the tables of 7.101 from the run files (local, no GPU):
  1. protocol of the 116 runs of the main analysis (config.json);
  2. steady tails (calm cycles at the end of each stage, steady-tail rule 7.54);
  3. closure |M(end) − M(start)| / M_peak of the steady-tail cycles, and against the half loop width M_peak sin δ
     (estimate for an elliptic loop: sin δ = w_loop / (2π h b));
  4. Δe / w_loop = (w_loop − w_dis) / w_loop per cycle, the share of consecutive cycles with the same sign of ΔM, and
     the events (criteria fixed before the data were looked at: "recurrent jumps" if |mean Δe / w_loop| ≤ max(2 SE,
     0.02) and the share ≤ 0.5 + 2 SE (binomial); "still relaxing" if the mean < −max(2 SE, 0.02) or the share
     > 0.5 + 2 SE; else not determined);
  5. the share against the references 1/3 (fixed state with scatter) and 1/2 (wandering state), SE by a bootstrap
     over whole stages (added after the data were looked at);
  6. C′ per frequency and particle size;
  7. K1_k at 9 mT: closure and w_dis / w_loop per cycle (the w_dis error at atol 1e-5, 7.88);
  8. B0 (atol 1e-5) against B0b (atol 1e-6): closure and W_loop per common complete cycle from samples_00_amp00.csv;
     reading of Chris: "physical" if the median closure with 1e-6 is ≥ 0.5 × the median with 1e-5, per seed.

    python eval_closure_state.py > results/closure_state.md      # also results/closure_state.json
"""
import json
import math
import pathlib
import sys

import numpy as np

import analyze
import analyze_design as AD
import eval_atol30 as E3

HERE = pathlib.Path(__file__).resolve().parent
R = HERE / "runs"
AMPS = list(AD.AMPS_MT)
MU0 = 4e-7 * math.pi
N_BOOT = 2000
CPRIME = (("15 MHz", "large", "K1A10_%d_fhalf"), ("15 MHz", "small", "K1_%d_fhalf"),
          ("30 MHz", "large", "K1A10_%d"), ("30 MHz", "small", "K1_%d_f1"),
          ("60 MHz", "large", "K1A10_%d_f2"), ("60 MHz", "small", "K1_%d_f2"))
B0_SEEDS = (941, 942, 943, 944)


def same_sign_pairs(signs):
    """(pairs with the same sign, pairs) of consecutive values."""
    return sum(signs[j] == signs[j + 1] for j in range(len(signs) - 1)), max(len(signs) - 1, 0)


def boot_share(stages, rng, n_boot=N_BOOT):
    """Same-sign share over all pairs of the stages and its SE by a bootstrap over whole stages."""
    st = np.array(stages, dtype=float)
    p = st[:, 0].sum() / st[:, 1].sum()
    bs = []
    for _ in range(n_boot):
        x = st[rng.integers(0, len(st), len(st))]
        bs.append(x[:, 0].sum() / x[:, 1].sum())
    return float(p), float(np.std(bs, ddof=1)), math.sqrt(p * (1 - p) / st[:, 1].sum()), int(st[:, 1].sum())


def run_stages(d):
    """Per measured stage of run d: the cycle records with the signed ΔM / M_peak from the samples file."""
    summ = json.loads((d / "summary.json").read_text())
    cfg = json.loads((d / "config.json").read_text())
    ms = [s for s in summ["stages"] if s.get("measure")]
    files = sorted(d.glob("samples_*_amp*.csv"))
    out = []
    for i, st in enumerate(ms):
        a = np.loadtxt(files[i], delimiter=",", comments="#")
        for c in st["cycles"]:
            r = a[a[:, 0] == c["cycle"]]
            c["dM"] = float((r[-1, 3] - r[0, 3]) / c["M_peak"])
            assert abs(abs(c["dM"]) - c["closure"]) < 1e-6 * max(1.0, c["closure"])
        meas = E3.measured(st)
        tail = analyze.steady_tail(meas, AD.BAL_EVENT)
        out.append({"h": st["h_amp"], "cycles": st["cycles"], "meas": meas, "tail": tail,
                    "max_cycles": cfg.get("max_cycles_per_amp")})
    return out


def is_event(c):
    return c["w_loop"] <= 0 or abs(c["w_dis"] / c["w_loop"] - 1.0) > AD.BAL_EVENT


def de(c):
    return (c["w_loop"] - c["w_dis"]) / c["w_loop"]


def stage_measures(s):
    tail = s["tail"]
    ev = [c for c in s["meas"] if is_event(c)]
    w = np.array([c["w_loop"] for c in tail])
    return {"n_tail": len(tail), "n_meas": len(s["meas"]), "n_cycles": len(s["cycles"]),
            "all_calm": len(ev) == 0, "ev_release": sum(c["w_loop"] <= 0 or de(c) < 0 for c in ev),
            "ev_uptake": sum(c["w_loop"] > 0 and de(c) >= 0 for c in ev),
            "clo_mean": float(np.mean([c["closure"] for c in tail])) if tail else float("nan"),
            "clo": [c["closure"] for c in tail],
            "half": [c["closure"] / (c["w_loop"] / (2 * math.pi * s["h"] * c["b_peak"])) for c in tail if c["w_loop"] > 0],
            "de": float(np.mean([de(c) for c in tail])) if tail else float("nan"),
            "pairs": same_sign_pairs([np.sign(c["dM"]) for c in tail]),
            "w_sd": float(w.std(ddof=1) / w.mean()) if len(w) > 1 else float("nan"),
            "at_cap_without_tail": len(s["cycles"]) >= (s["max_cycles"] or 10 ** 9) and len(tail) < 2}


def reading_drift(m, se, q, qse):
    lim = max(2 * se, 0.02)
    if abs(m) <= lim and q <= 0.5 + 2 * qse:
        return "recurrent jumps"
    if m < -lim or q > 0.5 + 2 * qse:
        return "still relaxing"
    return "not determined"


def b0_cycles(d):
    """{cycle: (closure, W_loop)} of the complete cycles (257 samples) in samples_00_amp00.csv, or None."""
    f = d / "samples_00_amp00.csv"
    if not f.exists():
        return None
    a = np.loadtxt(f, delimiter=",", comments="#")
    out = {}
    for c in np.unique(a[:, 0]).astype(int):
        r = a[a[:, 0] == c]
        if len(r) < 257:
            continue
        M, H = r[:, 3], r[:, 2]
        trapz = getattr(np, "trapezoid", None) or np.trapz
        out[int(c)] = (float(abs(M[-1] - M[0]) / (0.5 * (M.max() - M.min()))), float(trapz(H, MU0 * (H + M))))
    return out


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    rng = np.random.default_rng(1)
    res = {}
    recs, _ = AD.load_runs(R, reruns=True)
    names = [r["name"] for r in recs]
    L = ["# Loop closure and settling of the main run (PROTOKOLL 7.101, eval_closure_state.py)", ""]

    # 1. protocol
    cfgs = {}
    for n in names:
        a = json.loads((R / n / "config.json").read_text())
        key = (a["protocol"], a["cycles_per_amp"], a["max_cycles_per_amp"], a.get("min_calm") or 0)
        cfgs.setdefault(key, []).append(n)
    L += ["## 1. Protocol (%d runs)" % len(names), "", "| protocol | cycles_per_amp | max_cycles_per_amp | min_calm | runs |",
          "|---|---|---|---|---|"]
    L += ["| %s | %d | %d | %d | %d |" % (k + (len(v),)) for k, v in sorted(cfgs.items())]
    S = {n: run_stages(R / n) for n in names}
    used = [len(s["cycles"]) for k, v in cfgs.items() if k[3] > 0 for n in v for s in S[n]]
    L += ["", "Cycles per stage in the runs with min_calm > 0: %d … %d." % (min(used), max(used)), ""]
    res["protocol"] = {"|".join(map(str, k)): len(v) for k, v in cfgs.items()}

    M = {b: [stage_measures(S[n][i]) for n in names] for i, b in enumerate(AMPS)}

    # 2. steady tails
    L += ["## 2. Steady tails (812 stages)", "", "| B (mT) | stages with tail ≥ 2 | smallest tail | all measured cycles calm |",
          "|---|---|---|---|"]
    for b in AMPS:
        t = [x["n_tail"] for x in M[b]]
        L.append("| %g | %d / %d | %d | %d |" % (b, sum(v >= 2 for v in t), len(t), min(t), sum(x["all_calm"] for x in M[b])))

    # 3. closure
    L += ["", "## 3. Closure (mean of the steady-tail cycles per stage; stages above a limit)", "",
          "| B (mT) | > 2 % | > 0.5 % | > 0.2 % | > 0.1 % | median | median of closure / half width | cycles > 25 % "
          "of the half width | > 50 % |", "|---|---|---|---|---|---|---|---|---|"]
    tot = np.zeros(4, int)
    for b in AMPS:
        v = np.array([x["clo_mean"] for x in M[b]])
        cnt = [int((v > t).sum()) for t in (0.02, 0.005, 0.002, 0.001)]
        tot += cnt
        hw = np.array([h for x in M[b] for h in x["half"]])
        L.append("| %g | %d | %d | %d | %d | %.2f %% | %.3f | %.1f %% | %.1f %% |" % (
            b, *cnt, 100 * np.median(v), np.median(hw), 100 * (hw > 0.25).mean(), 100 * (hw > 0.5).mean()))
    L.append("| sum | %d | %d | %d | %d | | | | |" % tuple(tot))

    # 4. drift or recurrent jumps (pre-registered)
    L += ["", "## 4. Monotonic change or recurrent jumps (criteria fixed before the data were looked at)", "",
          "| B (mT) | Δe / w_loop mean ± SE | same-sign share (pairs) | binomial SE | stages with an event | events "
          "release / uptake | reading |", "|---|---|---|---|---|---|---|"]
    res["drift"] = {}
    for b in AMPS:
        d_ = np.array([x["de"] for x in M[b]])
        m, se = d_.mean(), d_.std(ddof=1) / math.sqrt(len(d_))
        same = sum(x["pairs"][0] for x in M[b])
        n = sum(x["pairs"][1] for x in M[b])
        q, qse = same / n, math.sqrt(0.25 / n)
        rd = reading_drift(m, se, q, qse)
        rel, upt = sum(x["ev_release"] for x in M[b]), sum(x["ev_uptake"] for x in M[b])
        nev = sum(not x["all_calm"] for x in M[b])
        res["drift"]["%g" % b] = {"de_mean": m, "de_se": se, "share": q, "pairs": n, "stages_with_event": nev,
                                  "release": rel, "uptake": upt, "reading": rd}
        L.append("| %g | %+.2f ± %.2f %% | %.1f %% (%d) | %.1f %% | %d | %d / %d | %s |" % (
            b, 100 * m, 100 * se, 100 * q, n, 100 * qse, nev, rel, upt, rd))

    # 5. share against 1/3 and 1/2, stage bootstrap
    L += ["", "## 5. Same-sign share against 1/3 (fixed state) and 1/2 (wandering state); SE: bootstrap over whole "
          "stages, %d draws (added after the data were looked at)" % N_BOOT, "",
          "| group | share | SE bootstrap | SE binomial at the share | pairs | distance to 1/2 (SE) | distance to 1/3 (SE) |",
          "|---|---|---|---|---|---|---|"]
    res["share"] = {}
    groups = [("%g mT" % b, [x["pairs"] for x in M[b]]) for b in AMPS]
    groups.append(("all B", [x["pairs"] for b in AMPS for x in M[b]]))
    MC = {}
    for f, size, pat in CPRIME:
        MC[(f, size)] = [stage_measures(s) for k in range(3) for s in run_stages(R / (pat % k))]
        groups.append(("C′ %s %s" % (f, size), [x["pairs"] for x in MC[(f, size)]]))
    for g, st in groups:
        p, se, sb, n = boot_share(st, rng)
        res["share"][g] = {"share": p, "se_boot": se, "se_binom": sb, "pairs": n}
        L.append("| %s | %.1f %% | %.1f %% | %.1f %% | %d | %+.1f | %+.1f |" % (
            g, 100 * p, 100 * se, 100 * sb, n, (p - 0.5) / se, (p - 1 / 3) / se))

    # 6. C′ per frequency
    L += ["", "## 6. C′ per frequency (corner K1, 3 seeds × 7 amplitudes)", "",
          "| f | particle | stages with tail ≥ 2 | largest number of cycles | at the limit without tail | events release / "
          "uptake | median closure | median cycle SD of w |", "|---|---|---|---|---|---|---|---|"]
    for (f, size), xs in MC.items():
        L.append("| %s | %s | %d / %d | %d | %d | %d / %d | %.3f %% | %.1f %% |" % (
            f, size, sum(x["n_tail"] >= 2 for x in xs), len(xs), max(x["n_cycles"] for x in xs),
            sum(x["at_cap_without_tail"] for x in xs), sum(x["ev_release"] for x in xs), sum(x["ev_uptake"] for x in xs),
            100 * np.median([c for x in xs for c in x["clo"]]), 100 * np.median([x["w_sd"] for x in xs])))

    # 7. K1_k at 9 mT
    L += ["", "## 7. K1_k at 9 mT (main run, atol 1e-5): closure and w_dis / w_loop per cycle", ""]
    for k in range(3):
        s = S["K1_%d" % k][0] if "K1_%d" % k in S else run_stages(R / ("K1_%d" % k))[0]
        L.append("- K1_%d: closure %s; w_dis / w_loop %s; c (eval_atol30) %.3f %%" % (
            k, " / ".join("%.5f" % c["closure"] for c in s["cycles"]),
            " / ".join("%.3f" % (c["w_dis"] / c["w_loop"]) for c in s["cycles"]),
            100 * E3.closure_median(E3.load(R / ("K1_%d" % k))[0])))

    # 8. B0 against B0b
    L += ["", "## 8. B0 (atol 1e-5) against B0b (atol 1e-6): common complete cycles", ""]
    res["b0"] = {}
    rd = []
    for s in B0_SEEDS:
        a, b = b0_cycles(R / ("B0_F_s%d" % s)), b0_cycles(R / ("B0b_F_s%d" % s))
        if a is None or b is None:
            L.append("- s%d: sample files not local" % s)
            continue
        com = sorted(set(a) & set(b))
        r = float(np.median([b[c][0] for c in com]) / np.median([a[c][0] for c in com]))
        rd.append(r >= 0.5)
        res["b0"][s] = {"cycles": com, "ratio_closure": [b[c][0] / a[c][0] for c in com],
                        "ratio_wloop": [b[c][1] / a[c][1] for c in com], "median_ratio": r}
        L.append("- s%d, cycles %s: closure 1e-6 / 1e-5 %s; W_loop 1e-6 / 1e-5 %s; median ratio %.4f" % (
            s, com, " / ".join("%.4f" % (b[c][0] / a[c][0]) for c in com),
            " / ".join("%.6f" % (b[c][1] / a[c][1]) for c in com), r))
    if rd:
        L += ["", "**Reading: physical in %d of %d seeds, numerical in %d.**" % (sum(rd), len(rd), len(rd) - sum(rd))]
    print("\n".join(L))
    (HERE / "results" / "closure_state.json").write_text(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    sys.exit(main())
