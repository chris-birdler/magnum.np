#!/usr/bin/env python
"""
Damping test at small B (PLAN 4.7, rules fixed before the runs): centre point, seeds C0 … C7, α 0.01. Per seed and
amplitude (9, 20 mT) one run cond_C<k>_a001_b<mT>: the design history up to the tested amplitude, then 30 cycles.
- Value at α 0.01: mean w_loop of the calm tail (steady rule, 7.54: calm = w_loop > 0 and |w_dis / w_loop − 1| ≤ 0.20)
  of the measured cycles of the tested stage; fewer than 2 calm cycles: missing.
- Reference at α 0.1: the design value of the same seed (steady rule ≥ 2 with reruns, analyze_design.run_record).
- Ratio per seed: ln P(α 0.01) − ln P(α 0.1); per amplitude mean ± SE, 90 % CI with t(n − 1); fewer than 3 seeds: no
  value. Also the design runs C<k>_a001 (9 cycles) for comparison, the length of the calm tail, the event cycles, and
  R1 (share of the cell loss at f, acc.json) of the tested stage at α 0.01 (window with or without event cycles).

    python eval_alpha_cond.py > results/alpha_cond.md
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import deck_data as DD

R = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = range(8)
AMPS = (9, 20)
MIN_N = 3


def tested(name):
    """Last stage of a run (the tested amplitude) and its measured cycles (as analyze.py: from 2 after a drive
    correction, else from 1)."""
    st = json.loads((R / name / "summary.json").read_text())["stages"][-1]
    return st, st["cycles"][2:] if "h_amp_initial" in st else st["cycles"][1:]


def calm_value(name):
    st, meas = tested(name)
    tail = analyze.steady_tail(meas, AD.BAL_EVENT)
    ev = [c["cycle"] for c in meas if c["w_loop"] <= 0 or abs(c["w_dis"] / c["w_loop"] - 1.0) > AD.BAL_EVENT]
    return (float(np.mean([c["w_loop"] for c in tail])) if len(tail) >= 2 else None), len(tail), len(meas), ev


def summary(v):
    v = np.asarray(v, float)
    if len(v) < MIN_N:
        return None
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.95, len(v) - 1) * se
    return v.mean(), se, v.std(ddof=1), v.mean() - h, v.mean() + h, len(v)


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2          # steady rule ≥ 2 calm cycles (7.54), as deck_data
    L = ["# Damping test at small B: α 0.01 against α 0.1 (PLAN 4.7)", "",
         "Centre point, seeds C0 … C7. α 0.01: calm tail (≥ 2 calm cycles) of 30 cycles at the tested amplitude after "
         "the design history. α 0.1: design value of the same seed (steady rule ≥ 2, reruns). Ratio = P(α 0.01) / "
         "P(α 0.1), mean of ln over the seeds, SE, SD, 90 %% CI with t; fewer than %d seeds: no value." % MIN_N, ""]
    L += ["## Per seed", "", "| seed | B (mT) | calm tail (cycles) | event cycles | P(α 0.01) / P(α 0.1), long run | "
          "the same from the design run C<k>_a001 (9 cycles) | R1 at α 0.01 (acc window calm?) |",
          "|---|---|---|---|---|---|---|"]
    ratios = {b: [] for b in AMPS}
    design = {b: [] for b in AMPS}
    r1 = {b: [] for b in AMPS}
    for k in SEEDS:
        ref = AD.run_record(AD.rerun_dir(R / ("C%d" % k)) or R / ("C%d" % k))
        des = AD.run_record(R / ("C%d_a001" % k)) if (R / ("C%d_a001" % k) / "DONE").exists() else None
        for b in AMPS:
            name = "cond_C%d_a001_b%d" % (k, b)
            if not (R / name / "DONE").exists():
                L.append("| C%d | %d | – | – | not run | – | – |" % (k, b))
                continue
            w, nt, nm, ev = calm_value(name)
            # ref: ln w_loop at the nominal B (analyze_design.run_record); w: w_loop at the measured B (drive corrected
            # to the nominal b after cycle 0)
            lr = math.log(w) - ref["lnw@%gmT" % b] if (w is not None and ref is not None) else None
            if lr is not None:
                ratios[b].append(lr)
            ld = (des["lnw@%gmT" % b] - ref["lnw@%gmT" % b]) if (des is not None and ref is not None
                                                                and np.isfinite(des["lnw@%gmT" % b])) else None
            if ld is not None:
                design[b].append(ld)
            st = DD.acc_stages(R / name)[-1] if (R / name / "acc.json").exists() else None
            if st is not None:
                r1[b].append((st["R1"], st["calm"]))
            L.append("| C%d | %d | %d of %d | %s | %s | %s | %s |" % (
                k, b, nt, nm, ", ".join(map(str, ev)) if ev else "none",
                "× %.2f" % math.exp(lr) if lr is not None else "no calm tail",
                "× %.2f" % math.exp(ld) if ld is not None else "no steady value",
                ("%.3f (%s)" % (st["R1"], "calm" if st["calm"] else "with events")) if st is not None else "–"))
    L += ["", "## Per amplitude", "", "| B (mT) | P(α 0.01) / P(α 0.1), long runs | SD of ln | n | "
          "design runs (9 cycles), n | R1 at α 0.01, median (calm windows / all) |", "|---|---|---|---|---|---|"]
    for b in AMPS:
        s = summary(ratios[b])
        sd_ = summary(design[b])
        cal = [x for x, c in r1[b] if c]
        L.append("| %d | %s | %s | %d | %s | %s / %s |" % (
            b, ("× %.2f (ln %+.2f ± %.2f, 90 %% CI × %.2f … %.2f)" % (math.exp(s[0]), s[0], s[1], math.exp(s[3]),
                                                                      math.exp(s[4]))) if s else "no value",
            ("%.2f" % s[2]) if s else "–", len(ratios[b]),
            ("× %.2f, n %d" % (math.exp(sd_[0]), sd_[5])) if sd_ else "n %d: no value" % len(design[b]),
            ("%.3f (n %d)" % (np.median(cal), len(cal))) if cal else "–",
            ("%.3f (n %d)" % (np.median([x for x, _ in r1[b]]), len(r1[b]))) if r1[b] else "–"))
    L += ["", "Note: a loss proportional to α would give × 0.10."]
    print("\n".join(L))


if __name__ == "__main__":
    main()
