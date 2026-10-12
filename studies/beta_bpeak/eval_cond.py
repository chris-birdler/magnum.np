#!/usr/bin/env python
"""
Conditioning test (PLAN 4.5 D, rules fixed before the runs): centre point, seeds C0 … C5 (n 6). Per seed 3 runs from
the init.pt of the design run, with the design history up to the tested stage, then 30 cycles (0 … 29) at the tested
amplitude: cond_Ck_b9 (9 mT), cond_Ck_b20 (9 → 20 mT), cond_Ck_b35 (9 → 20 → 35 mT).
- Design estimator: mean w_loop of the calm tail (steady rule, 7.54) of the measured cycles up to cycle 6 (the cycles
  of a design stage); fewer than 2 calm cycles: missing.
- Late value: mean w_loop of the cycles 28 and 29.
- Drift d = ln(late) − ln(design) per seed and amplitude (ln P = ln w + const).
- Decision per amplitude, mean of d over the seeds, t(n − 1): "design value stands" if the 90 % CI lies inside
  ±0.05; "systematic error" if the 95 % CI excludes 0 and |mean| ≥ 0.05; else "not determined".
- β 9 → 20 and 20 → 35 mT from the design estimators and from the late values (per seed, then mean ± SE).
- Events per seed: cycles 2 … 29 of the tested stage with w_loop ≤ 0 or |w_dis / w_loop − 1| > 0.20.

    python eval_cond.py > results/conditioning.md
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

from analyze import steady_tail
from analyze_design import BAL_EVENT

R = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = range(6)
AMPS = (9, 20, 35)
MARGIN = 0.05
DESIGN_LAST, LATE = 6, (28, 29)


def tested_stage(name):
    st = json.loads((R / name / "summary.json").read_text())["stages"]
    return [s for s in st if s.get("measure")][-1]


def is_event(c):
    return c["w_loop"] <= 0 or abs(c["w_dis"] / c["w_loop"] - 1.0) > BAL_EVENT


def values(stage):
    """(design w or None, late w or None, B in mT, event cycles) of one tested stage."""
    cy = stage["cycles"]
    first = 2 if "h_amp_initial" in stage else 1
    tail = steady_tail(cy[first:DESIGN_LAST + 1], BAL_EVENT)
    design = np.mean([c["w_loop"] for c in tail]) if len(tail) >= 2 else None
    late_c = [c for c in cy if c["cycle"] in LATE]
    late = np.mean([c["w_loop"] for c in late_c]) if len(late_c) == len(LATE) else None
    B = 1500.0 * np.mean([c["b_peak"] for c in tail or cy[first:]])
    events = [c["cycle"] for c in cy[first:] if is_event(c)]
    return design, late, B, events


def decide(d):
    d = np.asarray(d, float)
    n = len(d)
    if n < 3:
        return n, None, None, None, None, "not determined (n < 3)"
    m, se = d.mean(), d.std(ddof=1) / math.sqrt(n)
    h90, h95 = stats.t.ppf(0.95, n - 1) * se, stats.t.ppf(0.975, n - 1) * se
    if -MARGIN < m - h90 and m + h90 < MARGIN:
        cls = "design value stands"
    elif (m - h95 > 0 or m + h95 < 0) and abs(m) >= MARGIN:
        cls = "systematic error"
    else:
        cls = "not determined"
    return n, m, se, (m - h90, m + h90), (m - h95, m + h95), cls


def main():
    V = {}
    for s in SEEDS:
        for b in AMPS:
            name = "cond_C%d_b%d" % (s, b)
            V[(s, b)] = values(tested_stage(name)) if (R / name / "DONE").exists() else None
    print("# Conditioning test (PLAN 4.5 D): drift of the loss density over 30 cycles\n")
    print("Centre point, seeds C0 … C5. Drift d = ln P(cycles 28 … 29) − ln P(design estimator: calm tail of the "
          "cycles ≤ 6). Decision: 90 %% CI of the mean inside ±%.2f → design value stands; 95 %% CI excludes 0 and "
          "|mean| ≥ %.2f → systematic error; else not determined (t, n − 1).\n" % (MARGIN, MARGIN))
    print("## Drift per seed and amplitude (events = cycles 2 … 29 with |w_dis/w_loop − 1| > 0.20)\n")
    print("| seed | " + " | ".join("%d mT: d (events)" % b for b in AMPS) + " |")
    print("|---|---|---|---|")
    D = {b: [] for b in AMPS}
    for s in SEEDS:
        cells = []
        for b in AMPS:
            v = V[(s, b)]
            if v is None:
                cells.append("not run")
                continue
            design, late, _, ev = v
            if design is None or late is None:
                cells.append("missing (%s)" % ("< 2 calm cycles ≤ 6" if design is None else "no cycle 28/29"))
                continue
            d = math.log(late / design)
            D[b].append(d)
            cells.append("%+.3f (%d%s)" % (d, len(ev), ": " + ", ".join(map(str, ev)) if ev else ""))
        print("| C%d | %s |" % (s, " | ".join(cells)))
    print("\n## Decision per amplitude\n")
    print("| B (mT) | n | mean d ± SE | 90 % CI | 95 % CI | factor P(late)/P(design) | decision |")
    print("|---|---|---|---|---|---|---|")
    for b in AMPS:
        n, m, se, c90, c95, cls = decide(D[b])
        if m is None:
            print("| %d | %d | – | – | – | – | %s |" % (b, n, cls))
        else:
            print("| %d | %d | %+.3f ± %.3f | %+.3f … %+.3f | %+.3f … %+.3f | × %.3f | %s |"
                  % (b, n, m, se, c90[0], c90[1], c95[0], c95[1], math.exp(m), cls))
    print("\n## β from the design estimators and from the late values (per seed, mean ± SE)\n")
    print("| window | n | β design | β late | late − design |")
    print("|---|---|---|---|---|")
    for b1, b2 in ((9, 20), (20, 35)):
        rows = []
        for s in SEEDS:
            v1, v2 = V[(s, b1)], V[(s, b2)]
            if v1 is None or v2 is None or None in (v1[0], v1[1], v2[0], v2[1]):
                continue
            lb = math.log(v2[2] / v1[2])
            rows.append((math.log(v2[0] / v1[0]) / lb, math.log(v2[1] / v1[1]) / lb))
        if len(rows) < 3:
            print("| %d → %d mT | %d | – | – | not determined (n < 3) |" % (b1, b2, len(rows)))
            continue
        a = np.array(rows)
        f = lambda x: "%.3f ± %.3f" % (x.mean(), x.std(ddof=1) / math.sqrt(len(x)))
        print("| %d → %d mT | %d | %s | %s | %s |" % (b1, b2, len(a), f(a[:, 0]), f(a[:, 1]), f(a[:, 1] - a[:, 0])))
    print("\nDeterminism: the earlier stages of cond_Ck_b20 / _b35 repeat the design stages of Ck; the gate "
          "check_reruns.py compared them cycle by cycle (exit code 0, runs/mc2b_logs/queue_mc2b.log).")


if __name__ == "__main__":
    main()
