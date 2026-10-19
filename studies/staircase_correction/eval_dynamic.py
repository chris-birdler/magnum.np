#!/usr/bin/env python
"""
Evaluation of the dynamic arms of PROTOKOLL 3.5 (D1, D2): loss density P at 9 mT, d/l_ex 96, φ 0.55, seeds 941 … 944.

Value of a run: ln P of the calm tail (eval_conv.run_record of the β(B_peak) study: at least 2 calm cycles).
Arms: D3 (dx 3 voxel), D3f (dx 3, volume fraction alone, q = 2), D3P / D3L (dx 3, pairwise / cell-local correction,
q = 2, p = 3), D3P4 / D3L4 (the same with q = 4), D15 (dx 1.5 voxel), D15g (dx 1.5, geometry of dx 3 locked).

Reading of D1 (fixed before the results), paired by seed:
    Δ_c = ln P(corrected arm) − ln P(D3)        geometry effect of the correction at dx 3
    Δ_g = ln P(D15) − ln P(D15g)                 geometry effect of the fine staircase (q = 2) at dx 1.5
    d   = Δ_c − Δ_g                              per seed; mean ± SE
  "correction reproduces the geometry effect"  |mean d| ≤ 2 SE and ≤ 0.15
  "not reproduced"                             |mean d| > 2 SE and > 0.15
  "not determined"                             else, or fewer than 3 seeds with all four values

    python eval_dynamic.py --runs /root/runs > results/dynamic_D.md
"""
import argparse
import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
import eval_conv as EC  # noqa: E402

SEEDS = (941, 942, 943, 944)
ARMS = [("D3", "dx 3, voxel model (production)"), ("D3f", "dx 3, volume fraction alone, q 2"),
        ("D3P", "dx 3, pairwise correction, q 2, p 3"), ("D3L", "dx 3, cell-local correction, q 2, p 3"),
        ("D3P4", "dx 3, pairwise correction, q 4, p 3"), ("D3L4", "dx 3, cell-local correction, q 4, p 3"),
        ("D15", "dx 1.5, voxel model of dx 1.5"), ("D15g", "dx 1.5, voxel geometry of dx 3 locked")]
CORRECTED = ["D3P", "D3L", "D3P4", "D3L4", "D3f"]


def stats(v):
    v = np.asarray(v, float)
    n = len(v)
    if n < 2:
        return {"n": n, "mean": float(v.mean()) if n else float("nan"), "se": float("nan")}
    return {"n": n, "mean": float(v.mean()), "se": float(v.std(ddof=1) / math.sqrt(n)), "sd": float(v.std(ddof=1))}


def bloch_end(d):
    """Bloch points of the prepared state from stdout.log (line '[predemag] after the relaxation ... Bloch points +a -b')."""
    try:
        lines = [l for l in (d / "stdout.log").read_text().splitlines() if "after the relaxation" in l and "Bloch" in l]
        return lines[-1].split("Bloch points")[-1].strip() if lines else "?"
    except OSError:
        return "?"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="/root/runs")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    root = pathlib.Path(a.runs)
    rec = {(arm, s): EC.run_record(root / ("%s_s%d" % (arm, s))) for arm, _ in ARMS for s in SEEDS}
    lnP = {k: r["lnP"] for k, r in rec.items()}
    print("# Dynamic arms at d/l_ex 96, φ 0.55, 9 mT (PROTOKOLL 3.5, D1 and D2)\n")
    print("| arm | meaning | P (W/cm³) per seed 941 / 942 / 943 / 944 | geometric mean | SD of ln P | status / Bloch points of the prepared state |")
    print("|---|---|---|---|---|---|")
    for arm, txt in ARMS:
        vals = [lnP[(arm, s)] for s in SEEDS]
        ok = [v for v in vals if np.isfinite(v)]
        per = " / ".join("%.3g" % math.exp(v) if np.isfinite(v) else "–" for v in vals)
        gm = "%.3g" % math.exp(np.mean(ok)) if ok else "–"
        sd = "%.2f" % np.std(ok, ddof=1) if len(ok) > 1 else "–"
        st = "; ".join("%s %s" % (rec[(arm, s)]["status"], bloch_end(root / ("%s_s%d" % (arm, s)))) for s in SEEDS)
        print("| %s | %s | %s | %s | %s | %s |" % (arm, txt, per, gm, sd, st))
    print()
    out = {"lnP": {"%s_s%d" % k: v for k, v in lnP.items()}}

    def paired(x, y):
        return stats([lnP[(x, s)] - lnP[(y, s)] for s in SEEDS
                      if np.isfinite(lnP[(x, s)]) and np.isfinite(lnP[(y, s)])])

    dg = paired("D15", "D15g")
    dm = paired("D15g", "D3")
    print("Geometry effect of the fine staircase at dx 1.5, Δ_g = ln P(D15) − ln P(D15g): %s" % fmt(dg))
    print("Resolution of m at the dx 3 geometry, ln P(D15g) − ln P(D3): %s\n" % fmt(dm))
    print("| corrected arm | Δ_c = ln P(arm) − ln P(D3) | d = Δ_c − Δ_g (per seed) | reading |")
    print("|---|---|---|---|")
    out["delta_g"], out["delta_m"], out["arms"] = dg, dm, {}
    for arm in CORRECTED:
        dc = paired(arm, "D3")
        d = stats([(lnP[(arm, s)] - lnP[("D3", s)]) - (lnP[("D15", s)] - lnP[("D15g", s)]) for s in SEEDS
                   if all(np.isfinite(lnP[(x, s)]) for x in (arm, "D3", "D15", "D15g"))])
        if d["n"] < 3:
            r = "not determined (n %d)" % d["n"]
        elif abs(d["mean"]) <= 2 * d["se"] and abs(d["mean"]) <= 0.15:
            r = "correction reproduces the geometry effect"
        elif abs(d["mean"]) > 2 * d["se"] and abs(d["mean"]) > 0.15:
            r = "not reproduced"
        else:
            r = "not determined"
        if arm in ("D3P4", "D3L4", "D3f"):
            r += " (information: %s)" % ("q 4 target, the reference is q 2" if arm != "D3f" else "no near-field correction")
        print("| %s | %s | %s | %s |" % (arm, fmt(dc), fmt(d), r))
        out["arms"][arm] = {"delta_c": dc, "d": d, "reading": r}
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


def fmt(s):
    if s["n"] < 2 or not np.isfinite(s.get("se", np.nan)):
        return "n %d" % s["n"]
    return "%+.3f ± %.3f (SE, n %d)" % (s["mean"], s["se"], s["n"])


if __name__ == "__main__":
    main()
