#!/usr/bin/env python
"""
Evaluation of the dynamic arms of PROTOKOLL 3.5 (D1, D2): loss density P at 9 mT, d/l_ex 96, φ 0.55, seeds 941 … 944.

Value of a run: ln P of the calm tail (eval_conv.run_record of the β(B_peak) study: at least 2 calm cycles).
Arms: see ARMS. The main corrected arm is D3Pc (pairwise form with the centroid term, q = 2, p = 3; PROTOKOLL 3.7,
D4, D8); the other corrected arms are information.

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
        ("D3Pc", "dx 3, pairwise correction with centroid term, q 2, p 3 (main arm, D4/D8)"),
        ("D3Pc4", "dx 3, pairwise correction with centroid term, q 4, p 3"),
        ("D3Lc", "dx 3, cell-local correction with centroid term, q 2, p 3 (not usable by R3; information)"),
        ("D3P", "dx 3, pairwise correction without centroid term, q 2, p 3 (information)"),
        ("D15", "dx 1.5, voxel model of dx 1.5"), ("D15g", "dx 1.5, voxel geometry of dx 3 locked")]
CORRECTED = ["D3Pc", "D3Pc4", "D3Lc", "D3P", "D3f"]


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
    ap.add_argument("--design", choices=["A", "B"], default="A",
                    help="A: arms with own preparation (D…); B: transferred start from D15 (T…, PROTOKOLL 3.11)")
    a = ap.parse_args()
    root = pathlib.Path(a.runs)
    global ARMS, CORRECTED
    if a.design == "B":
        ARMS = [("T3", "dx 3, voxel model, start from D15"),
                ("T3Pc", "dx 3, pairwise correction with centroid term, q 2, p 3, start from D15 (main arm)"),
                ("D15", "dx 1.5, voxel model of dx 1.5 (the source of the transfer)"),
                ("T15g", "dx 1.5, voxel geometry of dx 3 locked, start from D15")]
        CORRECTED = ["T3Pc"]
    base, fine, lock = ("T3", "D15", "T15g") if a.design == "B" else ("D3", "D15", "D15g")
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

    dg = paired(fine, lock)
    dm = paired(lock, base)
    print("Design %s. Geometry effect of the fine staircase at dx 1.5, Δ_g = ln P(%s) − ln P(%s): %s" % (a.design, fine, lock, fmt(dg)))
    print("Resolution of m at the dx 3 geometry, ln P(%s) − ln P(%s): %s\n" % (lock, base, fmt(dm)))
    print("| corrected arm | Δ_c = ln P(arm) − ln P(%s) | d = Δ_c − Δ_g (per seed) | reading |" % base)
    print("|---|---|---|---|")
    out["delta_g"], out["delta_m"], out["arms"] = dg, dm, {}
    for arm in CORRECTED:
        dc = paired(arm, base)
        d = stats([(lnP[(arm, s)] - lnP[(base, s)]) - (lnP[(fine, s)] - lnP[(lock, s)]) for s in SEEDS
                   if all(np.isfinite(lnP[(x, s)]) for x in (arm, base, fine, lock))])
        if d["n"] < 3:
            r = "not determined (n %d)" % d["n"]
        elif abs(d["mean"]) <= 2 * d["se"] and abs(d["mean"]) <= 0.15:
            r = "correction reproduces the geometry effect"
        elif abs(d["mean"]) > 2 * d["se"] and abs(d["mean"]) > 0.15:
            r = "not reproduced"
        else:
            r = "not determined"
        if arm not in ("D3Pc", "T3Pc"):
            r += " (information)"
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
