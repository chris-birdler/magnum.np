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
  "correction reproduces the geometry effect"  |mean d| ≤ 2 SE and ≤ 0.15, and 2 SE ≤ 0.15
  "not reproduced"                             |mean d| > 2 SE and > 0.15
  "not determined"                             else, or fewer than 3 seeds with all four values
The condition 2 SE ≤ 0.15 (power of the test) was added by decision of Chris 2026-10-20 (PROTOKOLL 3.14), after
design A met the first two conditions only because its SE was large.

Design B (--design B, PROTOKOLL 3.11 … 3.19) also gives the proofs of Chris (2026-10-20) and the means of all arms:
    N1  ln P(T15g) − ln P(T3)    resolution of m at the fixed morphology of dx 3 (mesh dependence)
    N3  ln P(T3Pc) − ln P(X3)    the correction against the exact stray field of the fine staircase (same m)
    m   ln P(D15) − ln P(X3)     resolution of m on the fine staircase (not a property of the stray field)
  reading of N1 and N3 (decision of Chris 2026-10-20, K4): "within ±0.15" if the 95 % interval of the mean (t distribution)
  lies inside ±0.15; "outside ±0.15" if it lies outside; else "not determined".

    python eval_dynamic.py --runs /root/runs > results/dynamic_D.md
    python eval_dynamic.py --design B --runs RUNS1 RUNS2 --seeds 941 … 952 --plot results/proofs_B.png

Information, not a reading (PROTOKOLL 3.26, 3.28; exclusion after the data): --exclude_events EVENTS.json (output of
check_events.py) leaves out the runs with changed cells, so each effect uses only the seed pairs without an event in
either arm.
"""
import argparse
import json
import math
import pathlib
import sys

import numpy as np
from scipy.stats import t as t_dist

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
    global ARMS, CORRECTED, SEEDS
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", default=["/root/runs"], help="run roots, searched in this order")
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    ap.add_argument("--json", default=None)
    ap.add_argument("--plot", default=None, help="design B: figure of the effects with their 95 % intervals")
    ap.add_argument("--design", choices=["A", "B"], default="A",
                    help="A: arms with own preparation (D…); B: transferred start from D15 (T…, PROTOKOLL 3.11)")
    ap.add_argument("--exclude_events", default=None,
                    help="events JSON of check_events.py: leave out the runs with an event (information only, 3.28)")
    a = ap.parse_args()
    roots = [pathlib.Path(r) for r in a.runs]
    SEEDS = tuple(a.seeds)

    def folder(name):
        for r in roots:
            if (r / name).exists():
                return r / name
        return roots[0] / name
    if a.design == "B":
        ARMS = [("T3", "dx 3, voxel model, start from D15"),
                ("T3Pc", "dx 3, pairwise correction with centroid term, q 2, p 3, start from D15 (main arm)"),
                ("D15", "dx 1.5, voxel model of dx 1.5 (the source of the transfer)"),
                ("T15g", "dx 1.5, voxel geometry of dx 3 locked, start from D15"),
                ("X3", "dx 3, exact stray field of the fine staircase (reference of the correction, 3.18), start from D15")]
        CORRECTED = ["T3Pc"]
    base, fine, lock = ("T3", "D15", "T15g") if a.design == "B" else ("D3", "D15", "D15g")
    rec = {(arm, s): EC.run_record(folder("%s_s%d" % (arm, s))) for arm, _ in ARMS for s in SEEDS}
    lnP = {k: r["lnP"] for k, r in rec.items()}
    print("# Dynamic arms at d/l_ex 96, φ 0.55, 9 mT (PROTOKOLL 3.5, D1 and D2)\n")
    if a.exclude_events:
        ev = json.loads(pathlib.Path(a.exclude_events).read_text())
        missing = ["%s_s%d" % k for k in lnP if "%s_s%d" % k not in ev]
        assert not missing, "no event record for %s" % missing
        left = sorted(k for k in lnP if ev["%s_s%d" % k]["changed"])
        for k in left:
            lnP[k] = float("nan")
        print("**Information, not a reading (exclusion after the data, PROTOKOLL 3.28):** the runs with an event "
              "(changed cells, %s) are left out: %s.\n" % (a.exclude_events, ", ".join("%s_s%d" % k for k in left)))
    print("| arm | meaning | P (W/cm³) per seed %s | n | geometric mean (W/cm³) | SD of ln P | status / Bloch points of the prepared state |"
          % " / ".join(str(s) for s in SEEDS))
    print("|---|---|---|---|---|---|---|")
    for arm, txt in ARMS:
        vals = [lnP[(arm, s)] for s in SEEDS]
        ok = [v for v in vals if np.isfinite(v)]
        per = " / ".join("%.3g" % math.exp(v) if np.isfinite(v) else "–" for v in vals)
        gm = "%.3g" % math.exp(np.mean(ok)) if ok else "–"
        sd = "%.2f" % np.std(ok, ddof=1) if len(ok) > 1 else "–"
        st = "; ".join("%s %s" % (rec[(arm, s)]["status"], bloch_end(folder("%s_s%d" % (arm, s)))) for s in SEEDS)
        print("| %s | %s | %s | %d | %s | %s | %s |" % (arm, txt, per, len(ok), gm, sd, st))
    print()
    out = {"lnP": {"%s_s%d" % k: v for k, v in lnP.items()}}

    def paired(x, y):
        return stats([lnP[(x, s)] - lnP[(y, s)] for s in SEEDS
                      if np.isfinite(lnP[(x, s)]) and np.isfinite(lnP[(y, s)])])

    dg = paired(fine, lock)
    dm = paired(lock, base)
    print("Design %s. Geometry effect of the fine staircase at dx 1.5, Δ_g = ln P(%s) − ln P(%s): %s" % (a.design, fine, lock, fmt(dg)))
    print("Resolution of m at the dx 3 geometry, ln P(%s) − ln P(%s): %s\n" % (lock, base, fmt(dm)))
    print("| corrected arm | Δ_c = ln P(arm) − ln P(%s) | d = Δ_c − Δ_g (per seed) | 95 %% interval of mean d (t) | reading |" % base)
    print("|---|---|---|---|---|")
    out["delta_g"], out["delta_m"], out["arms"] = dg, dm, {}
    for arm in CORRECTED:
        dc = paired(arm, base)
        d = stats([(lnP[(arm, s)] - lnP[(base, s)]) - (lnP[(fine, s)] - lnP[(lock, s)]) for s in SEEDS
                   if all(np.isfinite(lnP[(x, s)]) for x in (arm, base, fine, lock))])
        if d["n"] < 3:
            r = "not determined (n %d)" % d["n"]
        elif abs(d["mean"]) <= 2 * d["se"] and abs(d["mean"]) <= 0.15 and 2 * d["se"] <= 0.15:
            r = "correction reproduces the geometry effect"
        elif abs(d["mean"]) > 2 * d["se"] and abs(d["mean"]) > 0.15:
            r = "not reproduced"
        else:
            r = "not determined"
        if arm not in ("D3Pc", "T3Pc") or a.exclude_events:
            r += " (information)"
        ci = "–"
        if d["n"] >= 2 and np.isfinite(d["se"]):
            h = t_dist.ppf(0.975, d["n"] - 1) * d["se"]
            ci = "%+.2f … %+.2f" % (d["mean"] - h, d["mean"] + h)
            d["ci95"] = [d["mean"] - h, d["mean"] + h]
        print("| %s | %s | %s | %s | %s |" % (arm, fmt(dc), fmt(d), ci, r))
        out["arms"][arm] = {"delta_c": dc, "d": d, "reading": r}
    if a.design == "B":
        out["proofs"] = proofs(lnP, a.plot, info=bool(a.exclude_events))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


def interval(s):
    """95 % interval of the mean (t distribution) of a stats() dict, or None."""
    if s["n"] < 2 or not np.isfinite(s.get("se", np.nan)):
        return None
    h = t_dist.ppf(0.975, s["n"] - 1) * s["se"]
    return [s["mean"] - h, s["mean"] + h]


def proofs(lnP, plot=None, info=False):
    """Design B: the proofs N1 and N3, the resolution part m and the geometry effects, paired by seed. info: runs are
    left out after the data, so the readings of N1 and N3 are information only."""
    rows = [("N1", "resolution of m at fixed morphology (dx 3 geometry): ln P(T15g) − ln P(T3)", "T15g", "T3"),
            ("N3", "correction against the exact fine stray field (same m): ln P(T3Pc) − ln P(X3)", "T3Pc", "X3"),
            ("Δ_g", "geometry effect of the fine staircase at dx 1.5: ln P(D15) − ln P(T15g)", "D15", "T15g"),
            ("Δ_c", "geometry effect of the correction at dx 3: ln P(T3Pc) − ln P(T3)", "T3Pc", "T3"),
            ("m", "resolution of m on the fine staircase: ln P(D15) − ln P(X3)", "D15", "X3"),
            ("total", "corrected dx 3 against dx 1.5 with its own staircase: ln P(T3Pc) − ln P(D15)", "T3Pc", "D15")]
    print("\n## Proofs and mean effects (paired by seed; ln P; 95 % interval of the mean, t distribution)\n")
    print("| effect | meaning | per seed %s | mean ± SE (n) | 95 %% interval | factor on P (mean) | reading |"
          % " / ".join(str(s) for s in SEEDS))
    print("|---|---|---|---|---|---|---|")
    res = {}
    for key, txt, x, y in rows:
        per = [lnP[(x, s)] - lnP[(y, s)] if np.isfinite(lnP[(x, s)]) and np.isfinite(lnP[(y, s)]) else float("nan")
               for s in SEEDS]
        st = stats([v for v in per if np.isfinite(v)])
        ci = interval(st)
        if key in ("N1", "N3"):
            if ci is None or st["n"] < 3:
                r = "not determined (n %d)" % st["n"]
            elif -0.15 <= ci[0] and ci[1] <= 0.15:
                r = "within ±0.15"
            elif ci[0] > 0.15 or ci[1] < -0.15:
                r = "outside ±0.15"
            else:
                r = "not determined"
            if info:
                r += " (information)"
        else:
            r = "information"
        f4 = "%+.4f" if key == "N3" else "%+.3f"                 # N3 is about 10⁻⁴: more digits
        ms = ((f4 + " ± %.4f (SE, n %d)") % (st["mean"], st["se"], st["n"]) if st["n"] >= 2 and np.isfinite(st["se"])
              else "n %d" % st["n"])
        print("| %s | %s | %s | %s | %s | %s | %s |" % (
            key, txt, " / ".join(f4 % v if np.isfinite(v) else "–" for v in per), ms,
            (f4 + " … " + f4) % tuple(ci) if ci else "–", "%.4f" % math.exp(st["mean"]) if st["n"] else "–", r))
        res[key] = dict(st, per_seed=per, ci95=ci, reading=r)
    if plot and any(res[k]["n"] >= 2 for k in res):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        keys = [k for k, *_ in rows]
        fig, ax = plt.subplots(figsize=(6.4, 3.9))
        ax.axvspan(-0.15, 0.15, color="0.9", zorder=0, label="±0.15")
        ax.axvline(0.0, color="0.5", lw=0.8)
        for i, k in enumerate(keys):
            r = res[k]
            if r["n"] < 1:
                continue
            ax.plot(r["per_seed"], [i] * len(r["per_seed"]), "o", ms=3.5, color="0.6", zorder=2)
            if r["ci95"]:
                ax.plot(r["ci95"], [i, i], color="k", lw=2.0, zorder=3)
            ax.plot([r["mean"]], [i], "D", color="tab:red", ms=6, zorder=4)
        ax.set_yticks(range(len(keys)))
        labels = [r"$N_1$", r"$N_3$", r"$\Delta_\mathrm{g}$", r"$\Delta_\mathrm{c}$", r"$m$", "total"]
        ax.set_yticklabels(["%s  ($n = %d$)" % (t, res[k]["n"]) for t, k in zip(labels, keys)])
        ax.invert_yaxis()
        ax.set_xlabel(r"difference of $\ln P$")
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=3.5, color="0.6", label="seed"),
                           Line2D([], [], marker="D", ls="", ms=6, color="tab:red", label="mean"),
                           Line2D([], [], color="k", lw=2.0, label="95 % interval"),
                           Patch(color="0.9", label="±0.15")], loc="upper center", bbox_to_anchor=(0.5, -0.2),
                  ncol=4, fontsize=7, frameon=False)            # below the axes: no seed point is hidden
        fig.tight_layout()
        fig.savefig(plot, dpi=160)
    return res


def fmt(s):
    if s["n"] < 2 or not np.isfinite(s.get("se", np.nan)):
        return "n %d" % s["n"]
    return "%+.3f ± %.3f (SE, n %d)" % (s["mean"], s["se"], s["n"])


if __name__ == "__main__":
    main()
