#!/usr/bin/env python
"""
Evaluate runs of run_loops.py: W(B_peak), P(B_peak) and beta_eff(B_peak).

    python analyze.py runs/*                   # table + all_runs.csv + plots
    python analyze.py --pair runs/X_f30 runs/X_f15   # frequency separation

Per run and measured amplitude, W and B_peak are the MEAN over all cycles
except the first one (the first cycle after an amplitude change is a
transient). W_err is the standard error of that mean (cycle-to-cycle scatter;
it is not small for 4 particles, see HANDOVER.md).

beta_eff is the local slope d ln W / d ln B_peak (central differences on the
log-log data, one-sided at the ends). With a constant frequency, the slope of
the loss power P = f W is the same.

Frequency separation (--pair): with W(f) = W_h + c f (quasi-static part plus
a part linear in f, e.g. LLG damping, relaxation), two runs that differ ONLY
in the frequency give
    W_h = (f2 W(f1) - f1 W(f2)) / (f2 - f1).
The stages are matched by name (same drive amplitudes).
"""
import argparse
import csv
import json
import math
import pathlib
import sys

import numpy as np


def load_run(d):
    d = pathlib.Path(d)
    cfg = json.loads((d / "config.json").read_text())
    summ = json.loads((d / "summary.json").read_text())
    rows = []
    for st in summ["stages"]:
        if not st.get("measure") or not st.get("complete") or not st["cycles"]:
            continue
        cyc = st["cycles"][1:] if len(st["cycles"]) > 1 else st["cycles"]   # drop the transient cycle
        W = np.array([c["W_loop"] for c in cyc])
        rows.append({"stage": st["name"], "H_amp": st["H_amp"],
                     "B_peak": float(np.mean([c["B_peak"] for c in cyc])),
                     "W_loop": float(W.mean()),
                     "W_err": float(W.std(ddof=1) / math.sqrt(len(W))) if len(W) > 1 else float("nan"),
                     "W_dis": float(np.mean([c["W_dis"] for c in cyc])),
                     "closure": cyc[-1]["closure"], "dW_rel": cyc[-1].get("dW_rel", float("nan")),
                     "steady": st.get("steady", False), "n_cycles": len(st["cycles"]),
                     "mean_dt": float(np.mean([c["mean_dt"] for c in cyc])),
                     "wall_s": float(sum(c["wall_s"] for c in st["cycles"]))})
    rows.sort(key=lambda r: r["B_peak"])
    add_beta(rows, "W_loop", "beta")
    for r in rows:
        r["P_kW_m3"] = cfg["freq"] * r["W_loop"] / 1e3
        r["balance"] = r["W_dis"] / r["W_loop"] if r["W_loop"] != 0 else float("nan")
    return cfg, rows


def add_beta(rows, key, out_key):
    n = len(rows)
    lb = [math.log(r["B_peak"]) for r in rows]
    lw = [math.log(r[key]) if r[key] > 0 else float("nan") for r in rows]
    for i in range(n):
        if n < 2:
            rows[i][out_key] = float("nan"); continue
        j0, j1 = max(i - 1, 0), min(i + 1, n - 1)
        rows[i][out_key] = (lw[j1] - lw[j0]) / (lb[j1] - lb[j0])


def label(cfg):
    return ("d=%.2fum phi=%.3f sig=%.0fMPa xi=%.0fnm K1=%.0fk f=%.0fMHz %s seed=%d/%d" %
            (cfg["d"] * 1e6, cfg["phi_vox"], cfg["sigma_rms"] / 1e6, cfg["xi"] * 1e9, cfg["K1"] / 1e3,
             cfg["freq"] / 1e6, cfg["precision"], cfg["seed_axes"], cfg["seed_stress"]))


def print_table(name, cfg, rows):
    print("\n== %s\n   %s" % (name, label(cfg)))
    print("   %-6s %10s %10s %11s %8s %11s %7s %6s %8s %6s %3s" %
          ("stage", "H [A/m]", "B_pk [T]", "W [J/m3]", "W_err", "P [kW/m3]", "beta", "Wd/W", "closure", "dW", "n"))
    for r in rows:
        print("   %-6s %10.4g %10.4g %11.4g %8.2g %11.4g %7.3f %6.3f %8.2g %6.2g %3d%s" %
              (r["stage"], r["H_amp"], r["B_peak"], r["W_loop"], r["W_err"], r["P_kW_m3"], r["beta"],
               r["balance"], r["closure"], r["dW_rel"], r["n_cycles"], "" if r["steady"] else " !"))


def pair(d1, d2):
    c1, r1 = load_run(d1)
    c2, r2 = load_run(d2)
    f1, f2 = c1["freq"], c2["freq"]
    ignore = {"freq", "out", "max_cycles", "stages", "git_commit", "torch", "device", "gpu",
              "no_resume", "timer", "vti"}
    diff = [k for k in c1 if k not in ignore and c1.get(k) != c2.get(k)]
    if diff:
        print("WARNING: runs differ in more than the frequency: %s" % diff)
    m2 = {r["stage"]: r for r in r2}
    out = []
    for a in r1:
        b = m2.get(a["stage"])
        if b is None:
            continue
        W_h = (f2 * a["W_loop"] - f1 * b["W_loop"]) / (f2 - f1)
        out.append({"stage": a["stage"], "B_peak": 0.5 * (a["B_peak"] + b["B_peak"]),
                    "dB_rel": abs(a["B_peak"] - b["B_peak"]) / a["B_peak"],
                    "W_f1": a["W_loop"], "W_f2": b["W_loop"], "W_h": W_h,
                    "dyn_share_f1": 1.0 - W_h / a["W_loop"]})
    out.sort(key=lambda r: r["B_peak"])
    add_beta(out, "W_h", "beta_h")
    add_beta(out, "W_f1", "beta_f1")
    print("\n== frequency separation  f1=%.3g Hz  f2=%.3g Hz" % (f1, f2))
    print("   %-6s %10s %8s %11s %11s %11s %8s %7s %7s" %
          ("stage", "B_pk [T]", "dB/B", "W(f1)", "W(f2)", "W_h", "dyn(f1)", "beta_h", "beta_f1"))
    for r in out:
        print("   %-6s %10.4g %8.2g %11.4g %11.4g %11.4g %8.2f %7.3f %7.3f" %
              (r["stage"], r["B_peak"], r["dB_rel"], r["W_f1"], r["W_f2"], r["W_h"],
               r["dyn_share_f1"], r["beta_h"], r["beta_f1"]))
    return out


def plot(results, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: no plot"); return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for name, cfg, rows in results:
        if not rows:
            continue
        B = [r["B_peak"] for r in rows]
        ax1.errorbar(B, [r["W_loop"] for r in rows], yerr=[r["W_err"] for r in rows],
                     fmt="o-", ms=3, capsize=2, label=label(cfg))
        ax1.set_xscale("log"); ax1.set_yscale("log")
        ax2.semilogx(B, [r["beta"] for r in rows], "o-", ms=3)
    ax1.set_xlabel("B_peak [T]"); ax1.set_ylabel("W per cycle [J/m$^3$]")
    ax2.set_xlabel("B_peak [T]"); ax2.set_ylabel(r"$\beta_\mathrm{eff}$ = d ln W / d ln B")
    for b in (2.0, 3.0):
        ax2.axhline(b, color="0.7", lw=0.8, ls="--")
    ax1.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print("plot: %s" % path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--pair", nargs=2, metavar=("RUN_F1", "RUN_F2"))
    ap.add_argument("--csv", default="all_runs.csv")
    ap.add_argument("--plot", default="beta_vs_B.png")
    a = ap.parse_args()

    if a.pair:
        pair(*a.pair)
        return

    results = []
    for d in a.runs:
        p = pathlib.Path(d)
        if not (p / "summary.json").exists():
            continue
        cfg, rows = load_run(p)
        results.append((p.name, cfg, rows))
        print_table(p.name, cfg, rows)

    if not results:
        sys.exit("no runs with summary.json found")
    keys = ["run", "d", "phi_vox", "sigma_rms", "xi", "K1", "freq", "precision", "dx",
            "seed_axes", "seed_stress"]
    with open(a.csv, "w", newline="") as f:
        w = csv.writer(f)
        rk = list(results[0][2][0].keys()) if results[0][2] else []
        w.writerow(keys + rk)
        for name, cfg, rows in results:
            for r in rows:
                w.writerow([name] + [cfg.get(k) for k in keys[1:]] + [r[k] for k in rk])
    print("\ncsv: %s" % a.csv)
    plot(results, a.plot)


if __name__ == "__main__":
    main()
