#!/usr/bin/env python
"""
Evaluate runs of run_loops.py in reduced units:
w = W/K_d per cycle versus b_peak = B_peak/Js, and beta_eff = d ln w / d ln b.

    python analyze.py runs/*                          # tables, features, all_runs.csv, plot
    python analyze.py --pair runs/P_kappa00800 runs/C_fhalf   # frequency separation

Per run and amplitude, w and b_peak are the MEAN over all cycles except the
first one (transient). w_err is the standard error of that mean.

Flags per amplitude:
  !   stage not steady (max_cycles_per_amp reached)
  *   energy balance |w_dis/w_loop - 1| > 2 %
  #   jump in the number of large-angle neighbour pairs during the measured
      cycles (> 20 % change): vortex-core / Bloch-point events -> mesh dependent

Features per run (for the ranking of the parameters): beta at fixed
b = 0.01, 0.03, 0.1, 0.3 (log interpolation), beta_max and b at beta_max.

Frequency separation (--pair): W(f) = W_h + c f  ->  W_h = (f2 W1 - f1 W2)/(f2 - f1).
"""
import argparse
import csv
import json
import math
import pathlib
import sys

import numpy as np

B_REF = (0.01, 0.03, 0.1, 0.3)


def load_run(d):
    d = pathlib.Path(d)
    cfg = json.loads((d / "config.json").read_text())
    summ = json.loads((d / "summary.json").read_text())
    rows = []
    for st in summ["stages"]:
        if not st.get("measure") or not st.get("complete") or not st["cycles"]:
            continue
        cyc = st["cycles"][1:] if len(st["cycles"]) > 1 else st["cycles"]
        w = np.array([c["w_loop"] for c in cyc])
        n60 = [c.get("n_pairs_gt60", 0) for c in cyc]
        rows.append({"stage": st["name"], "h_amp": st["h_amp"],
                     "b_peak": float(np.mean([c["b_peak"] for c in cyc])),
                     "w_loop": float(w.mean()),
                     "w_err": float(w.std(ddof=1) / math.sqrt(len(w))) if len(w) > 1 else float("nan"),
                     "w_dis": float(np.mean([c["w_dis"] for c in cyc])),
                     "closure": cyc[-1]["closure"], "dW_rel": cyc[-1].get("dW_rel", float("nan")),
                     "steady": st.get("steady", False), "n_cycles": len(st["cycles"]),
                     "frac_gt60": float(max(c.get("frac_pairs_gt60", 0.0) for c in cyc)),
                     "n60_jump": bool(max(n60) > 1.2 * max(min(n60), 1)),
                     "mean_dt": float(np.mean([c["mean_dt"] for c in cyc]))})
    rows.sort(key=lambda r: r["b_peak"])
    add_beta(rows, "w_loop", "beta")
    for r in rows:
        r["balance"] = r["w_dis"] / r["w_loop"] if r["w_loop"] != 0 else float("nan")
    return cfg, rows


def add_beta(rows, key, out_key):
    n = len(rows)
    lb = [math.log(r["b_peak"]) for r in rows]
    lw = [math.log(r[key]) if r[key] > 0 else float("nan") for r in rows]
    for i in range(n):
        if n < 2:
            rows[i][out_key] = float("nan"); continue
        j0, j1 = max(i - 1, 0), min(i + 1, n - 1)
        rows[i][out_key] = (lw[j1] - lw[j0]) / (lb[j1] - lb[j0])


def features(rows):
    f = {}
    if len(rows) < 2:
        return f
    lb = np.log([r["b_peak"] for r in rows])
    be = np.array([r["beta"] for r in rows])
    for b in B_REF:
        lbr = math.log(b)
        f["beta@b=%g" % b] = float(np.interp(lbr, lb, be)) if lb[0] <= lbr <= lb[-1] else float("nan")
    i = int(np.nanargmax(be))
    f["beta_max"], f["b_at_beta_max"] = float(be[i]), float(rows[i]["b_peak"])
    return f


def label(cfg):
    return ("d/lex=%g phi=%.3f Q=%g(%s) kappa=%g xi/lex=%g dx/lex=%g f/fM=%.3g seed=%d" %
            (cfg["d_lex"], cfg["phi_vox"], cfg["Q"], cfg["aniso"][:3], cfg["kappa"], cfg["xi_lex"],
             cfg["dx_lex"], cfg["f_rel"], cfg["seed"]))


def flags(r):
    s = "" if r["steady"] else "!"
    if abs(r["balance"] - 1.0) > 0.02:
        s += "*"
    if r["n60_jump"]:
        s += "#"
    return s


def print_table(name, cfg, rows):
    print("\n== %s\n   %s" % (name, label(cfg)))
    print("   %-6s %9s %9s %10s %8s %7s %6s %8s %7s %3s %s" %
          ("stage", "h", "b_peak", "w_loop", "w_err", "beta", "wd/w", "closure", "f>60", "n", "flags"))
    for r in rows:
        print("   %-6s %9.4g %9.4g %10.4g %8.2g %7.3f %6.3f %8.2g %7.2g %3d %s" %
              (r["stage"], r["h_amp"], r["b_peak"], r["w_loop"], r["w_err"], r["beta"], r["balance"],
               r["closure"], r["frac_gt60"], r["n_cycles"], flags(r)))
    f = features(rows)
    if f:
        print("   features: " + "  ".join("%s=%.3g" % kv for kv in f.items()))


def pair(d1, d2):
    c1, r1 = load_run(d1)
    c2, r2 = load_run(d2)
    f1, f2 = c1["f_rel"], c2["f_rel"]
    ignore = {"f_rel", "freq", "out", "max_cycles", "stages", "git_commit", "torch", "device", "gpu",
              "no_resume", "timer", "vti", "keep_stage_ckpt", "init_from"}
    diff = [k for k in c1 if k not in ignore and c1.get(k) != c2.get(k)]
    if diff:
        print("WARNING: runs differ in more than the frequency: %s" % diff)
    m2 = {r["stage"]: r for r in r2}
    out = []
    for a in r1:
        b = m2.get(a["stage"])
        if b is None:
            continue
        w_h = (f2 * a["w_loop"] - f1 * b["w_loop"]) / (f2 - f1)
        out.append({"stage": a["stage"], "b_peak": 0.5 * (a["b_peak"] + b["b_peak"]),
                    "db_rel": abs(a["b_peak"] - b["b_peak"]) / a["b_peak"],
                    "w_f1": a["w_loop"], "w_f2": b["w_loop"], "w_h": w_h,
                    "dyn_share_f1": 1.0 - w_h / a["w_loop"]})
    out.sort(key=lambda r: r["b_peak"])
    add_beta(out, "w_h", "beta_h")
    add_beta(out, "w_f1", "beta_f1")
    print("\n== frequency separation  f1/fM=%.3g  f2/fM=%.3g" % (f1, f2))
    print("   %-6s %9s %7s %10s %10s %10s %8s %7s %7s" %
          ("stage", "b_peak", "db/b", "w(f1)", "w(f2)", "w_h", "dyn(f1)", "beta_h", "beta_f1"))
    for r in out:
        print("   %-6s %9.4g %7.2g %10.4g %10.4g %10.4g %8.2f %7.3f %7.3f" %
              (r["stage"], r["b_peak"], r["db_rel"], r["w_f1"], r["w_f2"], r["w_h"],
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
        b = [r["b_peak"] for r in rows]
        ax1.errorbar(b, [r["w_loop"] for r in rows], yerr=[r["w_err"] for r in rows],
                     fmt="o-", ms=3, capsize=2, label=name)
        ax2.semilogx(b, [r["beta"] for r in rows], "o-", ms=3, label=name)
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("b = B_peak / Js"); ax1.set_ylabel("w = W / K_d per cycle")
    ax2.set_xlabel("b = B_peak / Js"); ax2.set_ylabel(r"$\beta_\mathrm{eff}$ = d ln w / d ln b")
    for v in (2.0, 3.0):
        ax2.axhline(v, color="0.7", lw=0.8, ls="--")
    ax1.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print("plot: %s" % path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--pair", nargs=2, metavar=("RUN_F1", "RUN_F2"))
    ap.add_argument("--csv", default="all_runs.csv")
    ap.add_argument("--plot", default="beta_vs_b.png")
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

    keys = ["d_lex", "phi_vox", "Q", "aniso", "kappa", "xi_lex", "dx_lex", "f_rel", "seed"]
    rk = [k for k in results[0][2][0].keys()] if results[0][2] else []
    with open(a.csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["run"] + keys + rk)
        for name, cfg, rows in results:
            for r in rows:
                w.writerow([name] + [cfg.get(k) for k in keys] + [r[k] for k in rk])
    print("\ncsv: %s" % a.csv)
    plot(results, a.plot)


if __name__ == "__main__":
    main()
