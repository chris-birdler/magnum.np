#!/usr/bin/env python
"""
Fine mesh at the corner K3 over 6 amplitudes (PROTOKOLL 7.51 option 2, 7.57): dx 3 vs dx 1.5 l_ex, seeds 2031 … 2036,
9 / 20 / 35 / 50 / 70 / 100 mT (9 / 50 / 100 mT: mesh check K3; 20 / 35 / 70 mT: fill-in with accumulators). One job
per amplitude from the relaxed state (init.pt of the 9 mT job), 5 cycles. A finished rerun X_mc2 (PLAN 4.5) replaces X.
  P      loss density at 30 MHz (W/cm³), steady-tail rule with at least 2 calm cycles (PROTOKOLL 7.54), else missing
  β      local loss exponent between neighbour amplitudes, per seed from the realised flux densities
  RF/RN  non-smooth / non-repeating share of the dissipation (eval_jumps_time.py), only with 3 calm cycles
  R1     share of the loss density in the local fundamental (acc.json; 20 / 35 / 70 mT only; the accumulators sum
         all measured cycles, events included). Reference: the Sobol corner runs K3_0 … K3_2 at dx 3 (ascending
         protocol, the same seeds 2031 … 2033)
Errors: SE over the seeds; paired differences dx 3 − dx 1.5 per seed.

    python eval_k3fill.py > results/k3_fine.md          # also results/k3_fine.png
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import eval_jumps_time as EJ

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE / "runs"
SEEDS = range(2031, 2037)
AMPS = (9, 20, 35, 50, 70, 100)
MESH = (("3", "dx 3 (14 nm)"), ("15", "dx 1.5 (7 nm)"))
F_HZ = 30e6


def folder(dx, s, b):
    d = RUNS / ("MK3_dx%s_s%d_b%d" % (dx, s, b))
    return AD.rerun_dir(d) or d


def loss(d):
    """(w, B in mT) of a job with the steady-tail rule (≥ 2 calm cycles), or None."""
    if not (d / "DONE").exists():
        return None
    _, rows = analyze.load_run(d)
    r = rows[0]
    return (r["w_loop"], 1500.0 * r["b_peak"]) if r["w_loop"] > 0 else None


def ms(v, f="%.2f"):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)])
    if len(v) < 2:
        return "n %d" % len(v)
    return (f + " ± " + f + " (%d)") % (v.mean(), v.std(ddof=1) / math.sqrt(len(v)), len(v))


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    kd = json.loads((folder("3", 2031, 9) / "config.json").read_text())["units"]["Kd"]
    to_P = kd * F_HZ / 1e6
    L = {(dx, s, b): loss(folder(dx, s, b)) for dx, _ in MESH for s in SEEDS for b in AMPS}
    reruns = sorted(d.name for d in RUNS.glob("MK3_*_mc2") if (d / "DONE").exists())
    print("# Fine mesh at K3 over 6 amplitudes (dx 3 vs dx 1.5)\n")
    print("K3: A 20 pJ/m at d = 1 µm (d/l_ex 212), L_eff 18, r_p 1.5, φ 0.69; 6 seeds. Steady-tail rule with at least 2 "
          "calm cycles; reruns used: %s. Missing jobs are left out (n per cell). Errors: SE over the seeds.\n"
          % (", ".join(reruns) or "none"))

    print("## Loss density P (W/cm³, 30 MHz)\n")
    print("| B (mT) | P dx 3, geometric mean | P dx 1.5 | P(dx 3) / P(dx 1.5), paired | missing (dx 3 / dx 1.5) |")
    print("|---|---|---|---|---|")
    curves = {}
    for b in AMPS:
        lp = {dx: [math.log(L[(dx, s, b)][0] * to_P) for s in SEEDS if L[(dx, s, b)]] for dx, _ in MESH}
        d = [math.log(L[("3", s, b)][0] / L[("15", s, b)][0]) for s in SEEDS if L[("3", s, b)] and L[("15", s, b)]]
        for dx, _ in MESH:
            curves.setdefault(dx, []).append(np.exp(np.mean(lp[dx])) if lp[dx] else np.nan)
        dd = np.array(d)
        rat = "× %.2f (± %.0f %%, n %d)" % (math.exp(dd.mean()), 100 * dd.std(ddof=1) / math.sqrt(len(dd)), len(dd)) \
            if len(dd) > 1 else "n %d" % len(dd)
        print("| %d | %.1f | %.1f | %s | %d / %d |" % (b, curves["3"][-1], curves["15"][-1], rat,
                                                      6 - len(lp["3"]), 6 - len(lp["15"])))

    print("\n## Local loss exponent β between neighbour amplitudes (per seed, realised flux densities)\n")
    print("| segment | β dx 3 | β dx 1.5 | difference dx 3 − dx 1.5, paired |")
    print("|---|---|---|---|")
    beta = {}
    for b0, b1 in zip(AMPS[:-1], AMPS[1:]):
        per = {}
        for dx, _ in MESH:
            per[dx] = {}
            for s in SEEDS:
                p0, p1 = L[(dx, s, b0)], L[(dx, s, b1)]
                if p0 and p1:
                    per[dx][s] = math.log(p1[0] / p0[0]) / math.log(p1[1] / p0[1])
        diff = [per["3"][s] - per["15"][s] for s in SEEDS if s in per["3"] and s in per["15"]]
        beta[(b0, b1)] = per
        print("| %d → %d mT | %s | %s | %s |" % (b0, b1, ms(per["3"].values()), ms(per["15"].values()), ms(diff, "%+.2f")))

    print("\n## Jumps (only jobs with 3 calm cycles) and R1\n")
    print("| B (mT) | RF dx 3 | RF dx 1.5 | seeds RF(dx 3) > RF(dx 1.5) (sign test p) | RN dx 3 | RN dx 1.5 | R1 dx 3 | "
          "R1 dx 1.5 | R1 Sobol K3_0 … 2 (dx 3, ascending) |")
    print("|---|---|---|---|---|---|---|---|---|")
    k3 = {}
    for n in ("K3_0", "K3_1", "K3_2"):
        st = json.loads((RUNS / n / "acc.json").read_text())["stages"]
        for k, s in enumerate(sorted(st, key=lambda x: x["B_peak_mT"])):
            k3.setdefault(AD.AMPS_MT[k], []).append(s["R1"])
    rf = {}
    for b in AMPS:
        M = {}
        for dx, _ in MESH:
            M[dx] = {}
            for s in SEEDS:
                d = folder(dx, s, b)
                if (d / "DONE").exists():
                    m = EJ.run_measures(d, 1)[0]
                    if m is not None and m["steady"]:
                        M[dx][s] = m
        pair = [s for s in SEEDS if s in M["3"] and s in M["15"]]
        gt = sum(M["3"][s]["RF"] > M["15"][s]["RF"] for s in pair)
        r1 = {}
        for dx, _ in MESH:
            v = []
            for s in SEEDS:
                f = folder(dx, s, b) / "acc.json"
                if f.exists():
                    v.append(json.loads(f.read_text())["stages"][0]["R1"])
            r1[dx] = ms(v) if v else "–"
        rf[b] = {dx: [M[dx][s]["RF"] for s in M[dx]] for dx, _ in MESH}
        print("| %d | %s | %s | %d/%d (%s) | %s | %s | %s | %s | %s |" % (
            b, ms(rf[b]["3"], "%.3f"), ms(rf[b]["15"], "%.3f"), gt, len(pair),
            "%.3f" % stats.binomtest(gt, len(pair), 0.5).pvalue if pair else "–",
            ms([M["3"][s]["RN"] for s in M["3"]], "%.3f"), ms([M["15"][s]["RN"] for s in M["15"]], "%.3f"),
            r1["3"], r1["15"], ms(k3.get(float(b), []))))

    fig_out(curves, beta, rf, to_P, L)


def fig_out(curves, beta, rf, to_P, L):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in pathlib.Path("/usr/share/fonts/truetype/crosextra").glob("Carlito*.ttf"):
        font_manager.fontManager.addfont(str(f))
    plt.rcParams.update({"font.family": "Carlito", "font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "legend.frameon": False})
    col = {"3": "#9B111E", "15": "#1F4E79"}
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    ax = axs[0]
    for dx, lab in MESH:
        for s in SEEDS:
            pts = [(L[(dx, s, b)][1], L[(dx, s, b)][0] * to_P) for b in AMPS if L[(dx, s, b)]]
            ax.loglog(*zip(*pts), "-", color=col[dx], alpha=0.25, lw=1)
        ax.loglog(AMPS, curves[dx], "o-", color=col[dx], lw=2.2, label=lab)
    ax.set_xticks(AMPS)
    ax.set_xticklabels([str(b) for b in AMPS])
    ax.set_xlabel("peak flux density B (mT)")
    ax.set_ylabel("loss density P at 30 MHz (W/cm³)")
    ax.set_title("Loss density at K3, thin lines = seeds")
    ax.legend()
    ax = axs[1]
    mids = [math.sqrt(b0 * b1) for b0, b1 in beta]
    for dx, lab in MESH:
        m = [np.mean(list(beta[k][dx].values())) if beta[k][dx] else np.nan for k in beta]
        e = [np.std(list(beta[k][dx].values()), ddof=1) / math.sqrt(len(beta[k][dx])) if len(beta[k][dx]) > 1 else np.nan
             for k in beta]
        ax.errorbar(mids, m, yerr=2 * np.array(e), fmt="o-", color=col[dx], lw=2, capsize=3, label=lab)
    ax.axhline(2.0, color="#666666", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xticks(AMPS)
    ax.set_xticklabels([str(b) for b in AMPS])
    ax.set_xlabel("peak flux density B (mT), segment centre")
    ax.set_ylabel("local loss exponent β")
    ax.set_title("Local β, mean ± 2 SE over seeds")
    ax.legend()
    ax = axs[2]
    for dx, lab in MESH:
        for b in AMPS:
            v = rf[b][dx]
            ax.scatter(np.full(len(v), b * (0.97 if dx == "3" else 1.03)), v, color=col[dx], s=18, alpha=0.6)
        ax.plot(AMPS, [np.median(rf[b][dx]) if rf[b][dx] else np.nan for b in AMPS], "-", color=col[dx], lw=2, label=lab)
    ax.set_xscale("log")
    ax.set_xticks(AMPS)
    ax.set_xticklabels([str(b) for b in AMPS])
    ax.set_xlabel("peak flux density B (mT)")
    ax.set_ylabel("non-smooth share RF of the dissipation")
    ax.set_title("Jumps: dots = seeds, line = median")
    ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / "results" / "k3_fine.png", dpi=200)


if __name__ == "__main__":
    main()
