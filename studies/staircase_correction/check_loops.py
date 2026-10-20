#!/usr/bin/env python
"""
Design B of the dynamic test (PROTOKOLL 3.16): where does the loss of each arm form? Loop energy and phase lag of each
particle in the last cycle (samples_00_amp00.csv of run_loops: H, M_par and the mean M_par of each particle).

- w_k: the share of particle k in the loop energy (∮ H dM_par,k, normalised so that the shares sum to 1), times the
  loop energy w of the box. The particle means of the arm with the fraction surface count each cell fully, so only
  the shares are used.
- Δ_c,k = ln w_k(T3Pc) − ln w_k(T3), Δ_g,k = ln w_k(D15) − ln w_k(T15g): the effects of 3.15 for each particle.
- Fundamental of M_par,k(t): amplitude, phase lag behind H, and the RMS of the rest relative to the amplitude.

    python check_loops.py --runs RUNS --out results   (writes loops_B.md, loops_B.json, loops_B.png)
"""
import argparse
import json
import pathlib

import numpy as np

MU0 = 4e-7 * np.pi
SEEDS = (941, 942, 943, 944)
ARMS = ("T3", "T3Pc", "D15", "T15g")


def last_cycle(path):
    a = np.loadtxt(path, delimiter=",", comments="#")
    return a[a[:, 0] == a[:, 0].max()]


def analyse(b):
    """loop energy w (J/m³), particle loop energies w_k (shares times w), fundamental of each particle."""
    H, M, P = b[:, 2], b[:, 3], b[:, 9:13]
    Hc, Mc, Pc = np.append(H, H[0]), np.append(M, M[0]), np.vstack([P, P[:1]])
    Hm = 0.5 * (Hc[1:] + Hc[:-1])
    w = MU0 * np.sum(Hm * np.diff(Mc))
    wk = MU0 * np.sum(Hm[:, None] * np.diff(Pc, axis=0), axis=0)
    wk = wk / wk.sum() * w
    n = len(b)
    fH = np.fft.rfft(H)[1]
    fund = []
    for k in range(4):
        F = np.fft.rfft(P[:, k])
        one = np.zeros_like(F)
        one[1] = F[1]
        amp = 2 * abs(F[1]) / n
        rest = P[:, k] - P[:, k].mean() - np.fft.irfft(one, n)
        fund.append({"amp": float(amp), "lag": float(np.angle(fH) - np.angle(F[1])), "rest": float(rest.std() / amp)})
    T = (b[1, 1] - b[0, 1]) * n
    return {"w": float(w), "wk": [float(x) for x in wk], "fund": fund, "H0": float(2 * abs(fH) / n), "f": float(1 / T)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--out", default=str(pathlib.Path(__file__).resolve().parent / "results"))
    a = ap.parse_args()
    root, out = pathlib.Path(a.runs), pathlib.Path(a.out)
    res = {"%s_s%d" % (arm, s): analyse(last_cycle(root / ("%s_s%d" % (arm, s)) / "samples_00_amp00.csv"))
           for s in SEEDS for arm in ARMS}
    lines = ["# Design B: loss of each particle in the last cycle (PROTOKOLL 3.16)\n",
             "Drive: H0 %.0f … %.0f A/m, f %.1f MHz. w_k in 10⁻³ J/m³; lag: phase lag of the fundamental of M_par,k "
             "behind H in mrad; rest: RMS of M_par,k without its mean and fundamental, relative to the amplitude.\n"
             % (min(r["H0"] for r in res.values()), max(r["H0"] for r in res.values()),
                np.mean([r["f"] for r in res.values()]) / 1e6),
             "| seed | arm | w_k, particles 0 / 1 / 2 / 3 | lag (mrad) 0 / 1 / 2 / 3 | rest (max) |", "|---|---|---|---|---|"]
    dc, dg, sd = [], [], []
    for s in SEEDS:
        for arm in ARMS:
            r = res["%s_s%d" % (arm, s)]
            lines.append("| %d | %s | %s | %s | %.3f |" % (
                s, arm, " / ".join("%.2f" % (1e3 * x) for x in r["wk"]),
                " / ".join("%.2f" % (1e3 * f["lag"]) for f in r["fund"]), max(f["rest"] for f in r["fund"])))
        g = {arm: np.array(res["%s_s%d" % (arm, s)]["wk"]) for arm in ARMS}
        dc += list(np.log(g["T3Pc"] / g["T3"]))
        dg += list(np.log(g["D15"] / g["T15g"]))
        sd += [s] * 4
    dc, dg = np.array(dc), np.array(dg)
    lines += ["", "| seed | particle | Δ_c,k | Δ_g,k | d_k = Δ_c,k − Δ_g,k |", "|---|---|---|---|---|"]
    for i in range(len(dc)):
        lines.append("| %d | %d | %+.2f | %+.2f | %+.2f |" % (sd[i], i % 4, dc[i], dg[i], dc[i] - dg[i]))
    summ = {"corr": float(np.corrcoef(dc, dg)[0, 1]), "rms_d": float(np.sqrt(np.mean((dc - dg) ** 2))),
            "rms_g": float(np.sqrt(np.mean(dg ** 2))), "rms_c_none": float(np.sqrt(np.mean(dg ** 2)))}
    lines += ["", "Over the 16 particle-seed pairs (not independent: the particles of one box couple through the stray "
              "field): correlation of Δ_c,k and Δ_g,k %.2f; RMS of d_k %.2f; RMS of Δ_g,k %.2f (= RMS of d_k without "
              "a correction)." % (summ["corr"], summ["rms_d"], summ["rms_g"])]
    (out / "loops_B.md").write_text("\n".join(lines) + "\n")
    (out / "loops_B.json").write_text(json.dumps({"runs": res, "summary": summ}, indent=1))
    print("\n".join(lines))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    lim = 1.05 * max(np.abs(dc).max(), np.abs(dg).max())
    ax.plot([-lim, lim], [-lim, lim], color="0.4", lw=1.0, label="correction exact")
    ax.axhline(0.0, color="0.75", lw=1.0, ls="--", label="no correction")
    for s, c in zip(SEEDS, ("tab:blue", "tab:orange", "tab:green", "tab:red")):
        m = np.array(sd) == s
        ax.scatter(dg[m], dc[m], color=c, s=36, label="seed %d" % s, zorder=3)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\Delta_{\mathrm{g},k}$: fine staircase at dx 1.5 (ln $w_k$)")
    ax.set_ylabel(r"$\Delta_{\mathrm{c},k}$: correction at dx 3 (ln $w_k$)")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "loops_B.png", dpi=160)


if __name__ == "__main__":
    main()
