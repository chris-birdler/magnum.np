#!/usr/bin/env python
"""
Absolute loss against β: do the high-β configurations have the lower loss, and where do the loss curves cross?
Exploratory (hints), α 0.1 Sobol runs (PROTOKOLL 7.47).
  (A) data: the runs split into thirds by their β all; mean loss P = f E[w] (W/cm³, core volume) per nominal
      amplitude per third, 95 % CI by bootstrap over the runs; ratio high-β / low-β third per amplitude
  (B) model: main fit (Gamma GLM per amplitude, 15 terms, as analyze_design) → range effect of each factor on
      ln E[w] at all 7 amplitudes, 95 % CI by bootstrap; the amplitude where the sign changes is the crossing
  figure results/crossover.png

    python eval_crossover.py > results/crossover.md
"""
import json
import math
import pathlib

import numpy as np

import analyze_design as AD

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE / "runs"
F_HZ = 30e6
N_BOOT = 1000


def main():
    recs, _ = AD.load_runs(RUNS)
    kd = json.loads((RUNS / recs[0]["name"] / "config.json").read_text())["units"]["Kd"]
    W = np.exp(np.array([[r["lnw@%gmT" % m] for m in AD.AMPS_MT] for r in recs]))     # reduced (w / K_d)
    P = W * kd * F_HZ / 1e6                                                              # W/cm³
    b = np.array([r["beta_all"] for r in recs])
    rng = np.random.default_rng(20261010)
    n = len(recs)
    print("# Absolute loss against β (exploratory, n = %d α 0.1 runs)\n" % n)

    # (A) thirds by β all
    q1, q2 = np.quantile(b, [1 / 3, 2 / 3])
    groups = {"low β": b <= q1, "mid β": (b > q1) & (b <= q2), "high β": b > q2}
    print("## A. Runs split into thirds by β all\n")
    print("β all: low ≤ %.2f < mid ≤ %.2f < high; group means %s.\n" % (q1, q2, ", ".join(
        "%s %.2f (n %d)" % (k, b[g].mean(), g.sum()) for k, g in groups.items())))
    lo, hi = groups["low β"], groups["high β"]
    print("| B̂ (mT) | P low β (W/cm³) | P mid β | P high β | ratio high / low (95 % CI) |")
    print("|---|---|---|---|---|")
    ratios = []
    for k, m in enumerate(AD.AMPS_MT):
        bs = []
        for _ in range(N_BOOT):
            il = rng.choice(np.where(lo)[0], lo.sum())
            ih = rng.choice(np.where(hi)[0], hi.sum())
            bs.append(P[ih, k].mean() / P[il, k].mean())
        r = P[hi, k].mean() / P[lo, k].mean()
        ci = np.percentile(bs, [2.5, 97.5])
        ratios.append((m, r, ci))
        print("| %g | %.1f | %.1f | %.1f | %.2f (%.2f … %.2f) |" % (m, P[lo, k].mean(), P[groups["mid β"], k].mean(),
                                                               P[hi, k].mean(), r, *ci))
    cross = [(ratios[i][0], ratios[i + 1][0]) for i in range(len(ratios) - 1) if (ratios[i][1] - 1) * (ratios[i + 1][1] - 1) < 0]
    print("\nCrossing of the mean loss curves (ratio = 1) between: %s.\n" % (", ".join("%g and %g mT" % c for c in cross) or "none"))

    # (B) model: range effect of each factor on ln E[w] at every amplitude
    X = AD.coded(recs)
    T = AD.terms(X)
    est = AD.glm_coefs(T, W)

    def eff(coefs):
        out = {}
        for j, f in enumerate(AD.FACTORS):
            Th, Tl = AD.contrast(X, j)
            out[f] = (Th @ coefs.T - Tl @ coefs.T).mean(0)           # [7] range effect on ln E[w]
        return out
    e0 = eff(est)
    boots = []
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        try:
            boots.append(eff(AD.glm_coefs(T[i], W[i])))
        except Exception:                                          # noqa: BLE001
            pass
    print("## B. Main fit: range effect of each factor on ln E[w] per amplitude (95 %% CI, %d resamples)\n" % len(boots))
    print("| B̂ (mT) | " + " | ".join(AD.FACTORS) + " |")
    print("|---" * (len(AD.FACTORS) + 1) + "|")
    for k, m in enumerate(AD.AMPS_MT):
        cells = []
        for f in AD.FACTORS:
            bb = np.array([x[f][k] for x in boots])
            lo_, hi_ = np.percentile(bb, [2.5, 97.5])
            cells.append("%+.2f (%+.2f … %+.2f)%s" % (e0[f][k], lo_, hi_, " *" if lo_ > 0 or hi_ < 0 else ""))
        print("| %g | %s |" % (m, " | ".join(cells)))
    print("\n\\* CI excludes 0. φ contains the dilution ln(0.69/0.55) = 0.23. Range effect +0.30 ≈ +35 % loss.\n")
    for f in AD.FACTORS:
        s = e0[f]
        cr = [(AD.AMPS_MT[i], AD.AMPS_MT[i + 1]) for i in range(6) if s[i] * s[i + 1] < 0]
        if cr:
            # linear interpolation in ln B for the zero crossing
            i = AD.AMPS_MT.index(cr[0][0])
            x0, x1 = math.log(AD.AMPS_MT[i]), math.log(AD.AMPS_MT[i + 1])
            z = math.exp(x0 + (x1 - x0) * s[i] / (s[i] - s[i + 1]))
            print("- %s: sign change between %g and %g mT (zero at ≈ %.0f mT)" % (f, *cr[0], z))
        else:
            print("- %s: no sign change in 9 … 150 mT" % f)

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), dpi=150)
    for k_, g in groups.items():
        ax[0].loglog(AD.AMPS_MT, P[g].mean(0), "o-", label="%s (β all %.2f, n %d)" % (k_, b[g].mean(), g.sum()))
    ax[0].set_xlabel("B̂ (mT)")
    ax[0].set_ylabel("P (W/cm³) at 30 MHz")
    ax[0].set_title("Mean loss of the β thirds")
    ax[0].legend(fontsize=8)
    ax[0].grid(True, which="both", alpha=0.3)
    for f in AD.FACTORS:
        bb = np.array([x[f] for x in boots])
        lo_, hi_ = np.percentile(bb, [2.5, 97.5], axis=0)
        ax[1].plot(AD.AMPS_MT, e0[f], "o-", label=f)
        ax[1].fill_between(AD.AMPS_MT, lo_, hi_, alpha=0.15)
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("B̂ (mT)")
    ax[1].set_ylabel("range effect on ln E[w]")
    ax[1].set_title("Main fit: factor effect on the loss vs B̂")
    ax[1].legend(fontsize=8)
    ax[1].grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "results" / "crossover.png")


if __name__ == "__main__":
    main()
