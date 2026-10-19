#!/usr/bin/env python
"""
β from power-law fits P = k · B^β (least squares of ln P on ln B over all points of a range), corner K3 and base
point, dx 3 against dx 1.5 (PROTOKOLL 7.58). Data and rules as eval_k3fill.py: one job per amplitude, steady-tail rule
with at least 2 calm cycles, a finished rerun X_mc2 replaces X, realised flux densities.
Ranges: K3 9 … 100 mT (6 amplitudes), the windows of 3 neighbour amplitudes, and 9 / 50 / 100 mT (the amplitudes of
the base point); base point 9 / 50 / 100 mT. Per mesh the seeds with all points of the range on that mesh; the
difference between the meshes only on the seeds with all points on both meshes (paired).
  β of the mean   fit to the arithmetic mean of P over the seeds (the estimand of the design: β of E[P]); 95 % CI by
                  bootstrap over the seeds (2000 resamples, the same seeds for both meshes)
  β per seed      fit per seed, mean ± SE; difference dx 3 − dx 1.5 paired per seed (t-test)

    python eval_k3beta.py > results/k3_beta_fits.md        # also results/k3_beta_fits.png
"""
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import eval_k3fill as E

HERE = pathlib.Path(__file__).resolve().parent
N_BOOT = 2000
POINTS = {"K3 (A 20 pJ/m)": ("MK3_dx%s_s%d_b%d", range(2031, 2037), (9, 20, 35, 50, 70, 100)),
          "base (A 10 pJ/m)": ("M300_dx%s_s%d_b%d", range(921, 925), (9, 50, 100))}


def slope(B, P):
    return float(np.polyfit(np.log(B), np.log(P), 1)[0])


def compute(seed=20261011):
    """Rows of the fit table: list of dicts (point, range, per mesh: β of the mean, CI, per-seed values)."""
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    rng = np.random.default_rng(seed)
    rows = []
    for lab, (fmt, seeds, amps) in POINTS.items():
        L = {(dx, s, b): E.loss(AD.rerun_dir(E.RUNS / (fmt % (dx, s, b))) or E.RUNS / (fmt % (dx, s, b)))
             for dx in ("3", "15") for s in seeds for b in amps}
        ranges = [amps] if len(amps) < 4 else \
            [amps] + [amps[i:i + 3] for i in range(len(amps) - 2)] + [(9, 50, 100)]
        for rg in ranges:
            row = {"point": lab, "range": list(rg)}
            for dx in ("3", "15"):
                ok = [s for s in seeds if all(L[(dx, s, b)] for b in rg)]
                W = np.array([[L[(dx, s, b)][0] for b in rg] for s in ok])
                Bm = np.array([[L[(dx, s, b)][1] for b in rg] for s in ok])
                boots = []
                for _ in range(N_BOOT):
                    i = rng.integers(0, len(ok), len(ok))
                    boots.append(slope(Bm[i].mean(0), W[i].mean(0)))
                row[dx] = {"mean": slope(Bm.mean(0), W.mean(0)), "ci": [float(x) for x in np.percentile(boots, [2.5, 97.5])],
                           "per_s": {int(s): slope(Bm[k], W[k]) for k, s in enumerate(ok)}}
            both = [s for s in seeds if s in row["3"]["per_s"] and s in row["15"]["per_s"]]
            row["diff"] = [row["3"]["per_s"][s] - row["15"]["per_s"][s] for s in both]
            rows.append(row)
    return rows


def main():
    print("# β from power-law fits, dx 3 against dx 1.5\n")
    print("Fit P = k · B^β by least squares of ln P on ln B over all amplitudes of the range (realised B). Steady-tail "
          "rule (≥ 2 calm cycles). Per mesh the seeds with all points of the range on that mesh; the difference "
          "between the meshes only on the common seeds (paired). β of the mean: fit to the arithmetic mean of P over the "
          "seeds, 95 % CI by bootstrap over the seeds. β per seed: mean ± SE; difference dx 3 − dx 1.5 paired per seed "
          "(t-test).\n")
    print("| point | range (points) | seeds dx 3 / dx 1.5 / both | β of the mean, dx 3 | β of the mean, dx 1.5 | "
          "β per seed, dx 3 | β per seed, dx 1.5 | dx 3 − dx 1.5, paired (p) |")
    print("|---|---|---|---|---|---|---|---|")
    rows = compute()
    se = lambda v: np.std(v, ddof=1) / math.sqrt(len(v)) if len(v) > 1 else float("nan")
    res = {}
    for row in rows:
        rg = tuple(row["range"])
        d = np.array(row["diff"])
        p = stats.ttest_1samp(d, 0.0).pvalue if len(d) > 1 else float("nan")
        per = {dx: np.array(list(row[dx]["per_s"].values())) for dx in ("3", "15")}
        k3 = row["point"].startswith("K3")
        name = ("%d … %d" % (rg[0], rg[-1])) if rg != (9, 50, 100) or not k3 else "9 / 50 / 100"
        print("| %s | %s mT (%d) | %d / %d / %d | %.2f (%.2f … %.2f) | %.2f (%.2f … %.2f) | %.2f ± %.2f | %.2f ± %.2f | "
              "%+.2f ± %.2f (%.3f) |" % (row["point"], name, len(rg), len(per["3"]), len(per["15"]), len(d),
                                        row["3"]["mean"], *row["3"]["ci"], row["15"]["mean"], *row["15"]["ci"],
                                        per["3"].mean(), se(per["3"]), per["15"].mean(), se(per["15"]),
                                        d.mean() if len(d) else float("nan"), se(d), p))
        res[(row["point"], rg)] = {dx: {"mean": row[dx]["mean"], "ci": row[dx]["ci"]} for dx in ("3", "15")}
    print("\nA seed with a missing job (no 2 calm cycles; rerun pending) is left out of that range only. The windows "
          "overlap, so they are not independent. With 2 … 5 seeds the bootstrap interval is coarse (few distinct "
          "resamples); the SE of the per-seed values is the more honest error.")
    fig_out(res)


def fig_out(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in pathlib.Path("/usr/share/fonts/truetype/crosextra").glob("Carlito*.ttf"):
        font_manager.fontManager.addfont(str(f))
    plt.rcParams.update({"font.family": "Carlito", "font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "legend.frameon": False})
    col = {"3": "#9B111E", "15": "#1F4E79"}
    lab = {"3": "dx 3 (14 nm)", "15": "dx 1.5 (7 nm)"}
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    k3 = [(rg, r) for (p, rg), r in res.items() if p.startswith("K3")]
    full = [x for x in k3 if len(x[0]) == 6][0]
    k3 = [x for x in k3 if x[0] != (9, 50, 100)]
    for j, dx in enumerate(("3", "15")):
        off = 0.97 if dx == "3" else 1.03
        win = [(rg, r) for rg, r in k3 if len(rg) == 3]
        xc = [math.exp(np.mean(np.log(rg))) * off for rg, _ in win]
        m = [r[dx]["mean"] for _, r in win]
        lo = [r[dx]["mean"] - r[dx]["ci"][0] for _, r in win]
        hi = [r[dx]["ci"][1] - r[dx]["mean"] for _, r in win]
        ax.errorbar(xc, m, yerr=[lo, hi], fmt="o-", color=col[dx], lw=2, capsize=4, label="%s, 3-point windows" % lab[dx])
        f = full[1][dx]
        ax.axhspan(f["ci"][0], f["ci"][1], color=col[dx], alpha=0.10)
        ax.axhline(f["mean"], color=col[dx], lw=1.2, ls="--", label="%s, 9 … 100 mT: β = %.2f" % (lab[dx], f["mean"]))
    ax.axhline(2.0, color="#666666", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xticks([9, 20, 35, 50, 70, 100])
    ax.set_xticklabels(["9", "20", "35", "50", "70", "100"])
    ax.set_xlabel("peak flux density B (mT), window centre")
    ax.set_ylabel("loss exponent β (fit of the mean loss density)")
    ax.set_title("K3 (A 20 pJ/m): β from power-law fits, 95 % CI by bootstrap over seeds", fontsize=12)
    ax.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(HERE / "results" / "k3_beta_fits.png", dpi=200)


if __name__ == "__main__":
    main()
