#!/usr/bin/env python
"""
Mesh check at the design corner K3 (PLAN 4.0 rule A) against the base point (mesh300, PROTOKOLL 7.25):
  K3    d/l_ex 212 (A = 20 pJ/m at d = 1 µm), L_eff 18, r_p 1.5, φ 0.69, seeds 2031 … 2036
  base  d/l_ex 300 (A = 10 pJ/m), L_eff 12, r_p 0, φ 0.65, seeds 921 … 924
  each: dx 3 and dx 1.5 l_ex, 9 / 50 / 100 mT (one job per amplitude, 5 cycles, the 9 mT job relaxes)
Outputs: per seed Δ ln w = ln w(dx 3) − ln w(dx 1.5), mean ± SE, the contrast K3 − base (Welch), the mean loss curves
of both meshes, and the figure results/k3_mesh.png. Works with incomplete data (pairs that exist).

    python eval_k3.py > results/k3_mesh.md
"""
import json
import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
import numpy as np                                                 # noqa: E402
from matplotlib import font_manager                                # noqa: E402
from scipy import stats                                            # noqa: E402

import analyze                                                     # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE / "runs"
AMPS = (9, 50, 100)
SETS = {"K3 (A 20 pJ/m, small particles)": ("MK3_dx%s_s%d_b%d", range(2031, 2037)),
        "base point (A 10 pJ/m)": ("M300_dx%s_s%d_b%d", range(921, 925))}
F_HZ = 30e6
for f in pathlib.Path("/usr/share/fonts/truetype/crosextra").glob("Carlito*.ttf"):
    font_manager.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": "Carlito", "font.size": 12, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False})
RED, BLUE, GREY = "#9B111E", "#1F4E79", "#666666"


def w(name):
    d = RUNS / name
    if not (d / "DONE").exists():
        return None
    _, rows = analyze.load_run(d)
    return rows[0]["w_loop"] if rows and rows[0]["w_loop"] > 0 else None


def wb(name):
    d = RUNS / name
    if not (d / "DONE").exists():
        return None
    _, rows = analyze.load_run(d)
    r = rows[0]
    return (r["w_loop"], 1500.0 * r["b_peak"]) if r["w_loop"] > 0 else None


SEGS = ("9 → 50 mT", "50 → 100 mT", "9 → 100 mT (fit)")


def seed_betas(fmt, s, dx):
    """β of one seed and mesh in the 2 segments and over the 3 amplitudes (realised peak flux densities)."""
    p = {b: wb(fmt % (dx, s, b)) for b in AMPS}
    out = {}
    if p[9] and p[50]:
        out[SEGS[0]] = math.log(p[50][0] / p[9][0]) / math.log(p[50][1] / p[9][1])
    if p[50] and p[100]:
        out[SEGS[1]] = math.log(p[100][0] / p[50][0]) / math.log(p[100][1] / p[50][1])
    if all(p.values()):
        out[SEGS[2]] = float(np.polyfit(np.log([p[b][1] for b in AMPS]), np.log([p[b][0] for b in AMPS]), 1)[0])
    return out


def fig_beta():
    """β per seed, coarse vs fine mesh, K3 and base point, 3 panels (segments); lines join the same seed."""
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.6), sharey=True)
    xpos = {(0, "3"): 0, (0, "15"): 1, (1, "3"): 2.6, (1, "15"): 3.6}
    for ax, seg in zip(axs, SEGS):
        for j, (lab, (fmt, seeds)) in enumerate(SETS.items()):
            col = RED if j == 0 else GREY
            vals = {dx: [] for dx in ("3", "15")}
            for sd in seeds:
                b3, b15 = seed_betas(fmt, sd, "3").get(seg), seed_betas(fmt, sd, "15").get(seg)
                if b3 is not None:
                    vals["3"].append(b3)
                if b15 is not None:
                    vals["15"].append(b15)
                if b3 is not None and b15 is not None:
                    ax.plot([xpos[(j, "3")], xpos[(j, "15")]], [b3, b15], "-", color=col, alpha=0.35, lw=1)
                for dx, b in (("3", b3), ("15", b15)):
                    if b is not None:
                        ax.scatter(xpos[(j, dx)], b, color=col, s=26, alpha=0.75, zorder=3)
            for dx in ("3", "15"):
                v = np.array(vals[dx])
                if len(v) > 1:
                    ax.errorbar(xpos[(j, dx)] + 0.18, v.mean(), yerr=2 * v.std(ddof=1) / math.sqrt(len(v)), fmt="s",
                                color=col, ms=8, capsize=4, lw=1.6, zorder=4)
                    ax.text(xpos[(j, dx)] + 0.3, v.mean(), "%.2f" % v.mean(), fontsize=10, color=col, va="center")
        ax.set_xticks([0, 1, 2.6, 3.6])
        ax.set_xticklabels(["K3\ncoarse", "K3\nfine", "base\ncoarse", "base\nfine"])
        ax.set_xlim(-0.4, 4.3)
        ax.axhline(2.0, color=GREY, lw=0.8, ls="--")
        ax.set_title("β " + seg)
    axs[0].set_ylabel("loss exponent β")
    fig.text(0.5, 0.005, "coarse = dx 3 l_ex (14 / 10 nm), fine = dx 1.5 l_ex (7 / 5 nm). Dots = seeds, lines join the "
             "same seed, squares = mean ± 2 standard errors. K3: A 20 pJ/m (small particles), base: A 10 pJ/m.",
             ha="center", fontsize=10, color=GREY)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(HERE / "results" / "k3_beta.png", dpi=200)


def main():
    kd = None
    res = {}
    for lab, (fmt, seeds) in SETS.items():
        r = {"dlnw": {}, "lnw": {}}
        for b in AMPS:
            pairs, l3, l15 = [], [], []
            for s in seeds:
                a, c = w(fmt % ("3", s, b)), w(fmt % ("15", s, b))
                if a:
                    l3.append(math.log(a))
                if c:
                    l15.append(math.log(c))
                if a and c:
                    pairs.append(math.log(a / c))
            r["dlnw"][b] = pairs
            r["lnw"][b] = (l3, l15)
        res[lab] = r
        if kd is None:
            cfg = json.loads((RUNS / (fmt % ("3", list(seeds)[0], 9)) / "config.json").read_text())
            kd = cfg["units"]["Kd"]
    labs = list(SETS)
    print("# Mesh check K3 against the base point (dx 3 vs dx 1.5)\n")
    print("Δ ln w = ln w(dx 3) − ln w(dx 1.5) per seed, mean ± standard error; contrast K3 − base with Welch.\n")
    print("| B (mT) | K3: per seed | K3 mean ± SE (n) | base: per seed | base mean ± SE (n) | K3 − base (p) |")
    print("|---|---|---|---|---|---|")
    for b in AMPS:
        k, bb = np.array(res[labs[0]]["dlnw"][b]), np.array(res[labs[1]]["dlnw"][b])
        ms = lambda v: "%+.2f ± %.2f (%d)" % (v.mean(), v.std(ddof=1) / math.sqrt(len(v)), len(v)) if len(v) > 1 else "n %d" % len(v)
        con = ""
        if len(k) > 1 and len(bb) > 1:
            con = "%+.2f ± %.2f (%.3f)" % (k.mean() - bb.mean(), math.sqrt(k.var(ddof=1) / len(k) + bb.var(ddof=1) / len(bb)),
                                          stats.ttest_ind(k, bb, equal_var=False).pvalue)
        print("| %d | %s | %s | %s | %s | %s |" % (b, " ".join("%+.2f" % x for x in k), ms(k), " ".join("%+.2f" % x for x in bb),
                                               ms(bb), con))

    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6))
    ax = axs[0]
    to_P = kd * F_HZ / 1e6
    for lab, ls in zip(labs, ("-", "--")):
        for i, (dx, col) in enumerate((("dx 3 (coarse)", RED), ("dx 1.5 (fine)", BLUE))):
            m = [np.exp(np.mean(res[lab]["lnw"][b][i])) * to_P if res[lab]["lnw"][b][i] else np.nan for b in AMPS]
            ax.loglog(AMPS, m, "o" + ls, color=col, lw=2, ms=6, label="%s, %s" % (lab.split(" (")[0], dx))
    ax.set_xlabel("peak flux density B (mT)")
    ax.set_ylabel("loss density P at 30 MHz (W/cm³), geometric mean")
    ax.set_xticks(AMPS)
    ax.set_xticklabels([str(b) for b in AMPS])
    ax.legend(fontsize=10)
    ax.set_title("Loss of both meshes")
    ax = axs[1]
    for j, (lab, col, off) in enumerate(zip(labs, (RED, GREY), (-0.12, 0.12))):
        for i, b in enumerate(AMPS):
            v = np.array(res[lab]["dlnw"][b])
            x = i + off
            ax.scatter(np.full(len(v), x), v, color=col, s=28, alpha=0.7, zorder=3,
                       label=lab if i == 0 else None)
            if len(v) > 1:
                se = v.std(ddof=1) / math.sqrt(len(v))
                ax.errorbar(x + 0.06, v.mean(), yerr=2 * se, fmt="s", color=col, ms=7, capsize=4, lw=1.5)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(len(AMPS)))
    ax.set_xticklabels(["%d mT" % b for b in AMPS])
    ax.set_ylabel("ln w(dx 3) − ln w(dx 1.5)")
    ax.set_title("Mesh error: dots = seeds, squares = mean ± 2 SE")
    ax.text(0.98, 0.04, "below 0: the coarse mesh gives too little loss", transform=ax.transAxes, ha="right", fontsize=10,
            color=GREY)
    ax.legend(fontsize=10, loc="upper right")
    fig.tight_layout()
    fig.savefig(HERE / "results" / "k3_mesh.png", dpi=200)
    fig_beta()


if __name__ == "__main__":
    main()
