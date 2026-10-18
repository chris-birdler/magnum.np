#!/usr/bin/env python
"""
Plots of the corner test (PLAN 4.13, Chris 2026-10-18): loss density P over B_peak per corner, the single runs and the
mean over the 3 seeds with a shaded error band, old (DemagFieldPBC, atol 1e-5) against new (PLAN 4.13 setup).
Mean: geometric mean of P over the seeds; band: mean × exp(± SE of ln P) (n 3). Second figure: factor P_new / P_old per
amplitude, mean over the 12 runs (geometric) with its 90 % CI (t, n − 1) as a band, and the corner means as lines.
Values: the steady-tail rule (≥ 2 calm cycles), the nominal amplitudes of the main fit (analyze_design.run_record).

    python plot_corners_nx.py --runs <study runs folder> --out <folder for the PNG files>
"""
import argparse
import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib import font_manager                 # noqa: E402
import numpy as np                                  # noqa: E402
from scipy import stats                             # noqa: E402

import analyze                                      # noqa: E402
import analyze_design as AD                         # noqa: E402
import eval_corners_nx as EN                        # noqa: E402
import units as U                                   # noqa: E402

RED, BLUE = "#9B111E", "#1F4E79"
CORNERS = {"K1": "d/l_ex 212, φ 0.62, L_eff 30, r_p 1.5", "K2": "d/l_ex 300, φ 0.62, L_eff 12, r_p 3",
           "K3": "d/l_ex 212, φ 0.69, L_eff 18, r_p 1.5", "K4": "d/l_ex 300, φ 0.55, L_eff 18, r_p 1.5"}


def style():
    for f in pathlib.Path("/usr/share/fonts/truetype/crosextra").glob("Carlito*.ttf"):
        font_manager.fontManager.addfont(str(f))
    plt.rcParams.update({"font.family": "Carlito", "font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 200})


def lnP(runs_dir):
    """{run: (ln P old, ln P new)} for the 12 corner runs, P in W/cm³ at the 7 nominal amplitudes."""
    u = U.units()
    lnF = math.log(u["Kd"] * U.F_REL_DEFAULT * u["f_M"] / 1e6)
    out = {}
    for c in CORNERS:
        for s in range(3):
            n = "%s_%d" % (c, s)
            out[n] = (EN.lnw(EN.old_dir(runs_dir, n)) + lnF, EN.lnw(runs_dir / (n + "_nx")) + lnF)
    return out


def band(L):
    """(geometric mean, lower, upper) of P from the ln P rows L: mean × exp(± SE)."""
    m = np.nanmean(L, 0)
    se = np.nanstd(L, 0, ddof=1) / np.sqrt(np.sum(np.isfinite(L), 0))
    return np.exp(m), np.exp(m - se), np.exp(m + se)


def fig_corners(D, out):
    B = np.asarray(AD.AMPS_MT)
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, sharey=True)
    for ax, c in zip(axs.flat, CORNERS):
        for k, (col, lab) in enumerate(((RED, "old (DemagFieldPBC, atol 1e-5)"), (BLUE, "new (setup of PLAN 4.13)"))):
            L = np.array([D["%s_%d" % (c, s)][k] for s in range(3)])
            for row in L:
                ax.plot(B, np.exp(row), color=col, lw=0.8, alpha=0.45, marker="o", ms=3)
            m, lo, hi = band(L)
            ax.fill_between(B, lo, hi, color=col, alpha=0.2, lw=0)
            ax.plot(B, m, color=col, lw=2.4, label=lab + ": mean ± SE (n 3)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title("%s   (%s)" % (c, CORNERS[c]), fontsize=12)
        ax.set_xticks(B)
        ax.set_xticklabels(["%g" % b for b in B])
        ax.minorticks_off()
    for ax in axs[1]:
        ax.set_xlabel("$B_\\mathrm{peak}$ (mT)")
    for ax in axs[:, 0]:
        ax.set_ylabel("$P$ (W/cm³)")
    axs[0, 0].legend(loc="upper left", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "corners_nx_loss.png")
    plt.close(fig)


def fig_factor(D, out):
    B = np.asarray(AD.AMPS_MT)
    d = np.array([D[n][1] - D[n][0] for n in D])
    m = d.mean(0)
    h = stats.t.ppf(0.95, len(d) - 1) * d.std(0, ddof=1) / math.sqrt(len(d))
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.fill_between(B, np.exp(m - h), np.exp(m + h), color="0.7", alpha=0.5, lw=0, label="all 12 runs: 90 % CI")
    ax.plot(B, np.exp(m), color="k", lw=2.4, label="all 12 runs: geometric mean")
    for c, col in zip(CORNERS, ("#9B111E", "#1F4E79", "#D98E04", "#4D8C57")):
        dc = np.array([D["%s_%d" % (c, s)][1] - D["%s_%d" % (c, s)][0] for s in range(3)])
        ax.plot(B, np.exp(dc.mean(0)), color=col, lw=1.2, marker="o", ms=3, label="%s (n 3)" % c)
    ax.axhline(1.0, color="0.4", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(B)
    ax.set_xticklabels(["%g" % b for b in B])
    ax.set_yticks([0.4, 0.5, 0.7, 1.0, 1.4, 2.0])
    ax.set_yticklabels(["0.4", "0.5", "0.7", "1", "1.4", "2"])
    ax.minorticks_off()
    ax.set_xlabel("$B_\\mathrm{peak}$ (mT)")
    ax.set_ylabel("$P_\\mathrm{new} / P_\\mathrm{old}$")
    ax.legend(fontsize=9, ncol=2, loc="lower left")
    fig.tight_layout()
    fig.savefig(out / "corners_nx_factor.png")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(pathlib.Path(__file__).resolve().parent / "runs"))
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    style()
    D = lnP(pathlib.Path(a.runs))
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    fig_corners(D, out)
    fig_factor(D, out)


if __name__ == "__main__":
    main()
