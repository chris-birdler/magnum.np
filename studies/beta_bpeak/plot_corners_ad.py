#!/usr/bin/env python
"""
Plots of the corner evaluation with three arms (Chris 2026-10-19; eval_corners_ad.py, PROTOKOLL 7.119): old main study
(DemagFieldPBC), nx (setup of PLAN 4.13, relaxed random start with Bloch points), ad (setup of PLAN 4.13, option C,
no Bloch point in the start state). Input: the JSON of eval_corners_ad.py --json.
  corners_ad_loss.png    P over B_peak per corner: geometric mean over the 3 microstructures (line and markers), band
                         from the smallest to the largest of the 3 runs (scatter); open marker = at least 1 of the 3 runs has a Bloch point in this stage (ad: at
                         a sample; nx: at the start or the end of the stage; old: no stage information); dashed lines =
                         power law fitted to the mean ln P over 9 … 50 mT (grey band).
  corners_ad_beta.png    β per corner and arm (mean ± SE, n 3) in the windows 9 … 150, 9 … 50 and 70 … 150 mT.
  corners_ad_factor.png  factor on P per amplitude, exp(mean paired Δ ln P) with the band exp(± SE) (n 12).
Colors: categorical palette of the dataviz skill, validated (slots 1 … 3 for the arms, 4 … 6 for the comparisons);
each arm also has its own marker.

    python plot_corners_ad.py --data results/corners_ad_data.json --out results/figures_corners_ad
"""
import argparse
import json
import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                  # noqa: E402

from matplotlib.ticker import NullFormatter        # noqa: E402

from plot_corners_nx import CORNERS, style          # noqa: E402

ARMS = ("old", "nx", "ad")
COLOR = {"old": "#2a78d6", "nx": "#eb6834", "ad": "#1baf7a"}
MARK = {"old": "o", "nx": "s", "ad": "^"}
LABEL = {"old": "old main study (DemagFieldPBC)", "nx": "new solver, start with Bloch points",
         "ad": "new solver, option C (no Bloch point at start)"}
CMP = {"ad/old": ("#eda100", "D"), "ad/nx": ("#e87ba4", "v"), "nx/old": ("#008300", "P")}
WIN = {"full": "9 … 150 mT", "low": "9 … 50 mT", "high": "70 … 150 mT"}
INK = "#3d3d3a"
TITLE = {"K1": r"K1: $d/l_\mathrm{ex}$ 212, $\varphi$ 0.62, $L_\mathrm{eff}/l_\mathrm{ex}$ 30, $r_\mathrm{p}$ 1.5",
         "K2": r"K2: $d/l_\mathrm{ex}$ 300, $\varphi$ 0.62, $L_\mathrm{eff}/l_\mathrm{ex}$ 12, $r_\mathrm{p}$ 3",
         "K3": r"K3: $d/l_\mathrm{ex}$ 212, $\varphi$ 0.69, $L_\mathrm{eff}/l_\mathrm{ex}$ 18, $r_\mathrm{p}$ 1.5",
         "K4": r"K4: $d/l_\mathrm{ex}$ 300, $\varphi$ 0.55, $L_\mathrm{eff}/l_\mathrm{ex}$ 18, $r_\mathrm{p}$ 1.5"}


def bp_flag(data, run, arm, k):
    bp = data["bp"][run]
    return arm in ("ad", "nx") and k in (bp.get(arm) or [])


def loss_figure(data, out):
    B = np.array(data["amps_mT"])
    lo = B <= 50.0
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, sharey=True)
    for ax, c in zip(axs.flat, CORNERS):
        runs = [r for r in data["runs"] if r.startswith(c + "_")]
        ax.axvspan(8, 55, color="#e6e6e3", zorder=0, lw=0)
        for i, a in enumerate(ARMS):
            L = np.array([data["lnP"][r][a] for r in runs], float)
            m = np.nanmean(L, 0)
            x = B
            ax.fill_between(B, np.exp(np.nanmin(L, 0)), np.exp(np.nanmax(L, 0)), color=COLOR[a], alpha=0.18, lw=0,
                            zorder=1)
            ax.plot(B, np.exp(m), "-", color=COLOR[a], lw=1.2, zorder=2)
            for k in range(len(B)):
                open_ = any(bp_flag(data, r, a, k) for r in runs)
                ax.plot(x[k], math.exp(m[k]), MARK[a], ms=8, mec=COLOR[a], mew=1.8,
                        mfc="white" if open_ else COLOR[a], zorder=3, label=LABEL[a] if k == 0 else None)
            p = np.polyfit(np.log(B[lo]), m[lo], 1)
            xx = np.array([9.0, 50.0])
            ax.plot(xx, np.exp(np.polyval(p, np.log(xx))), "--", color=COLOR[a], lw=2, zorder=2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(TITLE[c], fontsize=12, color=INK, loc="left")
        ax.text(9.2, 0.97, "fit window", transform=ax.get_xaxis_transform(), fontsize=10, color="#6b6b67", va="top")
        ax.grid(True, which="major", color="#d9d9d6", lw=0.6)
        ax.set_xticks(B)
        ax.set_xticklabels(["%g" % b for b in B])
    for ax in axs[1]:
        ax.set_xlabel(r"$B_\mathrm{peak}$ (mT)")
    for ax in axs[:, 0]:
        ax.set_ylabel(r"$P$ (W/cm$^3$)")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=11, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out / "corners_ad_loss.png")
    plt.close(fig)


def beta_figure(data, out):
    cs = list(CORNERS)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.6), sharey=True)
    for ax, w in zip(axs, ("full", "low", "high")):
        for i, a in enumerate(ARMS):
            v = [data["corner_beta"][w][c][a] for c in cs]
            x = np.arange(len(cs)) + 0.2 * (i - 1)
            m = np.array([t[0] if t else np.nan for t in v])
            e = np.array([t[1] if t else np.nan for t in v])
            ax.errorbar(x, m, yerr=e, fmt=MARK[a], ms=8, color=COLOR[a], mec=COLOR[a], elinewidth=1.4, capsize=0,
                        label=LABEL[a])
        ax.set_xticks(np.arange(len(cs)))
        ax.set_xticklabels(cs)
        ax.set_title("window %s" % WIN[w], fontsize=12, color=INK, loc="left")
        ax.grid(True, axis="y", color="#d9d9d6", lw=0.6)
        ax.axhline(2.0, color="#8a8a86", lw=0.8, ls=":")
    axs[0].set_ylabel(r"$\beta$")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=11, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(out / "corners_ad_beta.png")
    plt.close(fig)


def factor_figure(data, out):
    B = np.array(data["amps_mT"])
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, (k, (col, mk)) in enumerate(CMP.items()):
        s = data["factor"][k]
        m = np.array([t[0] if t else np.nan for t in s])
        e = np.array([t[1] if t else np.nan for t in s])
        ax.fill_between(B, np.exp(m - e), np.exp(m + e), color=col, alpha=0.18, lw=0)
        ax.plot(B, np.exp(m), mk + "-", ms=8, color=col, lw=1.5, label="%s" % {
                        "ad/old": "option C / old main study", "ad/nx": "option C / start with Bloch points",
                        "nx/old": "start with Bloch points / old main study"}[k])
    ax.axhline(1.0, color=INK, lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(B)
    ax.set_xticklabels(["%g" % b for b in B])
    ax.set_yticks([0.5, 0.7, 1.0, 1.4, 2.0])
    ax.set_yticklabels(["0.5", "0.7", "1", "1.4", "2"])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$B_\mathrm{peak}$ (mT)")
    ax.set_ylabel(r"factor on $P$")
    ax.grid(True, which="major", color="#d9d9d6", lw=0.6)
    ax.legend(fontsize=10, loc="lower left")
    fig.tight_layout()
    fig.savefig(out / "corners_ad_factor.png")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    data = json.loads(pathlib.Path(a.data).read_text())
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    style()
    loss_figure(data, out)
    beta_figure(data, out)
    factor_figure(data, out)


if __name__ == "__main__":
    main()
