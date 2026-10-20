#!/usr/bin/env python
"""
Plots of the preparation tests at 9 mT (Chris 2026-10-20; PROTOKOLL 7.125 … 7.129). P = loss density in W/cm³
(P = f K_d w, as eval_corners_ad.py); w_l = mean w_loop of the last 5 cycles of the 9 mT stage; the error bar is the SD
of w_loop over these 5 cycles (cycle-to-cycle scatter of one state).
  satstatic_states.png   P of each prepared state of K4_1 and K2_0, grouped by preparation; filled marker = no BP in
                         the prepared state and in the 9 mT stage, open marker = BPs; grey band = geometric mean × exp(± SD
                         of ln P) over the states without BPs (the spread between states).
  history_pairs.png      K3 and K4 history test (7.126, 7.127): P of the nx state and of the option C state per
                         microstructure, 16 cycles each.
  satstatic_cycles.png   P per cycle of the 9 mT stage for the 9 runs of the static saturation test.
Colors: categorical palette of plot_corners_ad.py (validated with the dataviz skill), fixed order; each group also has
its own marker.

    python plot_satstatic.py --runs <runs folder> --out <figure folder>
"""
import argparse
import json
import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                  # noqa: E402

import units as U                                   # noqa: E402
from plot_corners_nx import style                   # noqa: E402

INK = "#3d3d3a"
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
# BPs of the relaxed random start of the nx runs (no Bloch data in these old summaries): PROTOKOLL 7.118, table
NX_START_BP = {"K3_0": 0, "K3_1": 2, "K3_2": 1, "K4_0": 2, "K4_1": 2, "K4_2": 2, "K2_0": 9}
# (run suffix, label, group); groups: 0 relaxed random start, 1 AC from a relaxed random start, 2 protocol S,
# 3 static saturation test
PREP = [("_nx", "random start (cell blocks)", 0), ("_nx_ip", "random start (per particle)", 0),
        ("_nx_pd150", "AC from 150 mT", 1), ("_nx_pd300", "AC from 300 mT (option C)", 1),
        ("_nx_pdsat", r"AC from $1-\varphi$", 1), ("_nx_ps", "protocol S", 2),
        ("_nx_sad", r"A, $1-\varphi$", 3), ("_nx_sad_rep", r"A, $1-\varphi$, rerun", 3),
        ("_nx_sbd", r"B, $1-\varphi$", 3), ("_nx_sa200", r"A, $2M_\mathrm{s}$", 3), ("_nx_sb200", r"B, $2M_\mathrm{s}$", 3)]
GROUP = ["relaxed random start (no AC)", "AC demagnetization from the relaxed random start", "protocol S",
         "static saturation test (A: uniform start, B: state after AC from 300 mT)"]
GMARK = ["s", "^", "D", "o"]


def stage9(run_dir):
    """w_loop per cycle of the 9 mT stage and the BP flag of the run (None if the run is missing)."""
    try:
        s = json.loads((run_dir / "summary.json").read_text())
    except (OSError, ValueError):
        return None
    st = [x for x in s.get("stages", []) if x.get("name") == "amp00" and x.get("cycles")]
    if not st:
        return None
    cyc = st[0]["cycles"]
    i = s.get("init", {})
    pre = i.get("predemag") or {}
    bp = len(pre.get("bloch_points", [])) if pre else 0
    bp += max([c.get("bloch", {}).get("max") or 0 for c in cyc], default=0)
    if not pre and i.get("bloch_virgin"):
        bp += sum(i["bloch_virgin"])
    return np.array([c["w_loop"] for c in cyc]), bp


def logticks(ax):
    """Plain numbers on a log axis (1, 2, 3, 4, 6, 10 …)."""
    fmt = matplotlib.ticker.FuncFormatter(lambda v, _: ("%g" % v) if ("%g" % v)[0] in "12346" else "")
    ax.yaxis.set_major_formatter(fmt)
    ax.yaxis.set_minor_formatter(fmt)


def p_of(w, F):
    t = w[-5:]
    return float(t.mean() * F), float(t.std(ddof=1) * F)


def plot_states(runs, F, out):
    fig, axes = plt.subplots(1, 2, figsize=(13, 7.2), sharey=True)
    for ax, ms in zip(axes, ("K4_1", "K2_0")):
        x, free = 0.0, []
        ticks, labels, last_g = [], [], None
        for suf, lab, g in PREP:
            d = stage9(runs / (ms + suf))
            if d is None:
                continue
            w, bp = d
            if suf == "_nx":
                bp += NX_START_BP[ms]
            if last_g is not None and g != last_g:
                x += 0.6
            last_g = g
            P, sd = p_of(w, F)
            ax.errorbar([x], [P], yerr=[sd], fmt=GMARK[g], ms=9, color=PAL[g], mfc=PAL[g] if bp == 0 else "white",
                        mew=2, capsize=3, zorder=3)
            if bp == 0:
                free.append(P)
            ticks.append(x)
            labels.append(lab)
            x += 1.0
        lp = np.log(free)
        gm, sdl = math.exp(lp.mean()), lp.std(ddof=1)
        ax.axhspan(gm * math.exp(-sdl), gm * math.exp(sdl), color="#c9c8c2", alpha=0.45, zorder=0, lw=0)
        ax.axhline(gm, color="#8a8984", lw=1, ls="--", zorder=1)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=10, rotation=40, ha="right", rotation_mode="anchor")
        ax.set_yscale("log")
        logticks(ax)
        ax.set_title("%s\nstates without BPs: n %d, SD of ln P %.2f, max/min %.1f" % (
            ms, len(free), sdl, max(free) / min(free)), color=INK, fontsize=12)
        ax.grid(axis="y", which="both", color="#e6e5e0", lw=0.6)
    axes[0].set_ylabel(r"$P$ at 9 mT (W/cm$^3$)")
    h = [plt.Line2D([], [], ls="", marker=GMARK[g], color=PAL[g], ms=9, mew=2, label=GROUP[g]) for g in range(4)]
    h += [plt.Line2D([], [], ls="", marker="o", color=INK, mfc="white", ms=9, mew=2, label="open: BPs in the state or "
                     "in the 9 mT stage"),
          matplotlib.patches.Patch(color="#c9c8c2", alpha=0.45, label="geometric mean × exp(± SD of ln P), states "
                                   "without BPs")]
    fig.legend(handles=h, loc="lower center", ncol=2, fontsize=10, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(out / "satstatic_states.png")
    plt.close(fig)


def plot_history(runs, F, out):
    fig, ax = plt.subplots(figsize=(9, 5))
    names = ["K3_0", "K3_1", "K3_2", "K4_0", "K4_1", "K4_2"]
    for k, ms in enumerate(names):
        vals = []
        for j, (suf, col, mk) in enumerate((("_nx_h16", PAL[0], "s"), ("_ad_h16", PAL[1], "^"))):
            d = stage9(runs / (ms + suf))
            if d is None:
                continue
            w, _ = d
            P, sd = p_of(w, F)
            bp = NX_START_BP[ms] if j == 0 else 0
            ax.errorbar([k + (j - 0.5) * 0.25], [P], yerr=[sd], fmt=mk, ms=9, color=col,
                        mfc=col if bp == 0 else "white", mew=2, capsize=3, zorder=3)
            vals.append((k + (j - 0.5) * 0.25, P))
        if len(vals) == 2:
            ax.plot(*zip(*vals), color="#8a8984", lw=1, zorder=1)
            ax.text(k, max(v[1] for v in vals) * 1.08, "× %.2f" % (vals[1][1] / vals[0][1]), ha="center", fontsize=10,
                    color=INK)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names)
    ax.set_yscale("log")
    ax.margins(y=0.08)
    logticks(ax)
    ax.set_ylabel(r"$P$ at 9 mT (W/cm$^3$)")
    ax.grid(axis="y", which="both", color="#e6e5e0", lw=0.6)
    h = [plt.Line2D([], [], ls="", marker="s", color=PAL[0], ms=9, mew=2, label="nx: relaxed random start"),
         plt.Line2D([], [], ls="", marker="^", color=PAL[1], ms=9, mew=2, label="option C: AC from 300 mT"),
         plt.Line2D([], [], ls="", marker="s", color=INK, mfc="white", ms=9, mew=2, label="open: BPs in the start state")]
    ax.legend(handles=h, loc="upper left", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "history_pairs.png")
    plt.close(fig)


def plot_cycles(runs, F, out):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    runs_s = [("_nx_sad", r"A, $1-\varphi$", 0, "o"), ("_nx_sbd", r"B, $1-\varphi$", 1, "s"),
              ("_nx_sa200", r"A, $2M_\mathrm{s}$", 2, "^"), ("_nx_sb200", r"B, $2M_\mathrm{s}$", 3, "D"),
              ("_nx_sad_rep", r"A, $1-\varphi$, rerun", 4, "v")]
    for ax, ms in zip(axes, ("K4_1", "K2_0")):
        for suf, lab, c, mk in runs_s:
            d = stage9(runs / (ms + suf))
            if d is None:
                continue
            w, _ = d
            ax.plot(np.arange(len(w)), w * F, marker=mk, color=PAL[c], lw=2, ms=7, label=lab)
        ax.set_yscale("log")
        logticks(ax)
        ax.set_xlabel("cycle of the 9 mT stage")
        ax.set_title(ms, color=INK)
        ax.grid(axis="y", which="both", color="#e6e5e0", lw=0.6)
    axes[0].set_ylabel(r"$P$ per cycle (W/cm$^3$)")
    axes[0].legend(fontsize=10, loc="lower right")
    fig.tight_layout()
    fig.savefig(out / "satstatic_cycles.png")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    runs, out = pathlib.Path(a.runs), pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    style()
    u = U.units()
    F = u["Kd"] * U.F_REL_DEFAULT * u["f_M"] / 1e6           # w -> P (W/cm³)
    plot_states(runs, F, out)
    plot_history(runs, F, out)
    plot_cycles(runs, F, out)


if __name__ == "__main__":
    main()
