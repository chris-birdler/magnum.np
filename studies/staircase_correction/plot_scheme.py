#!/usr/bin/env python
"""
Schematic of design B (PROTOKOLL 3.28): the five arms as cross-sections of one particle surface and the effects as
arrows between them. Each effect is ln P(arrow head) − ln P(arrow tail). Not to scale: the real sphere (d/l_ex 96,
dx 3 l_ex = 10 nm) has 16 cells of 10 nm on its radius; the schematic has 5.3.

    python plot_scheme.py --plot results/scheme_B.png
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Arc, FancyArrowPatch, Patch, Rectangle  # noqa: E402

R = 5.3          # sphere radius in cells of 10 nm (schematic)
NC = 6           # cells of 10 nm on each side of a panel
BLUE = np.array([0.30, 0.50, 0.80])


def cells(h):
    """Cells of size h (1: 10 nm, 0.5: 5 nm) with the centre inside the sphere (voxel model)."""
    n = int(round(NC / h))
    c = (np.arange(n) + 0.5) * h
    return {(i, j) for i in range(n) for j in range(n) if c[i] ** 2 + c[j] ** 2 < R ** 2}


def outline(F, h):
    """Surface segments of the staircase F (cells of size h); the cuts at x 0 and y 0 are not surface."""
    seg = []
    for i, j in F:
        if (i + 1, j) not in F:
            seg.append([((i + 1) * h, j * h), ((i + 1) * h, (j + 1) * h)])
        if (i, j + 1) not in F:
            seg.append([(i * h, (j + 1) * h), ((i + 1) * h, (j + 1) * h)])
    return seg


def panel(ax, x0, y0, name, text, m_h, fill, field_h):
    """m_h: cell size of m; fill: 'voxel <h>' or 'fraction' (cells of 10 nm weighted by v); field_h: staircase of the field."""
    for k in np.arange(0, NC + 1e-9, m_h):
        ax.plot([x0 + k, x0 + k], [y0, y0 + NC], color="0.75", lw=0.4, zorder=1)
        ax.plot([x0, x0 + NC], [y0 + k, y0 + k], color="0.75", lw=0.4, zorder=1)
    if fill == "fraction":
        fine = cells(0.5)
        for i, j in cells(1.0) | {(a // 2, b // 2) for a, b in fine}:
            v = sum((2 * i + a, 2 * j + b) in fine for a in (0, 1) for b in (0, 1)) / 4.0
            ax.add_patch(Rectangle((x0 + i, y0 + j), 1, 1, color=1 - v * (1 - BLUE), lw=0, zorder=0))
    else:
        h = float(fill.split()[1])
        for i, j in cells(h):
            ax.add_patch(Rectangle((x0 + i * h, y0 + j * h), h, h, color=BLUE, lw=0, zorder=0))
    seg = [[(x0 + a, y0 + b) for a, b in s] for s in outline(cells(field_h), field_h)]
    ax.add_collection(LineCollection(seg, colors="k", lw=2.2, zorder=3))
    ax.add_patch(Arc((x0, y0), 2 * R, 2 * R, theta1=0, theta2=90, color="tab:red", lw=1.4, ls="--", zorder=4))
    ax.add_patch(Rectangle((x0, y0), NC, NC, fill=False, color="0.4", lw=0.8, zorder=5))
    ax.text(x0 + NC - 0.2, y0 + NC - 0.25, name, ha="right", va="top", fontsize=12, weight="bold", zorder=6)
    ax.text(x0 + NC / 2, y0 - 0.3, text, ha="center", va="top", fontsize=8, linespacing=1.4)


def arrow(ax, p, q, sym, label, lx, ly, ha="center", rad=0.0):
    """Arrow from p to q; the effect symbol at (lx, ly), its meaning below it."""
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=14, lw=1.6, color="0.15",
                                 connectionstyle="arc3,rad=%g" % rad, zorder=7))
    ax.text(lx, ly, sym, ha=ha, va="bottom", fontsize=13, zorder=8)
    ax.text(lx, ly - 0.1, label, ha=ha, va="top", fontsize=8.5, linespacing=1.4, zorder=8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", required=True)
    a = ap.parse_args()
    fig, ax = plt.subplots(figsize=(9.6, 8.4))
    X, Y = (0, 11, 22), (0, 12)                                  # columns and rows of the panels
    M = r"$\mathbf{m}$"
    panel(ax, X[0], Y[1], "T3", M + " on 10-nm cells\nfield: Newell, 10-nm staircase", 1.0, "voxel 1", 1.0)
    panel(ax, X[2], Y[1], "T15g", M + " on 5-nm cells, geometry of T3\nfield: Newell, 10-nm staircase",
          0.5, "voxel 1", 1.0)
    panel(ax, X[0], Y[0], "T3Pc", M + " on 10-nm cells, weight $v$\nfield: Newell + correction,\n5-nm staircase",
          1.0, "fraction", 0.5)
    panel(ax, X[1], Y[0], "X3", M + " on 10-nm cells, weight $v$\nfield: exact (fine FFT),\n5-nm staircase",
          1.0, "fraction", 0.5)
    panel(ax, X[2], Y[0], "D15", M + " on 5-nm cells, own geometry\nfield: Newell, 5-nm staircase", 0.5, "voxel 0.5", 0.5)

    arrow(ax, (X[0] + 6.3, Y[1] + 3), (X[2] - 0.3, Y[1] + 3),
          r"$N_1$", "only the cells of " + M + ": 10 nm → 5 nm\n(same 10-nm geometry)", 14, Y[1] + 5.4)
    arrow(ax, (X[0] + 3, Y[1] - 1.9), (X[0] + 3, Y[0] + 6.3),
          r"$\Delta_\mathrm{c}$", "field of the 5-nm staircase\n(correction), " + M + " stays coarse",
          X[0] + 3.5, Y[0] + 8.6, ha="left")
    arrow(ax, (X[2] + 3, Y[1] - 1.9), (X[2] + 3, Y[0] + 6.3),
          r"$\Delta_\mathrm{g}$", "5-nm geometry: field\nand surface of " + M, X[2] + 2.5, Y[0] + 8.6, ha="right")
    arrow(ax, (X[1] - 0.3, Y[0] + 3), (X[0] + 6.3, Y[0] + 3),
          r"$N_3$", "correction\nagainst exact field", (X[0] + 6 + X[1]) / 2, Y[0] + 5.2)
    arrow(ax, (X[1] + 6.3, Y[0] + 3), (X[2] - 0.3, Y[0] + 3),
          r"$m$", "cells of " + M + ": 10 → 5 nm,\nsame field geometry", (X[1] + 6 + X[2]) / 2, Y[0] + 5.2)
    arrow(ax, (X[2] + 3, Y[0] - 2.4), (X[0] + 3, Y[0] - 2.4),
          "total", "corrected 10-nm model against the 5-nm model", 14, Y[0] - 6.3, rad=-0.18)

    ax.text(14, Y[1] + 8.4, r"each effect $=\ \ln P(\mathrm{arrow\ head}) - \ln P(\mathrm{arrow\ tail})$",
            ha="center", va="center", fontsize=10)
    ax.text(14, Y[1] + 7.1, r"$d = \Delta_\mathrm{c} - \Delta_\mathrm{g} = N_1 + N_3 - m$", ha="center", va="center",
            fontsize=11, bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.5"))
    ax.legend(handles=[Patch(color=BLUE, label="magnetic cell"),
                       Patch(color=1 - 0.5 * (1 - BLUE), label="cell partly filled (weight $v$)"),
                       Line2D([], [], color="0.75", lw=0.8, label="cells of " + M),
                       Line2D([], [], color="k", lw=2.2, label="surface seen by the stray field"),
                       Line2D([], [], color="tab:red", lw=1.4, ls="--", label="smooth sphere (goal)")],
              loc="center", bbox_to_anchor=(14, 9.0), bbox_transform=ax.transData, fontsize=8, frameon=False)
    ax.text(14, Y[0] - 7.6, "Cross-section of one particle surface, not to scale (real: 16 cells of 10 nm on the "
            "radius). All arms start from the same prepared state.", ha="center", va="center", fontsize=7.5,
            color="0.3")
    ax.set_xlim(-0.8, 28.8)
    ax.set_ylim(Y[0] - 8.3, Y[1] + 9.2)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(a.plot, dpi=170)


if __name__ == "__main__":
    main()
