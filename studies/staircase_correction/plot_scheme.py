#!/usr/bin/env python
"""
Schematic of design B (PROTOKOLL 3.28): the five arms as cross-sections of one particle surface and the effects as
arrows between them. Each effect is ln P(arrow head) − ln P(arrow tail). Not to scale: the real sphere (d/l_ex 96,
dx 3 l_ex = 10 nm) has 16 cells of 10 nm on its radius; the schematic has 5.3.

Names (Chris, 2026-10-21): coarse 10 nm, fine 5 nm; "roughness" is the step size of the voxelised surface (before:
"staircase"). Title: a short name, Coarse (run code T3), Fine m (T15g), Fine (D15), Coarse m corrected (T3Pc: Newell on the
10-nm cells + correction), Coarse m exact (X3: Newell on all 5-nm cells). All arms use the Newell tensor with FFT. The line below the title always tells what m and the roughness are.

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
ORANGE = "tab:orange"
W = 1            # half width of the correction window in the schematic (runs: p = 3 cells)
M = r"$\mathbf{m}$"


def cells(h):
    """Cells of size h (1: 10 nm, 0.5: 5 nm) with the centre inside the sphere (voxel model)."""
    n = int(round(NC / h))
    c = (np.arange(n) + 0.5) * h
    return {(i, j) for i in range(n) for j in range(n) if c[i] ** 2 + c[j] ** 2 < R ** 2}


def outline(F, h):
    """Surface segments of the voxel surface F (cells of size h); the cuts at x 0 and y 0 are not surface."""
    seg = []
    for i, j in F:
        if (i + 1, j) not in F:
            seg.append([((i + 1) * h, j * h), ((i + 1) * h, (j + 1) * h)])
        if (i, j + 1) not in F:
            seg.append([(i * h, (j + 1) * h), ((i + 1) * h, (j + 1) * h)])
    return seg


def target_cell():
    """The cut cell of 10 nm nearest to the diagonal: centre of the correction window in the schematic."""
    fine, full = cells(0.5), cells(1.0)
    cut = {(a // 2, b // 2) for a, b in fine} - {c for c in full if all((2 * c[0] + a, 2 * c[1] + b) in fine
                                                                       for a in (0, 1) for b in (0, 1))}
    return min(cut, key=lambda c: abs(c[0] - c[1]) + abs(np.hypot(c[0] + 0.5, c[1] + 0.5) - R))


def panel(ax, x0, y0, title, sub, text, code, m_h, fill, field_h, window=False):
    """m_h: cell size of m; fill: 'voxel <h>' or 'fraction' (cells of 10 nm weighted by v); field_h: roughness of the
    surface that the stray field sees. window: the correction (pairwise form) for one target cell: inside the window
    the fine surface, outside each cell of 10 nm as a moment at the centroid of its filled part."""
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
    seg = outline(cells(field_h), field_h)
    if window:
        ti, tj = target_cell()
        lo, hi = (ti - W, tj - W), (ti + W + 1, tj + W + 1)
        inside = lambda x, y: lo[0] <= x <= hi[0] and lo[1] <= y <= hi[1]  # noqa: E731
        seg = [s for s in seg if inside((s[0][0] + s[1][0]) / 2, (s[0][1] + s[1][1]) / 2)]
        fine = cells(0.5)
        for i, j in cells(1.0) | {(a // 2, b // 2) for a, b in fine}:
            if inside(i + 0.5, j + 0.5):
                continue
            sub_c = [((2 * i + a + 0.5) / 2, (2 * j + b + 0.5) / 2) for a in (0, 1) for b in (0, 1)
                     if (2 * i + a, 2 * j + b) in fine]
            if sub_c:
                cx, cy = np.mean(sub_c, axis=0)
                ax.plot([x0 + cx], [y0 + cy], "o", ms=2.6, color="k", zorder=3)
        x1, y1 = max(lo[0], 0), max(lo[1], 0)
        ax.add_patch(Rectangle((x0 + x1, y0 + y1), min(hi[0], NC) - x1, min(hi[1], NC) - y1, fill=False,
                               color=ORANGE, lw=1.6, ls="--", zorder=6))
        ax.add_patch(Rectangle((x0 + ti, y0 + tj), 1, 1, fill=False, color=ORANGE, lw=2.2, zorder=6))
    seg = [[(x0 + a, y0 + b) for a, b in s] for s in seg]
    ax.add_collection(LineCollection(seg, colors="k", lw=2.2, zorder=3))
    ax.add_patch(Arc((x0, y0), 2 * R, 2 * R, theta1=0, theta2=90, color="tab:red", lw=1.4, ls="--", zorder=4))
    ax.add_patch(Rectangle((x0, y0), NC, NC, fill=False, color="0.4", lw=0.8, zorder=5))
    ax.text(x0 + NC / 2, y0 + NC + 0.95, title, ha="center", va="bottom", fontsize=11.5, weight="bold")
    ax.text(x0 + NC / 2, y0 + NC + 0.2, sub, ha="center", va="bottom", fontsize=9)
    ax.text(x0 + NC / 2, y0 - 0.3, text, ha="center", va="top", fontsize=8, linespacing=1.4)
    ax.text(x0 + NC / 2, y0 - 0.5 - 0.72 * (text.count("\n") + 1), "run code " + code, ha="center", va="top",
            fontsize=7.5, color="0.45")


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
    fig, ax = plt.subplots(figsize=(10.4, 9.9))
    X, Y = (0, 12.5, 25), (0, 15)                                # columns and rows of the panels
    xc = X[1] + NC / 2
    panel(ax, X[0], Y[1], "Coarse", "coarse " + M + ", coarse roughness",
          M + " on 10-nm cells\nstray field: Newell on 10-nm cells", "T3", 1.0, "voxel 1", 1.0)
    panel(ax, X[2], Y[1], "Fine " + M, "fine " + M + ", coarse roughness",
          M + " on 5-nm cells, roughness of Coarse\nstray field: Newell on 5-nm cells", "T15g", 0.5, "voxel 1", 1.0)
    panel(ax, X[0], Y[0], "Coarse " + M + " corrected", "coarse " + M + ", fine roughness",
          M + " on 10-nm cells, weight $v$\nstray field: Newell on 10-nm cells\n+ correction in the window", "T3Pc",
          1.0, "fraction", 0.5, window=True)
    panel(ax, X[1], Y[0], "Coarse " + M + " exact", "coarse " + M + ", fine roughness",
          M + " on 10-nm cells, weight $v$\nstray field: Newell on all\n5-nm cells", "X3",
          1.0, "fraction", 0.5)
    panel(ax, X[2], Y[0], "Fine", "fine " + M + ", fine roughness",
          M + " on 5-nm cells\nstray field: Newell on 5-nm cells", "D15", 0.5, "voxel 0.5", 0.5)

    arrow(ax, (X[0] + 6.3, Y[1] + 3), (X[2] - 0.3, Y[1] + 3),
          r"$N_1$", "only " + M + " becomes fine\n(same coarse roughness)", xc, Y[1] + 5.0)
    arrow(ax, (X[0] + 3, Y[1] - 2.9), (X[0] + 3, Y[0] + 8.4),
          r"$\Delta_\mathrm{c}$", "the stray field gets the fine\nroughness (correction),\n" + M + " stays coarse",
          X[0] + 3.5, Y[0] + 11.2, ha="left")
    arrow(ax, (X[2] + 3, Y[1] - 2.9), (X[2] + 3, Y[0] + 8.4),
          r"$\Delta_\mathrm{g}$", "the roughness becomes fine:\nstray field and\nsurface of " + M,
          X[2] + 2.5, Y[0] + 11.2, ha="right")
    arrow(ax, (X[1] - 0.3, Y[0] + 3), (X[0] + 6.3, Y[0] + 3),
          r"$N_3$", "correction\nagainst exact field", (X[0] + 6 + X[1]) / 2, Y[0] + 5.2)
    arrow(ax, (X[1] + 6.3, Y[0] + 3), (X[2] - 0.3, Y[0] + 3),
          r"$m$", "only " + M + " becomes fine\n(same fine roughness)", (X[1] + 6 + X[2]) / 2, Y[0] + 5.2)
    arrow(ax, (X[2] + 3, Y[0] - 3.3), (X[0] + 3, Y[0] - 3.3),
          "total", "Coarse " + M + " corrected against Fine", xc, Y[0] - 7.3, rad=-0.16)

    ax.text(xc, Y[1] + 11.0, r"each effect $=\ \ln P(\mathrm{arrow\ head}) - \ln P(\mathrm{arrow\ tail})$",
            ha="center", va="center", fontsize=10)
    ax.text(xc, Y[1] + 9.7, r"$d = \Delta_\mathrm{c} - \Delta_\mathrm{g} = N_1 + N_3 - m$", ha="center", va="center",
            fontsize=11, bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.5"))
    ax.legend(handles=[Patch(color=BLUE, label="magnetic cell"),
                       Patch(color=1 - 0.5 * (1 - BLUE), label="cell partly filled (weight $v$)"),
                       Line2D([], [], color="0.75", lw=0.8, label="cells of " + M),
                       Line2D([], [], color="k", lw=2.2, label="surface seen by the stray field"),
                       Line2D([], [], color="tab:red", lw=1.4, ls="--", label="smooth sphere (goal)"),
                       Line2D([], [], color=ORANGE, lw=1.6, ls="--",
                              label="correction window around the marked cell\n(schematic ±1 cell, runs ±3 cells)"),
                       Line2D([], [], marker="o", ls="", ms=3, color="k",
                              label="outside the window: moment at the\ncentroid of the filled part")],
              loc="center", bbox_to_anchor=(xc, Y[0] - 10.2), bbox_transform=ax.transData, ncol=3, fontsize=8,
              frameon=False)
    ax.text(xc, Y[0] - 12.6, "Cross-section of one particle surface, not to scale (real: 16 cells of 10 nm on the "
            "radius). Coarse: 10 nm, fine: 5 nm.\nRoughness: step size of the surface. All arms start from the same "
            "prepared state.", ha="center", va="center", fontsize=7.5, color="0.3", linespacing=1.4)
    ax.set_xlim(-1.5, X[2] + NC + 1.5)
    ax.set_ylim(Y[0] - 13.5, Y[1] + 11.7)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(a.plot, dpi=170)


if __name__ == "__main__":
    main()
