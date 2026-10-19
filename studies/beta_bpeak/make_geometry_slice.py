#!/usr/bin/env python
"""
Slice through the simulation box (base point: d/l_ex 300 = 1 um, phi 0.65, dx 3 l_ex = 10 nm,
L_eff 12 l_ex = 40 nm): voxel particles, exact sphere outlines, Herzer cubes and the mesh.
Plane z = first cell layer: it cuts particle 0 (corner site, wraps around the box) and particle 3
(face centre). Right: zoom on the narrowest gap with every cell.

    python make_geometry_slice.py      # -> results/geometry_slice.png
"""
import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import ListedColormap
from matplotlib.patches import Circle, Rectangle
import numpy as np

from anisotropy_model import box_edge, cube_axes
from geometry import FCC_SITES, fcc_box

HERE = pathlib.Path(__file__).resolve().parent
LEX_NM = 3.342
D_LEX, PHI, DX, LEFF, SEED = 300.0, 0.65, 3.0, 12.0, 901
RED, BLUE, GREY = "#9B111E", "#1F4E79", "#D9D9D9"

for f in pathlib.Path("/usr/share/fonts/truetype/crosextra").glob("Carlito*.ttf"):
    font_manager.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": "Carlito", "font.size": 12})


def main():
    a_lex, d_lex = box_edge(D_LEX, PHI)
    N = int(round(a_lex / DX))
    a_lex = N * DX
    ids = fcc_box(N, d_lex * LEX_NM * 1e-9, DX * LEX_NM * 1e-9)["ids"]
    ax_c, _ = cube_axes(ids, DX, a_lex, d_lex, LEFF, SEED)
    k = 0
    sl = ids[:, :, k].T                                       # [y, x]
    axs = ax_c[:, :, k, :]
    # cube boundaries: a cube edge lies between two cells whose axes differ
    ex = np.zeros((N, N), bool)
    ey = np.zeros((N, N), bool)
    same_p = (sl >= 0)
    ex[:, :-1] = same_p[:, :-1] & (sl[:, :-1] == sl[:, 1:]) & \
        (np.abs((axs.transpose(1, 0, 2)[:, :-1] * axs.transpose(1, 0, 2)[:, 1:]).sum(-1)) < 1 - 1e-9)
    ey[:-1, :] = same_p[:-1, :] & (sl[:-1, :] == sl[1:, :]) & \
        (np.abs((axs.transpose(1, 0, 2)[:-1, :] * axs.transpose(1, 0, 2)[1:, :]).sum(-1)) < 1 - 1e-9)

    h = DX * LEX_NM                                           # cell size in nm
    A = a_lex * LEX_NM
    cmap = ListedColormap(["white", "#FCE4D6", "#DCE6F1", "#FCE4D6", "#DCE6F1"])
    z0 = (k + 0.5) * h

    fig, (p1, p2) = plt.subplots(1, 2, figsize=(13.3, 6.6), dpi=200,
                                 gridspec_kw={"width_ratios": [1, 1]})
    c, w = A / 4.0, 34 * h        # zoom: midpoint between particle 0 (corner) and 3 (face centre)
    zoom = (c - w / 2, c + w / 2, c - w / 2, c + w / 2)
    for pan, (x0, x1, y0, y1), cells in ((p1, (0, A, 0, A), False), (p2, zoom, True)):
        pan.imshow(sl + 1, origin="lower", extent=(0, A, 0, A), cmap=cmap, vmin=0, vmax=4,
                   interpolation="nearest")
        # Herzer cube boundaries
        yy, xx = np.nonzero(ex)
        for y, x in zip(yy, xx):
            pan.plot([(x + 1) * h] * 2, [y * h, (y + 1) * h], color=BLUE, lw=0.35 if not cells else 1.4)
        yy, xx = np.nonzero(ey)
        for y, x in zip(yy, xx):
            pan.plot([x * h, (x + 1) * h], [(y + 1) * h] * 2, color=BLUE, lw=0.35 if not cells else 1.4)
        # exact sphere outlines (cut at z0), all periodic images
        R = d_lex * LEX_NM / 2
        for s in FCC_SITES:
            for ix in (-1, 0, 1):
                for iy in (-1, 0, 1):
                    for iz in (-1, 0, 1):
                        cx, cy, cz = (s + np.array([ix, iy, iz])) * A
                        dz = z0 - cz
                        if abs(dz) < R:
                            pan.add_patch(Circle((cx, cy), math.sqrt(R * R - dz * dz), fill=False,
                                                 ec=RED, lw=1.2 if cells else 0.9))
        if cells:                                             # every cell edge
            for i in range(N + 1):
                pan.axvline(i * h, color="0.55", lw=0.3, zorder=0.5)
                pan.axhline(i * h, color="0.55", lw=0.3, zorder=0.5)
        pan.set_xlim(x0, x1)
        pan.set_ylim(y0, y1)
        pan.set_aspect("equal")
        pan.set_xlabel("x [nm]")
        pan.set_ylabel("y [nm]")
        for sp in ("top", "right"):
            pan.spines[sp].set_visible(False)
    p1.add_patch(Rectangle((zoom[0], zoom[2]), w, w, fill=False, ec="black", lw=1.0, ls="--"))
    p1.set_title("Box %.0f nm, cut at z = %.0f nm: particles 0 (corner) and 3 (centre)" % (A, z0), fontsize=12)
    p2.set_title("Zoom (dashed): cells %.0f nm (grey), Herzer cubes %.0f nm (blue)" % (h, LEFF * LEX_NM),
                 fontsize=12)
    fig.text(0.5, 0.01, "d = %.0f nm (%.0f l_ex), phi = %.2f, dx = %.0f nm, L_eff = %.0f nm, seed %d. "
             "Red: exact sphere outlines. Coloured: voxel particles (staircase surface). White: air."
             % (d_lex * LEX_NM, d_lex, PHI, h, LEFF * LEX_NM, SEED), ha="center", fontsize=11, color="#666666")
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    out = HERE / "results" / "geometry_slice.png"
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
