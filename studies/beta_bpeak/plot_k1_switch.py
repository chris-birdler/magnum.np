"""K1 investigation (PROTOKOLL 7.135, Chris 2026-10-21 option 1): figure of particle 0 in the last 50 mT cycle of
replay runs (for example K1_4_rep50 and K1_8_rep50, high loss; K1_5_rep50, normal loss; K3_0_rep50). Cut through the
particle centre, in the plane of the easy axis u_0 and the direction of the largest change (eval_k1_switch.py).
Columns: m.u_0 at H = 0 falling, m.u_0 at H = 0 rising, angle between the two (the colour scale ends at 45°; larger
angles have the colour of 45°). Runs on the instance, reads only.
    python3 plot_k1_switch.py <runs dir> <code dir> <out.png> <run>=<label> ...
The cut plane holds u_0 and the direction from the particle centre to the cells that turn by more than 30° (the
direction of K1_4 if fewer than 20 cells turn by more than 30°).
"""
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

runs, code, png = pathlib.Path(sys.argv[1]), sys.argv[2], sys.argv[3]
sys.path.insert(0, code)
from geometry import fcc_box, FCC_SITES          # noqa: E402
from anisotropy_model import particle_axes       # noqa: E402
import domains as DM                             # noqa: E402

P = 0
C_DIR = np.array([-0.154, -0.087, 0.984])        # direction of the largest change in K1_4, particle 0 (k1_switch.py)
RUNS = [tuple(a.split("=", 1)) for a in sys.argv[4:]]

fig, axs = plt.subplots(len(RUNS), 3, figsize=(10.5, 3.5 * len(RUNS)), constrained_layout=True)
for row, (name, label) in enumerate(RUNS):
    run = runs / name
    cfg = json.loads((run / "config.json").read_text())
    N, d, dx = cfg["N_used"], cfg["d"], cfg["dx"]
    ids = fcc_box(N, d, dx)["ids"]
    u = particle_axes(cfg.get("direction", [1.0, 0.0, 0.0]))[P]
    ma = DM.load_state(run / "m_amp03_c6_ph180.vti")["m"].astype(np.float64)
    mb = DM.load_state(run / "m_amp03_c6_ph360.vti")["m"].astype(np.float64)
    R = 0.5 * d
    x = (np.arange(N) + 0.5) * dx
    G = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1) - FCC_SITES[P] * N * dx
    G -= N * dx * np.round(G / (N * dx))
    turn = (ids == P) & (np.degrees(np.arccos(np.clip((ma * mb).sum(-1), -1.0, 1.0))) > 30.0)
    c = G[turn].mean(0) if turn.sum() >= 20 else C_DIR
    c = c / np.linalg.norm(c)
    e2 = c - (c @ u) * u
    e2 /= np.linalg.norm(e2)
    print(name, "cells > 30 deg:", int(turn.sum()), "cut direction", np.round(c, 3))
    s = np.linspace(-1.05, 1.05, 211)
    S1, S2 = np.meshgrid(s, s, indexing="ij")
    pos = FCC_SITES[P] * N * dx + R * (S1[..., None] * u + S2[..., None] * e2)
    idx = np.floor(pos / dx).astype(int) % N
    i, j, k = idx[..., 0], idx[..., 1], idx[..., 2]
    inside = ids[i, j, k] == P
    a, b = ma[i, j, k] @ u, mb[i, j, k] @ u
    ang = np.degrees(np.arccos(np.clip((ma[i, j, k] * mb[i, j, k]).sum(-1), -1.0, 1.0)))
    ext = [-1.05, 1.05, -1.05, 1.05]
    for col, (f, cmap, lim, title) in enumerate([
            (a, "RdBu_r", (-1, 1), r"$m\cdot u_0$, $H = 0$ falling"),
            (b, "RdBu_r", (-1, 1), r"$m\cdot u_0$, $H = 0$ rising"),
            (ang, "magma", (0, 45), r"angle between the two, °")]):
        ax = axs[row, col]
        im = ax.imshow(np.where(inside, f, np.nan).T, origin="lower", extent=ext, cmap=cmap, vmin=lim[0], vmax=lim[1])
        ax.set_xlabel(r"along $u_0$ [$R$]")
        ax.set_ylabel(r"across $u_0$ [$R$]")
        ax.set_title("%s: %s" % (label, title), fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.8)
fig.savefig(png, dpi=130)
print("written", png)
