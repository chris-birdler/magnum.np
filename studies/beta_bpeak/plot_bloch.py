#!/usr/bin/env python
"""
Picture of the particles and the Bloch points of a state (Chris 2026-10-18; K3_1_nx). Figure 1: 3D view of the four
particles (surface cells, one colour per particle, periodic box unwrapped around the particle with the Bloch points)
and the Bloch points of each stored state (census_bloch.py), joined by a line in the order of the stages (+: Q = +1,
−: Q = −1). Figure 2: per Bloch point of the chosen state a slice z = const through it: arrows of (m_x, m_y), colour
m_z, outline of the particle, the Bloch point marked.

    python plot_bloch.py --run runs/K3_1_nx --state init.pt --census runs/census_bloch_nx.json --out <folder>
"""
import argparse
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
import numpy as np                                    # noqa: E402

import analyze_defect as A                            # noqa: E402
from plot_corners_nx import style                     # noqa: E402

PCOL = ("#9B111E", "#1F4E79", "#D98E04", "#4D8C57")
STAGE_MT = {"init.pt": "start", "ckpt_00_amp00.pt": "9", "ckpt_01_amp01.pt": "20", "ckpt_02_amp02.pt": "35",
            "ckpt_03_amp03.pt": "50", "ckpt_04_amp04.pt": "70", "ckpt_05_amp05.pt": "100", "ckpt_06_amp06.pt": "150"}


def unwrap(p, ref, N):
    """Periodic image of the points p closest to ref."""
    return p - N * np.round((p - ref) / N)


def fig3d(ids, census, run, out, N):
    bp_all = [(f, b) for f, v in census.items() if f in STAGE_MT for b in v["bloch_points"]]
    pid = bp_all[0][1]["particle"] if bp_all else 0
    # reference: the centre of the particle with the Bloch points (circular mean over the periodic box)
    pc = np.argwhere(ids == pid).astype(float) * (2 * np.pi / N)
    ref = (np.angle(np.exp(1j * pc).mean(0)) % (2 * np.pi)) * N / (2 * np.pi)
    surf = np.zeros(ids.shape, bool)
    for a in range(3):
        for s in (1, -1):
            surf |= (ids >= 0) & (ids != np.roll(ids, s, axis=a))
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    rng = np.random.default_rng(0)
    for p in range(4):
        pos = np.argwhere(surf & (ids == p)).astype(float)
        pos = unwrap(pos, ref, N)
        keep = rng.random(len(pos)) < (0.6 if p == pid else 0.15)
        ax.scatter(*pos[keep].T, s=4 if p == pid else 2, color=PCOL[p], alpha=0.35 if p == pid else 0.12,
                   linewidths=0, depthshade=False)
        ax.plot([], [], ls="none", marker="o", ms=6, color=PCOL[p],
                label="particle %d%s" % (p, " (with the Bloch points)" if p == pid else ""))
    for q, col in ((1.0, "#000000"), (-1.0, "#7F7F7F")):
        tr = [(f, b) for f, b in bp_all if b["Q"] * q > 0]
        for f, b in tr:
            c = unwrap(np.asarray(b["centre"], float), ref, N)
            ax.scatter(*c, s=70, color=col, marker="P" if q > 0 else "X", depthshade=False)
            if f in ("init.pt", "ckpt_03_amp03.pt"):
                ax.text(*c + 3.0, "start" if f == "init.pt" else "after 50 mT", fontsize=9, color=col)
        if len(tr) > 1:
            pts = np.array([unwrap(np.asarray(b["centre"], float), ref, N) for _, b in tr])
            ax.plot(*pts.T, color=col, lw=1.5)
    ax.plot([], [], ls="none", ms=9, color="#000000", marker="P", label="Bloch point Q = +1 (path: start … 50 mT)")
    ax.plot([], [], ls="none", ms=9, color="#7F7F7F", marker="X", label="Bloch point Q = −1 (path: start … 50 mT)")
    ax.set_xlabel("x (cells of 3 $l_\\mathrm{ex}$)")
    ax.set_ylabel("y (cells)")
    ax.set_zlabel("z (cells)")
    ax.set_title("%s: particles (surface cells) and Bloch points" % run, fontsize=12)
    ax.legend(loc="upper left", fontsize=9)
    ax.view_init(elev=18, azim=-60)
    lo, hi = ref - N / 2, ref + N / 2
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    fig.tight_layout()
    fig.savefig(out / ("%s_bloch_3d.png" % run))
    plt.close(fig)


def fig_slices(m, ids, bps, run, state, out, N, half=16):
    if not bps:
        return
    fig, axs = plt.subplots(1, len(bps), figsize=(6.2 * len(bps), 5.8))
    axs = np.atleast_1d(axs)
    for ax, b in zip(axs, bps):
        cx, cy, cz = b["centre"]
        ix = np.arange(cx - half, cx + half + 1) % N
        iy = np.arange(cy - half, cy + half + 1) % N
        M = m[np.ix_(ix, iy, [cz % N])][:, :, 0]
        I = ids[np.ix_(ix, iy, [cz % N])][:, :, 0]
        X, Y = np.meshgrid(np.arange(-half, half + 1), np.arange(-half, half + 1), indexing="ij")
        inside = I == b["particle"]
        im = ax.pcolormesh(X, Y, np.where(inside, M[..., 2], np.nan), cmap="RdBu_r", vmin=-1, vmax=1, shading="nearest")
        ax.contour(X, Y, inside.astype(float), levels=[0.5], colors="k", linewidths=1.0)
        sk = (slice(None, None, 2), slice(None, None, 2))
        ax.quiver(X[sk], Y[sk], np.where(inside, M[..., 0], 0)[sk], np.where(inside, M[..., 1], 0)[sk], scale=22,
                  width=0.004, color="k")
        ax.plot(0, 0, marker="P" if b["Q"] > 0 else "X", ms=14, color="#D98E04", mec="k")
        ax.set_aspect("equal")
        ax.set_title("Bloch point Q = %+.0f, particle %d, slice z = %d" % (b["Q"], b["particle"], cz), fontsize=11)
        ax.set_xlabel("x − x$_\\mathrm{BP}$ (cells)")
        ax.set_ylabel("y − y$_\\mathrm{BP}$ (cells)")
        fig.colorbar(im, ax=ax, shrink=0.8, label="$m_z$")
    fig.suptitle("%s, %s: magnetization around the Bloch points (arrows: $m_x$, $m_y$; colour: $m_z$; line: particle "
                 "surface)" % (run, state), fontsize=11)
    fig.tight_layout()
    fig.savefig(out / ("%s_bloch_slices_%s.png" % (run, state.replace(".pt", ""))))
    plt.close(fig)


def main():
    import torch
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--state", default="init.pt")
    ap.add_argument("--state_dir", default=None, help="folder of the state file (default: --run)")
    ap.add_argument("--census", required=True)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    style()
    run = pathlib.Path(a.run)
    cfg = json.loads((run / "config.json").read_text())
    ids, _, N = A.geometry(cfg)
    census = json.loads(pathlib.Path(a.census).read_text())[run.name]
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    fig3d(ids, census, run.name, out, N)
    m = torch.load(str(pathlib.Path(a.state_dir or run) / a.state), map_location="cpu")["m"].double().numpy()
    m /= np.maximum(np.linalg.norm(m, axis=-1, keepdims=True), 1e-30)
    fig_slices(m, ids, census[a.state]["bloch_points"], run.name, a.state, out, N)


if __name__ == "__main__":
    main()
