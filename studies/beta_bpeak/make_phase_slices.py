#!/usr/bin/env python
"""
Cut images of m through the simulation box (plane z = first cell layer: particles 0 and 3) at 4 phases
of the last cycle, for 9 mT and 150 mT (slices.npz of phase_extra.py). Colour: m_x (drive direction);
arrows: in-plane m (every 6th cell); grey: air. Second row per amplitude: change of m_x against the
cycle mean (where m moves during the cycle).

    python make_phase_slices.py runs/PH_s921/slices.npz   # -> results/phase_slices_PH_s921.png
"""
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
for f in pathlib.Path("/usr/share/fonts/truetype/crosextra").glob("Carlito*.ttf"):
    font_manager.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": "Carlito", "font.size": 11})
DX_NM = 10.0


def main(path):
    path = pathlib.Path(path)
    z = np.load(path)
    stages = sorted(k for k in z.files if not k.endswith("_ids"))
    sel = [stages[0], stages[-1]]                       # smallest and largest amplitude (9 and 150 mT)
    labels = {stages[0]: "9 mT", stages[-1]: "150 mT"}
    phases = [3, 7, 11, 15]                             # 90, 180, 270, 360 deg of 16 phases
    names = ["H max (90°)", "H = 0 falling (180°)", "H min (270°)", "H = 0 rising (360°)"]
    fig, axs = plt.subplots(4, 4, figsize=(13.3, 13.3), dpi=150)
    for row, st in enumerate(sel):
        m = z[st].astype(np.float32)                    # [16, N, N, 3]
        ids = z[st + "_ids"]
        air = ids < 0
        mean = m.mean(0)
        N = ids.shape[0]
        ext = (0, N * DX_NM, 0, N * DX_NM)
        for col, (j, nm) in enumerate(zip(phases, names)):
            ax = axs[2 * row, col]
            mx = np.ma.masked_where(air, m[j, :, :, 0]).T
            ax.imshow(mx, origin="lower", cmap="RdBu_r", vmin=-1, vmax=1, extent=ext)
            ax.imshow(np.ma.masked_where(~air.T, np.ones_like(mx)), origin="lower", cmap="Greys", vmin=0, vmax=3,
                      extent=ext)
            s = 6
            X, Y = np.meshgrid((np.arange(0, N, s) + 0.5) * DX_NM, (np.arange(0, N, s) + 0.5) * DX_NM, indexing="ij")
            u, v = m[j, ::s, ::s, 0], m[j, ::s, ::s, 1]
            ok = ~air[::s, ::s]
            ax.quiver(X[ok], Y[ok], u[ok], v[ok], scale=40, width=0.002, color="black")
            ax.set_title("%s, %s: m_x" % (labels[st], nm), fontsize=11)
            ax.set_xticks([]), ax.set_yticks([])
            ax2 = axs[2 * row + 1, col]
            dm = np.ma.masked_where(air, m[j, :, :, 0] - mean[:, :, 0]).T
            lim = max(1e-3, float(np.abs(m[:, :, :, 0] - mean[None, :, :, 0])[:, ~air].max()))
            ax2.imshow(dm, origin="lower", cmap="PuOr", vmin=-lim, vmax=lim, extent=ext)
            ax2.set_title("%s, %s: m_x − cycle mean (±%.2g)" % (labels[st], nm, lim), fontsize=11)
            ax2.set_xticks([]), ax2.set_yticks([])
    fig.suptitle("%s: cut z = 5 nm through particles 0 (corners) and 3 (centre); 1 µm, base point, α 0.1" %
                 path.parent.name, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = HERE / "results" / ("phase_slices_%s.png" % path.parent.name)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    for p in sys.argv[1:]:
        main(p)
