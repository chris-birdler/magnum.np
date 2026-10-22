"""K1 investigation (PROTOKOLL 7.135, Chris 2026-10-21 option 1): where the magnetization switches between H = 0
falling (phase 180°) and H = 0 rising (phase 360°) in the last 50 mT cycle of a replay run. Per particle: switched
V_sw (angle weighted: the mean angle / 180°), the volume shares and cell counts with an angle above 15, 30, 45 and 90°,
the largest angle, the radial distribution of V_sw (r / R), and for the cells above 30°: their radial position, their
position along the particle easy axis u_p, their distance from the axis through the particle centre, m.u_p before and
after, and the azimuth concentration around u_p. Runs on the instance, reads only.
    python3 eval_k1_switch.py <runs dir> <code dir> K1_4_rep50 K1_5_rep50 ...
"""
import glob
import json
import pathlib
import sys

import numpy as np

runs, code = pathlib.Path(sys.argv[1]), sys.argv[2]
sys.path.insert(0, code)
from geometry import fcc_box, FCC_SITES          # noqa: E402
from anisotropy_model import particle_axes       # noqa: E402
import domains as DM                             # noqa: E402

def azimuth_conc(v, n):
    """Azimuth concentration: the mean resultant length of the azimuth of the vectors v around the axis n. 1: all at
    one azimuth (a patch); near 0: no single azimuth (several patches or a ring)."""
    e1 = np.cross(n, [1.0, 0.0, 0.0] if abs(n[0]) < 0.9 else [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    phi = np.arctan2(v @ e2, v @ e1)
    return float(abs(np.exp(1j * phi).mean()))


def centroid_dir(v):
    c = v.mean(0)
    return c / np.linalg.norm(c)


out = {}
for name in sys.argv[3:]:
    run = runs / name
    cfg = json.loads((run / "config.json").read_text())
    N, d, dx = cfg["N_used"], cfg["d"], cfg["dx"]
    ids = fcc_box(N, d, dx)["ids"]
    u = particle_axes(cfg.get("direction", [1.0, 0.0, 0.0]))
    files = sorted(glob.glob(str(run / "m_amp03_c*_ph*.vti")))
    cyc = max(int(pathlib.Path(f).name.split("_c")[1].split("_")[0]) for f in files)
    fa = run / ("m_amp03_c%d_ph180.vti" % cyc)       # H = 0 falling
    fb = run / ("m_amp03_c%d_ph360.vti" % cyc)       # H = 0 rising (end of the cycle)
    ma, mb = DM.load_state(fa)["m"].astype(np.float64), DM.load_state(fb)["m"].astype(np.float64)
    a = N * dx
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    ang = np.degrees(np.arccos(np.clip((ma * mb).sum(-1), -1.0, 1.0)))
    parts = []
    for p in range(4):
        s = FCC_SITES[p] * a
        dd = np.stack([v - c - a * np.round((v - c) / a) for v, c in zip((X, Y, Z), s)], -1) / (0.5 * d)
        sel = ids == p
        r = np.sqrt((dd**2).sum(-1))
        ax = dd @ u[p]                                # position along the easy axis (units of R)
        rho = np.sqrt(np.maximum(r**2 - ax**2, 0.0))  # distance from the axis through the centre (units of R)
        rec = {"cells": int(sel.sum()), "V_sw_share": float((ang[sel] / 180.0).mean()),
               "max_deg": float(ang[sel].max())}
        for t in (15.0, 30.0, 45.0, 90.0):
            rec["share_gt%d" % t] = float((sel & (ang > t)).sum() / sel.sum())
            rec["n_gt%d" % t] = int((sel & (ang > t)).sum())
        sw = sel & (ang > 30.0)
        for lo, hi in ((0.0, 0.5), (0.5, 0.8), (0.8, 2.0)):
            sh = sel & (r >= lo) & (r < hi)
            rec["Vsw_share_r%.1f_%.1f" % (lo, hi)] = float((ang[sh] / 180.0).sum() / max((ang[sel] / 180.0).sum(), 1e-30))
        if sw.any():
            rec.update({
                "r_mean": float(r[sw].mean()), "r_q10_q90": [float(np.quantile(r[sw], q)) for q in (0.1, 0.9)],
                "axial_q10_q90": [float(np.quantile(ax[sw], q)) for q in (0.1, 0.9)],
                "rho_q10_q90": [float(np.quantile(rho[sw], q)) for q in (0.1, 0.9)],
                "m_u_before_mean": float((ma[sw] @ u[p]).mean()), "m_u_after_mean": float((mb[sw] @ u[p]).mean()),
                "m_u_before_abs_mean": float(np.abs(ma[sw] @ u[p]).mean()),
                "azimuth_conc": azimuth_conc(dd[sw], u[p]),
                "centroid_dir": [float(v) for v in centroid_dir(dd[sw])],
            })
        parts.append(rec)
    out[name] = {"cycle": cyc, "files": [fa.name, fb.name], "particles": parts}
print(json.dumps(out))
