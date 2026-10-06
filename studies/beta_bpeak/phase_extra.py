#!/usr/bin/env python
"""
Further evaluation of the 16-phase snapshots (PROTOKOLL 7.27), written to <run>/phase_extra.json and
<run>/slices.npz (small; for the cut images, rendered locally by make_phase_slices.py).
Per stage (16 phases x 3 cycles; the periodic part = mean over the cycles):
  where   share of p_fund (local fundamental, Gilbert) and of the cells in radial shells r/R_p < 0.5,
          0.5 ... 0.8, > 0.8 and in the 10 % cells of largest texture; per class the weighted spread of the
          local phase lag of the fundamental along the drive
  texture T(2 L_eff) and F90 at each of the 16 phases of the last cycle (change within a cycle)
  slices  m in the plane z = first cell layer (cuts particles 0 and 3) at the 16 phases of the last cycle

    python phase_extra.py runs/PH_s921 [...]
"""
import json
import math
import pathlib
import re
import sys

import numpy as np

import domains
from dissmap_batch import radii

MU0, GAMMA = 4e-7 * math.pi, 2.21276157e5


def wspread(lag, w):
    m = np.angle(np.sum(w * np.exp(1j * lag)))
    return float(m), float(math.sqrt(np.sum(w * np.angle(np.exp(1j * (lag - m))) ** 2) / np.sum(w)))


def run(d):
    d = pathlib.Path(d)
    cfg = json.loads((d / "config.json").read_text())
    alpha, Ms, freq = cfg["alpha"], cfg["Js"] / MU0, cfg["freq"]
    w = 2 * math.pi * freq
    e = np.array(cfg["direction"], float)
    e /= np.linalg.norm(e)
    files = {}
    for f in d.glob("m_*_c*_ph*.vti"):
        m = re.match(r"m_(.+)_c(\d+)_ph(\d+)\.vti", f.name)
        files.setdefault(m.group(1), {}).setdefault(int(m.group(2)), {})[int(m.group(3))] = f
    out = {"run": d.name, "stages": []}
    slices = {}
    rr = None
    for stage, cycles in sorted(files.items()):
        cis = sorted(cycles)
        phases = sorted(cycles[cis[0]])
        N = len(phases)
        th = 2 * math.pi * np.array(phases, float) / 360.0
        st0 = domains.load_state(cycles[cis[0]][phases[0]])
        mask, ids = st0["mask"], st0["ids"]
        if rr is None:
            rr = radii(ids, st0["dx"], ids.shape[0])
        r = rr[mask]
        Mm = np.zeros((N, int(mask.sum()), 3))
        last = []
        for ci in cis:
            for j, ph in enumerate(phases):
                s = domains.load_state(cycles[ci][ph])
                Mm[j] += s["m"][mask] / len(cis)
                if ci == cis[-1]:
                    last.append(s)
        m1 = np.tensordot(np.exp(-1j * th), Mm, axes=(0, 0)) * (2.0 / N)          # [cells, 3] complex
        pf = alpha * MU0 * Ms / GAMMA * w * w * (np.abs(m1) ** 2).sum(-1) / 2.0
        a1 = (m1 * e).sum(-1)
        lag = np.angle(a1 / (-1j))
        wgt = np.abs(a1) ** 2
        # texture of the cycle-mean state (smoothed 0.5 L_eff), top 10 % cells
        mean_state = domains.make_state(np.zeros_like(st0["m"]), ids, st0["dx"], st0["Leff"], "mean")
        full = np.zeros_like(st0["m"])
        full[mask] = Mm.mean(0) / np.maximum(np.linalg.norm(Mm.mean(0), axis=-1, keepdims=True), 1e-12)
        mean_state["m"] = full.astype(np.float32)
        g = domains.grad_norm(domains.smooth_m(mean_state["m"], ids, 0.5 * st0["Leff"] / st0["dx"]), ids, st0["dx"])[mask]
        hi = g >= np.quantile(g, 0.9)
        classes = {"r<0.5": r < 0.5, "0.5<=r<0.8": (r >= 0.5) & (r < 0.8), "r>=0.8": r >= 0.8, "texture top 10%": hi,
                   "rest of texture": ~hi}
        where = {}
        for k, sel in classes.items():
            mean_lag, spread = wspread(lag[sel], wgt[sel])
            where[k] = {"vol_share": float(sel.mean()), "p_fund_share": float(pf[sel].sum() / pf.sum()),
                        "lag_mean_deg": math.degrees(mean_lag), "lag_spread_deg": math.degrees(spread)}
        T, F90 = [], []
        for s in last:                                   # only T(2 L_eff) and F90 (cheaper than domains.measures)
            gs = domains.grad_norm(domains.smooth_m(s["m"], ids, domains.W_WALL * s["Leff"] / s["dx"]), ids, s["dx"])
            T.append(float(gs.sum() * s["dx"] ** 3 / math.pi))
            F90.append(float((domains.rotation_angles(s, window_leff=domains.W_ROT, n_lines=1500) > 90.0).mean()))
        B = 1000.0 * np.abs(np.loadtxt(next(d.glob("samples_*_%s.csv" % stage)), delimiter=",", comments="#")[:, 7]).max()
        out["stages"].append({"stage": stage, "B_peak_mT": float(B), "where": where,
                              "T_2Leff_vs_phase": T, "F90_vs_phase": F90, "phases_deg": phases})
        slices[stage] = np.stack([s["m"][:, :, 0, :] for s in last]).astype(np.float16)     # [16, N, N, 3]
        slices[stage + "_ids"] = ids[:, :, 0]
        print("%s %s: p_fund share r<0.5 %.2f (vol %.2f), texture top10 %.2f; T range %.3g … %.3g" %
              (d.name, stage, where["r<0.5"]["p_fund_share"], where["r<0.5"]["vol_share"],
               where["texture top 10%"]["p_fund_share"], min(T), max(T)), flush=True)
    out["stages"].sort(key=lambda r: r["B_peak_mT"])
    (d / "phase_extra.json").write_text(json.dumps(out, indent=1))
    np.savez_compressed(d / "slices.npz", **slices)


if __name__ == "__main__":
    for d in sys.argv[1:]:
        run(d)
