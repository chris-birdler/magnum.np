#!/usr/bin/env python
"""
Reduce the per-cell accumulators of one run (run_loops --acc, audit 3 item 9b) to small results, written to
<run>/acc.json and <run>/acc_slices.npz. Per measured stage (acc_<i>_<stage>.npz, cycles >= 2, all samples):
  check      p_box / p_dis: box mean of the local dissipation against the measured p_dis of the same samples
             (must be 1; the accumulation is exact)
  R1         p_fund / p_dis, p_fund = alpha mu0 Ms / gamma (omega^2 |m1|^2 / 2) per cell (fundamental from all
             samples: no aliasing); 1 - R1 = loss in harmonics and non-periodic motion
  where      for the classes r/R < 0.5, 0.5 ... 0.8, >= 0.8 (r distance to the particle centre, R = d/2), texture top 10 % (of the cycle-mean state m0,
             smoothed 0.5 L_eff, as phase_extra.py), rest, and the 10 % cells with the largest p: volume share,
             share of the local loss p, share of p_fund, mean and spread of the local lag of the fundamental
             along the drive (weighted by |m1_par|^2, relative to the drive fundamental h1)
  slices     p and |m1|^2 in the plane z = first cell layer (float16)
    python acc_reduce.py runs/S001 [...]
"""
import json
import math
import pathlib
import re
import sys

import numpy as np

import domains
from dissmap_batch import radii
from phase_extra import wspread

MU0, GAMMA = 4e-7 * math.pi, 2.21276157e5
CHECK_TOL = 1e-3


def run(d):
    d = pathlib.Path(d)
    cfg = json.loads((d / "config.json").read_text())
    alpha, Ms = cfg["alpha"], cfg["Js"] / MU0
    e = np.array(cfg["direction"], float)
    e /= np.linalg.norm(e)
    ids = domains._ids_from_config(cfg).astype(np.int8)
    mask = ids >= 0
    N = ids.shape[0]
    dx, Leff = cfg["dx_lex"], cfg["Leff_lex"]
    r = radii(ids, dx, N)[mask]
    out = {"run": d.name, "stages": []}
    slices = {}
    for f in sorted(d.glob("acc_*_*.npz")):
        if f.name.endswith(".tmp.npz"):
            continue
        stage = re.match(r"acc_(\d+)_(.+)\.npz", f.name).group(2)
        z = np.load(f)
        p, m0, m1, h1 = z["p"].astype(np.float64), z["m0"], z["m1"].astype(np.complex128), complex(z["h1"])
        if len(p) != int(mask.sum()):
            raise ValueError("%s: %d cells, the geometry of config.json has %d particle cells" % (f.name, len(p),
                                                                                                 int(mask.sum())))
        w = 2 * math.pi * float(z["freq"])
        check = float(z["p_box"]) / float(z["p_dis_samples"])
        pf = alpha * MU0 * Ms / GAMMA * w * w * (np.abs(m1) ** 2).sum(-1) / 2.0
        p_dis_cellmean = p.sum() / ids.size                  # box mean (void cells dissipate nothing)
        R1 = float(pf.sum() / ids.size / p_dis_cellmean)
        a1 = (m1 * e).sum(-1)
        lag = np.angle(a1 / h1) if abs(h1) > 0 else np.zeros_like(r)
        wgt = np.abs(a1) ** 2
        full = np.zeros(ids.shape + (3,), np.float32)
        full[mask] = m0 / np.maximum(np.linalg.norm(m0, axis=-1, keepdims=True), 1e-12)
        g = domains.grad_norm(domains.smooth_m(full, ids, 0.5 * Leff / dx), ids, dx)[mask]
        hi = g >= np.quantile(g, 0.9)
        ptop = p >= np.quantile(p, 0.9)
        classes = {"r<0.5": r < 0.5, "0.5<=r<0.8": (r >= 0.5) & (r < 0.8), "r>=0.8": r >= 0.8,
                   "texture top 10%": hi, "rest of texture": ~hi, "p top 10%": ptop}
        where = {}
        for k, sel in classes.items():
            mean_lag, spread = wspread(lag[sel], wgt[sel])
            where[k] = {"vol_share": float(sel.mean()), "p_share": float(p[sel].sum() / p.sum()),
                        "p_fund_share": float(pf[sel].sum() / pf.sum()),
                        "lag_mean_deg": math.degrees(mean_lag), "lag_spread_deg": math.degrees(spread)}
        B = 1000.0 * np.abs(np.loadtxt(next(d.glob("samples_*_%s.csv" % stage)), delimiter=",", comments="#")[:, 7]).max()
        out["stages"].append({"stage": stage, "B_peak_mT": float(B), "n_samples": int(z["n_samples"]),
                              "cycles_from": int(z["cycles_from"]), "check_p_box_over_p_dis": check,
                              "check_ok": abs(check - 1.0) < CHECK_TOL, "R1": R1,
                              "lag_drive_deg": math.degrees(float(np.angle(a1.sum() / h1))) if abs(h1) > 0 else None,
                              "where": where})
        fp = np.zeros(ids.shape, np.float32)
        fp[mask] = p
        fm = np.zeros(ids.shape, np.float32)
        fm[mask] = (np.abs(m1) ** 2).sum(-1)
        slices[stage + "_p"] = fp[:, :, 0].astype(np.float16)
        slices[stage + "_m1sq"] = fm[:, :, 0].astype(np.float16)
        print("%s %s: check %.6f  R1 %.3f  p share r<0.5 %.2f  r>=0.8 %.2f  texture top10 %.2f  p top10 %.2f" %
              (d.name, stage, check, R1, where["r<0.5"]["p_share"], where["r>=0.8"]["p_share"],
               where["texture top 10%"]["p_share"], where["p top 10%"]["p_share"]), flush=True)
    slices["ids"] = ids[:, :, 0]
    out["stages"].sort(key=lambda s: s["B_peak_mT"])
    (d / "acc.json").write_text(json.dumps(out, indent=1))
    np.savez_compressed(d / "acc_slices.npz", **slices)


if __name__ == "__main__":
    for d in sys.argv[1:]:
        run(d)
