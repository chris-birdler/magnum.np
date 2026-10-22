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
  particles  per particle q (ids == q, PROTOKOLL 7.137): share of the box loss, R1 = p_fund / p of the particle,
             R_harm and R_rest of the particle (with mk2), and for the classes core (r/R < 0.5), rim (>= 0.8), texture
             top 10 % and p top 10 % (thresholds of the whole box): the volume share of the class in the particle and
             the concentration c = (share of the particle loss in the class) / (volume share)
  slices     in the plane z = first cell layer, for x = p, harm (Σ_{k=2..K} p_k), rest (p − Σ_k p_k, signed):
             <stage>_x_rel = x / <stage>_x_scale (float16, float32 if float16 overflows), <stage>_x_scale = mean |x|
             and <stage>_x_mean = mean x of the particle cells of the plane (W/m³, float64); read the values with
             slice_values(). For x ≥ 0 (p, harm) scale = mean, thus <stage>_p_rel = p / (mean p of the plane).
             <stage>_m1sq = |m1|^2 (float16). Before 2026-10-11 p was stored in W/m³ as float16 under <stage>_p
             (values above 65504 overflow to inf, PROTOKOLL 7.61). Before 2026-10-14 there was no _scale: x_rel =
             x / x_mean, and a plane with x_mean ≤ 0 was stored as zeros (the rest at 60 MHz, PROTOKOLL 7.90)
    python acc_reduce.py runs/S001 [...]
    python acc_reduce.py --out DIR runs/S001 [...]     (writes DIR/<run>/acc.json and acc_slices.npz instead)
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
F16_MAX = 6.0e4          # largest relative value stored as float16 (float16 maximum 65504)



def relative_slice(plane, inside):
    """A plane relative to a scale (float16 has a maximum of 65504; p in W/m³ reaches 1e9), the scale and the signed
    mean over the particle cells (float64). Scale = mean |value| over the particle cells: equal to the mean for a
    quantity ≥ 0, and > 0 also for a signed quantity whose mean is ≤ 0 or near 0 (PROTOKOLL 7.90). float32 when a
    value relative to the scale does not fit into float16."""
    mean = float(plane[inside].mean()) if inside.any() else 0.0
    scale = float(np.abs(plane[inside]).mean()) if inside.any() else 0.0
    rel = plane / scale if scale > 0 else np.zeros_like(plane)
    dtype = np.float16 if np.abs(rel).max(initial=0.0) < F16_MAX else np.float32
    return rel.astype(dtype), np.float64(mean), np.float64(scale)


def slice_values(z, stage, key):
    """Values of a stored plane (W/m³): rel × scale, or rel × mean in files without a scale (before 2026-10-14; there a
    plane with a mean ≤ 0 is zeros, use plane_known)."""
    sk = "%s_%s_scale" % (stage, key)
    f = z[sk] if sk in z.files else z["%s_%s_mean" % (stage, key)]
    return z["%s_%s_rel" % (stage, key)].astype(np.float64) * float(f)


def plane_known(z, stage, key):
    """False for a plane that an old file (no scale) stored as zeros because its mean was ≤ 0."""
    return "%s_%s_scale" % (stage, key) in z.files or float(z["%s_%s_mean" % (stage, key)]) > 0

def harmonic_loss(mk2, alpha, Ms, w):
    """Loss density of the local motion at k w per cell (W/m³), k = 1 … K: p_k = alpha mu0 Ms / gamma (k w)^2
    |m_k|^2 / 2 (the same coefficient as pf; Gilbert damping: p = alpha mu0 Ms / gamma |dm/dt|^2). mk2: (cells, K)."""
    k = np.arange(1, mk2.shape[1] + 1, dtype=np.float64)
    return alpha * MU0 * Ms / GAMMA * w * w * k[None, :] ** 2 * mk2.astype(np.float64) / 2.0


def run(d, out_dir=None):
    """Reduce the accumulators of run folder d; write acc.json and acc_slices.npz into out_dir (default: d)."""
    d = pathlib.Path(d)
    out_dir = pathlib.Path(out_dir) if out_dir else d
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
        harm = None
        if "mk2" in z.files:                                 # run_loops --acc_harm (PROTOKOLL 7.72), calm cycles only
            pk = harmonic_loss(z["mk2"], alpha, Ms, w)
            ph, prest = pk[:, 1:].sum(1), p - pk.sum(1)        # harmonics k = 2 … K, rest (> K w, non-repeating)
            ratio = pk.sum(1)[p > 0] / p[p > 0]               # Bessel: sum_k p_k <= p in each cell
            harm = {"K": int(pk.shape[1]), "R_k": (pk.sum(0) / p.sum()).tolist(),
                    "R_rest": float(prest.sum() / p.sum()), "bessel_max": float(ratio.max()),
                    "bessel_frac_over_1e-3": float((ratio > 1.001).mean()),
                    "cycles_kept": [int(c) for c in z["cycles_kept"]], "cycles_acc": int(z["cycles_acc"])}
        where = {}
        for k, sel in classes.items():
            mean_lag, spread = wspread(lag[sel], wgt[sel])
            where[k] = {"vol_share": float(sel.mean()), "p_share": float(p[sel].sum() / p.sum()),
                        "p_fund_share": float(pf[sel].sum() / pf.sum()),
                        "lag_mean_deg": math.degrees(mean_lag), "lag_spread_deg": math.degrees(spread)}
            if harm is not None:
                where[k]["p_harm_share"] = float(ph[sel].sum() / ph.sum())
                where[k]["p_rest_share"] = float(prest[sel].sum() / prest.sum())
        pid = ids[mask]
        particles = []
        for q in range(int(pid.max()) + 1):
            sq = pid == q
            rec = {"p_share": float(p[sq].sum() / p.sum()), "R1": float(pf[sq].sum() / p[sq].sum())}
            if harm is not None:
                rec.update({"R_harm": float(ph[sq].sum() / p[sq].sum()), "R_rest": float(prest[sq].sum() / p[sq].sum())})
            for k in ("r<0.5", "r>=0.8", "texture top 10%", "p top 10%"):
                v = float(classes[k][sq].mean())
                rec["vol " + k] = v
                rec["c " + k] = float(p[sq & classes[k]].sum() / p[sq].sum() / v) if v > 0 else float("nan")
            particles.append(rec)
        B = 1000.0 * np.abs(np.loadtxt(next(d.glob("samples_*_%s.csv" % stage)), delimiter=",", comments="#")[:, 7]).max()
        out["stages"].append({"stage": stage, "B_peak_mT": float(B), "n_samples": int(z["n_samples"]),
                              "cycles_from": int(z["cycles_from"]), "check_p_box_over_p_dis": check,
                              "check_ok": abs(check - 1.0) < CHECK_TOL, "R1": R1,
                              "lag_drive_deg": math.degrees(float(np.angle(a1.sum() / h1))) if abs(h1) > 0 else None,
                              "where": where, "particles": particles, **({"harm": harm} if harm is not None else {})})
        fp = np.zeros(ids.shape, np.float32)
        fp[mask] = p
        fm = np.zeros(ids.shape, np.float32)
        fm[mask] = (np.abs(m1) ** 2).sum(-1)
        slices[stage + "_p_rel"], slices[stage + "_p_mean"], slices[stage + "_p_scale"] = relative_slice(
            fp[:, :, 0], ids[:, :, 0] >= 0)
        slices[stage + "_m1sq"] = fm[:, :, 0].astype(np.float16)
        if harm is not None:
            for key, v in (("harm", ph), ("rest", prest)):
                fh = np.zeros(ids.shape, np.float64)
                fh[mask] = v
                (slices[stage + "_%s_rel" % key], slices[stage + "_%s_mean" % key],
                 slices[stage + "_%s_scale" % key]) = relative_slice(fh[:, :, 0], ids[:, :, 0] >= 0)
            print("%s %s: R_k %s  rest %.3f  Bessel max %.4f  calm cycles %s of %d" % (
                d.name, stage, " ".join("%.3f" % x for x in harm["R_k"]), harm["R_rest"], harm["bessel_max"],
                harm["cycles_kept"], harm["cycles_acc"]), flush=True)
        print("%s %s: check %.6f  R1 %.3f  p share r<0.5 %.2f  r>=0.8 %.2f  texture top10 %.2f  p top10 %.2f" %
              (d.name, stage, check, R1, where["r<0.5"]["p_share"], where["r>=0.8"]["p_share"],
               where["texture top 10%"]["p_share"], where["p top 10%"]["p_share"]), flush=True)
    slices["ids"] = ids[:, :, 0]
    out["stages"].sort(key=lambda s: s["B_peak_mT"])
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "acc.json").write_text(json.dumps(out, indent=1))
    np.savez_compressed(out_dir / "acc_slices.npz", **slices)


if __name__ == "__main__":
    args = sys.argv[1:]
    top = None
    if args[:1] == ["--out"]:
        top, args = pathlib.Path(args[1]), args[2:]
    for d in args:
        run(d, top / pathlib.Path(d).name if top else None)
