#!/usr/bin/env python
"""
Where in the particle is the loss? (PROTOKOLL 7.23) From the 4 snapshots of the last cycle of each
stage (H max, H = 0 falling, H min, H = 0 rising) per cell the fundamental of m at the drive frequency:
    a = (m(90) - m(270)) / 2,  b = (m(360) - m(180)) / 2,   |m1|^2 = |a|^2 + |b|^2
and its Gilbert dissipation (LLG: p = alpha mu0 Ms / gamma |dm/dt|^2):
    p_fund = alpha mu0 Ms / gamma * omega^2 |m1|^2 / 2.
Per stage, written to <run>/dissmap.json:
  R          = <p_fund> / <p_dis>: the share of the measured dissipation (samples csv, same cycle) that
               the local rotation at the drive frequency explains; R << 1 means internal motion at
               higher frequencies (precession, ringing) that 4 snapshots cannot see
  shell      share of p_fund in the outer shell r/R_p > 0.8 (volume share 0.488) and in r/R_p < 0.5
               (volume share 0.125)
  texture    share of p_fund in the 10 % of the cells with the largest |grad m| (smoothed 0.5 L_eff)
  top10      share of p_fund in the 10 % of the cells with the largest p_fund (localization)
Limits: 4 samples per cycle: the 3rd harmonic aliases into the fundamental; motion at higher
frequencies is not seen (this is what R measures).

    python dissmap_batch.py runs/B300_s921 [...]
"""
import json
import math
import pathlib
import re
import sys

import numpy as np

import domains
from geometry import FCC_SITES

MU0, GAMMA = 4e-7 * math.pi, 2.21276157e5


def radii(ids, dx, N):
    """r / R_p per cell (distance to the own particle centre, minimum image)."""
    a = N * dx
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.full(ids.shape, np.nan)
    for p, s in enumerate(FCC_SITES * a):
        sel = ids == p
        d2 = 0.0
        for c, sc in zip((X, Y, Z), s):
            dd = c[sel] - sc
            dd -= a * np.round(dd / a)
            d2 = d2 + dd ** 2
        r[sel] = np.sqrt(d2)
    rp = np.nanmax(r)
    return r / rp


def run(d):
    d = pathlib.Path(d)
    cfg = json.loads((d / "config.json").read_text())
    alpha, Js, freq = cfg["alpha"], cfg["Js"], cfg["freq"]
    Ms = Js / MU0
    w = 2 * math.pi * freq
    groups = {}
    for f in sorted(d.glob("m_*_c*_ph*.vti")):
        m = re.match(r"m_(.+)_c(\d+)_ph(\d+)\.vti", f.name)
        groups.setdefault((m.group(1), int(m.group(2))), {})[int(m.group(3))] = f
    out = {"run": d.name, "stages": []}
    rr = None
    for (stage, ci), ph in sorted(groups.items()):
        if sorted(ph) != [90, 180, 270, 360]:
            continue
        S = {k: domains.load_state(v) for k, v in ph.items()}
        st = S[90]
        ids, mask = st["ids"], st["mask"]
        if rr is None:
            rr = radii(ids, st["dx"], ids.shape[0])
        a = 0.5 * (S[90]["m"] - S[270]["m"])
        b = 0.5 * (S[360]["m"] - S[180]["m"])
        m1 = (a * a).sum(-1) + (b * b).sum(-1)
        pf = np.where(mask, alpha * MU0 * Ms / GAMMA * w * w * m1 / 2.0, 0.0)
        p_mean = pf.mean()                                    # per unit-cell volume (as p_dis)
        csv = next(d.glob("samples_*_%s.csv" % stage))
        s = np.loadtxt(csv, delimiter=",", comments="#")
        cyc = s[s[:, 0].astype(int) == ci]
        p0 = cyc[:-1, 8].mean()
        tot = pf[mask].sum()
        g = domains.grad_norm(domains.smooth_m(st["m"], ids, 0.5 * st["Leff"] / st["dx"]), ids, st["dx"])
        hi = mask & (g >= np.quantile(g[mask], 0.9))
        top = mask & (pf >= np.quantile(pf[mask], 0.9))
        B = 1000.0 * np.abs(cyc[:, 7]).max()
        out["stages"].append({"stage": stage, "cycle": ci, "B_peak_mT": float(B), "R": float(p_mean / p0),
                              "shell_outer": float(pf[mask & (rr > 0.8)].sum() / tot),
                              "vol_outer": float((mask & (rr > 0.8)).sum() / mask.sum()),
                              "shell_inner": float(pf[mask & (rr < 0.5)].sum() / tot),
                              "vol_inner": float((mask & (rr < 0.5)).sum() / mask.sum()),
                              "texture_top10": float(pf[hi].sum() / tot), "top10": float(pf[top].sum() / tot)})
    out["stages"].sort(key=lambda r: r["B_peak_mT"])
    (d / "dissmap.json").write_text(json.dumps(out, indent=1))
    print("%s: %d stages" % (d.name, len(out["stages"])), flush=True)


if __name__ == "__main__":
    for d in sys.argv[1:]:
        run(d)
