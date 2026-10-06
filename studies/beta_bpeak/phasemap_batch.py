#!/usr/bin/env python
"""
Phase-resolved local response from the snapshots of the last cycles (PROTOKOLL 7.26), written to
<run>/phasemap.json. Needs --snap_phases N >= 8 and --snap_cycles C >= 2 (phase16: N = 16, C = 3).

Per stage and cell, from the N snapshots of each of the C cycles:
  m_mean(phase)  mean over the C cycles; m_np = deviation of each cycle from it (non-periodic part:
                 jumps, drift, ring-down that does not repeat)
  harmonics      of m_mean at k = 1 … N/2 - 1 (no aliasing below N/2)
  Gilbert dissipation of each part, p = alpha mu0 Ms / gamma <|dm/dt|^2>:
    p_k      = alpha mu0 Ms / gamma (k omega)^2 |m_k|^2 / 2           (periodic, harmonic k)
    p_np     = alpha mu0 Ms / gamma <|d m_np / dt|^2> (finite differences of the phase series; a
               lower bound: motion faster than T/N is not resolved)
  all per unit-cell volume, compared with p_dis (samples csv, same cycles).
  phase: in-phase / quadrature part of the local fundamental along the drive (lag of each cell).
Output per stage: R1 = p_1 / p_dis, R_harm = sum_{k>=2} p_k / p_dis, R_np = p_np / p_dis, R_rest =
1 - R1 - R_harm - R_np (motion faster than T/N or not periodic within the sampling); spread of the
local phase lag (weighted by |m_1|^2).

    python phasemap_batch.py runs/PH_s921 [...]
"""
import json
import math
import pathlib
import re
import sys

import numpy as np

import domains

MU0, GAMMA = 4e-7 * math.pi, 2.21276157e5


def phase_angles(phases):
    """Drive phase (rad) of the N snapshots of one cycle, from their INDEX: theta_j = 2 pi (j + 1) / N.
    The file names hold the phase rounded to whole degrees (run_loops: round(360 (j + 1) / samples)), e.g. 22
    and 68 for 22.5 and 67.5 deg; angles taken from the names make the grid uneven and leak the static part
    of m (|m| ~ 1) into the harmonics (audit 2026-10-06: p_4/p_dis up to 44). The names are only checked."""
    N = len(phases)
    want = [round(360.0 * (j + 1) / N) for j in range(N)]
    if list(phases) != want:
        raise ValueError("snapshot phases %s are not the %d equidistant phases %s" % (list(phases), N, want))
    return 2 * math.pi * np.arange(1, N + 1) / N


def run(d):
    d = pathlib.Path(d)
    cfg = json.loads((d / "config.json").read_text())
    alpha, Ms, freq = cfg["alpha"], cfg["Js"] / MU0, cfg["freq"]
    w = 2 * math.pi * freq
    coef = alpha * MU0 * Ms / GAMMA
    e = np.array(cfg["direction"], float)
    e /= np.linalg.norm(e)
    files = {}
    for f in d.glob("m_*_c*_ph*.vti"):
        m = re.match(r"m_(.+)_c(\d+)_ph(\d+)\.vti", f.name)
        files.setdefault(m.group(1), {}).setdefault(int(m.group(2)), {})[int(m.group(3))] = f
    out = {"run": d.name, "stages": []}
    for stage, cycles in sorted(files.items()):
        cis = sorted(cycles)
        phases = sorted(cycles[cis[0]])
        N, C = len(phases), len(cis)
        st0 = domains.load_state(cycles[cis[0]][phases[0]])
        mask = st0["mask"]
        M = np.zeros((C, N, int(mask.sum()), 3), np.float32)          # cycles, phases, cells
        for i, ci in enumerate(cis):
            for j, ph in enumerate(phases):
                M[i, j] = domains.load_state(cycles[ci][ph])["m"][mask]
        Mm = M.mean(0)                                                  # periodic part [N, cells, 3]
        if any(sorted(cycles[ci]) != phases for ci in cis):
            raise ValueError("%s: the cycles hold different phases" % stage)
        th = phase_angles(phases)                                       # drive phase of each snapshot
        F = [None] + [np.tensordot(np.exp(-1j * k * th), Mm, axes=(0, 0)) * (2.0 / N) for k in range(1, N // 2)]
        pk = [float(coef * (k * w) ** 2 * (np.abs(F[k]) ** 2).sum(-1).sum() / 2.0) for k in range(1, N // 2)]
        dt = 1.0 / (freq * N)
        dnp = np.diff(np.concatenate([M - Mm[None], (M - Mm[None])[:, :1]], axis=1), axis=1) / dt
        p_np = float(coef * (dnp ** 2).sum(-1).mean(axis=(0, 1)).sum())
        n_all = mask.size
        pk = [p / n_all for p in pk]
        p_np /= n_all
        # local phase lag of the fundamental along the drive (phase 0 = H = 0 rising; H = sin)
        a1 = (F[1] * e).sum(-1)                                        # complex, per cell
        lag = np.angle(a1 / (-1j))                                     # relative to sin
        wgt = np.abs(a1) ** 2
        lag_mean = float(np.angle(np.sum(wgt * np.exp(1j * lag))))
        lag_spread = float(math.sqrt(np.sum(wgt * np.angle(np.exp(1j * (lag - lag_mean))) ** 2) / np.sum(wgt)))
        csv = next(d.glob("samples_*_%s.csv" % stage))
        s = np.loadtxt(csv, delimiter=",", comments="#")
        sel = np.isin(s[:, 0].astype(int), cis)
        p_dis = float(s[sel, 8].mean())
        B = 1000.0 * np.abs(s[sel, 7]).max()
        R1, Rh, Rn = pk[0] / p_dis, sum(pk[1:]) / p_dis, p_np / p_dis
        out["stages"].append({"stage": stage, "B_peak_mT": float(B), "N": N, "C": C, "p_dis": p_dis,
                              "R1": R1, "R_harm": Rh, "R_np": Rn, "R_rest": 1.0 - R1 - Rh - Rn,
                              "p_k_over_p_dis": [p / p_dis for p in pk],
                              "lag_mean_rad": lag_mean, "lag_spread_rad": lag_spread})
        print("%s %s: R1 %.3f  R_harm %.3f  R_np %.3f  rest %.3f  lag %.3f ± %.3f rad" %
              (d.name, stage, R1, Rh, Rn, 1 - R1 - Rh - Rn, lag_mean, lag_spread), flush=True)
    out["stages"].sort(key=lambda r: r["B_peak_mT"])
    (d / "phasemap.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    for d in sys.argv[1:]:
        run(d)
