#!/usr/bin/env python
"""
Check of the grid-pinning hypothesis with the data at hand (post hoc, request of Chris 2026-10-13, PROTOKOLL 7.78).
Hypothesis: on the coarse mesh (cell 10 nm) the vortex core is about one cell wide; at small B it moves less than a
cell, so the grid holds it (less motion, less loss), at larger B it jumps from cell to cell (extra loss).
Per mesh-check point (base 0.71 µm, base 1 µm, K3 0.71 µm), mesh, seed and amplitude, from the last 2 calm cycles:
- χ_p: amplitude of the mean magnetisation of each particle along the drive, ΔM_p / Ms (mean over the 4 particles);
  a held core gives a smaller ΔM_p on the coarse mesh;
- x: core shift estimate x ≈ ΔM_p / Ms · R (rigid vortex, order of magnitude only), in nm and in coarse cells;
- sin δ = W_loop / (π μ0 H_peak M_peak) (loss angle of the linear response);
- jump share = RF / 2 of the last 3 cycles (eval_jumps_time; energy above the smooth periodic curve).
Output: means over the seeds, coarse / fine ratios paired per seed.

    python eval_pinning.py > results/pinning_check.md
"""
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import eval_jumps_time as EJ
import eval_k3fill as KF
import eval_size as ES

R = pathlib.Path(__file__).resolve().parent / "runs"
MU0 = 4e-7 * math.pi
POINTS = (("base 0.71 µm", lambda dx, s, b: ES.run_dir("M212_dx%s_s%d_b%d" % (dx, s, b)), range(931, 935), (9, 50, 100)),
          ("base 1.0 µm", lambda dx, s, b: ES.run_dir("M300_dx%s_s%d_b%d" % (dx, s, b)), range(921, 925), (9, 50, 100)),
          ("K3 0.71 µm", lambda dx, s, b: KF.folder({"15": "15", "3": "3"}[dx], s, b), KF.SEEDS, (9, 20, 35, 50, 70, 100)))


def measures(d):
    """(χ_p, x in nm, sin δ, jump share, cell nm) of the tested (last) stage of run folder d, or None."""
    if not (d / "DONE").exists():
        return None
    cfg = json.loads((d / "config.json").read_text())
    st = json.loads((d / "summary.json").read_text())["stages"][-1]
    meas = st["cycles"][2:] if "h_amp_initial" in st else st["cycles"][1:]
    tail = analyze.steady_tail(meas, AD.BAL_EVENT)
    if len(tail) < 2:
        return None
    keep = {c["cycle"] for c in tail[-2:]}
    f = sorted(d.glob("samples_*.csv"))[-1]
    a = np.loadtxt(f, delimiter=",", comments="#")
    a = a[np.isin(a[:, 0].astype(int), list(keep))]
    Ms = cfg["units"]["Ms"]
    Mp = a[:, 9:13]
    chi = float(np.mean(0.5 * (Mp.max(0) - Mp.min(0)))) / Ms
    lex = cfg["units"]["l_ex"] * 1e9
    Rn = 0.5 * cfg.get("d_lex_used", cfg["d_lex"]) * lex
    c = tail[-1]
    sind = c["W_loop"] / (math.pi * MU0 * c["H_peak"] * c["M_peak"])
    s = EJ.series(f)
    jump = EJ.measures(*s)["RF"] / 2 if s is not None else float("nan")
    return chi, chi * Rn, sind, jump, cfg["dx_lex"] * lex


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    L = ["# Grid-pinning check (post hoc, PROTOKOLL 7.78)", "",
         "Per point and B: mean over the seeds; coarse / fine paired per seed (geometric mean, range over the seeds). "
         "x = core shift estimate ΔM_p / Ms · R (rigid vortex, order of magnitude). Cell: 10 nm (coarse), 5 nm (fine).", "",
         "| point | B (mT) | n | x fine (nm) | x coarse (nm) | x coarse / cell | ΔM_p coarse / fine | sin δ coarse / fine | "
         "jump share fine | jump share coarse |", "|---|---|---|---|---|---|---|---|---|---|"]
    for lab, folder, seeds, amps in POINTS:
        for b in amps:
            rows = []
            for s in seeds:
                mf, mc = measures(folder("15", s, b)), measures(folder("3", s, b))
                if mf and mc:
                    rows.append((mf, mc))
            if len(rows) < 2:
                L.append("| %s | %d | %d | – | – | – | – | – | – | – |" % (lab, b, len(rows)))
                continue
            gm = lambda v: math.exp(np.mean(np.log(v)))
            rc = [mc[0] / mf[0] for mf, mc in rows]
            rs = [mc[2] / mf[2] for mf, mc in rows]
            L.append("| %s | %d | %d | %.1f | %.1f | %.2f | × %.2f (%.2f … %.2f) | × %.2f (%.2f … %.2f) | %.3f | %.3f |" % (
                lab, b, len(rows), np.mean([mf[1] for mf, _ in rows]), np.mean([mc[1] for _, mc in rows]),
                np.mean([mc[1] / mc[4] for _, mc in rows]), gm(rc), min(rc), max(rc), gm(rs), min(rs), max(rs),
                np.nanmean([mf[3] for mf, _ in rows]), np.nanmean([mc[3] for _, mc in rows])))
    print("\n".join(L))


if __name__ == "__main__":
    main()
