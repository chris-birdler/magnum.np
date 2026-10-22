#!/usr/bin/env python
"""
K1 investigation (PROTOKOLL 7.135, Chris 2026-10-21): the loop area of each particle in the runs K<c>_<k>_sat.
Per run and stage: w_p = mu0 * closed integral of H dM_p over one cycle (M_p: the mean of M_s m.e of the cells of
particle p, column Mp<p> of samples_<i>_<stage>.csv; H: the applied field), in J/m³ of the particle volume, mean of
the last 2 cycles of the stage. Particle order of geometry.FCC_SITES; cos of the angle between the particle easy
axis and the drive from anisotropy_model.particle_axes.
Also dM0_p: M_p at H = 0 on the falling branch minus M_p at H = 0 on the rising branch of the last cycle (linear
interpolation at the sign changes of H; the cycle is closed), over M_s. q_ell_p: dM0_p over the value of an
elliptic loop with the same area, 2 w_p / (pi mu0 H_hat M_s) (H_hat: half the peak-to-peak value of H in the last
cycle); 1 for a linear response with a phase lag.
Check: phi_vox / 4 * sum_p w_p must be equal to the loop area of the box (mu0 * closed integral of H dM_par, same
cycles) and to W_loop of summary.json (the 4 particles have the same number of cells; ms_void is ignored).

    python eval_k1_loops.py [--runs DIR] [--stages 00_amp00 03_amp03] [--json FILE] > report.md
"""
import argparse
import json
import pathlib

import numpy as np

from anisotropy_model import particle_axes

MU0 = 4e-7 * np.pi
HERE = pathlib.Path(__file__).resolve().parent
LABEL = {"00_amp00": "9 mT", "01_amp01": "20 mT", "02_amp02": "35 mT", "03_amp03": "50 mT", "04_amp04": "70 mT",
         "05_amp05": "100 mT", "06_amp06": "150 mT"}


def loop(H, M):
    """mu0 * closed integral of H dM over the samples of one cycle (trapezoid rule, the loop is closed)."""
    Hc, Mc = np.append(H, H[0]), np.append(M, M[0])
    return MU0 * float(np.sum(0.5 * (Hc[1:] + Hc[:-1]) * np.diff(Mc)))


def zero_branches(H, M):
    """M at H = 0 on the falling and on the rising branch of one closed cycle (linear interpolation; the samples
    cover one period without its end point)."""
    out = {}
    for i in range(len(H)):
        j = (i + 1) % len(H)
        if H[i] > 0.0 >= H[j] or H[i] < 0.0 <= H[j]:
            out["falling" if H[i] > 0.0 else "rising"] = M[i] + (M[j] - M[i]) * H[i] / (H[i] - H[j])
    return out["falling"], out["rising"]


def run_record(run, stage):
    a = np.loadtxt(run / ("samples_%s.csv" % stage), delimiter=",")
    cfg = json.loads((run / "config.json").read_text())
    cyc = sorted(set(int(c) for c in a[:, 0]))[-2:]
    wp, wb = [], []
    for c in cyc:
        s = a[a[:, 0] == c]
        wp.append([loop(s[:, 2], s[:, 9 + p]) for p in range(4)])
        wb.append(loop(s[:, 2], s[:, 3]))
    summ = json.loads((run / "summary.json").read_text())
    st = [x for x in summ["stages"] if x["name"] == stage.split("_")[1]][0]
    W = {c["cycle"]: c["W_loop"] for c in st["cycles"]}
    wp = np.mean(wp, axis=0)
    ms = cfg["Js"] / MU0
    s = a[a[:, 0] == cyc[-1]][:cfg["samples"]]        # the last row (t = T) is the first row of the next cycle
    dm0 = [float(np.subtract(*zero_branches(s[:, 2], s[:, 9 + p])) / ms) for p in range(4)]
    h_hat = 0.5 * (s[:, 2].max() - s[:, 2].min())
    q_ell = [float(np.pi * MU0 * h_hat * ms * d / (2.0 * w)) for d, w in zip(dm0, wp)]
    return {"cycles": cyc, "w_p": [float(v) for v in wp], "dM0_p": dm0, "q_ell_p": q_ell, "w_box": float(np.mean(wb)),
            "w_sum": float(cfg["phi_vox"] / 4.0 * wp.sum()),
            "W_loop": float(np.mean([W[c] for c in cyc])) if all(c in W for c in cyc) else None,
            "cos_axis_drive": [float(abs(u @ np.asarray(cfg.get("direction", [1.0, 0.0, 0.0]))))
                               for u in particle_axes(cfg.get("direction", [1.0, 0.0, 0.0]))]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--stages", nargs="+", default=["00_amp00", "03_amp03"])
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    root = pathlib.Path(a.runs)
    names = sorted((d.name for d in root.iterdir() if d.is_dir() and d.name.startswith("K") and
                    d.name.endswith("_sat") and (d / "DONE").exists()),
                   key=lambda n: (n[:2], int(n.split("_")[1])))
    out = {n: {st: run_record(root / n, st) for st in a.stages} for n in names}
    print("# Loop area of each particle, runs K<c>_<k>_sat\n")
    print("w_p: mu0 * closed integral of H dM_p, J/m³ of the particle volume, mean of the last 2 cycles of the stage. "
          "dM0_p: M_p at H = 0 falling minus rising, last cycle, over M_s (p 0 … 3). q_ell_p: dM0_p over the value "
          "of an elliptic loop with the same area (1: linear response with a phase lag). "
          "Check: phi_vox / 4 * sum w_p over the box loop area (mu0 * closed integral of H dM_par) and over W_loop "
          "(summary.json).\n")
    cosd = out[names[0]][a.stages[0]]["cos_axis_drive"]
    print("cos of the angle between the particle easy axis and the drive: %s\n" %
          ", ".join("p%d %.3f" % (p, c) for p, c in enumerate(cosd)))
    for st in a.stages:
        print("## %s\n" % LABEL.get(st, st))
        print("| run | w_0 | w_1 | w_2 | w_3 | dM0_p, 1e-4 | q_ell_p | sum check / box loop | sum check / W_loop |")
        print("|---|---|---|---|---|---|---|---|---|")
        for n in names:
            r = out[n][st]
            q = r["w_sum"] / r["W_loop"] if r["W_loop"] else float("nan")
            print("| %s | %s | %s | %s | %.4f | %.4f |" % (
                n, " | ".join("%.2f" % v for v in r["w_p"]), " ".join("%.2f" % (1e4 * v) for v in r["dM0_p"]),
                " ".join("%.2f" % v for v in r["q_ell_p"]), r["w_sum"] / r["w_box"], q))
        print()
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
