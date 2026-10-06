#!/usr/bin/env python
"""
Sobol design v2 (PLAN 4.0, PROTOKOLL 7.13 ... 7.29). Writes the job files and the design table:
  jobs/sobol_main.txt     96 Sobol points + 8 centre replicates, alpha 0.1, 16-phase snapshots at 9/50/150 mT
  jobs/sobol_corners.txt  4 corners x 3 seeds, alpha 0.1, snapshots as main
  jobs/sobol_a001.txt     the 8 centre replicates at alpha 0.01, 9 cycles, snapshots as main (test |m1|^2 ~ 1/alpha),
                          started from the init.pt of the alpha 0.1 run of the same seed (Chris 2026-10-06:
                          the alpha 0.01 Sobol subset is dropped, --n_a001 0; a later own alpha study can use it)
  results/sobol_design.csv  point, seed, factors (nominal and realised d/l_ex, phi), checks, cost
The design is extensible: the mapping of a Sobol point does not depend on N (L_eff level = floor(4 u)),
seed = 1001 + point index. Every point is checked with the geometry of the study (contact rule over 26
neighbours at dx 3, local wall width >= 2 cells, d/L_eff >= 7).

    python sobol_design.py [--n 96]
"""
import argparse
import csv
import math
import pathlib

import numpy as np
from scipy.stats import qmc

import jobs
from anisotropy_model import box_edge
from geometry import fcc_box

HERE = pathlib.Path(__file__).resolve().parent
DESIGN_SEED = 20261006
L_LEVELS = (12.0, 18.0, 24.0, 30.0)
D_RANGE = (212.0, 300.0)
RP_RANGE = (0.0, 3.0)
PHI_RANGE = (0.55, 0.69)
CENTRE = dict(Leff_lex=18.0, r_p=1.5, phi=0.62, d_lex=252.0)
CORNERS = [("K1", dict(d_lex=212.0, Leff_lex=30.0, r_p=1.5, phi=0.62)),
           ("K2", dict(d_lex=300.0, Leff_lex=12.0, r_p=3.0, phi=0.62)),
           ("K3", dict(d_lex=212.0, Leff_lex=18.0, r_p=1.5, phi=0.69)),
           ("K4", dict(d_lex=300.0, Leff_lex=18.0, r_p=1.5, phi=0.55))]
SNAP = dict(snap_phases=16, snap_cycles=2, snap_mT=[9.0, 50.0, 150.0])   # 2 cycles: periodic vs non-periodic part (audit R1)
ACC = dict(acc=True)    # per-cell accumulators at all 7 amplitudes, all runs (audit 3 item 9b, Chris 2026-10-06)
LEX = 3.342e-9


def template():
    """The baseline job (PROTOKOLL 7.17) as the template: same protocol, relaxation and amplitudes."""
    name, a = jobs.matrix()["baseline"][0]
    a = dict(a)
    for k in ("snap_phases", "snap_cycles", "snap_mT", "seed"):
        a.pop(k, None)
    return a


def point(u):
    return dict(Leff_lex=L_LEVELS[min(int(math.floor(4 * u[0])), 3)],
                r_p=round(float(RP_RANGE[0] + (RP_RANGE[1] - RP_RANGE[0]) * u[1]), 6),
                phi=round(float(PHI_RANGE[0] + (PHI_RANGE[1] - PHI_RANGE[0]) * u[2]), 6),
                d_lex=round(float(math.exp(math.log(D_RANGE[0]) + (math.log(D_RANGE[1]) - math.log(D_RANGE[0])) * u[3])), 4))


def check(f, dx=3.0):
    a_lex, d_used = box_edge(f["d_lex"], f["phi"])
    N = int(round(a_lex / dx))
    try:
        g = fcc_box(N, d_used * LEX, dx * LEX)
        contact = "ok"
        gap = g["gap_nom"] / (dx * LEX)
    except ValueError as e:
        contact, gap = "FAIL: %s" % e, float("nan")
    Q = 1.0 / f["Leff_lex"] ** 2
    wall = 1.0 / math.sqrt(Q * (1.0 + f["r_p"])) / dx
    return {"d_used": d_used, "a_lex": a_lex, "N": N, "contact": contact, "gap_cells": gap,
            "wall_cells": wall, "d_over_L": d_used / f["Leff_lex"],
            "ok": contact == "ok" and wall >= 2.0 - 1e-9 and d_used / f["Leff_lex"] >= 7.0 - 1e-9}


def cost(a):
    return jobs.cost_h(a, 0.012, 0.75, 3.5, 1.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=96)
    ap.add_argument("--n_a001", type=int, default=0)
    a = ap.parse_args()
    T = template()
    U = qmc.Sobol(d=4, scramble=True, seed=DESIGN_SEED).random(a.n)
    rows, main_jobs, a001_jobs, corner_jobs = [], [], [], []

    def add(name, f, seed, kind, alpha=0.1, snap=True, parent=None):
        j = dict(T, **f, seed=seed, alpha=alpha, **ACC)
        if snap:
            j.update(SNAP)
        if alpha != 0.1:
            # init_wait_h 0 (audit R3): the chain starts the subset after the main queue; a missing parent init.pt
            # makes the job fail at once instead of blocking a GPU
            j.update(cycles_per_amp=9, max_cycles_per_amp=9, init_from="%s/init.pt" % parent, init_wait_h=0.0)
        c = check(f)
        if not c["ok"]:
            raise SystemExit("design point %s fails the checks: %s" % (name, c))
        rows.append(dict(name=name, kind=kind, seed=seed, alpha=alpha, **{k: round(v, 5) for k, v in f.items()},
                         d_used=round(c["d_used"], 3), N=c["N"], gap_cells=round(c["gap_cells"], 2),
                         wall_cells=round(c["wall_cells"], 2), d_over_L=round(c["d_over_L"], 2),
                         cost_V100h=round(cost(j), 2)))
        return name, j

    for i, u in enumerate(U):
        f = point(u)
        main_jobs.append(add("S%03d" % i, f, 1001 + i, "sobol"))
        if i < a.n_a001:
            a001_jobs.append(add("S%03d_a001" % i, f, 1001 + i, "sobol alpha 0.01", alpha=0.01, snap=False,
                                 parent="S%03d" % i))
    for r in range(8):
        main_jobs.append(add("C%d" % r, CENTRE, 1901 + r, "centre"))
        a001_jobs.append(add("C%d_a001" % r, CENTRE, 1901 + r, "centre alpha 0.01", alpha=0.01, snap=True,
                             parent="C%d" % r))
    for tag, f in CORNERS:
        for r in range(3):
            corner_jobs.append(add("%s_%d" % (tag, r), f, 2001 + 10 * int(tag[1]) + r, "corner %s" % tag))
    out = HERE / "jobs"
    for fname, js in (("sobol_main.txt", main_jobs), ("sobol_corners.txt", corner_jobs), ("sobol_a001.txt", a001_jobs)):
        with open(out / fname, "w") as fh:
            for name, j in js:
                fh.write("%s %s\n" % (name, jobs.to_cli(j)))
    with open(HERE / "results" / "sobol_design.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    tot = {"main alpha 0.1 (96 + 8)": sum(r["cost_V100h"] for r in rows if r["alpha"] == 0.1 and not r["kind"].startswith("corner")),
           "corners": sum(r["cost_V100h"] for r in rows if r["kind"].startswith("corner")),
           "alpha 0.01 (%d + 8, model; overestimates alpha 0.01 ~2x, PROTOKOLL)" % a.n_a001: sum(r["cost_V100h"] for r in rows if r["alpha"] != 0.1)}
    lv = np.bincount([L_LEVELS.index(r["Leff_lex"]) for r in rows if r["kind"] == "sobol"], minlength=4)
    print("main %d + centre 8 + corners %d + alpha 0.01 %d jobs" % (a.n, len(corner_jobs), len(a001_jobs)))
    print("points per L_eff level (main): %s" % lv.tolist())
    print("ranges realised: d_used %.1f … %.1f, r_p %.2f … %.2f, phi %.3f … %.3f; min gap %.2f cells, min wall %.2f cells, "
          "min d/L %.2f" % (min(r["d_used"] for r in rows), max(r["d_used"] for r in rows),
                            min(r["r_p"] for r in rows), max(r["r_p"] for r in rows),
                            min(r["phi"] for r in rows), max(r["phi"] for r in rows),
                            min(r["gap_cells"] for r in rows), min(r["wall_cells"] for r in rows),
                            min(r["d_over_L"] for r in rows)))
    print("cost model (V100-h): " + ", ".join("%s %.0f" % kv for kv in tot.items()) + ", total %.0f" % sum(tot.values()))


if __name__ == "__main__":
    main()
