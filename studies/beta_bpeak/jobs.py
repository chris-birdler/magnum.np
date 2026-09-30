#!/usr/bin/env python
"""
Job lists of the beta(B_peak) study (fictitious nanocrystalline powder,
reduced units) and a cost estimate.

    python jobs.py                  # writes jobs/*.txt, prints the cost table
    python jobs.py --speedup 2.2    # cost for a GPU 2.2x faster than a V100

Job file format: one job per line, "<name> <run_loops.py arguments>", sorted
by estimated cost (largest first).

Groups:
  bench     1 cycle of a d/l_ex = 300 run (timing -> calibrate the cost model)
  meshtest  d/l_ex = 96, dx = 3 / 2 / 1.5 l_ex, floor (Q_eff = 0) and L_eff = 12:
            does the grid pin the unresolved vortex cores? (ran 2026-09-30: cycles
            not steady at b <= 0.25 within 3-4 cycles -> steady test)
  steady    same geometry, 2 amplitudes (b ~ 0.25 and 0.06), 12 cycles each,
            dx = 3 vs 1.5, floor and L_eff = 12, alpha = 0.02 and 0.1 (L12):
            mesh test on steady cycles. Runs BEFORE the pilot.
  numfloor  same geometry, dx = 3, L_eff = 12, B_peak = 150 / 50 / 10 mT, 8 cycles:
            numerical floor of the small loops (atol 1e-5 vs 1e-6, fp32 vs fp64).

Amplitudes: the study range is B_peak = 10 ... 150 mT (Chris, 2026-09-30):
b = B_peak/Js with Js = 1.5 T. run_loops.py --b_list sets the drive per
stage from the measured b/h of the stage before.
  pilot   d/l_ex = 150, 4 amplitudes: pinning signal against the floor, mesh,
          damping, drive direction, frequency, realisation scatter. The
          production matrix is decided after the pilot (HANDOVER.md).

Cost model (calibrated on V100, 2026-09-30: bench_base and meshtest):
  steps/cycle = 2 pi / (f_rel * dt_tau),  dt_tau = dt * gamma Ms
                dt_tau = 0.75 at dx = 3 l_ex, scales with (dx/3)^0.5
                (measured: 0.76 at dx = 3, d/l_ex = 300; 0.56 at dx = 1.5)
  s/step      = max(s_min_step, s_per_step_Mcell * n_cells / 1e6) / speedup
                (0.012 s/Mcell on large grids; small grids are limited by the
                overhead s_min_step = 0.0075 s)
  cycles      = 1 + n_stages * cycles_per_amp   (fixed number of cycles), else
                1 + n_stages * cycles_avg
                cycles_avg = 3.5: with 10 % scatter per cycle (target, dW_tol = 0.10)
                about half of the stages need the 4th cycle
The unit "V100-h" is only a cost unit. Any fp32 GPU can run the jobs.
"""
import argparse
import math
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from anisotropy_model import box_edge   # noqa: E402
from units import F_REL_DEFAULT         # noqa: E402

JS_T = 1.5                                         # Js of the unit system (units.py)


def b_of_mT(*mT):
    return [round(x / 1000.0 / JS_T, 6) for x in mT]


B_STUDY = b_of_mT(150, 100, 70, 50, 35, 20, 9)     # production amplitudes (9 mT: 10 mT inside the range)
B_PILOT = b_of_mT(150, 100, 50, 25, 9)

BASE = dict(d_lex=300.0, dx_lex=3.0, phi=0.65, Leff_lex=18.0, r_p=0.0,
            alpha=0.02, b_list=B_STUDY, n_amp=len(B_STUDY), cycles_per_amp=7, max_cycles_per_amp=7,
            samples=256, precision="single", seed=1)
PILOT = dict(BASE, d_lex=150.0, b_list=B_PILOT, n_amp=len(B_PILOT))
# frozen protocol of the mesh test and the steady test (they ran with it)
_OLD = dict(d_lex=150.0, dx_lex=3.0, phi=0.65, Leff_lex=18.0, r_p=0.0,
            alpha=0.02, n_amp=4, h_factor=0.25, cycles_per_amp=3, max_cycles_per_amp=4, dW_tol=0.10,
            samples=256, precision="single", seed=1)
MESHTEST = dict(_OLD, d_lex=96.0)       # a = 144 l_ex: a, a/2 and L_eff on cell faces for dx = 3, 2, 1.5
# steady test: b ~ 0.25 and 0.06 (from the meshtest), 12 cycles each, no early stop
STEADY = dict(MESHTEST, h_list=[0.02917, 0.007292], n_amp=2, cycles_per_amp=12, max_cycles_per_amp=12)


def job(base, name, **kw):
    a = dict(base)
    a.update(kw)
    return name, a


def matrix():
    bench = [job(BASE, "bench_base", max_cycles=1, relax_tau=100.0)]

    meshtest = []
    for dx in (3.0, 2.0, 1.5):
        tag = ("%g" % dx).replace(".", "")
        meshtest.append(job(MESHTEST, "M_floor_dx%s" % tag, Q_eff=0.0, dx_lex=dx))
        meshtest.append(job(MESHTEST, "M_L12_dx%s" % tag, Leff_lex=12.0, dx_lex=dx))

    steady = []
    for dx in (3.0, 1.5):
        tag = ("%g" % dx).replace(".", "")
        steady.append(job(STEADY, "S_floor_dx%s" % tag, Q_eff=0.0, dx_lex=dx))
        steady.append(job(STEADY, "S_L12_dx%s" % tag, Leff_lex=12.0, dx_lex=dx))
        steady.append(job(STEADY, "S_L12_a010_dx%s" % tag, Leff_lex=12.0, dx_lex=dx, alpha=0.1))

    numfloor = []
    FLOOR = dict(MESHTEST, Leff_lex=12.0, b_list=b_of_mT(150, 50, 10), n_amp=3,
                 cycles_per_amp=8, max_cycles_per_amp=8)
    del FLOOR["h_factor"], FLOOR["dW_tol"]
    for prec in ("single", "double"):
        for atol in (1e-5, 1e-6):
            numfloor.append(job(FLOOR, "F_%s_atol%.0e" % (prec[:2], atol), precision=prec, atol=atol))

    pilot = []
    for L in (12.0, 18.0, 24.0, 30.0):                              # Q_eff = 6.9e-3 ... 1.1e-3
        pilot.append(job(PILOT, "T_L%02.0f" % L, Leff_lex=L))
    pilot.append(job(PILOT, "T_floor", Q_eff=0.0))                   # no anisotropy: floor
    for rp in (1.0, 3.0):                                           # particle-scale stress
        pilot.append(job(PILOT, "T_L18_rp%.0f" % rp, r_p=rp))
    pilot.append(job(PILOT, "T_L18_alpha010", alpha=0.1))           # damping loss
    pilot.append(job(PILOT, "T_L18_d111", direction=[1.0, 1.0, 1.0]))  # cube lattice vs drive
    pilot.append(job(PILOT, "T_L18_f2", f_rel=2.0 * F_REL_DEFAULT))  # dynamic share (pair with T_L18)
    pilot.append(job(PILOT, "T_L30_dx2", Leff_lex=30.0, dx_lex=2.0))  # mesh, same cubes as T_L30
    pilot.append(job(PILOT, "T_L30_seed2", Leff_lex=30.0, seed=2))    # realisation scatter
    return {"bench": bench, "meshtest": meshtest, "steady": steady, "numfloor": numfloor, "pilot": pilot}


def n_cells(a):
    a_lex, _ = box_edge(a["d_lex"], a["phi"])
    N = int(round(a_lex / a["dx_lex"]))
    return N, N**3


def cost_h(a, s_per_step_Mcell, dt_tau, cycles_avg, speedup, s_min_step=0.0075):
    N, n = n_cells(a)
    f_rel = a.get("f_rel", F_REL_DEFAULT)
    steps = 2.0 * math.pi / (f_rel * dt_tau * (a["dx_lex"] / 3.0) ** 0.5)
    steps *= (1e-5 / a.get("atol", 1e-5)) ** 0.2                      # RKF45: dt ~ atol^(1/5)
    t_step = max(s_min_step, s_per_step_Mcell * n / 1e6) / speedup * (2.0 if a["precision"] == "double" else 1.0)
    n_st = len(a.get("h_list") or a.get("b_list") or []) or a["n_amp"]
    per_amp = a["cycles_per_amp"] if a["cycles_per_amp"] == a["max_cycles_per_amp"] else cycles_avg
    cycles = a.get("max_cycles", 1 + n_st * per_amp)
    return cycles * steps * t_step / 3600.0


def to_cli(a):
    parts = []
    for k, v in a.items():
        if v is True:
            parts.append("--%s" % k)
        elif isinstance(v, (list, tuple)):
            parts.append("--%s %s" % (k, " ".join(repr(float(x)) for x in v)))
        else:
            parts.append("--%s %s" % (k, repr(v) if isinstance(v, float) else v))
    return " ".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s_per_step_Mcell", type=float, default=0.012)
    ap.add_argument("--dt_tau", type=float, default=0.75)
    ap.add_argument("--cycles_avg", type=float, default=3.5)
    ap.add_argument("--speedup", type=float, default=1.0)
    ap.add_argument("--outdir", default=str(HERE / "jobs"))
    a = ap.parse_args()

    outdir = pathlib.Path(a.outdir)
    outdir.mkdir(exist_ok=True)
    for f in outdir.glob("*.txt"):
        f.unlink()
    total = 0.0
    for group, jobs in matrix().items():
        costed = sorted(((cost_h(j, a.s_per_step_Mcell, a.dt_tau, a.cycles_avg, a.speedup), n, j) for n, j in jobs),
                        key=lambda x: -x[0])
        with open(outdir / ("%s.txt" % group), "w") as f:
            for c, n, j in costed:
                f.write("%s %s\n" % (n, to_cli(j)))
        g = sum(c for c, _, _ in costed)
        total += g
        print("\n%-11s %2d jobs  %6.1f V100-h" % (group, len(costed), g))
        for c, n, j in costed:
            print("   %-20s N=%3d  %6.2f h" % (n, n_cells(j)[0], c))
    print("\nTOTAL %.1f V100-h (speedup %.2g)" % (total, a.speedup))


if __name__ == "__main__":
    main()
