#!/usr/bin/env python
"""
Job lists of the beta(B_peak) study (fictitious material, reduced units) and
a cost estimate.

    python jobs.py                  # writes jobs/*.txt, prints the cost table
    python jobs.py --speedup 2.2    # cost for a GPU 2.2x faster than a V100

Job file format: one job per line, "<name> <run_loops.py arguments>", sorted
by estimated cost (largest first).

Cost model (CALIBRATE with jobs/bench.txt, then change the defaults):
  steps/cycle = 2 pi / (f_rel * dt_tau),  dt_tau = dt * gamma Ms
                dt_tau = 0.6 at dx = 3 l_ex (CPU test: 0.69), scales with (dx/l_ex)^2
  s/step      = s_per_step_Mcell * n_cells / 1e6 / speedup
                (0.011: V100 fp32, RKF45, from bench/fp32_validation)
  cycles      = 1 + n_amp * cycles_avg
The unit "V100-h" is only a cost unit. Any fp32 GPU can run the jobs.
"""
import argparse
import math
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from geometry import a_from          # noqa: E402
from units import F_REL_DEFAULT      # noqa: E402

BASE = dict(d_lex=300.0, dx_lex=3.0, phi=0.65, Q=0.02, aniso="uniaxial", kappa=0.008, xi_lex=15.0,
            alpha=0.1, n_amp=9, h_factor=0.5, cycles_per_amp=3, max_cycles_per_amp=4, samples=256,
            precision="single", seed=1)


def job(name, **kw):
    a = dict(BASE)
    a.update(kw)
    return name, a


def matrix():
    bench = [job("bench_base", max_cycles=1, relax_tau=100.0)]

    prod = []
    for k in (0.0, 0.0025, 0.008, 0.025):                       # pinning strength (base: 0.008)
        prod.append(job("P_kappa%05.0f" % (k * 1e5), kappa=k))
    for dl in (150.0, 225.0):                                   # particle size
        prod.append(job("P_d%03.0f" % dl, d_lex=dl))
    for ph in (0.47, 0.71):                                     # fill factor
        prod.append(job("P_phi%03.0f" % (ph * 100), phi=ph))
    prod.append(job("P_Q010", Q=0.01))                          # anisotropy
    for s in (2, 3):                                            # realisations (stress field, initial m)
        prod.append(job("P_seed%d" % s, seed=s))

    ctrl = []
    ctrl.append(job("C_d150_kappa0", d_lex=150.0, kappa=0.0))                   # floor, dx = 3 l_ex
    ctrl.append(job("C_d150_kappa0_dx2", d_lex=150.0, kappa=0.0, dx_lex=2.0))   # floor, dx = 2 l_ex
    ctrl.append(job("C_d150_dx2", d_lex=150.0, dx_lex=2.0))                     # signal, dx = 2 l_ex
    ctrl.append(job("C_fhalf", f_rel=F_REL_DEFAULT / 2.0))                      # quasi-static check
    return {"bench": bench, "production": prod, "control": ctrl}


def n_cells(a):
    N = int(round(a_from(a["d_lex"], a["phi"]) / a["dx_lex"]))
    return N, N**3


def cost_h(a, s_per_step_Mcell, dt_tau, cycles_avg, speedup):
    N, n = n_cells(a)
    f_rel = a.get("f_rel", F_REL_DEFAULT)
    steps = 2.0 * math.pi / (f_rel * dt_tau * (a["dx_lex"] / 3.0) ** 2)
    t_step = s_per_step_Mcell * n / 1e6 / speedup * (2.0 if a["precision"] == "double" else 1.0)
    cycles = a.get("max_cycles", 1 + a["n_amp"] * cycles_avg)
    return cycles * steps * t_step / 3600.0


def to_cli(a):
    parts = []
    for k, v in a.items():
        if v is True:
            parts.append("--%s" % k)
        else:
            parts.append("--%s %s" % (k, repr(v) if isinstance(v, float) else v))
    return " ".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s_per_step_Mcell", type=float, default=0.011)
    ap.add_argument("--dt_tau", type=float, default=0.6)
    ap.add_argument("--cycles_avg", type=float, default=3.3)
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
