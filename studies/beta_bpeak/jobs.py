#!/usr/bin/env python
"""
Generate the job lists of the beta(B_peak) study and estimate the cost.

    python jobs.py                 # writes jobs/*.txt and prints the cost table
    python jobs.py --speedup 2.2   # cost for a GPU that is 2.2x faster than a V100

Job file format: one job per line, "<name> <run_loops.py arguments>".
The files are sorted by estimated cost (largest first) for load balance.

Cost model (CALIBRATE with the benchmark jobs, then update the defaults):
  time/step  = s_per_step_Mcell * n_cells / 1e6 / speedup
               (0.011 s per Mcell: V100 fp32, RKF45, 13.8 Mcell run of
                bench/fp32_validation; fp64 ~2x)
  dt         = dt_ps * (dx / 10.2569 nm)^2      (exchange-limited step)
               (dt_ps = 2.0: CPU smoke test gave 2.6 ps at 10 nm; conservative)
  cycles     = n_sat + n_amp * cycles_avg       (cycles_avg: 3 ... 4, adaptive)
"""
import argparse
import pathlib

DX = 10.2569e-9          # common cell size: a(d=1um, phi=0.65) / 144
DX7 = 0.73849e-6 / 105   # 7.03 nm, convergence test at d = 0.5 um

# base settings shared by all jobs (the defaults of run_loops.py are the
# material assumptions; they are repeated here so that job files are explicit)
BASE = dict(Js=1.80, A=15e-12, K1=30e3, lambda_s=5e-6, alpha=0.1, xi=50e-9,
            freq=30e6, n_amp=9, H_factor=0.5, cycles_per_amp=3, max_cycles_per_amp=4,
            samples=256, precision="single")


# Orientation sets of the 4 crystals, selected from seeds 1..20000 so that the
# statistics along the drive (x) are close to an isotropic powder:
# E_a/K1 = mean(a1^2 a2^2 + a2^2 a3^2 + a3^2 a1^2) = 0.20 (std 0.087),
# nearest-easy-axis cosine = 0.831 (std 0.100).  Seed 1 was a 2.7-sigma outlier
# (all 4 crystals with an easy axis within 22 deg of the field: E_a/K1 = 0.083).
AXES_SEEDS = {1: 16295, 2: 13903, 3: 6982}


def job(name, d, N, dx=DX, sigma=0.0, seed=1, **kw):
    a = dict(BASE)
    a.update(d=d, N=N, dx=dx, sigma_rms=sigma, seed_axes=AXES_SEEDS[seed], seed_stress=seed, seed_m=seed)
    a.update(kw)
    return name, a


def matrix():
    bench, val, val64, prod = [], [], [], []

    # --- benchmark: 1 cycle, timing only -------------------------------------
    bench.append(job("bench_N072", 0.5e-6, 72, sigma=150e6, max_cycles=1, relax_time=0.5e-9))
    bench.append(job("bench_N144", 1.0e-6, 144, sigma=150e6, max_cycles=1, relax_time=0.5e-9))

    # --- production -----------------------------------------------------------
    S = [0.0, 150e6, 400e6]
    for d, N in [(0.5e-6, 72), (0.75e-6, 108), (1.0e-6, 144)]:            # size series, phi = 0.65
        for s in S:
            kw = {"keep_stage_ckpt": True} if N == 72 else {}               # fp32 twins of the fp64 checks
            prod.append(job("P_d%03d_N%03d_s%03d" % (d * 1e8, N, s / 1e6), d, N, sigma=s, **kw))
    for N in (160, 140):                                                   # phi = 0.47 / 0.71 at d = 1 um
        for s in S:
            prod.append(job("P_d100_N%03d_s%03d" % (N, s / 1e6), 1.0e-6, N, sigma=s))
    for seed in (2, 3):                                                    # more realizations
        for d, N, s in [(1.0e-6, 144, 150e6), (1.0e-6, 144, 400e6), (0.75e-6, 108, 150e6)]:
            prod.append(job("P_d%03d_N%03d_s%03d_seed%d" % (d * 1e8, N, s / 1e6, seed), d, N, sigma=s, seed=seed))
    prod.append(job("P_d100_N144_s150_K20", 1.0e-6, 144, sigma=150e6, K1=20e3))   # K1 sensitivity

    # --- validation (fp32, any GPU) --------------------------------------------
    val.append(job("V_dx7_d050_s000", 0.5e-6, 105, dx=DX7, sigma=0.0))           # grid pinning
    val.append(job("V_dx7_d050_s150", 0.5e-6, 105, dx=DX7, sigma=150e6))
    val.append(job("V_f15_d050_N072_s000", 0.5e-6, 72, sigma=0.0, freq=15e6))    # quasi-static check
    val.append(job("V_f15_d050_N072_s150", 0.5e-6, 72, sigma=150e6, freq=15e6))
    val.append(job("V_f15_d100_N144_s150", 1.0e-6, 144, sigma=150e6, freq=15e6))
    val.append(job("V_asc_d050_N072_s150", 0.5e-6, 72, sigma=150e6, order="ascending"))  # protocol

    # --- validation fp64 (ONLY on V100 / A100 / H100: consumer GPUs have 1/64 FP64) ---
    val64.append(job("V_fp64_d050_N072_s000", 0.5e-6, 72, sigma=0.0, precision="double"))
    val64.append(job("V_fp64_d050_N072_s150", 0.5e-6, 72, sigma=150e6, precision="double"))

    return {"bench": bench, "validation": val, "validation_fp64": val64, "production": prod}


def to_cli(a):
    parts = []
    for k, v in a.items():
        if v is True:
            parts.append("--%s" % k)
        else:
            parts.append("--%s %s" % (k, repr(v) if isinstance(v, float) else v))
    return " ".join(parts)


def cost_h(a, s_per_step_Mcell, dt_ps, cycles_avg, speedup):
    n = a["N"] ** 3
    dt = dt_ps * 1e-12 * (a["dx"] / DX) ** 2
    t_step = s_per_step_Mcell * n / 1e6 / speedup * (2.0 if a["precision"] == "double" else 1.0)
    period = 1.0 / a["freq"]
    if "max_cycles" in a:
        cycles = a["max_cycles"]
    else:
        n_amp = a["n_amp"]
        cycles = 1 + n_amp * cycles_avg + (n_amp if a.get("order") == "ascending" else 0)
    return cycles * period / dt * t_step / 3600.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s_per_step_Mcell", type=float, default=0.011, help="V100 fp32 seconds per RKF45 step per Mcell")
    ap.add_argument("--dt_ps", type=float, default=2.0, help="mean RKF45 step at dx = 10.26 nm [ps]")
    ap.add_argument("--cycles_avg", type=float, default=3.3, help="mean cycles per amplitude (3 ... 4)")
    ap.add_argument("--speedup", type=float, default=1.0, help="GPU speed relative to a V100 (fp32)")
    ap.add_argument("--outdir", default=str(pathlib.Path(__file__).resolve().parent / "jobs"))
    a = ap.parse_args()

    outdir = pathlib.Path(a.outdir)
    outdir.mkdir(exist_ok=True)
    total = 0.0
    for group, jobs in matrix().items():
        costed = sorted(((cost_h(j, a.s_per_step_Mcell, a.dt_ps, a.cycles_avg, a.speedup), n, j) for n, j in jobs),
                        key=lambda x: -x[0])
        with open(outdir / ("%s.txt" % group), "w") as f:
            for c, n, j in costed:
                f.write("%s %s\n" % (n, to_cli(j)))
        g = sum(c for c, _, _ in costed)
        total += g
        print("\n%-16s %2d jobs  %7.1f GPU-h" % (group, len(costed), g))
        for c, n, j in costed:
            print("   %-28s N=%3d  %8.2f h" % (n, j["N"], c))
    print("\nTOTAL %.1f GPU-h (speedup %.2g vs V100 fp32)" % (total, a.speedup))


if __name__ == "__main__":
    main()
