#!/usr/bin/env python
"""
AC BH-loops of a periodic FCC powder model for the beta(B_peak) study.

One run = one parameter set. The run does:
  1. random initial magnetization, short relaxation at H = 0,
  2. one (or more) saturating cycle(s)  -> resets the domain topology,
  3. descending drive amplitudes (geometric series). Each amplitude runs
     at least cycles_per_amp cycles, then more cycles (up to
     max_cycles_per_amp) until the cycle is steady (closure_tol, dW_tol).
     The descending series is at the same time an AC demagnetization.
     The first cycle of each amplitude is a transient; analyze.py averages
     the other cycles.
     (--order ascending: saturating cycle, AC demag down to the smallest
     amplitude with 1 cycle each, then measure with ascending amplitudes.)

For each cycle the script writes (see summary.json):
  B_peak   half peak-to-peak of the macroscopic flux density  [T]
  W_loop   loop energy  oint H dB = mu0 oint H dM_x           [J/m^3 per cycle]
  W_dis    LLG dissipation  int <p_dis> dt over the cycle     [J/m^3 per cycle]
  closure  |M(end) - M(start)| / M_peak   (M along the drive)  [-]
  dW_rel   |W_loop(c) - W_loop(c-1)| / W_loop(c)              [-]
All energies are per unit volume of the CORE (voids included).

Energy balance check: for a closed cycle W_loop == W_dis (the stray-field,
exchange and anisotropy energies are state functions). A difference shows
insufficient sampling, an open cycle or integrator errors.

Precision: set_precision() is called before the mesh is created. state.t is
kept in float64 on purpose (fp32 time would give ~1 % errors in the time
step accumulation, see HANDOVER.md). All diagnostic averages and the loop
integrals are accumulated in float64.

Example (benchmark, 1 cycle):
    CUDA_DEVICE=0 python run_loops.py --d 0.5e-6 --N 72 --dx 10.2569e-9 \
        --sigma_rms 150e6 --max_cycles 1 --out runs/bench_N72
"""
import argparse
import json
import math
import os
import pathlib
import subprocess
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

MU0 = 1.2566370614e-6
_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy >= 2.0 renamed trapz


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_argument_group("geometry")
    g.add_argument("--d", type=float, required=True, help="particle diameter [m]")
    g.add_argument("--dx", type=float, required=True, help="cell size [m]")
    g.add_argument("--N", type=int, default=None, help="cells per edge (a = N*dx); alternative: --phi")
    g.add_argument("--phi", type=float, default=None, help="target packing fraction (N = round(a/dx))")

    g = p.add_argument_group("material (defaults: FeSiCr assumption, see HANDOVER.md)")
    g.add_argument("--Js", type=float, default=1.80, help="mu0*Ms [T]")
    g.add_argument("--A", type=float, default=15e-12, help="exchange stiffness [J/m]")
    g.add_argument("--K1", type=float, default=30e3, help="cubic anisotropy K1 [J/m^3]")
    g.add_argument("--lambda_s", type=float, default=5e-6, help="isotropic saturation magnetostriction [-]")
    g.add_argument("--alpha", type=float, default=0.1, help="Gilbert damping (numerical choice, see HANDOVER.md)")
    g.add_argument("--Js_void", type=float, default=1e-8, help="mu0*Ms in the void [T] (must be > 0)")

    g = p.add_argument_group("stress field")
    g.add_argument("--sigma_rms", type=float, default=0.0, help="rms of each stress component [Pa]")
    g.add_argument("--xi", type=float, default=50e-9, help="stress correlation length [m]")

    g = p.add_argument_group("seeds")
    g.add_argument("--seed_axes", type=int, default=1)
    g.add_argument("--seed_stress", type=int, default=1)
    g.add_argument("--seed_m", type=int, default=1)

    g = p.add_argument_group("drive / protocol")
    g.add_argument("--freq", type=float, default=30e6, help="drive frequency [Hz]")
    g.add_argument("--direction", type=float, nargs=3, default=[1.0, 0.0, 0.0])
    g.add_argument("--H_sat", type=float, default=None, help="saturating amplitude [A/m] (default (1-phi)*Ms)")
    g.add_argument("--n_sat", type=int, default=1, help="number of saturating cycles")
    g.add_argument("--H_max", type=float, default=None, help="largest measured amplitude [A/m] (default (1-phi)*Ms/3)")
    g.add_argument("--H_factor", type=float, default=0.5, help="ratio of neighbouring amplitudes")
    g.add_argument("--n_amp", type=int, default=9, help="number of measured amplitudes")
    g.add_argument("--H_list", type=float, nargs="*", default=None, help="explicit amplitudes [A/m] (overrides H_max/H_factor/n_amp)")
    g.add_argument("--cycles_per_amp", type=int, default=3, help="minimum cycles per measured amplitude (1st = transient)")
    g.add_argument("--max_cycles_per_amp", type=int, default=4, help="maximum cycles per measured amplitude")
    g.add_argument("--closure_tol", type=float, default=0.02, help="steady if |M_end-M_start|/M_peak < tol ...")
    g.add_argument("--dW_tol", type=float, default=0.03, help="... and |W_c - W_(c-1)|/W_c < tol")
    g.add_argument("--order", choices=["descending", "ascending"], default="descending")
    g.add_argument("--samples", type=int, default=256, help="samples per cycle")
    g.add_argument("--relax_time", type=float, default=2e-9, help="relaxation at H=0 before the protocol [s]")

    g = p.add_argument_group("numerics")
    g.add_argument("--precision", choices=["single", "double"], default="single")
    g.add_argument("--atol", type=float, default=1e-5, help="RKF45 tolerance")

    g = p.add_argument_group("run control / output")
    g.add_argument("--out", type=str, required=True, help="output directory")
    g.add_argument("--max_cycles", type=int, default=None, help="stop after this many cycles (benchmark)")
    g.add_argument("--vti", action="store_true", help="write m as .vti after each stage")
    g.add_argument("--keep_stage_ckpt", action="store_true",
                   help="keep m after every stage as ckpt_<k>_<stage>.pt (for --init_from)")
    g.add_argument("--no_resume", action="store_true", help="ignore an existing checkpoint")
    g.add_argument("--init_from", type=str, default=None,
                   help="start from m of this checkpoint.pt (no relaxation); e.g. with --n_sat 0 "
                        "--H_list <small amplitudes> --precision double to repeat the smallest amplitudes in fp64")
    g.add_argument("--timer", action="store_true", help="print the magnum.np timer report at the end")
    return p.parse_args(argv)


def git_commit():
    try:
        return subprocess.check_output(["git", "-C", str(HERE), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def build_stages(args, Ms, phi):
    H_sat = args.H_sat if args.H_sat is not None else (1.0 - phi) * Ms
    if args.H_list:
        amps = sorted(args.H_list, reverse=True)
    else:
        H_max = args.H_max if args.H_max is not None else (1.0 - phi) * Ms / 3.0
        amps = [H_max * args.H_factor**k for k in range(args.n_amp)]
    stages = [{"name": "sat", "H_amp": H_sat, "n_cycles": args.n_sat, "measure": False}]
    if args.order == "descending":
        stages += [{"name": "amp%02d" % k, "H_amp": H, "n_cycles": args.cycles_per_amp, "measure": True}
                   for k, H in enumerate(amps)]
    else:
        stages += [{"name": "demag%02d" % k, "H_amp": H, "n_cycles": 1, "measure": False}
                   for k, H in enumerate(amps)]
        stages += [{"name": "amp%02d" % k, "H_amp": H, "n_cycles": args.cycles_per_amp, "measure": True}
                   for k, H in reversed(list(enumerate(amps)))]
    return stages


def main(argv=None):
    args = parse_args(argv)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # --- magnum.np import and precision (before any Mesh / State) ---------
    from magnumnp import (set_precision, set_log_level, Mesh, State, constants, Timer,
                          DemagFieldPBC, ExchangeField, CubicAnisotropyField2, LLGSolver,
                          accumulate_h, write_vti)
    import torch
    set_log_level(25)
    set_precision(args.precision)
    if args.timer:
        Timer.enable()

    from geometry import fcc_box, a_from, random_rotations
    from stress import random_stress_field
    from field_terms_extra import StressAnisotropyField, SinusoidalDrive

    # --- geometry ----------------------------------------------------------
    if (args.N is None) == (args.phi is None):
        raise SystemExit("Give exactly one of --N or --phi.")
    N = args.N if args.N is not None else int(round(a_from(args.d, args.phi) / args.dx))
    geo = fcc_box(N, args.d, args.dx)
    ids = geo["ids"]
    mag = ids >= 0
    phi = geo["phi_vox"]

    Ms = args.Js / MU0
    Ms_void = args.Js_void / MU0

    mesh = Mesh((N, N, N), (args.dx, args.dx, args.dx), pbc=(1, 1, 1))
    state = State(mesh)
    dt_ = state.dtype
    dev = state.device

    def T(a):  # numpy -> tensor with run dtype on the run device
        return torch.tensor(a, dtype=dt_, device=dev)

    mag4 = mag[..., None]
    state.material = {"alpha": args.alpha}
    state.material["Ms"] = T(np.where(mag4, Ms, Ms_void))
    state.material["A"] = T(np.where(mag4, args.A, 0.0))
    state.material["Kc1"] = T(np.where(mag4, args.K1, 0.0))
    state.material["Kc2"] = 0.0
    state.material["Kc3"] = 0.0

    R = random_rotations(4, args.seed_axes)          # [4,3,3], rows = cubic axes
    ax1 = np.zeros((N, N, N, 3)); ax1[...] = [1.0, 0.0, 0.0]
    ax2 = np.zeros((N, N, N, 3)); ax2[...] = [0.0, 1.0, 0.0]
    for p in range(4):
        ax1[ids == p] = R[p, 0]
        ax2[ids == p] = R[p, 1]
    state.material["Kc_axis1"] = T(ax1)
    state.material["Kc_axis2"] = T(ax2)

    # --- field terms -------------------------------------------------------
    terms = [DemagFieldPBC(), ExchangeField(), CubicAnisotropyField2()]
    sig_rms_real = 0.0
    if args.sigma_rms > 0.0 and args.lambda_s != 0.0:
        sig = random_stress_field(N, args.dx, args.xi, args.sigma_rms, args.seed_stress, mask=mag)
        sig[~mag] = 0.0
        sig_rms_real = float(sig[mag].std())
        terms.append(StressAnisotropyField(T(sig), args.lambda_s))
        del sig
    drive = SinusoidalDrive(state, tuple(args.direction))
    terms.append(drive)
    llg = LLGSolver(terms, atol=args.atol)

    # --- diagnostics helpers ----------------------------------------------
    Ms_t = state.material["Ms"]                                  # [N,N,N,1]
    ids_t = torch.tensor(ids.astype(np.int64) + 1, device=dev).reshape(-1)   # 0 = void
    counts = torch.bincount(ids_t, minlength=5).double()[1:]
    e_dir = torch.tensor(args.direction, dtype=torch.float64)
    e_dir = (e_dir / torch.linalg.norm(e_dir)).to(dev)
    p_coef = args.alpha * constants.gamma * constants.mu_0 / (1.0 + args.alpha**2)

    @torch.no_grad()
    def diagnostics():
        m = state.m
        Mvec = (Ms_t * m).double().mean(dim=(0, 1, 2))                     # [3] A/m, mesh average
        Mpar = float(torch.dot(Mvec, e_dir))
        m_par = (m.double() * e_dir).sum(-1).reshape(-1) * Ms_t.double().reshape(-1)
        Mp = (torch.bincount(ids_t, weights=m_par, minlength=5)[1:] / counts)  # per particle, A/m
        H = accumulate_h(term.h(state) for term in terms)
        mxH = torch.linalg.cross(m, H)
        p_dis = float((p_coef * Ms_t * (mxH * mxH).sum(-1, keepdim=True)).double().mean())
        return Mpar, Mvec.cpu().numpy(), Mp.cpu().numpy(), p_dis

    # --- config ------------------------------------------------------------
    stages = build_stages(args, Ms, phi)
    config = dict(vars(args))
    config.update({"N_used": N, "a": geo["a"], "phi_nom": geo["phi_nom"], "phi_vox": phi,
                   "gap_nom": geo["gap_nom"], "d_cells": geo["d_cells"], "Ms": Ms,
                   "sigma_rms_realised": sig_rms_real, "n_cells": N**3,
                   "l_ex": math.sqrt(2 * args.A / (MU0 * Ms**2)),
                   "delta_wall": math.sqrt(args.A / args.K1) if args.K1 > 0 else None,
                   "stages": stages, "git_commit": git_commit(),
                   "torch": torch.__version__, "device": str(dev),
                   "gpu": torch.cuda.get_device_name() if torch.cuda.is_available() else "cpu"})
    (out / "config.json").write_text(json.dumps(config, indent=1))

    # --- resume or initialise ---------------------------------------------
    ckpt = out / "checkpoint.pt"
    summ_file = out / "summary.json"
    summary = {"stages": []}
    first_stage = 0
    if ckpt.exists() and summ_file.exists() and not args.no_resume:
        c = torch.load(ckpt, map_location=dev)
        state.m = c["m"].to(dtype=dt_, device=dev)
        summary = json.loads(summ_file.read_text())
        first_stage = c["next_stage"]
        summary["stages"] = summary["stages"][:first_stage]   # drop an incomplete (benchmark) stage
        print("[resume] from stage %d" % first_stage, flush=True)
    elif args.init_from:
        c = torch.load(args.init_from, map_location=dev)
        state.m = c["m"].to(dtype=dt_, device=dev)
        print("[init] m from %s" % args.init_from, flush=True)
    else:
        rng = np.random.default_rng(args.seed_m)
        m0 = rng.standard_normal((N, N, N, 3))
        m0 /= np.linalg.norm(m0, axis=-1, keepdims=True)
        state.m = T(m0)
        state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)
        drive.H_amp = 0.0
        t_w = time.time()
        llg.step(state, args.relax_time)
        print("[relax] %.3g s done in %.1f s wall" % (args.relax_time, time.time() - t_w), flush=True)

    # --- protocol ----------------------------------------------------------
    period = 1.0 / args.freq
    dt_s = period / args.samples
    drive.freq = args.freq
    cycles_done = sum(len(s["cycles"]) for s in summary["stages"])

    for si in range(first_stage, len(stages)):
        st = stages[si]
        if st["n_cycles"] <= 0:          # e.g. --n_sat 0
            summary["stages"].append({"name": st["name"], "H_amp": st["H_amp"], "measure": st["measure"],
                                      "cycles": [], "complete": True})
            continue
        state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)   # small t -> small rounding
        drive.t0 = 0.0
        drive.H_amp = st["H_amp"]
        rec = {"name": st["name"], "H_amp": st["H_amp"], "measure": st["measure"], "cycles": []}
        fcsv = open(out / ("samples_%02d_%s.csv" % (si, st["name"])), "w")
        fcsv.write("# cycle,t,H,M_par,Mx,My,Mz,B,p_dis,Mp0,Mp1,Mp2,Mp3\n")

        ci = 0
        complete = False
        while True:
            if args.max_cycles is not None and cycles_done >= args.max_cycles:
                break
            t_w = time.time()
            steps0 = state._step
            rows = []
            t = float(state.t)
            Mpar, Mv, Mp, pd = diagnostics()
            rows.append((t, drive.value(t), Mpar, *Mv, pd, *Mp))
            for k in range(args.samples):
                llg.step(state, dt_s)
                t = float(state.t)
                Mpar, Mv, Mp, pd = diagnostics()
                rows.append((t, drive.value(t), Mpar, *Mv, pd, *Mp))
            wall = time.time() - t_w
            steps = state._step - steps0

            a = np.array(rows)                  # float64
            t_a, H_a, M_a, p_a = a[:, 0], a[:, 1], a[:, 2], a[:, 6]
            B_a = MU0 * (H_a + M_a)
            W_loop = float(_trapz(H_a, B_a))  # = mu0 oint H dM (H dH closes)
            W_dis = float(_trapz(p_a, t_a))
            M_peak = float(0.5 * (M_a.max() - M_a.min()))
            cyc = {"cycle": ci,
                   "B_peak": float(0.5 * (B_a.max() - B_a.min())),
                   "H_peak": float(0.5 * (H_a.max() - H_a.min())),
                   "M_peak": M_peak,
                   "W_loop": W_loop, "W_dis": W_dis,
                   "closure": float(abs(M_a[-1] - M_a[0]) / max(M_peak, 1e-30)),
                   "steps": int(steps), "wall_s": wall,
                   "mean_dt": period / max(steps, 1)}
            if ci > 0:
                W_prev = rec["cycles"][-1]["W_loop"]
                cyc["dW_rel"] = float(abs(W_loop - W_prev) / max(abs(W_loop), 1e-300))
            rec["cycles"].append(cyc)
            for r in rows:
                fcsv.write("%d," % ci + ",".join("%.9e" % v for v in (*r[:3], *r[3:6], MU0 * (r[1] + r[2]), *r[6:])) + "\n")
            fcsv.flush()
            cycles_done += 1
            print("[%s c%d] H=%.4g A/m  B_peak=%.4g T  W_loop=%.4g  W_dis=%.4g J/m3  closure=%.2g  dW=%s  "
                  "steps=%d  dt=%.3g ps  wall=%.0f s" %
                  (st["name"], ci, st["H_amp"], cyc["B_peak"], W_loop, W_dis, cyc["closure"],
                   ("%.2g" % cyc["dW_rel"]) if "dW_rel" in cyc else "-",
                   steps, 1e12 * cyc["mean_dt"], wall), flush=True)
            ci += 1

            # stop criterion: fixed count for non-measured stages; for measured
            # stages at least n_cycles, then until the cycle is steady
            if ci >= st["n_cycles"]:
                if not st["measure"]:
                    complete = True
                    break
                steady = (cyc["closure"] < args.closure_tol and
                          cyc.get("dW_rel", 1.0) < args.dW_tol)
                if steady or ci >= args.max_cycles_per_amp:
                    rec["steady"] = bool(steady)
                    complete = True
                    break
        fcsv.close()

        rec["complete"] = complete
        summary["stages"].append(rec)
        summ_file.write_text(json.dumps(summary, indent=1))
        if not complete:
            break
        torch.save({"m": state.m.detach().cpu(), "next_stage": si + 1}, ckpt)
        if args.keep_stage_ckpt:
            torch.save({"m": state.m.detach().cpu(), "next_stage": si + 1},
                       out / ("ckpt_%02d_%s.pt" % (si, st["name"])))
        if args.vti:
            write_vti({"m": state.m}, str(out / ("m_%02d_%s.vti" % (si, st["name"]))), state)

    if args.max_cycles is None or cycles_done < args.max_cycles:
        (out / "DONE").write_text("ok\n")
    if args.timer:
        Timer.print_report()


if __name__ == "__main__":
    main()
