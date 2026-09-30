#!/usr/bin/env python
"""
AC BH-loops of a periodic FCC powder model (fictitious nanocrystalline
material, reduced units).

Inputs are dimensionless (see units.py). Js and A only select the SI unit
system of the output (default Js = 1.5 T, A = 10 pJ/m -> l_ex = 3.34 nm,
f/f_M default -> 30 MHz).

Anisotropy (see anisotropy_model.py): only
  * Herzer residual anisotropy: cubes of edge L_eff with random easy axes,
    Q_eff = K_eff/K_d = (l_ex/L_eff)^2,
  * residual stress on the particle scale: K_p = r_p K_eff per particle,
    deterministic axes.
Box edge a is a multiple of 12 l_ex and L_eff a multiple of 6 l_ex, so that
dx = 3, 2 and 1.5 l_ex see the same cubes. d is adjusted (<= 2.5 %) to keep phi.

One run = one parameter set:
  1. random initial magnetization, short relaxation at H = 0,
  2. one saturating cycle -> resets the domain topology,
  3. descending drive amplitudes h = H/Ms (geometric series). Each amplitude
     runs at least cycles_per_amp cycles, then more (up to
     max_cycles_per_amp) until the cycle is steady (closure_tol, dW_tol).
     dW_tol = 0.10 is the target scatter of w from cycle to cycle.
     The descending series is also an AC demagnetization. The first cycle of
     each amplitude is a transient; analyze.py averages the other cycles.

Per cycle (summary.json), SI and reduced:
  B_peak, b_peak = B_peak/Js      half peak-to-peak flux density (core average)
  W_loop, w_loop = W_loop/K_d     oint H dB = mu0 oint H dM   (per cycle)
  W_dis,  w_dis                   LLG dissipation int <p_dis> dt (per cycle)
  closure                         |M(end)-M(start)| / M_peak
  dW_rel                          |W(c) - W(c-1)| / W(c)
  max_angle_deg                   largest angle between neighbour cells of the
                                  same particle during the cycle
  n_pairs_gt60, frac_pairs_gt60   largest number (fraction) of neighbour pairs with
                                  > 60 deg during the cycle; per sample in the csv.
                                  Jumps in the per-sample count show nucleation or
                                  annihilation of vortex cores / Bloch points
Energies are per unit volume of the CORE (voids included).

Validity diagnostics: W_loop == W_dis for a closed cycle (energy balance).
Pairs with > 60 deg mark structures that the grid does not resolve (vortex
cores, Bloch lines, Bloch points). At dx = 3 l_ex a vortex core is narrower
than one cell, so max_angle_deg is close to 180 deg whenever a vortex exists
(CPU test). Use frac_pairs_gt60 and its per-sample jumps, not max_angle_deg,
to find nucleation/annihilation events.

Precision: set_precision() before the mesh; state.t in float64 (fp32 time
gives ~1 % errors in the accumulated time step); diagnostics in float64.

Robustness: checkpoint.pt and summary.json are written atomically (tmp file
+ rename). A resume stops with an error if the physics or protocol arguments
differ from the run that wrote the checkpoint.

Example (benchmark, 1 cycle):
    CUDA_DEVICE=0 python run_loops.py --d_lex 300 --Leff_lex 18 --max_cycles 1 --out runs/bench
"""
import argparse
import json
import math
import os
import pathlib
import subprocess
import sys
import time
import warnings

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from units import units, F_REL_DEFAULT, JS_REF, A_REF, MU0   # noqa: E402

_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy >= 2.0 renamed trapz

# arguments that do not change the physics or the protocol (allowed to differ on resume)
RUN_CONTROL = {"out", "max_cycles", "vti", "keep_stage_ckpt", "no_resume", "init_from", "timer",
               "allow_new_commit"}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_argument_group("reduced geometry")
    g.add_argument("--d_lex", type=float, required=True, help="particle diameter / l_ex")
    g.add_argument("--dx_lex", type=float, default=3.0, help="cell size / l_ex")
    g.add_argument("--phi", type=float, default=0.65, help="packing fraction (exact; d is adjusted to the box quantum)")
    g.add_argument("--a_quant", type=float, default=12.0, help="box edge quantum / l_ex (0 = off)")

    g = p.add_argument_group("reduced material")
    g.add_argument("--Leff_lex", type=float, default=18.0,
                   help="Herzer exchange length L_eff / l_ex = cube edge (multiple of 6 for mesh-identical cubes)")
    g.add_argument("--Q_eff", type=float, default=None,
                   help="K_eff / K_d (default (l_ex/L_eff)^2 = Herzer; 0 = no Herzer anisotropy)")
    g.add_argument("--r_p", type=float, default=0.0, help="particle-scale stress anisotropy K_p / K_eff")
    g.add_argument("--alpha", type=float, default=0.02, help="Gilbert damping (numerical choice, see HANDOVER)")
    g.add_argument("--ms_void", type=float, default=1e-8, help="Ms_void / Ms (must be > 0)")

    g = p.add_argument_group("unit system (only for SI output)")
    g.add_argument("--Js", type=float, default=JS_REF, help="mu0 Ms [T]")
    g.add_argument("--A", type=float, default=A_REF, help="exchange stiffness [J/m]")

    g = p.add_argument_group("seeds")
    g.add_argument("--seed", type=int, default=1, help="realisation: cube axes and offsets, initial m")

    g = p.add_argument_group("drive / protocol (reduced)")
    g.add_argument("--f_rel", type=float, default=F_REL_DEFAULT, help="f / f_M (default: 30 MHz at Js = 1.5 T)")
    g.add_argument("--direction", type=float, nargs=3, default=[1.0, 0.0, 0.0])
    g.add_argument("--h_sat", type=float, default=None, help="saturating amplitude H/Ms (default 1-phi)")
    g.add_argument("--n_sat", type=int, default=1)
    g.add_argument("--h_max", type=float, default=None, help="largest measured H/Ms (default (1-phi)/3)")
    g.add_argument("--h_factor", type=float, default=0.5)
    g.add_argument("--n_amp", type=int, default=9)
    g.add_argument("--h_list", type=float, nargs="*", default=None, help="explicit amplitudes H/Ms")
    g.add_argument("--b_list", type=float, nargs="*", default=None,
                   help="target amplitudes b = B_peak/Js (descending). The drive amplitude of each stage is "
                        "h = b_target / mu, mu = b/h of the stage before (first stage: --mu_guess)")
    g.add_argument("--mu_guess", type=float, default=8.0, help="b/h for the first --b_list stage")
    g.add_argument("--cycles_per_amp", type=int, default=3, help="minimum cycles per amplitude (1st = transient)")
    g.add_argument("--max_cycles_per_amp", type=int, default=4)
    g.add_argument("--closure_tol", type=float, default=0.02)
    g.add_argument("--dW_tol", type=float, default=0.10,
                   help="target scatter of w per cycle (relative); the next cycle runs if |dW|/W is larger")
    g.add_argument("--samples", type=int, default=256, help="samples per cycle (>= 256 for the energy balance)")
    g.add_argument("--relax_tau", type=float, default=500.0, help="relaxation time at H = 0 in units 1/(gamma Ms)")

    g = p.add_argument_group("numerics")
    g.add_argument("--precision", choices=["single", "double"], default="single")
    g.add_argument("--atol", type=float, default=1e-5, help="RKF45 tolerance")

    g = p.add_argument_group("run control / output")
    g.add_argument("--out", type=str, required=True)
    g.add_argument("--max_cycles", type=int, default=None, help="stop after this many cycles (benchmark)")
    g.add_argument("--vti", action="store_true", help="write m as .vti after each stage")
    g.add_argument("--keep_stage_ckpt", action="store_true", help="keep m after every stage")
    g.add_argument("--no_resume", action="store_true", help="start again: delete the old outputs in --out")
    g.add_argument("--allow_new_commit", action="store_true",
                   help="resume although the code (git commit) changed since the checkpoint")
    g.add_argument("--init_from", type=str, default=None, help="start from m of this checkpoint (no relaxation)")
    g.add_argument("--timer", action="store_true")
    return p.parse_args(argv)


def git_commit():
    """HEAD commit; '-dirty' if tracked files of the study or the library have local changes."""
    try:
        h = subprocess.check_output(["git", "-C", str(HERE), "rev-parse", "HEAD"], text=True).strip()
        dirty = subprocess.check_output(["git", "-C", str(HERE), "status", "--porcelain", "--untracked-files=no",
                                         "--", str(HERE), str(HERE.parent.parent / "magnumnp")], text=True).strip()
        return h + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def write_atomic(path, text):
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def save_atomic(obj, path):
    import torch
    tmp = path.with_name(path.name + ".tmp")
    torch.save(obj, tmp)
    os.replace(tmp, path)


def is_multiple(x, q, tol=1e-9):
    return abs(x / q - round(x / q)) < tol


def build_stages(args, phi):
    h_sat = args.h_sat if args.h_sat is not None else (1.0 - phi)
    if args.b_list:
        stages = [{"name": "sat", "h_amp": h_sat, "n_cycles": args.n_sat, "measure": False}]
        stages += [{"name": "amp%02d" % k, "h_amp": None, "b_target": b, "n_cycles": args.cycles_per_amp,
                    "measure": True} for k, b in enumerate(sorted(args.b_list, reverse=True))]
        return stages
    if args.h_list:
        amps = sorted(args.h_list, reverse=True)
    else:
        h_max = args.h_max if args.h_max is not None else (1.0 - phi) / 3.0
        amps = [h_max * args.h_factor**k for k in range(args.n_amp)]
    stages = [{"name": "sat", "h_amp": h_sat, "n_cycles": args.n_sat, "measure": False}]
    stages += [{"name": "amp%02d" % k, "h_amp": h, "n_cycles": args.cycles_per_amp, "measure": True}
               for k, h in enumerate(amps)]
    return stages


def main(argv=None):
    args = parse_args(argv)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    from magnumnp import (set_precision, set_log_level, Mesh, State, constants, Timer,
                          DemagFieldPBC, ExchangeField, UniaxialAnisotropyField,
                          LLGSolver, accumulate_h, write_vti)
    import torch
    set_log_level(25)
    set_precision(args.precision)
    if args.timer:
        Timer.enable()

    from geometry import fcc_box
    from anisotropy_model import box_edge, cube_axes, particle_axes, CUBE_QUANT
    from field_terms_extra import SinusoidalDrive

    # --- units -------------------------------------------------------------
    U = units(args.Js, args.A)
    Ms, Kd, lex = U["Ms"], U["Kd"], U["l_ex"]
    freq = args.f_rel * U["f_M"]

    # --- geometry (box edge quantized, d adjusted to keep phi) --------------
    a_lex, d_lex = box_edge(args.d_lex, args.phi, args.a_quant)
    N = int(round(a_lex / args.dx_lex))
    if not is_multiple(a_lex, args.dx_lex):
        warnings.warn("a = %g l_ex is not a multiple of dx = %g l_ex: a -> %g l_ex" % (a_lex, args.dx_lex, N * args.dx_lex))
        a_lex = N * args.dx_lex
    d, dx = d_lex * lex, args.dx_lex * lex
    geo = fcc_box(N, d, dx)
    ids = geo["ids"]
    mag = ids >= 0
    phi = geo["phi_vox"]

    mesh = Mesh((N, N, N), (dx, dx, dx), pbc=(1, 1, 1))
    state = State(mesh)
    dt_ = state.dtype
    dev = state.device

    def T(a):
        return torch.tensor(a, dtype=dt_, device=dev)

    mag4 = mag[..., None]
    state.material = {"alpha": args.alpha}
    state.material["Ms"] = T(np.where(mag4, Ms, args.ms_void * Ms))
    state.material["A"] = T(np.where(mag4, args.A, 0.0))

    e_d = np.asarray(args.direction, float); e_d /= np.linalg.norm(e_d)
    L = args.Leff_lex
    Q_eff = args.Q_eff if args.Q_eff is not None else 1.0 / L**2
    K_eff = Q_eff * Kd
    K_p = args.r_p * K_eff
    if not (is_multiple(L, CUBE_QUANT) and is_multiple(L, args.dx_lex)):
        warnings.warn("L_eff = %g l_ex is not a multiple of %g and of dx: the cubes depend on the mesh" % (L, CUBE_QUANT))
    if L < 3.0 * args.dx_lex:
        warnings.warn("L_eff = %g l_ex is less than 3 cells: the anisotropy cubes are not resolved" % L)
    Q_max = Q_eff + args.r_p * Q_eff                     # cube and particle axes parallel
    delta_min = 1.0 / math.sqrt(Q_max) if Q_max > 0 else float("inf")
    if delta_min < 2.0 * args.dx_lex:
        warnings.warn("local wall width %.3g l_ex < 2 cells: walls are not resolved" % delta_min)

    terms = [DemagFieldPBC(), ExchangeField()]
    n_cubes = [0, 0, 0, 0]
    if K_eff > 0.0:
        ax_c, n_cubes = cube_axes(ids, args.dx_lex, a_lex, d_lex, L, args.seed)
        state.material["Ku"] = T(np.where(mag4, K_eff, 0.0))
        state.material["Ku_axis"] = T(ax_c)
        terms.append(UniaxialAnisotropyField())
        del ax_c
    ax_p = particle_axes(e_d)
    if K_p > 0.0:
        kp_ax = np.zeros((N, N, N, 3))
        for p in range(4):
            kp_ax[ids == p] = ax_p[p]
        state.material["Kp"] = T(np.where(mag4, K_p, 0.0))
        state.material["Kp_axis"] = T(kp_ax)
        terms.append(UniaxialAnisotropyField(Ku="Kp", Ku_axis="Kp_axis"))
        del kp_ax
    aniso_info = {"Q_eff": Q_eff, "K_eff": K_eff, "L_eff": L * lex, "L_over_dx": L / args.dx_lex,
                  "K_p": K_p, "Q_p": args.r_p * Q_eff, "cos_theta_p": np.abs(ax_p @ e_d).tolist(),
                  "delta_min_lex": delta_min, "delta_min_over_dx": delta_min / args.dx_lex,
                  "n_cubes_per_particle": [int(n) for n in n_cubes]}

    drive = SinusoidalDrive(state, tuple(e_d))
    terms.append(drive)
    llg = LLGSolver(terms, atol=args.atol)

    # --- diagnostics -------------------------------------------------------
    Ms_t = state.material["Ms"]
    ids_t = torch.tensor(ids.astype(np.int64) + 1, device=dev).reshape(-1)
    counts = torch.bincount(ids_t, minlength=5).double()[1:]
    e_t = torch.tensor(e_d, dtype=torch.float64, device=dev)
    p_coef = args.alpha * constants.gamma * constants.mu_0 / (1.0 + args.alpha**2)
    ids_n = torch.tensor(ids.astype(np.int64), device=dev)
    pair_masks = [(ids_n >= 0) & (ids_n == torch.roll(ids_n, 1, dims=ax)) for ax in range(3)]
    n_pairs = int(sum(int(pm.sum()) for pm in pair_masks))
    cos60 = 0.5

    @torch.no_grad()
    def diagnostics():
        m = state.m
        Mvec = (Ms_t * m).double().mean(dim=(0, 1, 2))
        Mpar = float(torch.dot(Mvec, e_t))
        m_par = (m.double() * e_t).sum(-1).reshape(-1) * Ms_t.double().reshape(-1)
        Mp = torch.bincount(ids_t, weights=m_par, minlength=5)[1:] / counts
        H = accumulate_h(term.h(state) for term in terms)
        mxH = torch.linalg.cross(m, H)
        p_dis = float((p_coef * Ms_t * (mxH * mxH).sum(-1, keepdim=True)).double().mean())
        cmin, n60 = 1.0, 0
        for ax in range(3):
            c = (m * torch.roll(m, 1, dims=ax)).sum(-1)[pair_masks[ax]]
            cmin = min(cmin, float(c.min()))
            n60 += int((c < cos60).sum())
        return Mpar, Mvec.cpu().numpy(), Mp.cpu().numpy(), p_dis, cmin, n60

    # --- config ------------------------------------------------------------
    stages = build_stages(args, args.phi)        # nominal phi: same amplitudes on every mesh
    config = dict(vars(args))
    config.update({"units": U, "N_used": N, "n_cells": N**3, "phi_vox": phi, "phi_nom": geo["phi_nom"],
                   "a_lex": a_lex, "d_lex_used": d_lex, "d": d, "dx": dx, "a": geo["a"],
                   "gap_nom_lex": geo["gap_nom"] / lex, "gap_cells": geo["gap_nom"] / dx,
                   "freq": freq, "anisotropy": aniso_info, "stages": stages, "git_commit": git_commit(),
                   "torch": torch.__version__, "device": str(dev),
                   "gpu": torch.cuda.get_device_name() if torch.cuda.is_available() else "cpu"})

    # --- resume or initialise ---------------------------------------------
    ckpt = out / "checkpoint.pt"
    summ_file = out / "summary.json"
    cfg_file = out / "config.json"
    resume = ckpt.exists() and summ_file.exists() and not args.no_resume
    if resume and cfg_file.exists():
        old = json.loads(cfg_file.read_text())
        diff = [k for k in vars(args) if k not in RUN_CONTROL and old.get(k) != config[k]]
        if diff:
            sys.exit("ERROR: resume with other arguments than the checkpoint run: %s. "
                     "Use --no_resume or another --out." % ", ".join("%s %r -> %r" % (k, old.get(k), config[k]) for k in diff))
        if old.get("git_commit") != config["git_commit"] and not args.allow_new_commit:
            sys.exit("ERROR: the code changed since the checkpoint (%s -> %s). One run would mix two code "
                     "versions. Use --no_resume, or --allow_new_commit if the change does not touch the physics."
                     % (old.get("git_commit"), config["git_commit"]))
    if args.no_resume:
        for f in [ckpt, summ_file, out / "DONE"] + list(out.glob("samples_*.csv")) + list(out.glob("ckpt_*.pt")):
            f.unlink(missing_ok=True)
    write_atomic(cfg_file, json.dumps(config, indent=1))

    summary = {"stages": []}
    first_stage = 0
    if resume:
        c = torch.load(ckpt, map_location=dev)
        state.m = c["m"].to(dtype=dt_, device=dev)
        summary = json.loads(summ_file.read_text())
        first_stage = c["next_stage"]
        summary["stages"] = summary["stages"][:first_stage]
        print("[resume] from stage %d" % first_stage, flush=True)
    elif args.init_from:
        c = torch.load(args.init_from, map_location=dev)
        state.m = c["m"].to(dtype=dt_, device=dev)
        print("[init] m from %s" % args.init_from, flush=True)
    else:
        rng = np.random.default_rng([args.seed, 0])      # own stream (cube axes use [seed, 1000 + p])
        m0 = rng.standard_normal((N, N, N, 3))
        m0 /= np.linalg.norm(m0, axis=-1, keepdims=True)
        state.m = T(m0)
        state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)
        drive.H_amp = 0.0
        t_w = time.time()
        llg.step(state, args.relax_tau * U["t_M"])
        print("[relax] %.3g s done in %.1f s wall" % (args.relax_tau * U["t_M"], time.time() - t_w), flush=True)

    # --- protocol ----------------------------------------------------------
    period = 1.0 / freq
    dt_s = period / args.samples
    drive.freq = freq
    cycles_done = sum(len(s["cycles"]) for s in summary["stages"])

    for si in range(first_stage, len(stages)):
        st = stages[si]
        if st["n_cycles"] <= 0:
            summary["stages"].append({"name": st["name"], "h_amp": st["h_amp"], "measure": st["measure"],
                                      "cycles": [], "complete": True})
            continue
        state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)
        drive.t0 = 0.0
        if st.get("b_target") is not None:
            # h from the measured b/h of the stage before (deterministic, also on resume)
            prev = [r for r in summary["stages"] if r.get("measure") and r.get("cycles")]
            if prev:
                mu = float(np.mean([c["b_peak"] for c in prev[-1]["cycles"]])) / prev[-1]["h_amp"]
            else:
                mu = args.mu_guess
            st = dict(st, h_amp=st["b_target"] / mu)
        drive.H_amp = st["h_amp"] * Ms
        rec = {"name": st["name"], "h_amp": st["h_amp"], "H_amp": st["h_amp"] * Ms,
               "measure": st["measure"], "git_commit": config["git_commit"], "cycles": []}
        if st.get("b_target") is not None:
            rec["b_target"] = st["b_target"]
        fcsv = open(out / ("samples_%02d_%s.csv" % (si, st["name"])), "w")
        fcsv.write("# cycle,t,H,M_par,Mx,My,Mz,B,p_dis,Mp0,Mp1,Mp2,Mp3,n_pairs_gt60  (SI)\n")

        ci = 0
        complete = False
        while True:
            if args.max_cycles is not None and cycles_done >= args.max_cycles:
                break
            t_w = time.time()
            steps0 = state._step
            rows = []
            cmin_c, n60_c = 1.0, 0

            def sample():
                nonlocal cmin_c, n60_c
                t = float(state.t)
                Mpar, Mv, Mp, pd, cmin, n60 = diagnostics()
                cmin_c, n60_c = min(cmin_c, cmin), max(n60_c, n60)
                rows.append((t, drive.value(t), Mpar, *Mv, pd, *Mp, n60))

            sample()
            for k in range(args.samples):
                llg.step(state, dt_s)
                sample()
            wall = time.time() - t_w
            steps = state._step - steps0

            a = np.array(rows)
            t_a, H_a, M_a, p_a = a[:, 0], a[:, 1], a[:, 2], a[:, 6]
            B_a = MU0 * (H_a + M_a)
            W_loop = float(_trapz(H_a, B_a))
            W_dis = float(_trapz(p_a, t_a))
            M_peak = float(0.5 * (M_a.max() - M_a.min()))
            B_peak = float(0.5 * (B_a.max() - B_a.min()))
            cyc = {"cycle": ci,
                   "B_peak": B_peak, "b_peak": B_peak / args.Js,
                   "H_peak": float(0.5 * (H_a.max() - H_a.min())), "M_peak": M_peak,
                   "W_loop": W_loop, "W_dis": W_dis, "w_loop": W_loop / Kd, "w_dis": W_dis / Kd,
                   "closure": float(abs(M_a[-1] - M_a[0]) / max(M_peak, 1e-30)),
                   "max_angle_deg": float(math.degrees(math.acos(max(-1.0, min(1.0, cmin_c))))),
                   "n_pairs_gt60": int(n60_c), "frac_pairs_gt60": n60_c / max(n_pairs, 1),
                   "steps": int(steps), "wall_s": wall, "mean_dt": period / max(steps, 1)}
            if ci > 0:
                W_prev = rec["cycles"][-1]["W_loop"]
                cyc["dW_rel"] = float(abs(W_loop - W_prev) / max(abs(W_loop), 1e-300))
            rec["cycles"].append(cyc)
            for r in rows:
                fcsv.write("%d," % ci + ",".join("%.9e" % v for v in (*r[:3], *r[3:6], MU0 * (r[1] + r[2]), *r[6:])) + "\n")
            fcsv.flush()
            cycles_done += 1
            print("[%s c%d] h=%.4g  b_peak=%.4g  w_loop=%.4g  w_dis=%.4g  closure=%.2g  dW=%s  "
                  "max_angle=%.0f deg  steps=%d  dt=%.3g ps  wall=%.0f s" %
                  (st["name"], ci, st["h_amp"], cyc["b_peak"], cyc["w_loop"], cyc["w_dis"], cyc["closure"],
                   ("%.2g" % cyc["dW_rel"]) if "dW_rel" in cyc else "-", cyc["max_angle_deg"],
                   steps, 1e12 * cyc["mean_dt"], wall), flush=True)
            ci += 1

            if ci >= st["n_cycles"]:
                if not st["measure"]:
                    complete = True
                    break
                closed = cyc["closure"] < args.closure_tol
                in_scatter = cyc.get("dW_rel", 1.0) < args.dW_tol
                steady = closed and in_scatter
                if steady or ci >= args.max_cycles_per_amp:
                    rec["steady"] = bool(steady)
                    rec["closed"] = bool(closed)
                    rec["in_scatter"] = bool(in_scatter)
                    complete = True
                    break
        fcsv.close()

        rec["complete"] = complete
        summary["stages"].append(rec)
        write_atomic(summ_file, json.dumps(summary, indent=1))
        if not complete:
            break
        save_atomic({"m": state.m.detach().cpu(), "next_stage": si + 1}, ckpt)
        if args.keep_stage_ckpt:
            torch.save({"m": state.m.detach().cpu(), "next_stage": si + 1},
                       out / ("ckpt_%02d_%s.pt" % (si, st["name"])))
        if args.vti:
            write_vti({"m": state.m}, str(out / ("m_%02d_%s.vti" % (si, st["name"]))), state)

    if args.max_cycles is None or cycles_done < args.max_cycles:
        write_atomic(out / "DONE", "ok\n")
    if args.timer:
        Timer.print_report()


if __name__ == "__main__":
    main()
