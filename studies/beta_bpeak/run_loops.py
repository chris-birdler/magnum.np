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
    deterministic axes (an isotropic tetrahedral set, PROTOKOLL 7.37).
Box edge a is a multiple of 12 l_ex and L_eff a multiple of 6 l_ex, so that
dx = 3, 2 and 1.5 l_ex see the same cubes. d is adjusted (<= 2.5 %) to keep phi.

One run = one parameter set:
  1. random initial magnetization, short relaxation at H = 0,
  2. one saturating cycle -> resets the domain topology,
  3. descending drive amplitudes h = H/Ms (geometric series). Each amplitude
     runs at least cycles_per_amp cycles, then more (up to
     max_cycles_per_amp) until the cycle is steady (closure_tol, dW_tol).
     dW_tol = 0.10 is the target scatter of w from cycle to cycle.
     With --min_calm N > 0 (PLAN 4.5) this stop rule is replaced: after
     cycles_per_amp cycles a stage continues until its calm tail (calm
     cycles at the end of the measured cycles; calm: w_loop > 0 and
     |w_dis / w_loop - 1| <= calm_tol) is >= N, at most max_cycles_per_amp.
     --cycles_last_amp M: the last measured stage runs M cycles.
     With --b_list the drive of a stage is set from b/h of the stage before
     and corrected once after cycle 0 with b/h of cycle 0 (analyze.py then
     discards cycles 0 and 1).
     The descending series is also an AC demagnetization. The first cycle of
     each amplitude is a transient; analyze.py averages the other cycles.

Per cycle (summary.json), SI and reduced:
  B_peak, b_peak = B_peak/Js      half peak-to-peak flux density (core average)
  W_loop, w_loop = W_loop/K_d     oint H dB = mu0 oint H dM   (per cycle)
  W_dis,  w_dis                   LLG dissipation int <p_dis> dt (per cycle)
  closure                         |M(end)-M(start)| / M_peak
  offset                          loop centre / amplitude = (M_max + M_min)/(M_max - M_min):
                                  a minor loop around a remanent state (not AC-demagnetized)
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
import re
import subprocess
import sys
import time
import warnings

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))   # magnumnp of this checkout (git_commit records it), not an installed copy
#                                               (PROTOKOLL 7.108: a worktree run had loaded the magnumnp of /root/magnum.np)

from units import units, F_REL_DEFAULT, JS_REF, A_REF, MU0   # noqa: E402

_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy >= 2.0 renamed trapz

# arguments that do not change the physics or the protocol (allowed to differ on resume)
ACC_FROM_CYCLE = 2      # --acc: cycles kept from this index (analyze skips 1, or 2 after a drive correction; one more
#                         after a second correction, acc_first)
NEW_KEYS = {"min_calm": 0, "calm_tol": 0.20, "cycles_last_amp": None,   # added 2026-10-11 (PLAN 4.5): an old config
            "acc_harm": 0, "init_transfer": None,                  # without them is equal to these defaults
            "strayfield": "fd", "strayfield_shells": 4,            # (acc_harm: 2026-10-12, PROTOKOLL 7.72;
            "b_tol": 0.0,                                          # b_tol: 2026-10-15, PROTOKOLL 7.107)
            "atol_cycles": None, "rtol_cycles": None,              # 2026-10-16, PROTOKOLL 7.108
            "strayfield_cycles": None, "geom_dx_lex": None,       # geom_dx_lex: 2026-10-17, PLAN 4.12 V3g
            "precision_cycles": None,                             # 2026-10-17, decision Chris (fp64 test)
            "relax_tol_E_abs": None,                              # 2026-10-17, PLAN 4.13 S7 (D1)
            "init_stage": None,                                   # 2026-10-18, solver test at 150 mT (Chris)
            "predemag_h": None, "predemag_q": 0.8,                # 2026-10-18, Bloch points (Chris)
            "init_mode": "block",                                 # 2026-10-18, Bloch points, way 2 (Chris)
            "bloch_free_tries": 0}                                # 2026-10-18, Bloch points (Chris)
#                                                                    strayfield: 2026-10-15, PROTOKOLL 7.105)
RUN_CONTROL = {"out", "max_cycles", "vti", "keep_stage_ckpt", "no_resume", "init_from", "init_stage", "timer",
               "allow_new_commit", "allow_unconverged", "init_wait_h", "relax_maxiter"}
# relax_maxiter is only a limit: a relaxation that converged gives the same state for any larger limit
# (deterministic, audit 3, item 5); run_queue raises it after a gate failure (exit code 3)
# arguments that set the relaxed virgin state; an --init_from state must agree in all of them
# (the drive, the amplitudes and alpha can differ: the relaxation runs at alpha = 1 and H = 0)
INIT_KEYS = ["d_lex", "dx_lex", "phi", "a_quant", "surface", "Leff_lex", "Q_eff", "r_p", "ms_void", "Js", "A",
             "seed", "init_seed", "init_block", "protocol", "relax_maxiter", "relax_dm_tol", "relax_mode",
             "relax_chunk", "relax_tol_E", "relax_confirm", "precision", "atol", "strayfield", "strayfield_shells",
             "geom_dx_lex", "relax_tol_E_abs", "predemag_h", "predemag_q", "init_mode",
             "bloch_free_tries"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_argument_group("reduced geometry")
    g.add_argument("--d_lex", type=float, required=True, help="particle diameter / l_ex")
    g.add_argument("--dx_lex", type=float, default=3.0, help="cell size / l_ex")
    g.add_argument("--phi", type=float, default=0.65, help="packing fraction (exact; d is adjusted to the box quantum)")
    g.add_argument("--a_quant", type=float, default=12.0, help="box edge quantum / l_ex (0 = off)")
    g.add_argument("--geom_dx_lex", type=float, default=None,
                   help="lock the voxel geometry, the particle index and the cube axes at this coarser cell size (a "
                        "multiple of dx; voxel surface only): each fine cell repeats its coarse voxel, so only dx "
                        "changes between meshes (PLAN 4.12 V3g)")
    g.add_argument("--surface", choices=["voxel", "fraction"], default="voxel",
                   help="voxel: staircase surface; fraction: surface cells carry their magnetic volume "
                        "fraction f (8^3 samples); Ms, A, K_eff and K_p are scaled with f")

    g = p.add_argument_group("reduced material")
    g.add_argument("--Leff_lex", type=float, default=18.0,
                   help="Herzer exchange length L_eff / l_ex = cube edge (multiple of 6 for mesh-identical cubes)")
    g.add_argument("--Q_eff", type=float, default=None,
                   help="K_eff / K_d (default (l_ex/L_eff)^2 = Herzer; 0 = no Herzer anisotropy)")
    g.add_argument("--r_p", type=float, default=0.0, help="particle-scale stress anisotropy K_p / K_eff")
    g.add_argument("--alpha", type=float, default=0.02, help="Gilbert damping (numerical choice, see PROTOKOLL)")
    g.add_argument("--ms_void", type=float, default=1e-8, help="Ms_void / Ms (must be > 0)")

    g = p.add_argument_group("unit system (only for SI output)")
    g.add_argument("--Js", type=float, default=JS_REF, help="mu0 Ms [T]")
    g.add_argument("--A", type=float, default=A_REF, help="exchange stiffness [J/m]")

    g = p.add_argument_group("seeds")
    g.add_argument("--seed", type=int, default=1, help="realisation: cube axes and offsets, initial m")
    g.add_argument("--init_seed", type=int, default=None, help="seed of the random initial m (default --seed)")
    g.add_argument("--init_mode", choices=["block", "particle"], default="block",
                   help="virgin start (Chris 2026-10-18, Bloch points): block = random m on --init_block blocks; "
                        "particle = each particle uniform along its own random direction")
    g.add_argument("--bloch_free_tries", type=int, default=0,
                   help="virgin_asc converge (Chris 2026-10-18): if the prepared state (after --predemag_h if set) has a "
                        "Bloch point, a new random start (other stream of init_seed; the microstructure stays) and the "
                        "preparation again, at most this many attempts in all; 0 = one attempt, no check")
    g.add_argument("--init_block", type=float, default=0.0,
                   help="random initial m constant on blocks of this edge / l_ex in physical coordinates "
                        "(multiple of every dx: the same start on every mesh); 0 = one direction per cell")

    g = p.add_argument_group("drive / protocol (reduced)")
    g.add_argument("--f_rel", type=float, default=F_REL_DEFAULT, help="f / f_M (default: 30 MHz at Js = 1.5 T)")
    g.add_argument("--direction", type=float, nargs=3, default=[1.0, 0.0, 0.0])
    g.add_argument("--protocol", choices=["reset_desc", "virgin_asc", "acdemag_asc"], default="reset_desc",
                   help="reset_desc: random m, short relaxation, 1 saturating cycle, descending amplitudes; "
                        "virgin_asc: random m, full relaxation at H = 0 (T = 0 analogue of the thermally "
                        "demagnetized state after the anneal), ascending amplitudes; "
                        "acdemag_asc: saturating cycle with a decaying amplitude (AC demagnetization), "
                        "ascending amplitudes")
    g.add_argument("--h_sat", type=float, default=None, help="saturating amplitude H/Ms (default 1-phi)")
    g.add_argument("--n_sat", type=int, default=1)
    g.add_argument("--demag_q", type=float, default=0.8, help="acdemag: amplitude factor per half cycle")
    g.add_argument("--demag_f_factor", type=float, default=1.0,
                   help="acdemag: frequency of the demagnetizing cycles / f (only the end state is used)")
    g.add_argument("--demag_h_end", type=float, default=None,
                   help="acdemag: end amplitude H/Ms (default 0.3 x the first stage estimate)")
    g.add_argument("--relax_maxiter", type=int, default=5000, help="virgin: iterations of 1e-11 s at alpha = 1")
    g.add_argument("--relax_dm_tol", type=float, default=100.0, help="virgin: max |dm/dt|/gamma [A/m] to stop")
    g.add_argument("--relax_mode", choices=["fixed", "converge"], default="fixed",
                   help="virgin: fixed = magnum.np relax (max torque < relax_dm_tol or relax_maxiter); "
                        "converge = chunks of relax_chunk iterations until the relative energy change per "
                        "chunk < relax_tol_E (the max torque stays large at unresolved structures)")
    g.add_argument("--relax_chunk", type=int, default=500)
    g.add_argument("--relax_tol_E", type=float, default=1e-6)
    g.add_argument("--relax_tol_E_abs", type=float, default=None,
                   help="converge (PLAN 4.13 S7, D1): a chunk is quiet if |dE| / (K_d V) over the chunk is at most "
                        "this value (an absolute threshold; a rise below it counts as quiet). Replaces relax_tol_E")
    g.add_argument("--predemag_h", type=float, default=None,
                   help="virgin_asc (Bloch points, Chris 2026-10-18): after the virgin relaxation an AC demagnetization "
                        "at the drive frequency from this amplitude H/Ms (0: --h_sat, default 1 - phi), the amplitude "
                        "times --predemag_q at each zero crossing of the field down to 0.3 x the first b target / "
                        "mu_guess, then the relaxation again. With --init_from a virgin state without it: applied "
                        "after loading")
    g.add_argument("--predemag_q", type=float, default=0.8, help="--predemag_h: amplitude factor per half cycle")
    g.add_argument("--relax_confirm", type=int, default=2,
                   help="converge: number of consecutive chunks with a relative energy change < relax_tol_E "
                        "(one quiet chunk can come before an avalanche)")
    g.add_argument("--h_max", type=float, default=None, help="largest measured H/Ms (default (1-phi)/3)")
    g.add_argument("--h_factor", type=float, default=0.5)
    g.add_argument("--n_amp", type=int, default=9)
    g.add_argument("--h_list", type=float, nargs="*", default=None, help="explicit amplitudes H/Ms")
    g.add_argument("--b_list", type=float, nargs="*", default=None,
                   help="target amplitudes b = B_peak/Js (descending). The drive amplitude of each stage is "
                        "h = b_target / mu, mu = b/h of the stage before (first stage: --mu_guess)")
    g.add_argument("--mu_guess", type=float, default=8.0, help="b/h for the first --b_list stage")
    g.add_argument("--cycles_per_amp", type=int, default=3, help="minimum cycles per amplitude (1st = transient)")
    g.add_argument("--b_tol", type=float, default=0.0,
                   help="second drive correction with b of cycle 1 if |b / b_target - 1| > b_tol (0: off); then "
                        "cycle 2 is discarded too and the stage gets one more cycle (PROTOKOLL 7.107)")
    g.add_argument("--max_cycles_per_amp", type=int, default=4)
    g.add_argument("--min_calm", type=int, default=0,
                   help="PLAN 4.5: > 0 replaces the closure / dW stop rule: a measured stage continues after "
                        "cycles_per_amp cycles until its calm tail is >= min_calm (at most max_cycles_per_amp)")
    g.add_argument("--calm_tol", type=float, default=0.20, help="calm cycle: w_loop > 0, |w_dis/w_loop - 1| <= calm_tol")
    g.add_argument("--cycles_last_amp", type=int, default=None, help="cycles of the last measured stage (PLAN 4.5 D)")
    g.add_argument("--closure_tol", type=float, default=0.02)
    g.add_argument("--dW_tol", type=float, default=0.10,
                   help="target scatter of w per cycle (relative); the next cycle runs if |dW|/W is larger")
    g.add_argument("--samples", type=int, default=256, help="samples per cycle (>= 256 for the energy balance)")
    g.add_argument("--relax_tau", type=float, default=500.0, help="relaxation time at H = 0 in units 1/(gamma Ms)")

    g = p.add_argument_group("numerics")
    g.add_argument("--precision", choices=["single", "double"], default="single")
    g.add_argument("--precision_cycles", choices=["single", "double"], default=None,
                   help="precision of the drive cycles (default: --precision). Only with --init_from: the run loads a "
                        "state relaxed with --precision and runs all of its own work in this precision; not an "
                        "INIT_KEY (decision Chris 2026-10-17, fp64 test)")
    g.add_argument("--atol", type=float, default=1e-5, help="RKF45 tolerance")
    g.add_argument("--atol_cycles", type=float, default=None,
                   help="RKF45 atol of the drive cycles (default: --atol). Not an INIT_KEY: the relaxation uses --atol, "
                        "so a state of --init_from can run its cycles with another tolerance (PROTOKOLL 7.108)")
    g.add_argument("--rtol_cycles", type=float, default=None,
                   help="RKF45 rtol of the drive cycles (default: the RKF45 default 1e-5; PROTOKOLL 7.108)")
    g.add_argument("--strayfield", choices=["fd", "newell", "newell_ext"], default="fd",
                   help="fd: DemagFieldPBC (finite differences, Bruckner 2021; error at surfaces and cell-scale null "
                        "modes, PROTOKOLL 7.105); newell: DemagField with image shells and compensated outer surface "
                        "charges (true_pbc); newell_ext: the same with the kernel extrapolated in the number of shells "
                        "(1 / (2k + 1)^2 law of the image sum, PROTOKOLL 7.108)")
    g.add_argument("--strayfield_shells", type=int, default=4,
                   help="image shells per side for --strayfield newell")
    g.add_argument("--strayfield_cycles", choices=["fd", "newell", "newell_ext"], default=None,
                   help="stray-field form of the drive cycles (default: --strayfield). Not an INIT_KEY: an --init_from "
                        "state relaxed with newell can run its cycles with newell_ext (PROTOKOLL 7.108); fd: "
                        "DemagFieldPBC in the cycles of a run with a Newell --strayfield (solver test, 2026-10-18)")

    g = p.add_argument_group("run control / output")
    g.add_argument("--out", type=str, required=True)
    g.add_argument("--max_cycles", type=int, default=None, help="stop after this many cycles (benchmark)")
    g.add_argument("--vti", action="store_true", help="write m as .vti after each stage")
    g.add_argument("--snap_phases", type=int, default=0,
                   help="write m as .vti at N equal phase points of the LAST cycle of each measured stage "
                        "(N = 4: H max, H = 0 falling, H min, H = 0 rising) and geometry.vti once")
    g.add_argument("--snap_mT", type=float, nargs="*", default=None,
                   help="with --snap_phases: snapshots only in the stages whose b target is within 10 %% of these "
                        "B_peak values [mT] (default: all measured stages)")
    g.add_argument("--snap_cycles", type=int, default=1,
                   help="with --snap_phases: snapshots in the last N cycles of each measured stage")
    g.add_argument("--ringdown_periods", type=float, default=0.0,
                   help="after each measured stage: a quarter period more to H max, then H held constant at H max "
                        "for this many periods; p_dis, M and n60 at --ringdown_samples per period go to "
                        "ringdown_<i>_<stage>.csv; m and t are restored afterwards (the next stage is not changed)")
    g.add_argument("--ringdown_samples", type=int, default=512, help="samples per period in the ring-down")
    g.add_argument("--keep_stage_ckpt", action="store_true", help="keep m after every stage")
    g.add_argument("--acc", action="store_true",
                   help="per-cell accumulators over the cycles >= %d (>= %d after a second drive correction) of each "
                        "measured stage (all samples of a "
                        "cycle): mean local dissipation p, mean m (m0) and the fundamental m1 = 2/n sum m e^{-i w t}; "
                        "written to acc_<i>_<stage>.npz (particle cells only; audit 3, item 9b)"
                        % (ACC_FROM_CYCLE, ACC_FROM_CYCLE + 1))
    g.add_argument("--acc_harm", type=int, default=0,
                   help="with --acc (PROTOKOLL 7.72): K >= 2 adds per cell |m_k|^2 of the harmonics k = 1 … K "
                        "(per cycle, summed over the components, averaged over the kept cycles), and only calm "
                        "cycles (w_loop > 0, |w_dis / w_loop - 1| <= --calm_tol) enter all accumulators; per cycle "
                        "the box sums go to the cycle record. 0 (default): the accumulators as before")
    g.add_argument("--no_resume", action="store_true", help="start again: delete the old outputs in --out")
    g.add_argument("--allow_new_commit", action="store_true",
                   help="resume although the code (git commit) changed since the checkpoint")
    g.add_argument("--init_from", type=str, default=None,
                   help="start from the relaxed virgin state init.pt of another run (no relaxation); a "
                        "relative path is relative to the parent of --out (runs/); the arguments in "
                        "INIT_KEYS must agree")
    g.add_argument("--init_stage", type=str, default=None,
                   help="solver test (Chris 2026-10-18): start from a stage state file ckpt_<i>_<stage>.pt of another run "
                        "(--keep_stage_ckpt) at its next stage; the earlier stages of that run are copied from its "
                        "summary.json (marked copied_from), so the drive of the next stage is the same; all INIT_KEYS and "
                        "b_list must agree; a relative path is relative to the parent of --out")
    g.add_argument("--init_transfer", type=str, default=None,
                   help="virgin_asc (PLAN 4.10): start the relaxation from the relaxed state init.pt of a run on a finer "
                        "mesh (integer ratio of the cells) instead of a random state: block average of m over the "
                        "particle sub-cells, normalised; then the production relaxation. All INIT_KEYS except dx_lex "
                        "and relax_maxiter must agree; waits like --init_from (--init_wait_h)")
    g.add_argument("--init_wait_h", type=float, default=0.0,
                   help="wait up to this many hours for the --init_from file (it is written by a job that "
                        "still relaxes)")
    g.add_argument("--allow_unconverged", action="store_true",
                   help="continue although the virgin relaxation did not converge (default: stop, exit code 3)")
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


def replace_durable(tmp, path):
    """os.replace with the data of tmp on the disk before the rename and the directory entry on the disk after it
    (PROTOKOLL 7.102: after a hard restart of the instance a renamed init.pt had 0 bytes, because its data were still
    in the page cache)."""
    fd = os.open(tmp, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)
    os.replace(tmp, path)
    fd = os.open(pathlib.Path(path).parent, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def write_atomic(path, text):
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    replace_durable(tmp, path)


def save_atomic(obj, path):
    import torch
    tmp = path.with_name(path.name + ".tmp")
    torch.save(obj, tmp)
    replace_durable(tmp, path)


def strayfield_setup(args):
    """(mesh pbc, stray-field term) for --strayfield (PROTOKOLL 7.105). fd: DemagFieldPBC, true periodic by its k = 0 rule (the
    pbc count only switches on the periodic exchange); newell: DemagField with --strayfield_shells image shells per side and
    the compensation of the outer surface charges of the image block (true_pbc), the same zero mean field."""
    from magnumnp import DemagField, DemagFieldPBC
    if args.strayfield in ("newell", "newell_ext"):
        ext = args.strayfield == "newell_ext"
        if args.strayfield_shells < (2 if ext else 1):
            raise ValueError("--strayfield %s needs --strayfield_shells >= %d" % (args.strayfield, 2 if ext else 1))
        return (args.strayfield_shells,) * 3, DemagField(true_pbc=True, extrapolate=ext)
    return (1, 1, 1), DemagFieldPBC()


def cycle_strayfield(args):
    """The stray-field term of the drive cycles if --strayfield_cycles differs from --strayfield, else None
    (PROTOKOLL 7.108). The mesh keeps the pbc of --strayfield_shells: thus --strayfield must be a Newell form. fd
    (solver test 2026-10-18): DemagFieldPBC does not use the pbc count (its k = 0 rule), and the exchange is periodic
    for any pbc > 0, so the cycles have the physics of a --strayfield fd run."""
    import types
    if args.strayfield_cycles is None or args.strayfield_cycles == args.strayfield:
        return None
    if args.strayfield == "fd":
        raise ValueError("--strayfield_cycles needs a Newell form in --strayfield (the mesh pbc stays)")
    return strayfield_setup(types.SimpleNamespace(strayfield=args.strayfield_cycles,
                                                  strayfield_shells=args.strayfield_shells))[1]


def geometry_rep(args):
    """Number of fine cells per edge of a locked coarse voxel (1 without --geom_dx_lex)."""
    if not args.geom_dx_lex:
        return 1
    rep = int(round(args.geom_dx_lex / args.dx_lex))
    if rep < 1 or not is_multiple(args.geom_dx_lex, args.dx_lex):
        raise ValueError("--geom_dx_lex %g is not a multiple of --dx_lex %g" % (args.geom_dx_lex, args.dx_lex))
    if rep > 1 and args.surface != "voxel":
        raise ValueError("--geom_dx_lex needs --surface voxel")
    return rep


def upsample(x, rep):
    """Repeat each cell of x [n, n, n, ...] rep times along the three grid axes."""
    for ax in range(3):
        x = np.repeat(x, rep, axis=ax)
    return x


def locked_geometry(N, d, dx, rep, sub=0):
    """fcc_box on the mesh; with rep > 1 the voxel geometry of the cell size rep dx, repeated onto the mesh (ids and
    frac of each coarse voxel in its rep^3 fine cells; PLAN 4.12 V3g)."""
    from geometry import fcc_box
    if rep == 1:
        return fcc_box(N, d, dx, sub=sub)
    if N % rep:
        raise ValueError("N = %d is not a multiple of %d" % (N, rep))
    g = fcc_box(N // rep, d, dx * rep, sub=0)
    return dict(g, ids=upsample(g["ids"], rep), frac=upsample(g["frac"], rep))


def locked_cube_axes(ids, rep, dx_lex, a_lex, d_lex, L, seed):
    """cube_axes on the mesh; with rep > 1 on the coarse voxels (ids[::rep, ::rep, ::rep]) and repeated."""
    from anisotropy_model import cube_axes
    if rep == 1:
        return cube_axes(ids, dx_lex, a_lex, d_lex, L, seed)
    ax, n = cube_axes(ids[::rep, ::rep, ::rep], dx_lex * rep, a_lex, d_lex, L, seed)
    return upsample(ax, rep), n


def run_precision(args):
    """Precision of this run: --precision_cycles (only with --init_from, the run then does no relaxation) or
    --precision."""
    if args.precision_cycles is None or args.precision_cycles == args.precision:
        return args.precision
    if not args.init_from:
        raise ValueError("--precision_cycles %s needs --init_from (the relaxation keeps --precision %s)"
                         % (args.precision_cycles, args.precision))
    return args.precision_cycles


def cycle_tolerances(args):
    """Keyword arguments of the RKF45 solver for the drive cycles (PROTOKOLL 7.108), or None if the cycles use the
    solver of the relaxation (--atol, RKF45 default rtol). The error norm of RKF45 is err / (atol + rtol |m_i|)."""
    if args.atol_cycles is None and args.rtol_cycles is None:
        return None
    kw = {"atol": args.atol if args.atol_cycles is None else args.atol_cycles}
    if args.rtol_cycles is not None:
        kw["rtol"] = args.rtol_cycles
    return kw


def transfer_m(m_f, ids_f, ids_c):
    """Block average of the fine state m_f [Nf,Nf,Nf,3] onto the coarse grid of ids_c [Nc,Nc,Nc] (Nf = k Nc), normalised
    (PLAN 4.10). Each coarse cell takes the mean over the fine sub-cells with the same index (particle 0 … 3 or void -1);
    a coarse cell without such a sub-cell (surface mismatch) takes the plain block mean. Returns (m_c, stats)."""
    Nf, Nc = ids_f.shape[0], ids_c.shape[0]
    k = Nf // Nc
    if k * Nc != Nf or k < 1:
        raise ValueError("the fine mesh (%d cells) is not an integer multiple of the coarse mesh (%d cells)" % (Nf, Nc))
    m = np.asarray(m_f, np.float64).reshape(Nc, k, Nc, k, Nc, k, 3)
    same = ids_f.reshape(Nc, k, Nc, k, Nc, k) == ids_c[:, None, :, None, :, None]
    cnt = same.sum(axis=(1, 3, 5))
    v = (m * same[..., None]).sum(axis=(1, 3, 5))
    fall = cnt == 0
    v[fall] = m.sum(axis=(1, 3, 5))[fall]
    nrm = np.linalg.norm(v, axis=-1)
    zero = nrm < 1e-12
    v = np.where(zero[..., None], np.array([1.0, 0.0, 0.0]), v / np.maximum(nrm, 1e-300)[..., None])
    mag_c = ids_c >= 0
    stats = {"factor": int(k), "coarse_particle_cells": int(mag_c.sum()),
             "no_same_index_subcell": int((mag_c & fall).sum()),
             "partial_blocks": int((mag_c & ~fall & (cnt < k ** 3)).sum()),
             "zero_vector": int((mag_c & zero).sum()),
             "mean_norm_before": float(np.mean(nrm[mag_c & ~fall] / cnt[mag_c & ~fall])) if (mag_c & ~fall).any() else None}
    return v, stats


def virgin_m0(N, a_lex, dx_lex, init_block, init_seed, mode="block", ids=None, attempt=0):
    """Random initial m of the virgin state: isotropic unit vectors (normalised Gaussian), constant on blocks of
    edge init_block l_ex (0: per cell); own random stream [init_seed, 0] (the cube axes use [seed, 1000 + p]).
    mode "particle" (Chris 2026-10-18, Bloch points): each particle (ids) uniform along its own isotropic random
    direction (stream [init_seed, 2]); the void cells keep the block field (they carry no magnetization).
    attempt k > 0 (--bloch_free_tries): other streams [init_seed, 10 + k] (blocks) and [init_seed, 20 + k] (particles)."""
    if mode == "particle":
        m0 = virgin_m0(N, a_lex, dx_lex, init_block, init_seed, attempt=attempt)
        g = np.random.default_rng([init_seed, 2 if attempt == 0 else 20 + attempt]).standard_normal((4, 3))
        g /= np.linalg.norm(g, axis=-1, keepdims=True)
        for p in range(4):
            m0[ids == p] = g[p]
        return m0
    rng = np.random.default_rng([init_seed, 0 if attempt == 0 else 10 + attempt])
    if init_block > 0.0:
        nb = int(round(a_lex / init_block))
        tab = rng.standard_normal((nb, nb, nb, 3))
        ib = np.minimum(np.floor((np.arange(N) + 0.5) * dx_lex / init_block).astype(int), nb - 1)
        m0 = tab[np.ix_(ib, ib, ib)]
    else:
        m0 = rng.standard_normal((N, N, N, 3))
    return m0 / np.linalg.norm(m0, axis=-1, keepdims=True)


def cube_charge_t(m):
    """Topological charge of each elementary cube of 2 x 2 x 2 cell centres (torch, on the device of m): the sum of the
    oriented solid angles of its 6 faces / 4π; ±1 at a Bloch point (Chris 2026-10-18; analyze_defect.cube_charges)."""
    import torch

    def sa(a, b, c):
        num = (a * torch.linalg.cross(b, c, dim=-1)).sum(-1)
        den = 1.0 + (a * b).sum(-1) + (b * c).sum(-1) + (c * a).sum(-1)
        return 2.0 * torch.atan2(num, den)

    q = torch.zeros(m.shape[:3], dtype=m.dtype, device=m.device)
    for ax, (u, v) in enumerate(((1, 2), (2, 0), (0, 1))):
        b = torch.roll(m, -1, dims=u)
        c = torch.roll(b, -1, dims=v)
        d = torch.roll(m, -1, dims=v)
        om = sa(m, b, c) + sa(m, c, d)
        q += torch.roll(om, -1, dims=ax) - om
    return q / (4.0 * math.pi)


def bloch_mask_t(ids_t):
    """Cubes with all 8 corners in one particle (torch bool)."""
    import torch
    ok = ids_t >= 0
    for s in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)):
        ok &= torch.roll(ids_t, shifts=(-s[0], -s[1], -s[2]), dims=(0, 1, 2)) == ids_t
    return ok


def stage_start(src, config):
    """--init_stage (solver test, Chris 2026-10-18): (checkpoint dict, summary, first stage). src is a stage state file
    ckpt_<i>_<stage>.pt of another run; that run (the folder of src) must agree in all INIT_KEYS (relax_maxiter only a
    limit); b_list must be its b_list or the start of it with more stages than the next stage (then the run ends with
    the last stage of its b_list). Its stages before the next stage are copied (marked copied_from), so that the drive
    amplitude of the next stage comes from the same b/h."""
    import torch
    src = pathlib.Path(src)
    if not src.exists():
        sys.exit("ERROR: --init_stage %s does not exist" % src)
    scfg = json.loads((src.parent / "config.json").read_text())
    diff = [k for k in INIT_KEYS if k != "relax_maxiter" and scfg.get(k, NEW_KEYS.get(k)) != config.get(k)]
    if diff:
        sys.exit("ERROR: %s comes from a run with other arguments: %s" % (src, ", ".join(
            "%s %r -> %r" % (k, scfg.get(k, NEW_KEYS.get(k)), config.get(k)) for k in diff)))
    c = torch.load(str(src), map_location="cpu")
    ssum = json.loads((src.parent / "summary.json").read_text())
    k = int(c["next_stage"])
    # b_list: the same list, or its first stages (the run then ends after stage len(b_list) - 1; solver scan)
    bl, sbl = list(config.get("b_list") or []), list(scfg.get("b_list") or [])
    if not (bl == sbl[:len(bl)] and len(bl) > k):
        sys.exit("ERROR: --b_list %r is not the start of the b_list %r of %s with more than %d stages" % (bl, sbl, src, k))
    stages = [dict(s, copied_from=str(src.parent)) for s in ssum["stages"][:k]]
    summary = {"stages": stages, "init": dict(ssum.get("init", {}), source=str(src)),
               "init_stage": {"source": str(src), "first_stage": k}}
    return c, summary, k


def quiet_chunk(E, E_prev, norm, tol_rel, tol_abs=None):
    """Relaxation test of one chunk (PLAN 4.13 S7): (quiet, relative change, |dE| / norm). With tol_abs the chunk is
    quiet if |dE| / norm <= tol_abs, else if the relative change < tol_rel. The first chunk (E_prev None) is never
    quiet."""
    if E_prev is None:
        return False, float("nan"), float("nan")
    rel = abs(E - E_prev) / abs(E) if E != 0 else float("nan")
    d_abs = abs(E - E_prev) / norm
    return (d_abs <= tol_abs) if tol_abs is not None else (rel < tol_rel), rel, d_abs


SNAP_RE = re.compile(r"m_(.+)_c(\d+)_ph\d+\.vti$")


def clean_snapshots(out, stage, keep_from):
    """PLAN 4.13 S6: delete the snapshots of a stage with a cycle index < keep_from (keep_from None: all of them).
    Returns the deleted names."""
    gone = []
    for f in pathlib.Path(out).glob("m_%s_c*_ph*.vti" % stage):
        m = SNAP_RE.match(f.name)
        if m and m.group(1) == stage and (keep_from is None or int(m.group(2)) < keep_from):
            f.unlink(missing_ok=True)
            gone.append(f.name)
    return sorted(gone)


def is_calm(c, tol):
    """Calm cycle: w_loop > 0 and |w_dis / w_loop - 1| <= tol (PLAN 4.5; the opposite of an event, PROTOKOLL 7.53)."""
    return c["w_loop"] > 0 and abs(c["w_dis"] / c["w_loop"] - 1.0) <= tol


def calm_tail(cycles, tol):
    """Number of calm cycles at the end of cycles (calm: w_loop > 0 and |w_dis / w_loop - 1| <= tol; PLAN 4.5).
    Equal to the number of cycles that the steady-tail rule of analyze.py keeps (PROTOKOLL 7.54)."""
    n = 0
    for c in reversed(cycles):
        if not is_calm(c, tol):
            break
        n += 1
    return n


class HarmAcc:
    """Per-cell accumulators with harmonics (--acc --acc_harm K, PROTOKOLL 7.72). All tensors are allocated once at
    the start of the job (the GPU memory is claimed before the queue places the next job).
    Per cycle: sums of the local dissipation p, of m, of m e^{-i k w t} (k = 1 in float64 as before, k = 2 … K in
    float32), of h e^{-i w t} and of p_dis. At the end of a cycle: |m_k|^2 = |2/n sum m e^{-i k w t}|^2 summed over
    the components (k = 1 … K) and the box sums over the particle cells (to the cycle record). Only a calm cycle is
    added to the stage sums; the stage sums of p, m0, m1, h1, p_dis are then those of acc_reduce.py over the calm
    cycles, and mk2 is the mean of |m_k|^2 over the calm cycles."""

    def __init__(self, shape, K, device, mask):
        import torch
        assert K >= 2
        self.K, self.mask = K, mask
        f64, f32 = dict(dtype=torch.float64, device=device), dict(dtype=torch.float32, device=device)
        self.c = {"p": torch.zeros(shape[:3], **f64), "m0": torch.zeros(shape, **f64),
                  "re": torch.zeros(shape, **f64), "im": torch.zeros(shape, **f64)}
        self.c_re = torch.zeros((K - 1,) + tuple(shape), **f32)
        self.c_im = torch.zeros((K - 1,) + tuple(shape), **f32)
        self.s = {"p": torch.zeros(shape[:3], **f64), "m0": torch.zeros(shape, **f64),
                  "re": torch.zeros(shape, **f64), "im": torch.zeros(shape, **f64)}
        self.s_mk2 = torch.zeros((K,) + tuple(shape[:3]), **f32)
        self.new_cycle()
        self.new_stage()

    def new_cycle(self):
        for v in list(self.c.values()) + [self.c_re, self.c_im]:
            v.zero_()
        self.c_h1, self.c_psum, self.c_n = 0j, 0.0, 0

    def new_stage(self):
        for v in list(self.s.values()) + [self.s_mk2]:
            v.zero_()
        self.s_h1, self.s_psum, self.s_n, self.kept, self.n_cycles = 0j, 0.0, 0, [], 0

    def add(self, m, p_cell, th, h, p_dis):
        """One sample: m (float64), p_cell (float64), drive phase th = w (t - t0), drive value h, box p_dis."""
        c, s = math.cos(th), math.sin(th)
        self.c["p"] += p_cell
        self.c["m0"] += m
        self.c["re"].add_(m, alpha=c)
        self.c["im"].add_(m, alpha=-s)
        m32 = m.float()
        for j, k in enumerate(range(2, self.K + 1)):
            self.c_re[j].add_(m32, alpha=math.cos(k * th))
            self.c_im[j].add_(m32, alpha=-math.sin(k * th))
        self.c_h1 += h * complex(c, -s)
        self.c_psum += p_dis
        self.c_n += 1

    def cycle_mk2(self):
        """|m_k|^2 of the current cycle, k = 1 … K, shape (K, N, N, N), float32."""
        import torch
        g = (2.0 / self.c_n) ** 2
        out = torch.empty_like(self.s_mk2)
        out[0] = ((self.c["re"] ** 2 + self.c["im"] ** 2).sum(-1) * g).float()
        for j in range(self.K - 1):
            out[j + 1] = (self.c_re[j] ** 2 + self.c_im[j] ** 2).sum(-1) * g
        return out

    def end_cycle(self, calm, ci):
        """Close the cycle: box sums for the cycle record; a calm cycle goes into the stage sums."""
        mk2 = self.cycle_mk2()
        n = self.c_n
        box = {"calm": bool(calm), "n": n, "p_cells": float(self.c["p"][self.mask].sum()) / n,
               "mk2_cells": [float(mk2[k][self.mask].double().sum()) for k in range(self.K)]}
        self.n_cycles += 1
        if calm:
            for key in self.s:
                self.s[key] += self.c[key]
            self.s_mk2 += mk2
            self.s_h1 += self.c_h1
            self.s_psum += self.c_psum
            self.s_n += n
            self.kept.append(int(ci))
        self.new_cycle()
        return box


def measured_cycles(rec):
    """The cycles that analyze.py measures: from 1 on; from 2 on after a drive correction; from 3 on after a second
    correction (analyze.first_measured)."""
    from analyze import first_measured
    return rec["cycles"][first_measured(rec):]


def acc_first(rec):
    """First cycle that the accumulators keep: ACC_FROM_CYCLE, one later after a second drive correction."""
    return ACC_FROM_CYCLE + (1 if "h_amp_second" in rec else 0)


def second_correction(h, b_target, b_peak, tol):
    """New drive amplitude after cycle 1 (the first cycle at the drive corrected with cycle 0), or None if
    |b_peak / b_target - 1| <= tol or tol <= 0. Factor limited to 0.67 ... 1.5 as the first correction."""
    if tol <= 0 or abs(b_peak / b_target - 1.0) <= tol:
        return None
    return h * min(1.5, max(1.0 / 1.5, b_target / b_peak))


def stop_stage(rec, ci, n_cycles, args):
    """Stop decision of a measured stage after ci cycles (ci >= n_cycles). Returns (stop, info)."""
    cyc = rec["cycles"][-1]
    closed = cyc["closure"] < args.closure_tol
    in_scatter = cyc.get("dW_rel", 1.0) < args.dW_tol
    info = {"steady": bool(closed and in_scatter), "closed": bool(closed), "in_scatter": bool(in_scatter)}
    extra = 1 if "h_amp_second" in rec else 0                 # one more discarded cycle after a second correction
    max_cyc = max(args.max_cycles_per_amp + extra, n_cycles)
    if args.min_calm > 0:
        ct = calm_tail(measured_cycles(rec), args.calm_tol)
        info.update(calm_tail=ct, calm_ok=bool(ct >= args.min_calm))
        return ct >= args.min_calm or ci >= max_cyc, info
    return info["steady"] or ci >= args.max_cycles_per_amp + extra, info


def is_multiple(x, q, tol=1e-9):
    return abs(x / q - round(x / q)) < tol


def build_stages(args, phi):
    h_sat = args.h_sat if args.h_sat is not None else (1.0 - phi)
    asc = args.protocol in ("virgin_asc", "acdemag_asc")
    if args.protocol == "reset_desc":
        stages = [{"name": "sat", "h_amp": h_sat, "n_cycles": args.n_sat, "measure": False}]
    elif args.protocol == "acdemag_asc":
        if args.demag_h_end is not None:
            h_end = args.demag_h_end
        elif args.b_list:
            h_end = 0.3 * min(args.b_list) / args.mu_guess
        else:
            h_end = 0.3 * min(args.h_list)
        n_half = int(math.ceil(math.log(h_sat / h_end) / math.log(1.0 / args.demag_q)))
        stages = [{"name": "demag", "h_amp": h_sat, "decay": args.demag_q, "n_cycles": (n_half + 1) // 2,
                   "f_factor": args.demag_f_factor, "measure": False}]
    else:
        stages = []
    if args.b_list:
        bl = sorted(args.b_list, reverse=not asc)
        stages += [{"name": "amp%02d" % k, "h_amp": None, "b_target": b, "n_cycles": args.cycles_per_amp,
                    "measure": True} for k, b in enumerate(bl)]
        return stages
    if args.h_list:
        amps = sorted(args.h_list, reverse=not asc)
    else:
        h_max = args.h_max if args.h_max is not None else (1.0 - phi) / 3.0
        amps = [h_max * args.h_factor**k for k in range(args.n_amp)]
        if asc:
            amps = amps[::-1]
    stages += [{"name": "amp%02d" % k, "h_amp": h, "n_cycles": args.cycles_per_amp, "measure": True}
               for k, h in enumerate(amps)]
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
                          ExchangeField, UniaxialAnisotropyField,
                          LLGSolver, accumulate_h, write_vti, normalize)
    import torch
    set_log_level(25)
    set_precision(run_precision(args))           # --precision_cycles with --init_from (no relaxation in this run)
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
    rep = geometry_rep(args)
    geo = locked_geometry(N, d, dx, rep, sub=8 if args.surface == "fraction" else 0)
    ids = geo["ids"]
    mag = ids >= 0
    phi = geo["phi_vox"]

    pbc, strayfield = strayfield_setup(args)
    mesh = Mesh((N, N, N), (dx, dx, dx), pbc=pbc)
    state = State(mesh)
    dt_ = state.dtype
    dev = state.device

    def T(a):
        return torch.tensor(a, dtype=dt_, device=dev)

    mag4 = mag[..., None]
    frac4 = geo["frac"].astype(np.float64)[..., None]     # 1 inside, f at the surface, 0 in the void
    state.material = {"alpha": args.alpha}
    state.material["Ms"] = T(np.where(mag4, np.maximum(frac4 * Ms, args.ms_void * Ms), args.ms_void * Ms))
    state.material["A"] = T(np.where(mag4, frac4 * args.A, 0.0))

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

    terms = [strayfield, ExchangeField()]
    term_names = ["demag", "exchange"]                    # PLAN 4.13 S5: names of the energy terms, same order
    n_cubes = [0, 0, 0, 0]
    if K_eff > 0.0:
        ax_c, n_cubes = locked_cube_axes(ids, rep, args.dx_lex, a_lex, d_lex, L, args.seed)
        state.material["Ku"] = T(np.where(mag4, frac4 * K_eff, 0.0))
        state.material["Ku_axis"] = T(ax_c)
        terms.append(UniaxialAnisotropyField())
        term_names.append("aniso_cube")
        del ax_c
    ax_p = particle_axes(e_d)
    if K_p > 0.0:
        kp_ax = np.zeros((N, N, N, 3))
        for p in range(4):
            kp_ax[ids == p] = ax_p[p]
        state.material["Kp"] = T(np.where(mag4, frac4 * K_p, 0.0))
        state.material["Kp_axis"] = T(kp_ax)
        terms.append(UniaxialAnisotropyField(Ku="Kp", Ku_axis="Kp_axis"))
        term_names.append("aniso_particle")
        del kp_ax
    aniso_info = {"Q_eff": Q_eff, "K_eff": K_eff, "L_eff": L * lex, "L_over_dx": L / args.dx_lex,
                  "K_p": K_p, "Q_p": args.r_p * Q_eff, "cos_theta_p": np.abs(ax_p @ e_d).tolist(),
                  "delta_min_lex": delta_min, "delta_min_over_dx": delta_min / args.dx_lex,
                  "n_cubes_per_particle": [int(n) for n in n_cubes]}

    if args.snap_phases > 0:
        geo_f = {"Ms": state.material["Ms"], "particle": T(ids.astype(np.float64)[..., None])}
        if "Ku_axis" in state.material:
            geo_f["cube_axis"] = state.material["Ku_axis"]
        write_vti(geo_f, str(out / "geometry.vti"), state)
    drive = SinusoidalDrive(state, tuple(e_d))
    terms.append(drive)
    term_names.append("drive")
    llg = LLGSolver(terms, atol=args.atol)

    def term_energies():
        """PLAN 4.13 S5: energy of each field term / (K_d V), fp64 sums (the swapped stray-field term included)."""
        with torch.no_grad():
            return {nm: float(tm.E(state)) / (Kd * N**3 * dx**3) for nm, tm in zip(term_names, terms)}

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

    last_p = {}                                  # per-cell dissipation of the last diagnostics() call (--acc)
    harm = (HarmAcc((N, N, N, 3), args.acc_harm, dev, torch.as_tensor(ids >= 0, device=dev))
            if args.acc and args.acc_harm > 0 else None)     # --acc_harm (PROTOKOLL 7.72): allocated once

    @torch.no_grad()
    def diagnostics():
        m = state.m
        Mvec = (Ms_t * m).double().mean(dim=(0, 1, 2))
        Mpar = float(torch.dot(Mvec, e_t))
        m_par = (m.double() * e_t).sum(-1).reshape(-1) * Ms_t.double().reshape(-1)
        Mp = torch.bincount(ids_t, weights=m_par, minlength=5)[1:] / counts
        H = accumulate_h(term.h(state) for term in terms)
        mxH = torch.linalg.cross(m, H)
        p_cell = p_coef * Ms_t * (mxH * mxH).sum(-1, keepdim=True)
        last_p["p"] = p_cell
        p_dis = float(p_cell.double().mean())
        cmin, n60 = 1.0, 0
        for ax in range(3):
            c = (m * torch.roll(m, 1, dims=ax)).sum(-1)[pair_masks[ax]]
            cmin = min(cmin, float(c.min()))
            n60 += int((c < cos60).sum())
        return Mpar, Mvec.cpu().numpy(), Mp.cpu().numpy(), p_dis, cmin, n60

    cube_ok = bloch_mask_t(ids_n)

    @torch.no_grad()
    def bloch(positions=False):
        """Bloch points of the present state (Chris 2026-10-18): (n with Q = +1, n with Q = -1) of the cubes inside one
        particle; with positions also [[i, j, k, Q, particle], ...]."""
        q = cube_charge_t(state.m)
        sel = cube_ok & (q.abs() > 0.5)
        npos, nneg = int((sel & (q > 0)).sum()), int((sel & (q < 0)).sum())
        if not positions:
            return npos, nneg
        pts = [[int(i), int(j), int(k), int(round(float(q[i, j, k]))), int(ids[i, j, k])]
               for i, j, k in torch.nonzero(sel).cpu().numpy()]
        return npos, nneg, pts

    def relax_converge(tag):
        """The converge relaxation of virgin_asc (alpha 1, chunks, quiet_chunk; PLAN 4.13 S7). Returns (history,
        iterations, converged); the history has the Bloch points per chunk."""
        magf = torch.tensor(mag.reshape(-1), device=dev)
        hist, E_prev, it, converged, quiet = [], None, 0, False, 0
        while it < args.relax_maxiter:
            for _ in range(args.relax_chunk):
                state.t, state.m = llg._solver.step(state.t, state.m, 1e-11, state=state, alpha=1.0)
            normalize(state.m)
            it += args.relax_chunk
            with torch.no_grad():
                tq = (llg.dm(state.t, state.m, state=state, alpha=1.0).double().norm(dim=-1).reshape(-1)[magf]
                      / constants.gamma)
                E = float(llg.E(state))
            q, rel, d_abs = quiet_chunk(E, E_prev, Kd * N**3 * dx**3, args.relax_tol_E, args.relax_tol_E_abs)
            nb = bloch()
            hist.append({"it": it, "E_over_Kd": E / (Kd * N**3 * dx**3), "rel_dE": rel, "dE_over_KdV": d_abs,
                         "torque_rms": float(tq.pow(2).mean().sqrt()), "torque_max": float(tq.max()),
                         "n_bloch": list(nb)})
            print("[%s] it %6d  E/K_d %.10g  rel dE %.2e  |dE|/(K_d V) %.2e  torque rms %.3g  max %.3g A/m  "
                  "Bloch points +%d -%d" % (tag, it, hist[-1]["E_over_Kd"], rel, d_abs, hist[-1]["torque_rms"],
                                            hist[-1]["torque_max"], nb[0], nb[1]), flush=True)
            quiet = quiet + 1 if q else 0
            if quiet >= args.relax_confirm:
                converged = True
                break
            E_prev = E
        return hist, it, converged

    def predemag():
        """AC demagnetization after the virgin relaxation (Chris 2026-10-18), then the relaxation again. Returns the
        record for summary["init"]["predemag"]; writes bloch_predemag.csv (Bloch points at each sample)."""
        h0 = args.predemag_h if args.predemag_h > 0 else (args.h_sat if args.h_sat is not None else 1.0 - args.phi)
        h_end = 0.3 * min(args.b_list) / args.mu_guess if args.b_list else 0.01 * h0
        n_half = max(1, int(math.ceil(math.log(h0 / h_end) / math.log(1.0 / args.predemag_q))))
        dt_d = 1.0 / freq / args.samples
        state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)
        drive.t0, drive.freq, drive.H_amp = 0.0, freq, h0 * Ms
        t_w, halves = time.time(), []
        with open(out / "bloch_predemag.csv", "w") as fb:
            fb.write("# half_cycle,sample,h,M_par_over_Ms,n_bloch_pos,n_bloch_neg\n")
            for hh in range(n_half):
                mx, nbmax = 0.0, 0
                for k in range(args.samples // 2):
                    llg.step(state, dt_d)
                    with torch.no_grad():
                        mp = float(torch.dot((Ms_t * state.m).double().mean(dim=(0, 1, 2)), e_t)) / Ms
                    nb = bloch()
                    mx, nbmax = max(mx, abs(mp)), max(nbmax, nb[0] + nb[1])
                    fb.write("%d,%d,%.6e,%.6e,%d,%d\n" % (hh, k, drive.H_amp / Ms, mp, nb[0], nb[1]))
                halves.append({"h": drive.H_amp / Ms, "M_par_peak_over_Ms": mx, "n_bloch_max": nbmax,
                               "n_bloch_end": list(nb)})
                print("[predemag] half cycle %d: h %.4g  |M_par|/Ms peak %.4f  Bloch points max %d, end +%d -%d" %
                      (hh, drive.H_amp / Ms, mx, nbmax, nb[0], nb[1]), flush=True)
                drive.H_amp *= args.predemag_q
        drive.H_amp = 0.0
        rec = {"h0": h0, "q": args.predemag_q, "h_end": h_end, "n_half": n_half, "halves": halves,
               "wall_s": time.time() - t_w}
        hist, it, conv = relax_converge("relax2")
        Mpar0, Mv0, _, _, cmin0, n60_0 = diagnostics()
        npos, nneg, pts = bloch(positions=True)
        rec.update({"relax_history": hist, "relax_iterations": it, "relax_converged": conv,
                    "E_over_Kd": hist[-1]["E_over_Kd"], "M_over_Ms": [float(x) / Ms for x in Mv0],
                    "n_pairs_gt60": int(n60_0), "bloch_points": pts})
        print("[predemag] after the relaxation: converged=%s after %d iterations, M/Ms (%.4f, %.4f, %.4f), Bloch "
              "points +%d -%d" % (conv, it, *rec["M_over_Ms"], npos, nneg), flush=True)
        state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)
        return rec

    # --- config ------------------------------------------------------------
    stages = build_stages(args, args.phi)        # nominal phi: same amplitudes on every mesh
    if args.cycles_last_amp:
        [s for s in stages if s["measure"]][-1]["n_cycles"] = args.cycles_last_amp
    import magnumnp
    config = dict(vars(args))
    config.update({"magnumnp_path": str(pathlib.Path(magnumnp.__file__).resolve().parent), "units": U, "N_used": N, "n_cells": N**3, "phi_vox": phi, "phi_nom": geo["phi_nom"],
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
        diff = [k for k in vars(args) if k not in RUN_CONTROL and old.get(k, NEW_KEYS.get(k)) != config[k]]
        if diff:
            sys.exit("ERROR: resume with other arguments than the checkpoint run: %s. "
                     "Use --no_resume or another --out." % ", ".join("%s %r -> %r" % (k, old.get(k), config[k]) for k in diff))
        if old.get("git_commit") != config["git_commit"] and not args.allow_new_commit:
            sys.exit("ERROR: the code changed since the checkpoint (%s -> %s). One run would mix two code "
                     "versions. Use --no_resume, or --allow_new_commit if the change does not touch the physics."
                     % (old.get("git_commit"), config["git_commit"]))
    init_file = out / "init.pt"
    if args.no_resume:
        for f in [ckpt, summ_file, init_file, out / "DONE"] + list(out.glob("samples_*.csv")) + list(out.glob("ckpt_*.pt")):
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
    elif args.init_stage:
        if args.init_from:
            sys.exit("ERROR: --init_stage and --init_from exclude each other")
        src = pathlib.Path(args.init_stage)
        c, summary, first_stage = stage_start(src if src.is_absolute() else out.parent / src, config)
        state.m = c["m"].to(dtype=dt_, device=dev)
        state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)
        print("[init] m from %s, start at stage %d" % (args.init_stage, first_stage), flush=True)
    elif args.init_from or init_file.exists():
        # relaxed virgin state of another run (--init_from) or of this run before an interruption
        src = init_file
        if args.init_from and not init_file.exists():
            src = pathlib.Path(args.init_from)
            src = src if src.is_absolute() else out.parent / src
            t_end = time.time() + 3600.0 * args.init_wait_h
            while not src.exists() and time.time() < t_end:
                time.sleep(60.0)
            if not src.exists():
                sys.exit("ERROR: --init_from %s does not exist" % src)
        src = str(src)
        c = torch.load(src, map_location=dev)
        late_predemag = bool(args.predemag_h is not None and "init" in c
                             and c["config"].get("predemag_h", NEW_KEYS["predemag_h"]) is None)
        if "init" in c:
            diff = [k for k in INIT_KEYS if c["config"].get(k, NEW_KEYS.get(k)) != config.get(k)
                    and not (k == "relax_maxiter" and c["init"].get("relax_converged"))
                    and not (late_predemag and k in ("predemag_h", "predemag_q"))]
            if diff:
                sys.exit("ERROR: %s was relaxed with other arguments: %s" % (src, ", ".join(
                    "%s %r -> %r" % (k, c["config"].get(k, NEW_KEYS.get(k)), config.get(k)) for k in diff)))
            summary["init"] = dict(c["init"], source=src)
        elif args.protocol == "virgin_asc":
            sys.exit("ERROR: %s holds no relaxation record (old checkpoint format)" % src)
        state.m = c["m"].to(dtype=dt_, device=dev)
        state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)
        print("[init] m from %s" % src, flush=True)
        if late_predemag:
            # Chris 2026-10-18: the virgin state of another run, then the AC demagnetization and the relaxation here
            summary["init"]["bloch_virgin"] = list(bloch())
            summary["init"]["predemag"] = predemag()
            summary["init"]["relax_converged"] = bool(c["init"].get("relax_converged")
                                                      and summary["init"]["predemag"]["relax_converged"])
            save_atomic({"m": state.m.detach().cpu(), "init": summary["init"],
                         "config": {k: config.get(k) for k in INIT_KEYS}}, init_file)
    else:
        init_seed = args.seed if args.init_seed is None else args.init_seed
        if args.init_block > 0.0 and not (is_multiple(a_lex, args.init_block) and is_multiple(args.init_block, args.dx_lex)):
            warnings.warn("init_block %g l_ex is not a multiple of dx and a divisor of a: "
                          "the initial state depends on the mesh" % args.init_block)
        transfer = None
        if args.init_transfer:
            # PLAN 4.10: the relaxed state of a finer mesh, block averaged onto this mesh, as the start of the relaxation
            src = pathlib.Path(args.init_transfer)
            src = src if src.is_absolute() else out.parent / src
            t_end = time.time() + 3600.0 * args.init_wait_h
            while not src.exists() and time.time() < t_end:
                time.sleep(60.0)
            if not src.exists():
                sys.exit("ERROR: --init_transfer %s does not exist" % src)
            cf = torch.load(str(src), map_location="cpu")
            if args.surface != "voxel":
                sys.exit("ERROR: --init_transfer needs --surface voxel")
            diff = [k for k in INIT_KEYS if k not in ("dx_lex", "relax_maxiter")
                    and cf["config"].get(k, NEW_KEYS.get(k)) != config.get(k)]
            if diff:
                sys.exit("ERROR: %s was relaxed with other arguments: %s" % (src, ", ".join(
                    "%s %r -> %r" % (k, cf["config"].get(k, NEW_KEYS.get(k)), config.get(k)) for k in diff)))
            if not cf.get("init", {}).get("relax_converged"):
                sys.exit("ERROR: the relaxation of %s did not converge" % src)
            dxf = float(cf["config"]["dx_lex"])
            if not (is_multiple(a_lex, dxf) and is_multiple(args.dx_lex, dxf)):
                sys.exit("ERROR: dx %g l_ex of %s is not a divisor of dx %g and of a %g l_ex" % (dxf, src, args.dx_lex, a_lex))
            Nf = int(round(a_lex / dxf))
            rep_f = int(round(args.geom_dx_lex / dxf)) if args.geom_dx_lex else 1      # the same locked geometry
            ids_f = locked_geometry(Nf, d, dxf * lex, rep_f)["ids"]
            m_c, tstats = transfer_m(cf["m"].numpy(), ids_f, ids)
            state.m = T(m_c)
            transfer = dict(tstats, source=str(src), dx_lex_fine=dxf)
            print("[init] transferred from %s (dx %g -> %g l_ex): %s" % (src, dxf, args.dx_lex, tstats), flush=True)
        else:
            state.m = T(virgin_m0(N, a_lex, args.dx_lex, args.init_block, init_seed, args.init_mode, ids))
        state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)
        drive.H_amp = 0.0
        t_w = time.time()
        if args.protocol == "virgin_asc" and args.relax_mode == "converge":
            # T = 0 analogue of cooling in zero field: random m, relaxed at alpha = 1 in chunks until the
            # energy does not change any more (relative change per chunk < relax_tol_E, or with relax_tol_E_abs an
            # absolute change per chunk <= relax_tol_E_abs, PLAN 4.13 S7)
            tries = []
            for attempt in range(max(1, args.bloch_free_tries)):
                if attempt > 0:
                    # Chris 2026-10-18: a new random start of the same microstructure, the preparation again
                    state.m = T(virgin_m0(N, a_lex, args.dx_lex, args.init_block, init_seed, args.init_mode, ids,
                                          attempt))
                    state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)
                    drive.H_amp = 0.0
                hist, it, converged = relax_converge("relax" if attempt == 0 else "relax_try%d" % attempt)
                Mpar0, Mv0, _, _, cmin0, n60_0 = diagnostics()
                summary["init"] = {"protocol": args.protocol, "relax_mode": "converge", "relax_converged": converged,
                                   "relax_iterations": it, "relax_history": hist,
                                   "E_over_Kd": hist[-1]["E_over_Kd"],
                                   "M_over_Ms": [float(x) / Ms for x in Mv0], "n_pairs_gt60": int(n60_0),
                                   "max_angle_deg": float(math.degrees(math.acos(max(-1.0, min(1.0, cmin0))))),
                                   "wall_s": time.time() - t_w}
                print("[relax] virgin state: converged=%s after %d iterations, %.0f s wall" %
                      (converged, it, summary["init"]["wall_s"]), flush=True)
                summary["init"]["bloch_virgin"] = list(bloch())
                n_end = sum(summary["init"]["bloch_virgin"])
                if args.predemag_h is not None:
                    summary["init"]["predemag"] = predemag()
                    summary["init"]["relax_converged"] = bool(converged
                                                              and summary["init"]["predemag"]["relax_converged"])
                    n_end = len(summary["init"]["predemag"]["bloch_points"])
                tries.append({"attempt": attempt, "relax_iterations": it, "converged": converged,
                              "n_bloch_prepared": n_end})
                if args.bloch_free_tries > 0:
                    print("[bloch] attempt %d: %d Bloch points in the prepared state" % (attempt, n_end), flush=True)
                if args.bloch_free_tries <= 0 or n_end == 0 or not summary["init"]["relax_converged"]:
                    break
            summary["init"]["bloch_tries"] = tries
            summary["init"]["bloch_free"] = bool(tries[-1]["n_bloch_prepared"] == 0)
        elif args.protocol == "virgin_asc":
            # T = 0 analogue of cooling in zero field: random m, relaxed at alpha = 1 until |dm/dt| is small
            converged = bool(llg.relax(state, maxiter=args.relax_maxiter, dm_tol=args.relax_dm_tol, dt=1e-11))
            Mpar0, Mv0, _, _, cmin0, n60_0 = diagnostics()
            summary["init"] = {"protocol": args.protocol, "relax_converged": converged,
                               "E_over_Kd": float(llg.E(state)) / (Kd * N**3 * dx**3),
                               "M_over_Ms": [float(x) / Ms for x in Mv0], "n_pairs_gt60": int(n60_0),
                               "max_angle_deg": float(math.degrees(math.acos(max(-1.0, min(1.0, cmin0))))),
                               "wall_s": time.time() - t_w}
            print("[relax] virgin state: converged=%s  E/K_d=%.5g  M/Ms=(%.3g, %.3g, %.3g)  n60=%d  %.0f s wall" %
                  (converged, summary["init"]["E_over_Kd"], *summary["init"]["M_over_Ms"], n60_0,
                   summary["init"]["wall_s"]), flush=True)
        else:
            llg.step(state, args.relax_tau * U["t_M"])
            print("[relax] %.3g s done in %.1f s wall" % (args.relax_tau * U["t_M"], time.time() - t_w), flush=True)
        if transfer is not None and "init" in summary:
            summary["init"]["transfer"] = transfer
        if args.protocol == "virgin_asc":
            save_atomic({"m": state.m.detach().cpu(), "init": summary["init"],
                         "config": {k: config.get(k) for k in INIT_KEYS}}, init_file)

    # hard gate: the protocol requires a converged virgin state
    if (args.protocol == "virgin_asc" and not summary.get("init", {}).get("relax_converged", True)
            and not args.allow_unconverged):
        write_atomic(summ_file, json.dumps(summary, indent=1))
        print("ERROR: the virgin relaxation did not converge (%s). Raise --relax_maxiter and start "
              "again with --no_resume." % out, flush=True)
        sys.exit(3)

    # --- protocol ----------------------------------------------------------
    kw_cyc = cycle_tolerances(args)
    sf_cyc = cycle_strayfield(args)
    if sf_cyc is not None:
        terms[terms.index(strayfield)] = sf_cyc   # the diagnostics (p_dis) use the same list
        strayfield = sf_cyc
        print("[cycles] strayfield %s" % args.strayfield_cycles, flush=True)
    if kw_cyc is not None or sf_cyc is not None:
        llg = LLGSolver(terms, **(kw_cyc or {"atol": args.atol}))   # p_dis and the diagnostics use the terms
        if kw_cyc is not None:
            print("[cycles] RKF45 %s" % " ".join("%s %g" % kv for kv in sorted(kw_cyc.items())), flush=True)
    cycles_done = sum(len(s["cycles"]) for s in summary["stages"])

    for si in range(first_stage, len(stages)):
        st = stages[si]
        if st["n_cycles"] <= 0:
            summary["stages"].append({"name": st["name"], "h_amp": st["h_amp"], "measure": st["measure"],
                                      "cycles": [], "complete": True})
            continue
        state.t = torch.tensor(0.0, dtype=torch.float64, device=dev)
        drive.t0 = 0.0
        f_stage = freq * st.get("f_factor", 1.0)          # demagnetizing cycles may run faster
        drive.freq = f_stage
        period = 1.0 / f_stage
        dt_s = period / args.samples
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
        if args.snap_phases > 0:
            clean_snapshots(out, st["name"], None)        # PLAN 4.13 S6: files of an earlier attempt of this stage
        fcsv = open(out / ("samples_%02d_%s.csv" % (si, st["name"])), "w")
        if harm is not None and st["measure"]:
            harm.new_stage()
            harm.new_cycle()
        elif args.acc and st["measure"]:
            acc = {"p": torch.zeros(state.m.shape[:3], dtype=torch.float64, device=dev),
                   "m0": torch.zeros(state.m.shape, dtype=torch.float64, device=dev),
                   "re": torch.zeros(state.m.shape, dtype=torch.float64, device=dev),
                   "im": torch.zeros(state.m.shape, dtype=torch.float64, device=dev),
                   "h1": 0j, "p_sum": 0.0, "n": 0}
        fcsv.write("# cycle,t,H,M_par,Mx,My,Mz,B,p_dis,Mp0,Mp1,Mp2,Mp3,n_pairs_gt60  (SI)\n")
        fbl = open(out / ("bloch_%02d_%s.csv" % (si, st["name"])), "w")
        fbl.write("# cycle,sample,n_bloch_pos,n_bloch_neg  (Bloch points at each sample, Chris 2026-10-18)\n")

        ci = 0
        complete = False
        while True:
            if args.max_cycles is not None and cycles_done >= args.max_cycles:
                break
            t_w = time.time()
            steps0 = state._step
            rows, brows = [], []
            cmin_c, n60_c = 1.0, 0

            def sample():
                nonlocal cmin_c, n60_c
                t = float(state.t)
                Mpar, Mv, Mp, pd, cmin, n60 = diagnostics()
                cmin_c, n60_c = min(cmin_c, cmin), max(n60_c, n60)
                rows.append((t, drive.value(t), Mpar, *Mv, pd, *Mp, n60))
                brows.append(bloch())                    # Chris 2026-10-18: Bloch points at each sample

            e_start = term_energies()                     # PLAN 4.13 S5; the drive field is 0 at the cycle start
            sample()
            acc_on = args.acc and st["measure"] and ci >= acc_first(rec)
            snap = args.snap_phases > 0 and st["measure"] and ci >= st["n_cycles"] - args.snap_cycles
            if snap and args.snap_mT and st.get("b_target") is not None:
                bmT = 1000.0 * st["b_target"] * args.Js
                snap = any(abs(bmT / x - 1.0) < 0.1 for x in args.snap_mT)
            for k in range(args.samples):
                llg.step(state, dt_s)
                if st.get("decay") and (k + 1 == args.samples // 2 or k + 1 == args.samples):
                    drive.H_amp *= st["decay"]          # AC demagnetization: change at a zero crossing
                sample()
                if acc_on and harm is not None:
                    th = 2.0 * math.pi * f_stage * (float(state.t) - drive.t0)
                    harm.add(state.m.double(), last_p["p"][..., 0].double(), th, drive.value(float(state.t)),
                             rows[-1][6])
                elif acc_on:
                    # audit 3, item 9b: the samples after each step cover one cycle at equal phases
                    th = 2.0 * math.pi * f_stage * (float(state.t) - drive.t0)
                    c, sn = math.cos(th), math.sin(th)
                    md = state.m.double()
                    acc["p"] += last_p["p"][..., 0].double()
                    acc["m0"] += md
                    acc["re"] += md * c
                    acc["im"] -= md * sn
                    acc["h1"] += drive.value(float(state.t)) * complex(c, -sn)
                    acc["p_sum"] += rows[-1][6]
                    acc["n"] += 1
                if snap and (k + 1) % (args.samples // args.snap_phases) == 0:
                    deg = round(360.0 * (k + 1) / args.samples)
                    write_vti({"m": state.m}, str(out / ("m_%s_c%d_ph%03d.vti" % (st["name"], ci, deg))), state)
            wall = time.time() - t_w
            steps = state._step - steps0
            if not bool(torch.isfinite(state.m).all()):
                # audit 3, item 4: a NaN state must not block a GPU; exit code 4 (run_queue: FAIL, one retry)
                fcsv.close()
                print("ERROR: [%s] m is not finite after cycle %d" % (st["name"], ci), flush=True)
                sys.exit(4)

            a = np.array(rows)
            t_a, H_a, M_a, p_a = a[:, 0], a[:, 1], a[:, 2], a[:, 6]
            B_a = MU0 * (H_a + M_a)
            W_loop = float(_trapz(H_a, B_a))
            W_dis = float(_trapz(p_a, t_a))
            M_peak = float(0.5 * (M_a.max() - M_a.min()))
            M_mid = float(0.5 * (M_a.max() + M_a.min()))
            B_peak = float(0.5 * (B_a.max() - B_a.min()))
            cyc = {"cycle": ci,
                   "B_peak": B_peak, "b_peak": B_peak / args.Js,
                   "H_peak": float(0.5 * (H_a.max() - H_a.min())), "M_peak": M_peak,
                   "W_loop": W_loop, "W_dis": W_dis, "w_loop": W_loop / Kd, "w_dis": W_dis / Kd,
                   "closure": float(abs(M_a[-1] - M_a[0]) / max(M_peak, 1e-30)),
                   "M_mid": M_mid, "offset": M_mid / max(M_peak, 1e-30),
                   "max_angle_deg": float(math.degrees(math.acos(max(-1.0, min(1.0, cmin_c))))),
                   "n_pairs_gt60": int(n60_c), "frac_pairs_gt60": n60_c / max(n_pairs, 1),
                   "steps": int(steps), "wall_s": wall, "mean_dt": period / max(steps, 1), "E_start": e_start}
            if ci > 0:
                W_prev = rec["cycles"][-1]["W_loop"]
                cyc["dW_rel"] = float(abs(W_loop - W_prev) / max(abs(W_loop), 1e-300))
            rec["cycles"].append(cyc)
            if acc_on and harm is not None:
                cyc["acc"] = harm.end_cycle(is_calm(cyc, args.calm_tol), ci)
            for r in rows:
                fcsv.write("%d," % ci + ",".join("%.9e" % v for v in (*r[:3], *r[3:6], MU0 * (r[1] + r[2]), *r[6:])) + "\n")
            fcsv.flush()
            nbs = [a + b for a, b in brows]
            cyc["bloch"] = {"start": list(brows[0]), "end": list(brows[-1]), "max": max(nbs),
                            "mean": float(np.mean(nbs))}
            for k, (a, b) in enumerate(brows):
                fbl.write("%d,%d,%d,%d\n" % (ci, k, a, b))
            fbl.flush()
            if snap:
                clean_snapshots(out, st["name"], ci - args.snap_cycles + 1)   # PLAN 4.13 S6: keep the last cycles
            cycles_done += 1
            print("[%s c%d] h=%.4g  b_peak=%.4g  w_loop=%.4g  w_dis=%.4g  closure=%.2g  dW=%s  "
                  "max_angle=%.0f deg  steps=%d  dt=%.3g ps  wall=%.0f s" %
                  (st["name"], ci, st["h_amp"], cyc["b_peak"], cyc["w_loop"], cyc["w_dis"], cyc["closure"],
                   ("%.2g" % cyc["dW_rel"]) if "dW_rel" in cyc else "-", cyc["max_angle_deg"],
                   steps, 1e12 * cyc["mean_dt"], wall), flush=True)
            ci += 1

            if ci == 1 and st.get("b_target") is not None and st["n_cycles"] > 1:
                # correct the drive with b/h of the first cycle (it is discarded anyway);
                # the change happens at a zero crossing of the field
                # factor limited to 0.67 ... 1.5: cycle 0 is a transient (V100 data: b of cycle 0
                # within +-4 % of the later cycles), an outlier must not spoil the stage
                h_new = st["h_amp"] * min(1.5, max(1.0 / 1.5, st["b_target"] / cyc["b_peak"]))
                rec["h_amp_initial"] = st["h_amp"]
                st = dict(st, h_amp=h_new)
                drive.H_amp = h_new * Ms
                rec["h_amp"], rec["H_amp"] = h_new, h_new * Ms
                print("[%s] drive corrected: h %.4g -> %.4g (b target %.4g, b in cycle 0 %.4g)" %
                      (st["name"], rec["h_amp_initial"], h_new, st["b_target"], cyc["b_peak"]), flush=True)

            if ci == 2 and "h_amp_initial" in rec and "h_amp_second" not in rec and st.get("b_target") is not None:
                # cycle 0 is a transient from the start state, so its b can mislead the first correction (N0b_I15_s941:
                # b of cycle 0 6 % above the later cycles, the stage then 6 % below the target; PROTOKOLL 7.107).
                # Correct again with b of cycle 1, at a zero crossing; cycle 2 is then discarded too, and the stage
                # gets one more cycle
                h_new = second_correction(st["h_amp"], st["b_target"], cyc["b_peak"], args.b_tol)
                if h_new is not None:
                    rec["h_amp_second"] = st["h_amp"]
                    st = dict(st, h_amp=h_new, n_cycles=st["n_cycles"] + 1)
                    drive.H_amp = h_new * Ms
                    rec["h_amp"], rec["H_amp"] = h_new, h_new * Ms
                    print("[%s] drive corrected again: h %.4g -> %.4g (b target %.4g, b in cycle 1 %.4g)" %
                          (st["name"], rec["h_amp_second"], h_new, st["b_target"], cyc["b_peak"]), flush=True)

            if ci >= st["n_cycles"]:
                if not st["measure"]:
                    complete = True
                    break
                stop, info = stop_stage(rec, ci, st["n_cycles"], args)
                if stop:
                    rec.update(info)
                    complete = True
                    break
        fcsv.close()
        fbl.close()
        if complete:
            rec["E_end"] = term_energies()                # PLAN 4.13 S5: the end of the last cycle
            rec["bloch_points_end"] = bloch(positions=True)[2]

        if complete and st["measure"] and args.ringdown_periods > 0:
            # ring-down test: drive on for a quarter period (to H max), then H constant; record the decay of the
            # dissipation; m and t are restored, the solver step size adapts again
            m_save, t_save = state.m.detach().clone(), state.t.detach().clone()
            for _ in range(args.samples // 4):
                llg.step(state, dt_s)
            drive.hold = 1.0
            t_hold = float(state.t)
            dt_r = period / args.ringdown_samples
            with open(out / ("ringdown_%02d_%s.csv" % (si, st["name"])), "w") as frd:
                frd.write("# t_since_hold,H,M_par,p_dis,n_pairs_gt60  (SI)\n")
                for k in range(int(round(args.ringdown_periods * args.ringdown_samples)) + 1):
                    if k > 0:
                        llg.step(state, dt_r)
                    Mpar, _, _, pd, _, n60 = diagnostics()
                    frd.write("%.9e,%.9e,%.9e,%.9e,%d\n" % (float(state.t) - t_hold, drive.value(float(state.t)),
                                                          Mpar, pd, n60))
            drive.hold = None
            state.m, state.t = m_save, t_save
            rec["ringdown"] = "ringdown_%02d_%s.csv" % (si, st["name"])
            print("[%s] ring-down %.1f periods written" % (st["name"], args.ringdown_periods), flush=True)

        if complete and harm is not None and st["measure"]:
            n = harm.s_n
            info = {"n_samples": n, "cycles_kept": harm.kept, "cycles_acc": harm.n_cycles, "harm_K": harm.K}
            if n > 0:
                mask_t = harm.mask
                p_box = float(harm.s["p"].sum()) / n / harm.s["p"].numel()
                np_tmp = out / ("acc_%02d_%s.tmp.npz" % (si, st["name"]))
                np.savez(np_tmp, p=(harm.s["p"][mask_t] / n).float().cpu().numpy(),
                         m0=(harm.s["m0"][mask_t] / n).float().cpu().numpy(),
                         m1=((harm.s["re"][mask_t] + 1j * harm.s["im"][mask_t]) * (2.0 / n)).to(torch.complex64)
                         .cpu().numpy(),
                         h1=np.complex128(harm.s_h1 * 2.0 / n), n_samples=n, cycles_from=acc_first(rec),
                         p_box=p_box, p_dis_samples=harm.s_psum / n, freq=f_stage,
                         mk2=(harm.s_mk2[:, mask_t] / len(harm.kept)).T.contiguous().cpu().numpy(),
                         cycles_kept=np.array(harm.kept), cycles_acc=harm.n_cycles)
                replace_durable(np_tmp, out / ("acc_%02d_%s.npz" % (si, st["name"])))
                info.update(file="acc_%02d_%s.npz" % (si, st["name"]),
                            p_box_over_p_dis=p_box / max(harm.s_psum / n, 1e-300))
            else:
                info["file"] = None                   # no calm cycle: no accumulator file for this stage
            rec["acc"] = info
            print("[%s] accumulators: %d calm of %d cycles, %d samples%s" % (
                st["name"], len(harm.kept), harm.n_cycles, n,
                (", sum check p_box/p_dis = %.6f" % info["p_box_over_p_dis"]) if n > 0 else ""), flush=True)
        elif complete and args.acc and st["measure"] and acc["n"] > 0:
            n = acc["n"]
            mask_t = torch.as_tensor(ids >= 0, device=dev)
            p_box = float(acc["p"].sum()) / n / acc["p"].numel()
            np_tmp = out / ("acc_%02d_%s.tmp.npz" % (si, st["name"]))
            np.savez(np_tmp, p=(acc["p"][mask_t] / n).float().cpu().numpy(),
                     m0=(acc["m0"][mask_t] / n).float().cpu().numpy(),
                     m1=((acc["re"][mask_t] + 1j * acc["im"][mask_t]) * (2.0 / n)).to(torch.complex64).cpu().numpy(),
                     h1=np.complex128(acc["h1"] * 2.0 / n), n_samples=n, cycles_from=acc_first(rec),
                     p_box=p_box, p_dis_samples=acc["p_sum"] / n, freq=f_stage)
            replace_durable(np_tmp, out / ("acc_%02d_%s.npz" % (si, st["name"])))
            rec["acc"] = {"file": "acc_%02d_%s.npz" % (si, st["name"]), "n_samples": n,
                          "p_box_over_p_dis": p_box / max(acc["p_sum"] / n, 1e-300)}
            print("[%s] accumulators: %d samples, sum check p_box/p_dis = %.6f" %
                  (st["name"], n, rec["acc"]["p_box_over_p_dis"]), flush=True)
            del acc
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
