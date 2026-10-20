#!/usr/bin/env python
"""
Job lines of B0, step 1 of the mesh convergence series (PLAN 4.10, variant B): base microstructure (L_eff 12 l_ex,
r_p 0, φ 0.65, α 0.1, converged virgin relaxation) at d/l_ex 96, 9 mT only, seeds 941 … 944. Per seed five runs:
  F    dx 0.75 l_ex, virgin relaxation (the reference state of the seed);
  T15  dx 1.5,  start from the relaxed F state (block average, --init_transfer), then the production relaxation;
  T3   dx 3,    the same from F;
  I15  dx 1.5,  independent virgin relaxation from the same seed (as all mesh checks so far);
  I3   dx 3,    the same.
Template: the mesh212 line M212_dx15_s931_b9 (single amplitude, ≥ 5 cycles, ≤ 12 cycles, --min_calm 2, 256 samples,
no accumulators, no snapshots); only the name, --d_lex, --dx_lex, --seed, --init_transfer, --init_wait_h and --atol
change (atol 1e-6 since PROTOKOLL 7.88; names B0b_*).
Order: the four F runs first (longest, one per GPU), then I, then T (the T runs need the init.pt of their F run;
--init_wait_h covers a slow F relaxation).

    python conv_jobs.py > jobs/conv_b0.txt

Variant N0b (PLAN 4.11, PROTOKOLL 7.105): the same lines with the reference stray-field solver (--strayfield newell
--strayfield_shells 4) and the prefix N0b; the T arms transfer from N0b_F. Stage 1: I3 and I15; stage 2: F, T15, T3.

    python conv_jobs.py --newell 1 > jobs/newell_b0_s1.txt
    python conv_jobs.py --newell 2 > jobs/newell_b0_s2.txt
    python conv_jobs.py --newell 1b > jobs/newell_b0_s1b.txt     # the two low stage 1 runs again with --b_tol
    python conv_jobs.py --newell eps > jobs/newell_eps075.txt    # energy balance on dx 0.75 with 1024 samples
    python conv_jobs.py --newell tol > jobs/newell_tol075.txt    # the same with a 10 x tighter RKF45 tolerance
    python conv_jobs.py --newell kext > jobs/newell_kext15.txt   # dx 1.5 runs with the extrapolated kernel
    python conv_jobs.py --newell final > jobs/final_n1b.txt     # PLAN 4.12 V3, prefix N1b (replaced by V3g)
    python conv_jobs.py --newell locked > jobs/locked_n1g.txt   # PLAN 4.12 V3g, prefix N1g, geometry locked at dx 3
    python conv_jobs.py --newell tolbench > jobs/tolbench.txt   # cost of the cycle tolerance on the main-study mesh
    python conv_jobs.py --newell corners > jobs/corners_nx.txt  # the 12 corner runs, production setup of PLAN 4.13
    python conv_jobs.py --newell precision > jobs/precision_nx.txt   # PLAN 4.13 precision test (fp64 against fp32)
    python conv_jobs.py --newell solvertest > jobs/solvertest.txt    # FD against Newell at 150 mT from the same state
    python conv_jobs.py --newell solverscan > jobs/solverscan.txt    # FD at all 7 stages of one seed per corner
    python -c "import conv_jobs as CJ; [print(*x) for x in CJ.predemag_lines(runs_dir=<study runs>)]" > jobs/predemag.txt
    python conv_jobs.py --newell initmode > jobs/initmode.txt        # start: each particle uniform (Bloch points)
"""
import shlex

from rerun_jobs import get_arg, job_lines, set_arg

SEEDS = range(941, 945)
D_LEX = 96.0
ARMS = (("F", 0.75, False), ("I15", 1.5, False), ("I3", 3.0, False), ("T15", 1.5, True), ("T3", 3.0, True))
WAIT_H = 40.0
CHANGED = ("d_lex", "dx_lex", "seed", "init_transfer", "init_wait_h", "atol")
ATOL = "1e-06"          # PROTOKOLL 7.88: with atol 1e-5 the energy check fails on dx 0.75 (w_loop is not affected)
PREFIX = "B0b"          # the rerun of B0 with atol 1e-6 (the stopped B0 runs keep their names)
PREFIX_NEWELL = "N0b"   # B0b with the reference stray-field solver (PLAN 4.11)
NEWELL_ARGS = (("strayfield", ["newell"]), ("strayfield_shells", ["4"]))
NEWELL_STAGES = {1: ("I3", "I15"), 2: ("F", "T15", "T3")}
B_TOL = "0.01"         # stage 2: drive corrected again after cycle 1 if b is more than 1 % off the target
S1B = ("N0b_I15_s941", "N0b_I3_s944")   # stage 1 runs 5 … 6 % below the target b: again with --b_tol after stage 2
#                                         (decision Chris 2026-10-15, PROTOKOLL 7.107)


def name(arm, s, prefix=PREFIX):
    return "%s_%s_s%d" % (prefix, arm, s)


def lines(src="jobs/mesh212.txt"):
    tok0 = shlex.split(dict(job_lines(src))["M212_dx15_s931_b9"])
    out = []
    for arm, dx, transfer in ARMS:
        for s in SEEDS:
            t = set_arg(tok0, "d_lex", ["%.1f" % D_LEX])
            t = set_arg(t, "dx_lex", [repr(dx)])
            t = set_arg(t, "seed", [str(s)])
            t = set_arg(t, "atol", [ATOL])
            if transfer:
                t = set_arg(t, "init_transfer", ["%s/init.pt" % name("F", s)])
                t = set_arg(t, "init_wait_h", ["%g" % WAIT_H])
            out.append((name(arm, s), " ".join(t)))
    return out


def newell_lines(stage, src="jobs/mesh212.txt"):
    """The B0b lines of the arms of the stage with the reference stray-field solver (PLAN 4.11): prefix N0b,
    --strayfield newell --strayfield_shells 4, the T arms transfer from N0b_F of the same seed."""
    out = []
    for n, a in lines(src):
        arm, s = n.split("_")[1], int(n.split("_s")[1])
        if arm not in NEWELL_STAGES[stage]:
            continue
        t = shlex.split(a)
        for k, v in NEWELL_ARGS:
            t = set_arg(t, k, v)
        if arm.startswith("T"):
            t = set_arg(t, "init_transfer", ["%s/init.pt" % name("F", s, PREFIX_NEWELL)])
        if stage == 2:
            t = set_arg(t, "b_tol", [B_TOL])              # second drive correction (Chris 2026-10-15, 7.107)
        out.append((name(arm, s, PREFIX_NEWELL), " ".join(t)))
    return out


def newell_rerun_lines(src="jobs/mesh212.txt"):
    """The stage 1 lines of S1B with --b_tol B_TOL (the same drive rule as stage 2)."""
    return [(n, " ".join(set_arg(shlex.split(a), "b_tol", [B_TOL]))) for n, a in newell_lines(1, src) if n in S1B]


EPS_SEEDS = (941, 942)    # the F runs with the largest energy-balance error (+0.170, +0.091; PROTOKOLL 7.108)
EPS_SAMPLES = "1024"


def newell_eps_lines(src="jobs/mesh212.txt"):
    """Check of the energy balance on dx 0.75 (decision Chris 2026-10-16, PROTOKOLL 7.108): the stage 2 F line of
    EPS_SEEDS with --samples 1024 (instead of 256), started from the relaxed state of the same F run (--init_from;
    no new relaxation). Name <F run>_smp1024."""
    out = []
    for n, a in newell_lines(2, src):
        if not n.startswith(PREFIX_NEWELL + "_F_") or int(n.split("_s")[1]) not in EPS_SEEDS:
            continue
        t = set_arg(shlex.split(a), "samples", [EPS_SAMPLES])
        t = set_arg(t, "init_from", ["%s/init.pt" % n])
        out.append((n + "_smp" + EPS_SAMPLES, " ".join(t)))
    return out


TOL_CYCLES = (("atol_cycles", "1e-07"), ("rtol_cycles", "1e-06"))    # 10 x tighter than atol 1e-6 / rtol 1e-5
TOL_SUFFIX = "_atol1e-07_rtol1e-06"


def newell_tol_lines(src="jobs/mesh212.txt"):
    """Integrator test on dx 0.75 (decision Chris 2026-10-16, option C2, PROTOKOLL 7.108): the stage 2 F line of
    EPS_SEEDS with the RKF45 tolerance of the drive cycles 10 x tighter (--atol_cycles 1e-07 --rtol_cycles 1e-06),
    started from the relaxed state of the same F run (--init_from; no new relaxation; 256 samples as in stage 2)."""
    out = []
    for n, a in newell_lines(2, src):
        if not n.startswith(PREFIX_NEWELL + "_F_") or int(n.split("_s")[1]) not in EPS_SEEDS:
            continue
        t = shlex.split(a)
        for k, v in TOL_CYCLES:
            t = set_arg(t, k, [v])
        t = set_arg(t, "init_from", ["%s/init.pt" % n])
        out.append((n + TOL_SUFFIX, " ".join(t)))
    return out


KEXT_ARMS = ("T15", "I15")
KEXT_SUFFIX = "_kext"


def newell_kext_lines(src="jobs/mesh212.txt"):
    """Test of the shell extrapolation on the loss (decision Chris 2026-10-16, option 1', PROTOKOLL 7.108): the dx 1.5
    runs of N0b (T15 from stage 2, I15 from stage 1 with the stage 1b line of s941) again from their own relaxed state
    (--init_from, no new relaxation; --init_transfer removed) with the extrapolated kernel in the drive cycles
    (--strayfield_cycles newell_ext), all else equal. Name <run>_kext."""
    cur = dict(newell_lines(1, src))
    cur.update(newell_rerun_lines(src))
    cur.update(newell_lines(2, src))
    out = []
    for arm in KEXT_ARMS:
        for s in SEEDS:
            n = name(arm, s, PREFIX_NEWELL)
            t = set_arg(shlex.split(cur[n]), "init_transfer", None)
            t = set_arg(t, "strayfield_cycles", ["newell_ext"])
            t = set_arg(t, "init_from", ["%s/init.pt" % n])
            out.append((n + KEXT_SUFFIX, " ".join(t)))
    return out


PREFIX_FINAL = "N1b"     # PLAN 4.12 V3: the final stray-field setting (PROTOKOLL 7.109: extrapolated from 4 shells)
FINAL_ARGS = (("strayfield", ["newell_ext"]), ("strayfield_shells", ["4"]), ("b_tol", [B_TOL]),
              ("atol_cycles", ["1e-07"]), ("rtol_cycles", ["1e-06"]))
FINAL_ORDER = ("F", "I3", "I15", "T15", "T3")    # the long F runs first; the T arms wait for N1b_F/init.pt


def final_lines(src="jobs/mesh212.txt"):
    """PLAN 4.12 V3: stage 1 and 2 of N0b again with the final setting (prefix N1b): --strayfield newell_ext
    --strayfield_shells 4, --b_tol 0.01 on all lines, the tighter tolerance of the drive cycles, the T arms transfer
    from N1b_F of the same seed. All other arguments as N0b."""
    cur = dict(newell_lines(1, src))
    cur.update(newell_lines(2, src))
    out = []
    for arm in FINAL_ORDER:
        for s in SEEDS:
            t = shlex.split(cur[name(arm, s, PREFIX_NEWELL)])
            for k, v in FINAL_ARGS:
                t = set_arg(t, k, v)
            if arm.startswith("T"):
                t = set_arg(t, "init_transfer", ["%s/init.pt" % name("F", s, PREFIX_FINAL)])
            out.append((name(arm, s, PREFIX_FINAL), " ".join(t)))
    return out


PREFIX_LOCKED = "N1g"    # PLAN 4.12 V3g (decision Chris 2026-10-17): the final setting with the geometry locked at dx 3
GEOM_DX = "3.0"


def locked_lines(src="jobs/mesh212.txt"):
    """PLAN 4.12 V3g: the N1b lines (final_lines) with --geom_dx_lex 3.0 (voxel geometry, particle index and cube axes
    of dx 3 on all meshes) and the prefix N1g; the T arms transfer from N1g_F."""
    out = []
    for n, a in final_lines(src):
        t = set_arg(shlex.split(a), "geom_dx_lex", [GEOM_DX])
        n = n.replace(PREFIX_FINAL, PREFIX_LOCKED)
        if "_T" in n:
            t = set_arg(t, "init_transfer", ["%s/init.pt" % name("F", int(n.split("_s")[1]), PREFIX_LOCKED)])
        out.append((n, " ".join(t)))
    return out


TB_CORNERS = ("K1_0", "K2_0")      # d/l_ex 212 and 300 at dx 3 (the mesh of the main study)
TB_VARIANTS = (("a5", ()), ("a6", (("atol_cycles", ["1e-06"]),)),
               ("tight", (("atol_cycles", ["1e-07"]), ("rtol_cycles", ["1e-06"]))))


def tolbench_lines(src="jobs/sobol_corners.txt"):
    """Cost of the cycle tolerance on the mesh of the main study (decision Chris 2026-10-17: measure, then decide): the
    corner lines K1_0 and K2_0 with the final stray-field setting (newell_ext, 4 shells), one amplitude (9 mT), 3
    cycles, no snapshots; variant a5 relaxes with the tolerance of the main study (atol 1e-5) and runs its cycles with
    it; a6 (atol 1e-6) and tight (atol 1e-7, rtol 1e-6) run the cycles from the same relaxed state (--init_from)."""
    lines = dict(job_lines(src))
    out = []
    for c in TB_CORNERS:
        t = shlex.split(lines[c])
        for k in ("snap_phases", "snap_cycles", "snap_mT"):
            t = set_arg(t, k, None)
        for k, v in (("strayfield", ["newell_ext"]), ("strayfield_shells", ["4"]), ("n_amp", ["1"]),
                     ("b_list", ["0.006"]), ("cycles_per_amp", ["3"]), ("max_cycles_per_amp", ["3"])):
            t = set_arg(t, k, v)
        for v, extra in TB_VARIANTS:
            u = list(t)
            for k, val in extra:
                u = set_arg(u, k, val)
            if v != "a5":
                u = set_arg(u, "init_from", ["TB_%s_a5/init.pt" % c])
                u = set_arg(u, "init_wait_h", ["3"])
            out.append(("TB_%s_%s" % (c, v), " ".join(u)))
    return out


CORNERS_NX = ("K1", "K2", "K3", "K4")
NX_SUFFIX = "_nx"


def corner_newell_lines(src="jobs/sobol_corners.txt"):
    """Corner test with the production setup of PLAN 4.13 (approved by Chris 2026-10-17): the 12 corner runs of the
    old main study (same seeds and design values, name <run>_nx) with the final stray-field setting (newell_ext,
    4 shells), the tight cycle tolerance, the relaxation with atol 1e-6 (D2) and the absolute stop threshold 3e-10
    (D1, S7), the online stop (S1: >= 7 cycles, until 2 calm cycles, at most 15), the second drive correction (S2),
    --acc_harm 7 (S3), --keep_stage_ckpt (S4) and the snapshots of the old main study (16 phases, 2 cycles, 9 / 50 /
    150 mT). Order: by seed, then corner, so that d/l_ex 212 (K1, K3) and 300 (K2, K4) alternate (S11)."""
    out = []
    for n, a in job_lines(src):
        t = shlex.split(a)
        for k, v in NX_SET:
            t = set_arg(t, k, v)
        out.append((n + NX_SUFFIX, " ".join(t)))
    return sorted(out, key=lambda x: (x[0].split("_")[1], x[0].split("_")[0]))


NX_SET = (("strayfield", ["newell_ext"]), ("strayfield_shells", ["4"]), ("atol", ["1e-06"]),
          ("relax_tol_E_abs", ["3e-10"]), ("atol_cycles", ["1e-07"]), ("rtol_cycles", ["1e-06"]),
          ("cycles_per_amp", ["7"]), ("min_calm", ["2"]), ("max_cycles_per_amp", ["15"]), ("b_tol", ["0.01"]),
          ("acc", []), ("acc_harm", ["7"]), ("keep_stage_ckpt", []),
          ("snap_phases", ["16"]), ("snap_cycles", ["2"]), ("snap_mT", ["9.0", "50.0", "150.0"]))
PRECISION_RUNS = ("K1_0_nx", "K2_0_nx", "K3_0_nx", "K4_0_nx")


def precision_lines(src="jobs/sobol_corners.txt"):
    """Precision test of PLAN 4.13 (part 5, step 4): from the init.pt of one corner run per corner, the 9 mT stage only,
    all other settings as the corner run, without snapshots and accumulators; <run>_p64 with --precision_cycles
    double, <run>_p32 in fp32. Both start from the same init.pt."""
    corners = dict(corner_newell_lines(src))
    out = []
    for n in PRECISION_RUNS:
        t = shlex.split(corners[n])
        for k in ("snap_phases", "snap_cycles", "snap_mT", "acc", "acc_harm", "keep_stage_ckpt"):
            t = set_arg(t, k, None)
        t = set_arg(set_arg(t, "b_list", ["0.006"]), "n_amp", ["1"])
        t = set_arg(set_arg(t, "init_from", ["%s/init.pt" % n]), "init_wait_h", ["0.0"])
        out.append((n + "_p64", " ".join(set_arg(t, "precision_cycles", ["double"]))))
        out.append((n + "_p32", " ".join(t)))
    return out


SOLVER_RUNS = ("K1_0_nx", "K1_1_nx", "K1_2_nx")
SOLVER_CKPT = "ckpt_05_amp05.pt"           # the state after the 100 mT stage (--keep_stage_ckpt)


def solvertest_lines(src="jobs/sobol_corners.txt"):
    """Solver test at 150 mT (Chris 2026-10-18): from the state of K1_s_nx after the 100 mT stage, the 150 mT stage again
    with the same settings as the corner run, once with DemagFieldPBC in the cycles (<run>_150fd, --strayfield_cycles
    fd) and once with the final setting (<run>_150nx). Without snapshots and accumulators."""
    corners = dict(corner_newell_lines(src))
    out = []
    for n in SOLVER_RUNS:
        t = shlex.split(corners[n])
        for k in ("snap_phases", "snap_cycles", "snap_mT", "acc", "acc_harm", "keep_stage_ckpt"):
            t = set_arg(t, k, None)
        t = set_arg(t, "init_stage", ["%s/%s" % (n, SOLVER_CKPT)])
        out.append((n + "_150fd", " ".join(set_arg(t, "strayfield_cycles", ["fd"]))))
        out.append((n + "_150nx", " ".join(t)))
    return out


SCAN_RUNS = ("K1_0_nx", "K2_0_nx", "K3_0_nx", "K4_0_nx")


def scan_name(run, k):
    """Run name of the solver scan for stage k of corner run `run` (K1_0_nx at 150 mT: the solver test run)."""
    return run + "_150fd" if (run, k) == ("K1_0_nx", 6) else "%s_fd%02d" % (run, k)


def solverscan_lines(src="jobs/sobol_corners.txt"):
    """Solver scan (option (a), Chris 2026-10-18): one seed per corner, each of the 7 stages once more with DemagFieldPBC
    in the cycles, from the same state as the corner run (stage 0: its init.pt; stage k: its state after stage k − 1)
    and with the first k + 1 amplitudes, so that the run ends after stage k. The Newell arm is the corner run itself (the
    Newell restart of the solver test reproduced it exactly). K1_0_nx at 150 mT is K1_0_nx_150fd of the solver test."""
    corners = dict(corner_newell_lines(src))
    out = []
    for n in SCAN_RUNS:
        t0 = shlex.split(corners[n])
        for k in ("snap_phases", "snap_cycles", "snap_mT", "acc", "acc_harm", "keep_stage_ckpt"):
            t0 = set_arg(t0, k, None)
        bl = get_arg(t0, "b_list")
        for k in range(len(bl)):
            if scan_name(n, k).endswith("_150fd"):
                continue
            t = set_arg(set_arg(t0, "b_list", bl[:k + 1]), "n_amp", [str(k + 1)])
            if k == 0:
                t = set_arg(set_arg(t, "init_from", ["%s/init.pt" % n]), "init_wait_h", ["0.0"])
            else:
                t = set_arg(t, "init_stage", ["%s/ckpt_%02d_amp%02d.pt" % (n, k - 1, k - 1)])
            out.append((scan_name(n, k), " ".join(set_arg(t, "strayfield_cycles", ["fd"]))))
    return sorted(out, key=lambda x: (x[0][:2] not in ("K2", "K4"), x[0]))     # the longer d/l_ex 300 runs first


PD_RUNS = ("K2_0_nx", "K2_1_nx", "K4_1_nx", "K3_1_nx")


def predemag_lines(src="jobs/sobol_corners.txt", runs_dir="runs"):
    """Test of the AC demagnetization after the virgin relaxation (Chris 2026-10-18, Bloch points): the 4 corner runs
    with the most Bloch points, from their relaxed virgin state (init.pt), AC demagnetization from (a) the drive of
    their 150 mT stage (<run>_pd150), (b) twice that (<run>_pd300), (c) saturation 1 − φ (<run>_pdsat); then the
    relaxation again and the 9 mT stage. All other settings as the corner run, without snapshots and accumulators."""
    import json
    import pathlib
    corners = dict(corner_newell_lines(src))
    out = []
    for n in PD_RUNS:
        st = json.loads((pathlib.Path(runs_dir) / n / "summary.json").read_text())["stages"]
        h150 = [s for s in st if s.get("measure")][-1]["h_amp"]
        t = shlex.split(corners[n])
        for k in ("snap_phases", "snap_cycles", "snap_mT", "acc", "acc_harm", "keep_stage_ckpt"):
            t = set_arg(t, k, None)
        t = set_arg(set_arg(t, "b_list", ["0.006"]), "n_amp", ["1"])
        t = set_arg(set_arg(t, "init_from", ["%s/init.pt" % n]), "init_wait_h", ["0.0"])
        for tag, h in (("pd150", h150), ("pd300", 2.0 * h150), ("pdsat", 0.0)):
            out.append(("%s_%s" % (n, tag), " ".join(set_arg(t, "predemag_h", ["%.6g" % h]))))
    return sorted(out, key=lambda x: (x[0].split("_")[-1], x[0]))


def initmode_lines(src="jobs/sobol_corners.txt"):
    """Test of the start "each particle uniform along its own random direction" (Chris 2026-10-18, Bloch points, way 2):
    the 4 corner runs with the most Bloch points (PD_RUNS), virgin relaxation from this start (--init_mode particle,
    no --init_from), then the 9 mT stage (<run>_ip). All other settings as the corner run, without snapshots and
    accumulators."""
    corners = dict(corner_newell_lines(src))
    out = []
    for n in PD_RUNS:
        t = shlex.split(corners[n])
        for k in ("snap_phases", "snap_cycles", "snap_mT", "acc", "acc_harm", "keep_stage_ckpt"):
            t = set_arg(t, k, None)
        t = set_arg(set_arg(set_arg(t, "b_list", ["0.006"]), "n_amp", ["1"]), "init_mode", ["particle"])
        out.append((n + "_ip", " ".join(t)))
    return out


START_TRIES = 3


def start_lines(kind, src="jobs/sobol_corners.txt", runs_dir="runs"):
    """The 12 corner runs with a start state without Bloch points (Chris 2026-10-18; first runs of a new main study,
    the full setup of PLAN 4.13): kind "hs" = each particle uniform (--init_mode particle), names <run>_hs; kind "ad300"
    = AC demagnetization from twice the drive of the 150 mT stage of the corner run <run>_nx, "adsat" = from
    saturation (names <run>_ad and <run>_as). --bloch_free_tries 3: a prepared state with a Bloch point gets a new random start of
    the same microstructure (at most 3 attempts). The Bloch points are counted at each sample (run_loops)."""
    import json
    import pathlib
    out = []
    for n, a in corner_newell_lines(src):
        run = n[:-len(NX_SUFFIX)]
        t = set_arg(shlex.split(a), "bloch_free_tries", [str(START_TRIES)])
        if kind == "hs":
            t = set_arg(t, "init_mode", ["particle"])
            out.append((run + "_hs", " ".join(t)))
        else:
            if kind == "ad300":
                st = json.loads((pathlib.Path(runs_dir) / n / "summary.json").read_text())["stages"]
                h = 2.0 * [s for s in st if s.get("measure")][-1]["h_amp"]
            else:
                h = 0.0
            out.append((run + ("_ad" if kind == "ad300" else "_as"), " ".join(set_arg(t, "predemag_h", ["%.6g" % h]))))
    return out


PREC = (("P48_dx075_single", 0.75, {}), ("P48_dx075_double", 0.75, {"precision": ["double"]}),
        ("P48_dx075_atol6", 0.75, {"atol": ["1e-06"]}), ("P48_dx15_single", 1.5, {}))


def ad_continue_lines(src="jobs/sobol_corners.txt", runs_dir="runs"):
    """Option C (Chris 2026-10-18): the corner runs of start_lines("ad300") whose test run <run>_nx_pd300 ended without
    a Bloch point start from its prepared state (init.pt: the virgin relaxation of <run>_nx, AC from 300 mT, the
    relaxation again; the same INIT_KEYS) and run all 7 stages, the 9 mT stage again with snapshots and accumulators.
    --bloch_free_tries 0 as the test run: the state needs no new start."""
    import json
    import pathlib
    out = []
    for n, a in start_lines("ad300", src, runs_dir):
        pd = pathlib.Path(runs_dir) / ("%s%s_pd300" % (n[:-len("_ad")], NX_SUFFIX))
        if not (pd / "DONE").exists() or json.loads((pd / "summary.json").read_text())["init"]["predemag"]["bloch_points"]:
            continue
        t = set_arg(shlex.split(a), "bloch_free_tries", ["0"])
        t = set_arg(set_arg(t, "init_from", ["%s/init.pt" % pd.name]), "init_wait_h", ["0.0"])
        out.append((n, " ".join(t)))
    return out


def ad300_lines(src="jobs/sobol_corners.txt", runs_dir="runs"):
    """Option C, the 12 corner runs (jobs/corners_ad300.txt): start_lines("ad300"), the lines of ad_continue_lines in
    place of their random start."""
    cont = dict(ad_continue_lines(src, runs_dir))
    return [(n, cont.get(n, a)) for n, a in start_lines("ad300", src, runs_dir)]


def ad_rest_lines(src="jobs/sobol_corners.txt", runs_dir="runs"):
    """Option C: the 8 corner runs of start_lines("ad300") outside PD_RUNS (they start from the random start)."""
    skip = {r[:-len(NX_SUFFIX)] + "_ad" for r in PD_RUNS}
    return [(n, a) for n, a in start_lines("ad300", src, runs_dir) if n not in skip]


def satstart_lines(src="jobs/sobol_corners.txt", runs_dir="runs"):
    """Test of protocol S (Chris 2026-10-18, Bloch points): the 4 runs of PD_RUNS start uniform along the drive direction
    (--init_mode saturated, no relaxation before), AC demagnetization from saturation and the relaxation, AC from twice
    the drive of the 150 mT stage of <run>_nx (300 mT) and the relaxation; if a Bloch point remains, one further pass
    from 300 mT (--bloch_free_tries 2); then the 9 mT stage (<run>_ps). All other settings as the corner run, without
    snapshots and accumulators."""
    import json
    import pathlib
    corners = dict(corner_newell_lines(src))
    out = []
    for n in PD_RUNS:
        st = json.loads((pathlib.Path(runs_dir) / n / "summary.json").read_text())["stages"]
        h150 = [s for s in st if s.get("measure")][-1]["h_amp"]
        t = shlex.split(corners[n])
        for k in ("snap_phases", "snap_cycles", "snap_mT", "acc", "acc_harm", "keep_stage_ckpt"):
            t = set_arg(t, k, None)
        t = set_arg(set_arg(t, "b_list", ["0.006"]), "n_amp", ["1"])
        t = set_arg(set_arg(t, "init_mode", ["saturated"]), "bloch_free_tries", ["2"])
        out.append((n + "_ps", " ".join(set_arg(t, "predemag_h", ["0", "%.6g" % (2.0 * h150)]))))
    return out


SM_RUNS = ("K1_2_nx", "K2_0_nx", "K2_1_nx", "K2_2_nx", "K3_1_nx", "K3_2_nx", "K4_0_nx", "K4_1_nx", "K4_2_nx")
# the corner runs whose relaxed virgin state (init.pt) has Bloch points (census_bloch_nx_v2.json, PROTOKOLL 7.112)
SA_RUNS = ("K3_1_ad", "K4_0_ad", "K2_2_ad")
# the option C runs whose first attempt (AC from 300 mT and the relaxation) kept a Bloch point (K3_1_ad: PROTOKOLL
# 7.115; K4_0_ad, K2_2_ad: stdout.log "[bloch] attempt 0" on the box; PROTOKOLL 7.116)


def smooth_lines(src="jobs/sobol_corners.txt"):
    """Test of the Bloch-point cut (Chris 2026-10-19, --bloch_smooth 3): the relaxed virgin state of each run of SM_RUNS
    (init.pt, with Bloch points) is loaded, cut and relaxed (at most 3 steps), then the 9 mT stage runs (<run>_sm). All
    other settings as the corner run, without snapshots and accumulators."""
    corners = dict(corner_newell_lines(src))
    out = []
    for n in SM_RUNS:
        t = shlex.split(corners[n])
        for k in ("snap_phases", "snap_cycles", "snap_mT", "acc", "acc_harm", "keep_stage_ckpt"):
            t = set_arg(t, k, None)
        t = set_arg(set_arg(t, "b_list", ["0.006"]), "n_amp", ["1"])
        t = set_arg(set_arg(t, "init_from", ["%s/init.pt" % n]), "init_wait_h", ["0.0"])
        out.append((n + "_sm", " ".join(set_arg(t, "bloch_smooth", ["3"]))))
    return out


def smooth_ad_lines(src="jobs/sobol_corners.txt", runs_dir="runs"):
    """Test of the Bloch-point cut after one AC pass (Chris 2026-10-19; replaces the rest of the test of smooth_lines):
    the option C line of each run of SA_RUNS with one attempt (--bloch_free_tries 1: the random start, the relaxation,
    AC from 300 mT, the relaxation) and --bloch_smooth 3 (if a Bloch point remains: at most 3 cuts, each with the
    relaxation after it), then the 9 mT stage (<run>_sm). All other settings as the option C run, without snapshots and
    accumulators."""
    ad = dict(start_lines("ad300", src, runs_dir))
    out = []
    for n in SA_RUNS:
        t = shlex.split(ad[n])
        for k in ("snap_phases", "snap_cycles", "snap_mT", "acc", "acc_harm", "keep_stage_ckpt"):
            t = set_arg(t, k, None)
        t = set_arg(set_arg(t, "b_list", ["0.006"]), "n_amp", ["1"])
        t = set_arg(t, "bloch_free_tries", ["1"])
        out.append((n + "_sm", " ".join(set_arg(t, "bloch_smooth", ["3"]))))
    return out


SAT_H_CHECK = (None, 1.0, 2.0)     # --h_sat for the field check; None = the default 1 - phi


def _sat_base(n, corners):
    """The corner line of n (PD_RUNS) for a 9 mT test: no snapshots, no accumulators, one stage."""
    t = shlex.split(corners[n])
    for k in ("snap_phases", "snap_cycles", "snap_mT", "acc", "acc_harm", "keep_stage_ckpt"):
        t = set_arg(t, k, None)
    return set_arg(set_arg(t, "b_list", ["0.006"]), "n_amp", ["1"])


def satcheck_lines(src="jobs/sobol_corners.txt"):
    """Static saturation, field check (Chris 2026-10-19): for each run of PD_RUNS the uniform start along the drive
    (--init_mode saturated), the static saturation at --h_sat (SAT_H_CHECK; None: 1 - phi) with --sat_static, and the
    stop right after it (--sat_check: alignment of the particle cells). Names <run>_sc<h> (h in % of M_s; d: default)."""
    corners = dict(corner_newell_lines(src))
    out = []
    for n in PD_RUNS:
        for h in SAT_H_CHECK:
            t = _sat_base(n, corners)
            t = set_arg(set_arg(t, "init_mode", ["saturated"]), "predemag_h", ["0"])
            if h is not None:
                t = set_arg(t, "h_sat", ["%g" % h])
            t = set_arg(t, "bloch_free_tries", ["1"])
            out.append(("%s_sc%s" % (n, "d" if h is None else "%d" % round(100 * h)), " ".join(t) + " --sat_static --sat_check"))
    return out


SAT_TEST_RUNS = ("K4_1_nx", "K2_0_nx")    # history strong (K4_1) and Bloch points (K2_0)
SAT_TEST_H = (None, 2.0)                  # --h_sat: None = the default 1 - phi; 2 M_s


def satstatic_lines(runs=SAT_TEST_RUNS, levels=SAT_TEST_H, src="jobs/sobol_corners.txt"):
    """Static saturation test (Chris 2026-10-19): for each run and each --h_sat level (None: 1 - phi) two starts; each
    pass: the field up to h_sat (a smooth quarter-sine ramp, except the first pass of start A), held there and relaxed
    (--sat_static), the decaying AC field from its peak (factor --predemag_q per half cycle), the relaxation at H = 0;
    if a Bloch point remains, the next pass (at most 3, --bloch_free_tries 3); then the 9 mT stage.
      <run>_sa<h>  start A: uniform along the drive (--init_mode saturated); the first pass starts at h_sat at once
      <run>_sb<h>  start B: the prepared state after AC from 300 mT (<run>_pd300/init.pt, 7.114)
    Question 1 (does the saturation erase the history): sa against sb at the same level. Question 2 (does the level
    matter): the two levels at the same start. h: d = default 1 - phi, else in % of M_s."""
    corners = dict(corner_newell_lines(src))
    out = []
    for n in runs:
        for h in levels:
            tag = "d" if h is None else "%d" % round(100 * h)
            for kind in ("sa", "sb"):
                t = _sat_base(n, corners)
                t = set_arg(set_arg(t, "predemag_h", ["0"]), "bloch_free_tries", ["3"])
                if h is not None:
                    t = set_arg(t, "h_sat", ["%g" % h])
                if kind == "sa":
                    t = set_arg(t, "init_mode", ["saturated"])
                else:
                    t = set_arg(set_arg(t, "init_from", ["%s_pd300/init.pt" % n]), "init_wait_h", ["0.0"])
                out.append(("%s_%s%s" % (n, kind, tag), " ".join(t) + " --sat_static"))
    return out


CREP_RUNS = ("K4_1",)
CREP_SEEDS = (9001, 9002, 9003, 9004)


def crep_lines(runs=CREP_RUNS, seeds=CREP_SEEDS, src="jobs/sobol_corners.txt", ad_src="jobs/corners_ad300.txt"):
    """Option C repeated with new random starts (Chris 2026-10-20, PROTOKOLL 7.130): the same microstructure (--seed),
    a new random start (--init_seed s), the relaxation, AC from the field of the option C run <run>_ad (its
    --predemag_h in ad_src), the relaxation, at most 3 attempts (--bloch_free_tries 3); only the 9 mT stage, no
    snapshots and no accumulators (as the static saturation test). Names <run>_cr<s>. Question: can option C at the
    same microstructure also give a low 9 mT loss?"""
    corners = dict(corner_newell_lines(src))
    ad = dict(job_lines(ad_src))
    out = []
    for r in runs:
        h = get_arg(shlex.split(ad[r + "_ad"]), "predemag_h")
        for sd in seeds:
            t = _sat_base(r + NX_SUFFIX, corners)
            t = set_arg(set_arg(t, "predemag_h", h), "bloch_free_tries", [str(START_TRIES)])
            out.append(("%s_cr%d" % (r, sd), " ".join(set_arg(t, "init_seed", [str(sd)]))))
    return out


SAT_SET = (("init_mode", ["saturated"]), ("predemag_h", ["0"]), ("h_sat", ["2"]), ("bloch_free_tries", ["3"]))
SAT_SUFFIX = "_sat"
SAT_N = 4                      # PLAN 4.15 D1 (c): 4 seeds per corner in the first step
SAT_SMALL = ("K1", "K3")       # d/l_ex 212; K2 and K4 have d/l_ex 300


def sat_seed(c, k):
    """Seed of corner K<c>, seed index k (PLAN 4.15): k 0 ... 2 the seeds of the old corners (2011 ... 2013 for K1), k 3 ...
    8 the same scheme, then 3000 + 100 c + k."""
    return 2000 + 10 * c + 1 + k if k <= 8 else 3000 + 100 * c + k


def sat_corner_lines(n=SAT_N, kind="all", src="jobs/sobol_corners.txt", k0=0):
    """Corner study of PLAN 4.15 (approved by Chris 2026-10-20): the production line of corner_newell_lines (PLAN 4.13)
    with the preparation of PROTOKOLL 7.132: the uniform start along the drive (--init_mode saturated), the field 2 M_s
    at once (--h_sat 2), held and relaxed (--sat_static), the decaying AC field from its peak, the relaxation at H = 0,
    at most 3 passes (--bloch_free_tries 3); then all 7 amplitudes. Names K<c>_<k>_sat with the seed sat_seed(c, k),
    k = k0 ... n - 1 (k0 > 0: more seeds of the same corners, e.g. k 4 and 5 of K1 and K3, Chris 2026-10-20). kind:
    "small" (K1, K3), "large" (K2, K4) or "all". Order: by k, then corner."""
    corners = dict(corner_newell_lines(src))
    out = []
    for k in range(k0, n):
        for c in (1, 2, 3, 4):
            if kind != "all" and ("K%d" % c in SAT_SMALL) != (kind == "small"):
                continue
            t = shlex.split(corners["K%d_0%s" % (c, NX_SUFFIX)])    # the 3 old lines of a corner differ only in --seed
            for a, v in SAT_SET:
                t = set_arg(t, a, v)
            t = set_arg(t, "seed", [str(sat_seed(c, k))])
            out.append(("K%d_%d%s" % (c, k, SAT_SUFFIX), " ".join(t) + " --sat_static"))
    return out


K4H_RUNS = ("K4_0", "K4_1", "K4_2")
K4H_CYCLES = 16


def k4hist_lines(src="jobs/sobol_corners.txt", runs_dir="runs"):
    """K4 history test (Chris 2026-10-19), see hist_lines; with snapshots."""
    return hist_lines(K4H_RUNS, True, src, runs_dir)


K3H_RUNS = ("K3_0", "K3_1", "K3_2")


def k3hist_lines(src="jobs/sobol_corners.txt", runs_dir="runs"):
    """K3 history test (Chris 2026-10-19, option A): as the K4 test (hist_lines) for K3, without snapshots (at K3 the
    loss of option C at 9 ... 50 mT is 0.4 ... 0.6 times the loss of nx, with a large SE; PROTOKOLL 7.119)."""
    return hist_lines(K3H_RUNS, False, src, runs_dir)


def hist_lines(runs, snapshots, src="jobs/sobol_corners.txt", runs_dir="runs"):
    """History test of a corner (Chris 2026-10-19; first for K4, PROTOKOLL 7.119: at K4 the loss of option C at
    9 ... 50 mT is about 2 x the loss of the relaxed random start, and in the 9 mT stage it rises over the first 2 ... 3
    cycles). For each microstructure of runs the 9 mT stage with K4H_CYCLES cycles from both prepared states: <run>_ad_h16 from the option C
    state (<run>_ad/init.pt; K4_1: its AC-300 test state, as the run), <run>_nx_h16 from the relaxed random start of the
    corner test (<run>_nx/init.pt). With snapshots: m at 2 phases (180° and 360°: H = 0 falling and rising; drive
    H0 sin(ωt)) in every cycle (--snap_cycles
    K4H_CYCLES; 4 phases would need 25 GB of disk for K4); without: no snapshots. Accumulators as
    the runs; all other settings as the runs."""
    ad = dict(ad300_lines(src, runs_dir))
    nx = dict(corner_newell_lines(src))
    out = []
    for r in runs:
        for kind, line in (("ad", ad[r + "_ad"]), ("nx", nx[r + NX_SUFFIX])):
            t = shlex.split(line)
            t = set_arg(t, "keep_stage_ckpt", None)
            t = set_arg(set_arg(t, "b_list", ["0.006"]), "n_amp", ["1"])
            if kind == "nx" or get_arg(t, "init_from") is None:
                t = set_arg(t, "init_from", ["%s_%s/init.pt" % (r, kind)])
            t = set_arg(t, "init_wait_h", ["0.0"])
            t = set_arg(set_arg(t, "cycles_per_amp", [str(K4H_CYCLES)]), "max_cycles_per_amp", [str(K4H_CYCLES)])
            if snapshots:
                t = set_arg(set_arg(t, "snap_phases", ["2"]), "snap_cycles", [str(K4H_CYCLES)])
                t = set_arg(t, "snap_mT", ["9.0"])
            else:
                for k in ("snap_phases", "snap_cycles", "snap_mT"):
                    t = set_arg(t, k, None)
            out.append(("%s_%s_h16" % (r, kind), " ".join(t)))
    return out


def prec_lines(src="jobs/mesh212.txt"):
    """Precision test after B0 (PROTOKOLL 7.85, decision Chris 2026-10-14): the B0 line at d/l_ex 48 (8 × fewer cells
    than d/l_ex 96, the same dx / l_ex), seed 951, 9 mT, exactly 5 cycles (no steady rule: --min_calm 0,
    --max_cycles_per_amp 5). dx 0.75 in single and double precision and with atol 1e-6, dx 1.5 single as the reference.
    Question: why is w_dis / w_loop 2 … 5 in the B0 runs on dx 0.75?"""
    tok0 = shlex.split(dict(job_lines(src))["M212_dx15_s931_b9"])
    out = []
    for nm, dx, extra in PREC:
        t = set_arg(tok0, "d_lex", ["48.0"])
        t = set_arg(t, "dx_lex", [repr(dx)])
        t = set_arg(t, "seed", ["951"])
        t = set_arg(set_arg(t, "max_cycles_per_amp", ["5"]), "min_calm", ["0"])
        for k, v in extra.items():
            t = set_arg(t, k, v)
        out.append((nm, " ".join(t)))
    return out


if __name__ == "__main__":
    import sys
    if "--newell" in sys.argv:
        st = sys.argv[sys.argv.index("--newell") + 1]
        st, _, ks = st.partition(":")      # satsmall:4:6 -> the seed indices k 4 and 5 (sat_corner_lines k0, n)
        k0, n = (int(x) for x in ks.split(":")) if ks else (0, SAT_N)
        jobs = (newell_rerun_lines() if st == "1b" else newell_eps_lines() if st == "eps" else
                newell_tol_lines() if st == "tol" else newell_kext_lines() if st == "kext" else
                final_lines() if st == "final" else locked_lines() if st == "locked" else
                tolbench_lines() if st == "tolbench" else corner_newell_lines() if st == "corners" else
                precision_lines() if st == "precision" else solvertest_lines() if st == "solvertest" else
                solverscan_lines() if st == "solverscan" else predemag_lines() if st == "predemag" else
                initmode_lines() if st == "initmode" else
                satstart_lines() if st == "ps" else ad_continue_lines() if st == "adcont" else
                ad_rest_lines() if st == "adrest" else
                smooth_lines() if st == "smooth" else smooth_ad_lines() if st == "smoothad" else
                k4hist_lines() if st == "k4hist" else k3hist_lines() if st == "k3hist" else
                satcheck_lines() if st == "satcheck" else satstatic_lines() if st == "satstatic" else
                crep_lines() if st == "crep" else
                sat_corner_lines(n, "small", k0=k0) if st == "satsmall" else
                sat_corner_lines(n, "large", k0=k0) if st == "satlarge" else
                newell_lines(int(st)))
    else:
        jobs = prec_lines() if "--prec" in sys.argv else lines()
    for n, a in jobs:
        print(n, a)
