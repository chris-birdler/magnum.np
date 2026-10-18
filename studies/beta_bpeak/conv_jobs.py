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
"""
import shlex

from rerun_jobs import job_lines, set_arg

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


PREC = (("P48_dx075_single", 0.75, {}), ("P48_dx075_double", 0.75, {"precision": ["double"]}),
        ("P48_dx075_atol6", 0.75, {"atol": ["1e-06"]}), ("P48_dx15_single", 1.5, {}))


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
        jobs = (newell_rerun_lines() if st == "1b" else newell_eps_lines() if st == "eps" else
                newell_tol_lines() if st == "tol" else newell_kext_lines() if st == "kext" else
                final_lines() if st == "final" else locked_lines() if st == "locked" else
                tolbench_lines() if st == "tolbench" else corner_newell_lines() if st == "corners" else
                precision_lines() if st == "precision" else solvertest_lines() if st == "solvertest" else
                newell_lines(int(st)))
    else:
        jobs = prec_lines() if "--prec" in sys.argv else lines()
    for n, a in jobs:
        print(n, a)
