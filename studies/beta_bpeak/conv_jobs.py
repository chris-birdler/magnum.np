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


EPS_SEEDS = (941, 942)    # the F runs with the largest energy-balance error (+0.16, +0.09; PROTOKOLL 7.108)
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
        jobs = (newell_rerun_lines() if st == "1b" else newell_eps_lines() if st == "eps" else newell_lines(int(st)))
    else:
        jobs = prec_lines() if "--prec" in sys.argv else lines()
    for n, a in jobs:
        print(n, a)
