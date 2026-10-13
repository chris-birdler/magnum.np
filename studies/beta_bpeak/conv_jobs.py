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
no accumulators, no snapshots); only the name, --d_lex, --dx_lex, --seed, --init_transfer and --init_wait_h change.
Order: the four F runs first (longest, one per GPU), then I, then T (the T runs need the init.pt of their F run;
--init_wait_h covers a slow F relaxation).

    python conv_jobs.py > jobs/conv_b0.txt
"""
import shlex

from rerun_jobs import job_lines, set_arg

SEEDS = range(941, 945)
D_LEX = 96.0
ARMS = (("F", 0.75, False), ("I15", 1.5, False), ("I3", 3.0, False), ("T15", 1.5, True), ("T3", 3.0, True))
WAIT_H = 40.0
CHANGED = ("d_lex", "dx_lex", "seed", "init_transfer", "init_wait_h")


def name(arm, s):
    return "B0_%s_s%d" % (arm, s)


def lines(src="jobs/mesh212.txt"):
    tok0 = shlex.split(dict(job_lines(src))["M212_dx15_s931_b9"])
    out = []
    for arm, dx, transfer in ARMS:
        for s in SEEDS:
            t = set_arg(tok0, "d_lex", ["%.1f" % D_LEX])
            t = set_arg(t, "dx_lex", [repr(dx)])
            t = set_arg(t, "seed", [str(s)])
            if transfer:
                t = set_arg(t, "init_transfer", ["%s/init.pt" % name("F", s)])
                t = set_arg(t, "init_wait_h", ["%g" % WAIT_H])
            out.append((name(arm, s), " ".join(t)))
    return out


if __name__ == "__main__":
    for n, a in lines():
        print(n, a)
