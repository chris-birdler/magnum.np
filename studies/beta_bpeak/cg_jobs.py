#!/usr/bin/env python
"""
Job lines of the same-state mesh test (PROTOKOLL 7.63): dx 3 at K3 from the coarse-grained relaxed dx 1.5 state
(coarse_init.py -> runs/MK3cg_init/s<seed>.pt), 9 / 50 / 100 mT, the arguments of MK3_dx3_s<seed>_b9 with the
amplitude, --init_from, --min_calm 2 and --max_cycles_per_amp 12 (as the reruns of PLAN 4.5).

    python cg_jobs.py > jobs/meshK3cg.txt
"""
import shlex

import jobs
from rerun_jobs import get_arg, job_lines, set_arg

src = dict(job_lines("jobs/meshK3.txt"))
for s in range(2031, 2037):
    tok = shlex.split(src["MK3_dx3_s%d_b9" % s])
    assert get_arg(tok, "dx_lex") == ["3.0"] and get_arg(tok, "init_from") is None
    for mT in (9, 50, 100):
        t = set_arg(list(tok), "b_list", ["%r" % jobs.b_of_mT(mT)[0]])
        t = set_arg(t, "init_from", ["MK3cg_init/s%d.pt" % s])
        t = set_arg(t, "max_cycles_per_amp", ["12"])
        t = set_arg(t, "min_calm", ["2"])
        print("MK3cg_dx3_s%d_b%d %s" % (s, mT, " ".join(t)))
