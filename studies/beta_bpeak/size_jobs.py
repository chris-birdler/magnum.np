#!/usr/bin/env python
"""
Job lines of the size trend of the mesh error (PLAN 4.6, option 1 of Chris 2026-10-11): the base microstructure of
mesh300 (L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1, converged virgin relaxation) at d/l_ex 212 and 150, dx 1.5 and dx 3
with the same seed (paired), 9 / 50 / 100 mT. Template: the mesh300 lines of seed 921; only the name, --d_lex, --seed
and --init_from change, and --min_calm 2 --max_cycles_per_amp 12 are added (as the reruns of PLAN 4.5). The 50 and
100 mT jobs start from the relaxed state of the 9 mT job of the same seed and mesh (as mesh300).
Order: the dx 1.5 jobs at d/l_ex 212 first (longest), then dx 1.5 at 150, then dx 3.

    python size_jobs.py > jobs/meshsize.txt
"""
import shlex

from rerun_jobs import get_arg, job_lines, set_arg

SIZES = ((212, range(931, 935)), (150, range(941, 945)))
MESHES = ("dx15", "dx3")
AMPS = (9, 50, 100)
CHANGED = ("d_lex", "seed", "init_from", "init_wait_h", "min_calm", "max_cycles_per_amp")


def lines(src="jobs/mesh300.txt"):
    tmpl = dict(job_lines(src))
    out = []
    for mesh in MESHES:
        for d, seeds in SIZES:
            for mT in AMPS:
                for s in seeds:
                    tok = shlex.split(tmpl["M300_%s_s921_b%d" % (mesh, mT)])
                    t = set_arg(tok, "d_lex", ["%.1f" % d])
                    t = set_arg(t, "seed", [str(s)])
                    if get_arg(tok, "init_from") is not None:
                        t = set_arg(t, "init_from", ["M%d_%s_s%d_b9/init.pt" % (d, mesh, s)])
                    t = set_arg(t, "max_cycles_per_amp", ["12"])
                    t = set_arg(t, "min_calm", ["2"])
                    out.append(("M%d_%s_s%d_b%d" % (d, mesh, s, mT), " ".join(t)))
    return out


if __name__ == "__main__":
    for name, a in lines():
        print(name, a)
