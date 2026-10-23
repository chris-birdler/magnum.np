#!/usr/bin/env python
"""
Local test of K1 (PLAN 14, step 2): an arm of step_dump.py --k1_arm with an injected fault, to show that the
criteria see it in the window where it starts.

    python k1_inject.py --field TERM EPS --t_from T -- <step_dump.py arguments>
    python k1_inject.py --norm DELTA --t_from T -- <step_dump.py arguments>

--field: the field of the term class TERM (e.g. UniaxialAnisotropyField) times (1 + EPS) when state.t >= T.
--norm: the RKF45 error estimate times (1 + DELTA) when t >= T (the tree and the try record see the same error, so the
error norm of the tree changes by DELTA relative, and dt_opt with it).
"""
import argparse
import pathlib
import sys


def main(argv):
    i = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", nargs=2, metavar=("TERM", "EPS"))
    ap.add_argument("--norm", type=float)
    ap.add_argument("--t_from", type=float, required=True)
    a = ap.parse_args(argv[:i])
    rest = argv[i + 1:]
    j = rest.index("--")
    rl_dir = "."
    if "--rl_dir" in rest[:j]:
        rl_dir = rest[rest.index("--rl_dir") + 1]
    sys.path.insert(0, str(pathlib.Path(rl_dir).resolve()))
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import magnumnp
    if a.field:
        name, eps = a.field[0], float(a.field[1])
        cls = getattr(magnumnp, name)
        f_h = cls.h

        def h(self, state, *x, **k):
            r = f_h(self, state, *x, **k)
            return r * (1.0 + eps) if float(state.t) >= a.t_from else r

        cls.h = h
    if a.norm is not None:
        from magnumnp.solvers.ode_solvers.rkf45 import RKF45
        f_try = RKF45._try_step

        def try_step(self, t, x, **k):
            x1, t1, err = f_try(self, t, x, **k)
            return (x1, t1, err * (1.0 + a.norm)) if float(t) >= a.t_from else (x1, t1, err)

        RKF45._try_step = try_step
    import step_dump
    step_dump.main(rest)


if __name__ == "__main__":
    main(sys.argv[1:])
