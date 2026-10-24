#!/usr/bin/env python
"""
Local test of K1 (PLAN 14, step 2): an arm of step_dump.py --k1_arm with an injected fault, to show that the
criteria see it in the window where it starts.

    python k1_inject.py --field TERM EPS --t_from T -- <step_dump.py arguments>
    python k1_inject.py --norm DELTA --t_from T -- <step_dump.py arguments>

--field: the field of the term class TERM (e.g. UniaxialAnisotropyField) times (1 + EPS) when state.t >= T.
--norm: the RKF45 error estimate times (1 + DELTA) when t >= T (the tree and the try record see the same error, so the
error norm of the tree changes by DELTA relative, and dt_opt with it).
--spmv (detector control of PLAN 14b, D3; a tree with PairBlocks): in the calls N1 N2 ... of the product of PairBlocks
(counted only inside RKF45 tries, the second tries of step_dump --k1_redo included), the largest entry of the product
times (1 + EPS).

    python k1_inject.py --spmv EPS --calls N1 N2 ... -- <step_dump.py arguments>
"""
import argparse
import pathlib
import sys


def main(argv):
    i = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", nargs=2, metavar=("TERM", "EPS"))
    ap.add_argument("--norm", type=float)
    ap.add_argument("--spmv", type=float)
    ap.add_argument("--calls", type=int, nargs="+", default=[])
    ap.add_argument("--t_from", type=float)
    a = ap.parse_args(argv[:i])
    if (a.field or a.norm is not None) and a.t_from is None:
        ap.error("--field and --norm need --t_from")
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
    if a.spmv is not None:
        import magnumnp.field_terms.demag_staircase as ds
        from magnumnp.solvers.ode_solvers.rkf45 import RKF45
        if not hasattr(ds, "PairBlocks"):
            raise SystemExit("k1_inject: --spmv needs a tree with PairBlocks")
        f_mm, f_try_s = ds.PairBlocks.__matmul__, RKF45._try_step
        st, calls = {"in_try": 0, "n": 0}, set(a.calls)

        def try_step_s(self, t, x, **k):
            st["in_try"] += 1
            try:
                return f_try_s(self, t, x, **k)
            finally:
                st["in_try"] -= 1

        def matmul(self, x):
            r = f_mm(self, x)
            if st["in_try"]:
                st["n"] += 1
                if st["n"] in calls:
                    j = int(r.abs().argmax())
                    r = r.clone()
                    r[j] = r[j] * (1.0 + a.spmv)
                    print("[k1_inject] product call %d inside the tries: entry %d times (1 + %g)" % (st["n"], j, a.spmv),
                          flush=True)
            return r

        RKF45._try_step = try_step_s
        ds.PairBlocks.__matmul__ = matmul
    import step_dump
    step_dump.main(rest)


if __name__ == "__main__":
    main(sys.argv[1:])
