#!/usr/bin/env python
"""
K2 of PLAN 14 (single-precision kernels on the checkpoint states) and the data of K3: run run_loops.main with the
given arguments (the code path of production); at the first call of LLGSolver.step, for each checkpoint state of the
K1 reference file (step_dump.py --k1_ref), evaluate on m32 (the state cast to single precision; in a double-precision
run m32.double(), so that the rounding of the input is not in the error):
  - the field of each term, in the order of production (for the trees with B1: DemagField.h before the correction;
    the script records for each call of the correction whether it took the stored stray field, DemagField.last_h),
  - the effective field (the sum of the term fields in the order of LLGSolver.dm),
  - the LLG right-hand side (LLGSolver.dm of the tree),
  - one RKF45 try (_try_step with _dt and t of the reference, not step()): the six stage values (right-hand sides), the
    increment x1 - x (formed in double precision), the error estimate, its RMS and the error norm of the tree,
then stop. Without --against (the double-precision reference trees) it writes all quantities, one file for each
state, into --out. With --against (a folder of a reference tree) it writes for each quantity and state the error
e = RMS |q - q_ref| / RMS |q_ref| and the bias b = sum (q - q_ref) . q_ref / sum q_ref . q_ref (median over --reps
repetitions); a quantity with q_ref = 0 gets the RMS of q - q_ref (key "abs") and no ratio.

    CUDA_DEVICE=0 python s1_states.py --ref ref.pt --out k2/d538bde_double --kcache kc/B_double -- <run_loops arguments>
    CUDA_DEVICE=0 python s1_states.py --ref ref.pt --out k2/c56c0a5.pt --against k2/d538bde_double --reps 3 \
        --kcache kc/B_single -- <run_loops arguments>
"""
import argparse
import inspect
import pathlib
import statistics
import sys

import torch

OLD_NORM = "torch.linalg.norm((err / (atol + rtol*x.abs())).flatten(), torch.inf)"


def error_norm(solver, x, err):
    import magnumnp.solvers.ode_solvers.rkf45 as rk
    atol, rtol = solver._atol, solver._rtol
    if hasattr(rk, "_error_norm"):
        return float(rk._error_norm(err, x, atol, rtol))
    if OLD_NORM not in inspect.getsource(rk.RKF45._optimal_stepsize):
        raise SystemExit("s1_states: unknown error norm in RKF45._optimal_stepsize")
    return float(torch.linalg.norm((err / (atol + rtol*x.abs())).flatten(), torch.inf))


def evaluate(llg, state, c, last):
    """All quantities of one state c (m, t, dt of the reference) on m32."""
    solver = llg._solver
    x = c["m"].float().to(device=state.m.device, dtype=state.m.dtype)
    t = torch.tensor(c["t"], dtype=torch.float64, device=state.m.device)
    q = {}
    state.m, state.t = x.clone(), t.clone()
    last.clear()
    hs = [term.h(state) for term in llg._terms]
    for i, (term, h) in enumerate(zip(llg._terms, hs)):
        q["h%d_%s" % (i, type(term).__name__)] = h
    h = hs[0]
    for f in hs[1:]:
        h = h + f
    q["h_eff"] = h
    state.m, state.t = x.clone(), t.clone()
    q["rhs"] = llg.dm(state.t, state.m, state=state)
    state.m, state.t = x.clone(), t.clone()
    solver._dt = c["dt"]
    f0, stages = solver._f, []

    def f(tt, xx, **kw):
        r = f0(tt, xx, **kw)
        stages.append(r)
        return r

    solver._f = f
    try:
        x0 = state.m
        x1, t1, err = solver._try_step(state.t, x0, state=state)
    finally:
        solver._f = f0
    for i, s in enumerate(stages):
        q["stage%d" % (i + 1)] = s
    q["dx"] = x1.double() - x0.double()
    q["err"] = err
    sc = {"err_rms": float(err.double().pow(2).mean().sqrt()), "norm": error_norm(solver, x0, err)}
    return q, sc, list(last)


def compare(q, ref):
    out = {}
    for k, v in q.items():
        r = ref[k].to(device=v.device, dtype=torch.float64)
        d = v.double() - r
        rr = float(r.pow(2).sum())
        if rr == 0:
            out[k] = {"abs": float(d.pow(2).mean().sqrt())}
        else:
            out[k] = {"e": float(d.pow(2).sum().sqrt()) / rr**0.5, "b": float((d * r).sum()) / rr}
    return out


def main(argv):
    i = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="K1 reference file with the checkpoint states")
    ap.add_argument("--out", required=True, help="folder (reference tree) or file (with --against)")
    ap.add_argument("--against", help="folder of the double-precision reference tree")
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--kcache", help="kernel cache folder (one for each geometry and precision)")
    ap.add_argument("--rl_dir", default=".", help="folder of run_loops.py (default: the current folder)")
    a = ap.parse_args(argv[:i])
    sys.path.insert(0, str(pathlib.Path(a.rl_dir).resolve()))
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import magnumnp
    import run_loops
    if a.kcache:
        import kcache
        kcache.use(a.kcache)
    ref = torch.load(a.ref, map_location="cpu")
    f_step = magnumnp.LLGSolver.step

    def step(self, state, dt, **k):
        last = []
        dem = next((t for t in self._terms if type(t).__name__ == "DemagField"), None)
        for term in self._terms:
            if hasattr(term, "_demag"):                      # the staircase correction: did it take last_h?
                f_h = term.h

                def h_corr(st, f_h=f_h):
                    last.append(hasattr(dem, "last_h") and dem.last_h(st) is not None)
                    return f_h(st)

                term.h = h_corr
        names = [type(t).__name__ for t in self._terms]
        Ms = state.material["Ms"]
        meta = {"names": names, "dtype": str(state.m.dtype), "argv": argv[i + 1:], "ck": sorted(ref["ck"]),
                "reps": a.reps, "Ms": (torch.broadcast_to(Ms, state.m.shape[:3] + (1,))[..., 0].double().cpu()
                                       if torch.is_tensor(Ms) else float(Ms)), "Ms_full": float(state.material["Ms"].max())
                if torch.is_tensor(Ms) else float(Ms)}
        res = {"meta": meta, "states": {}}
        if a.against:
            ameta = torch.load(pathlib.Path(a.against) / "meta.pt", map_location="cpu")
            if ameta["names"] != names:
                raise SystemExit("s1_states: terms %s, reference %s" % (names, ameta["names"]))
        else:
            out = pathlib.Path(a.out)
            out.mkdir(parents=True, exist_ok=True)
            torch.save(meta, out / "meta.pt")
        for k0 in sorted(ref["ck"]):
            c = ref["ck"][k0]
            if a.against:
                rq = torch.load(pathlib.Path(a.against) / ("state_%03d.pt" % k0), map_location="cpu")
                reps = []
                for _ in range(a.reps):
                    q, sc, used = evaluate(self, state, c, last)
                    cmp = compare(q, rq["q"])
                    for kk, v in sc.items():
                        cmp[kk] = {"e": abs(v / rq["sc"][kk] - 1.0), "b": v / rq["sc"][kk] - 1.0}
                    reps.append((cmp, used))
                    del q
                med = {}
                for kk in reps[0][0]:
                    med[kk] = {f: statistics.median(r[0][kk][f] for r in reps) for f in reps[0][0][kk]}
                res["states"][k0] = {"q": med, "used_last": [r[1] for r in reps],
                                     "reps": [r[0] for r in reps]}
            else:
                q, sc, used = evaluate(self, state, c, last)
                torch.save({"q": {kk: v.detach().cpu() for kk, v in q.items()}, "sc": sc, "used_last": used},
                           pathlib.Path(a.out) / ("state_%03d.pt" % k0))
                res["states"][k0] = {"sc": sc, "used_last": used}
                del q
            print("[s1_states] state k0 %d done (used last_h: %s)" % (k0, res["states"][k0]["used_last"]), flush=True)
        if a.against:
            torch.save(res, a.out)
        print("[s1_states] %d states, terms %s, dtype %s; saved to %s" % (len(ref["ck"]), names, state.m.dtype, a.out),
              flush=True)
        raise SystemExit(0)

    magnumnp.LLGSolver.step = step
    run_loops.main(argv[i + 1:])
    print("[s1_states] the run ended before the first sample interval", flush=True)
    sys.exit(2)


if __name__ == "__main__":
    main(sys.argv[1:])
