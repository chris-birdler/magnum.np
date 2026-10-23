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
  - the energy of each term (term.E, as the diagnostics of run_loops call it),
  - for the trees with B1: the correction without the stored stray field (DemagField.last_h returns None), which must
    equal the correction of the old code (T6 of PLAN 14a),
  - one RKF45 try (_try_step at t of the reference, not step()) with a typical dt (the median dt of the accepted tries
    in interval k0 + 1 of the reference, without the tries clamped to the rest of the interval): the six stage values
    f_i (right-hand sides, not k_i = dt f_i), the increment x1 - x (formed in double precision), the error estimate,
    its RMS and the error norm of the tree,
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


def typical_dt(ref, k0):
    """dt of the try: the median dt of the accepted tries in interval k0 + 1 of the reference window, without the
    tries clamped to the rest of the interval (their dt is a sliver; the error estimate scales as dt^5, its round-off
    as dt). Falls back to the checkpoint dt."""
    w = ref["win"].get(k0, {})
    dts = sorted(r[1] for r in w.get("rows", []) if r[0] == k0 + 1 and r[5] and r[1] != r[4])
    return dts[len(dts) // 2] if dts else ref["ck"][k0]["dt"]


def evaluate(llg, state, c, last, null=False):
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
    sc = {}
    state.m, state.t = x.clone(), t.clone()
    for i, term in enumerate(llg._terms):                 # energies on one m (E of the correction may use last_h)
        sc["E%d_%s" % (i, type(term).__name__)] = float(term.E(state))
    from magnumnp.field_terms.demag import DemagField
    for i, term in enumerate(llg._terms):
        if hasattr(term, "_demag") and hasattr(DemagField, "last_h"):
            # B1 identity (T6): the correction without the stored stray field must be the old correction
            f_last = DemagField.last_h
            DemagField.last_h = lambda self, st: None
            try:
                state.m, state.t = x.clone(), t.clone()
                q["h%d_%s_noh1" % (i, type(term).__name__)] = term.h(state)
            finally:
                DemagField.last_h = f_last
    hh = hs[::-1] if null else hs                     # null run: the field sum in the reverse order
    h = hh[0]
    for f in hh[1:]:
        h = h + f
    q["h_eff"] = h
    state.m, state.t = x.clone(), t.clone()
    q["rhs"] = llg.dm(state.t, state.m, state=state)
    state.m, state.t = x.clone(), t.clone()
    solver._dt = c["dt_try"]
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
    sc.update(err_rms=float(err.double().pow(2).mean().sqrt()), norm=error_norm(solver, x0, err), dt_try=c["dt_try"])
    return q, sc, list(last)


def compare(q, ref):
    """e, b, and the RMS of the difference (abs) and of the reference (ref) in the units of q."""
    out = {}
    for k, v in q.items():
        r = ref[k].to(device=v.device, dtype=torch.float64)
        d = v.double() - r
        rr = float(r.pow(2).sum())
        out[k] = {"abs": float(d.pow(2).mean().sqrt()), "ref": float(r.pow(2).mean().sqrt())}
        if rr > 0:
            out[k].update(e=float(d.pow(2).sum().sqrt()) / rr**0.5, b=float((d * r).sum()) / rr)
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
    ap.add_argument("--null", action="store_true", help="null run: the field sum in the reverse order (trees without A1)")
    ap.add_argument("--k0_from", type=int, default=8, help="first checkpoint state (PLAN 14a: 8)")
    a = ap.parse_args(argv[:i])
    sys.path.insert(0, str(pathlib.Path(a.rl_dir).resolve()))
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import magnumnp
    import run_loops
    if a.kcache:
        import kcache
        kcache.use(a.kcache)
    if a.null:
        import step_dump
        step_dump.k1_null("reverse")
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
        meta = {"names": names, "dtype": str(state.m.dtype), "argv": argv[i + 1:], "ck": sorted(k for k in ref["ck"] if k >= a.k0_from),
                "reps": a.reps, "Ms": (torch.broadcast_to(Ms, state.m.shape[:3] + (1,))[..., 0].double().cpu()
                                       if torch.is_tensor(Ms) else float(Ms)), "Ms_full": float(state.material["Ms"].max())
                if torch.is_tensor(Ms) else float(Ms)}
        drive = [t for t in self._terms if hasattr(t, "H_amp")]
        meta["H_amp"] = float(drive[0].H_amp) if drive else None
        al = state.material["alpha"]
        meta["alpha"] = (torch.broadcast_to(al, state.m.shape[:3] + (1,))[..., 0].double().cpu()
                         if torch.is_tensor(al) else float(al))
        from magnumnp.common import constants
        meta["gamma"] = float(constants.gamma)
        meta["null"] = a.null
        cv = state.mesh.cell_volumes
        meta["volume"] = float(cv.sum()) if torch.is_tensor(cv) and cv.numel() > 1 else float(cv) * state.m[..., 0].numel()
        res = {"meta": meta, "states": {}}
        if a.against:
            ameta = torch.load(pathlib.Path(a.against) / "meta.pt", map_location="cpu")
            if ameta["names"] != names:
                raise SystemExit("s1_states: terms %s, reference %s" % (names, ameta["names"]))
        else:
            out = pathlib.Path(a.out)
            out.mkdir(parents=True, exist_ok=True)
            torch.save(meta, out / "meta.pt")
        for k0 in sorted(k for k in ref["ck"] if k >= a.k0_from):
            c = dict(ref["ck"][k0], dt_try=typical_dt(ref, k0))
            if a.against:
                rq = torch.load(pathlib.Path(a.against) / ("state_%03d.pt" % k0), map_location="cpu")
                reps = []
                for _ in range(a.reps):
                    q, sc, used = evaluate(self, state, c, last, a.null)
                    cmp = compare(q, rq["q"])
                    for kk, v in sc.items():
                        if kk == "dt_try":
                            continue
                        r0 = rq["sc"][kk]
                        cmp[kk] = ({"e": abs(v / r0 - 1.0), "b": v / r0 - 1.0, "abs": abs(v - r0), "ref": abs(r0)}
                                   if r0 != 0 else {"abs": abs(v - r0), "ref": 0.0})
                    reps.append((cmp, used))
                    del q
                med = {}
                for kk in reps[0][0]:
                    med[kk] = {f: statistics.median(r[0][kk][f] for r in reps) for f in reps[0][0][kk]}
                res["states"][k0] = {"q": med, "used_last": [r[1] for r in reps],
                                     "reps": [r[0] for r in reps]}
            else:
                q, sc, used = evaluate(self, state, c, last, a.null)
                torch.save({"q": {kk: v.detach().cpu() for kk, v in q.items()}, "sc": sc, "used_last": used,
                            "dt": c["dt"], "t": c["t"]},
                           pathlib.Path(a.out) / ("state_%03d.pt" % k0))
                res["states"][k0] = {"sc": sc, "used_last": used}
                del q
            print("[s1_states] state k0 %d done (used last_h: %s)" % (k0, res["states"][k0]["used_last"]), flush=True)
        if a.against:
            torch.save(res, a.out)
        print("[s1_states] %d states, terms %s, dtype %s; saved to %s" % (len(res["states"]), names, state.m.dtype, a.out),
              flush=True)
        raise SystemExit(0)

    magnumnp.LLGSolver.step = step
    run_loops.main(argv[i + 1:])
    print("[s1_states] the run ended before the first sample interval", flush=True)
    sys.exit(2)


if __name__ == "__main__":
    main(sys.argv[1:])
