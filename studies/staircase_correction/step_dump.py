#!/usr/bin/env python
"""
Check S2 of PLAN 11 (bit equality of a speed patch): run run_loops.main with the given arguments (the code path of
production) and save m after the given calls of LLGSolver.step (one call integrates one sample interval), then stop.
Two runs of the same code must give bitwise equal files (precondition); a patch that keeps the order of all
operations must give the same files as the old code.

    CUDA_DEVICE=0 python step_dump.py --at 1 256 --save m.pt -- <run_loops arguments>
    CUDA_DEVICE=0 python step_dump.py --at 256 --save m.pt --go_on -- <run_loops arguments with --max_cycles>
    CUDA_DEVICE=0 python step_dump.py --at 1 2 4 --trace --save m.pt --go_on -- <run_loops arguments>
    python step_dump.py --compare a.pt b.pt

--trace (precision test, PLAN 13) also saves after each call: the counts of RKF45 tries and of accepted steps, t,
and m in 2048 probe cells (magnetic cells, fixed choice with seed 0; the same cells in all arms of one geometry).

K1 of PLAN 14 (windows over the cycle, resync):

    CUDA_DEVICE=0 python step_dump.py --k1_ref ref.pt -- <run_loops arguments with --max_cycles 1>
    CUDA_DEVICE=0 python step_dump.py --k1_arm ref.pt --k1_out arm.pt [--k1_null] -- <run_loops arguments>

--k1_ref integrates the full cycle. Before sample interval k0 + 1 of each checkpoint k0 (--k1_ck) it saves m, t and
the RKF45 time step _dt; after each sample interval it saves t and _dt. --k1_arm reads that file: before interval
k0 + 1 it sets m (copy_ into state.m, so that the stray-field cache of B1, keyed by the tensor and its version, becomes
invalid), t and _dt to the values of the reference; it integrates only the windows k0 + 1 ... k0 + --k1_span and
outside them sets t to the value of the reference at the end of the interval. Both record in each window: the try
record (below), m in the probe cells after each accepted step, the full m after the intervals k0 + --k1_full, and t
and _dt after each interval. The reference records the tries of the full cycle (count of the A4 case, PLAN 14).
--k1_null sums the field terms in the reverse order (null run; only trees without A1, where LLGSolver.dm sums with
accumulate_h); the terms are evaluated in the order of production. --k1_null rotate sums them in the order rotated by
one (a second null run, for the scatter between two null runs). --k1_null kernel multiplies the Newell kernel by
(1 + 2^-52) (a coherent null run). In trees with B1 the run counts the calls of DemagField.last_h with and without the
stored field.

Try record, for each RKF45 try: the sample interval, dt, the error norm as the tree computes it (_error_norm of A3, else
the expression of _optimal_stepsize), dt_opt, and the branch tuple of rkf45.py: (decision, A4 clamp of dt to the rest
of the interval before the try, try past the end, norm class 0: below 0.5, 1: 0.5 ... 1.1, 2: above 1.1, clamp at
minscale, clamp at maxscale, clamp of r at 1, clamp at maxstep, retry 0: none, 1: dt_opt, 2: the rest). The wrappers
wrap only _try_step and the norm functions; the step control of the tree runs unchanged. The wrapper repeats the
control logic of rkf45.py to name the branch and checks it against the tree (the next dt, the _dt after the call, the
returned dt_opt); a difference counts in "mismatch" and makes the exit code 3.

Outlier test (finding of the K1 pilot, PLAN 14a): --k1_reps R integrates each window of an arm R times from the same
start state (m, t and _dt of the reference); the rows of the repetitions 1 ... R-1 go to win[k0]["reps"]. --k1_redo
adds to each row the RMS of the scaled error, err / (atol + rtol |x|), and the flat index of its largest entry; it
does each try with a norm above 1.1 once more with the same input (t, x, dt; _try_step has no side effects) and adds
the norm, the RMS and the index of the second try. A transient fault gives a normal second try. --k1_redo all does
each try once more and compares the second try with the first: the error estimate in tolerance units (largest
difference and its flat index) and x1 (largest difference). The output of each field term in each stage of the first
try is copied into a fixed buffer; at a difference above THR_ERR tolerance units or THR_X1 in x1 (an event), the run
records for each stage and term the largest difference of the second try to the first. The first stage and term with
a difference above the rounding noise names the part with the fault (later stages differ because their input differs).
"""
import argparse
import inspect
import pathlib
import sys

import torch

# branch tuple fields (see the docstring)
BRANCH = ("accept", "a4_clamp", "past_end", "norm_class", "c_minscale", "c_maxscale", "c_one", "c_maxstep", "retry")
REDO = ("err_rms", "imax", "redo_norm", "redo_err_rms", "redo_imax")
REDO_ALL = ("redo_derr", "redo_dx1", "redo_dimax")
THR_ERR, THR_X1 = 1e-6, 1e-12          # event thresholds of --k1_redo all (as glitch_stress.py)
OLD_NORM = "torch.linalg.norm((err / (atol + rtol*x.abs())).flatten(), torch.inf)"


def compare(a, b):
    A, B = torch.load(a, map_location="cpu"), torch.load(b, map_location="cpu")
    ok = True
    for k in sorted(A["m"]):
        eq = torch.equal(A["m"][k], B["m"][k])
        d = float((A["m"][k].double() - B["m"][k].double()).abs().max())
        print("after %4d sample intervals: bitwise equal %s, max |dm| %.3e, t %r / %r"
              % (k, eq, d, A["t"][k], B["t"][k]))
        ok &= eq and A["t"][k] == B["t"][k]
    print("RESULT bitwise equal: %s" % ok)
    return ok


def probe_cells(m):
    """2048 magnetic cells, fixed choice with seed 0 (the cells of section 13)."""
    mag = (m.detach().norm(dim=-1) > 0).reshape(-1).cpu().nonzero().squeeze(1)
    g = torch.Generator().manual_seed(0)
    # device: magnum.np sets the default device to cuda on the GPU, the generator is on the CPU
    return mag[torch.randperm(len(mag), generator=g, device="cpu")[:2048]].sort().values


class TryRecord:
    """Try record of RKF45 (PLAN 14). on: record the tries of the current LLGSolver.step call."""

    def __init__(self):
        import magnumnp.solvers.ode_solvers.rkf45 as rk
        self.rk, self.on, self.n, self.probe_on, self.probe_dev = rk, False, 0, False, None
        self.rows, self.probe_m, self.mismatch, self.accepted = [], [], 0, 0
        self.remaining = self.dt_expect = None
        self.redo, self.args = False, None
        self.cap, self.calls, self.buf, self.cmp, self.events, self.names = None, [], {}, [], [], []
        self.has_a3 = hasattr(rk, "_error_norm")
        f_try, f_opt = rk.RKF45._try_step, rk.RKF45._optimal_stepsize
        self.f_try = f_try
        if self.has_a3:
            f_norm = self.f_norm = rk._error_norm

            def error_norm(err, x, atol, rtol):
                r = f_norm(err, x, atol, rtol)
                if self.on:
                    self.norm = float(r)
                return r

            rk._error_norm = error_norm
        elif OLD_NORM not in inspect.getsource(f_opt):
            raise SystemExit("step_dump: unknown error norm in RKF45._optimal_stepsize")
        rec = self

        def try_step(solver, t, x, **kw):
            if not rec.on:
                return f_try(solver, t, x, **kw)
            dt = solver._dt
            a4 = dt != rec.dt_expect
            if a4 and not (rec.dt_expect > rec.remaining and dt == rec.remaining):
                rec.mismatch += 1
            if rec.redo == "all":
                rec.cap, rec.calls = "store", [0] * len(rec.names)
            x1, t1, err = f_try(solver, t, x, **kw)
            rec.cap = None
            rec.pending = (dt, a4, x1)
            if rec.redo:
                rec.args = (solver, t, x, kw)
            return x1, t1, err

        def optimal_stepsize(solver, x, err, atol, rtol):
            dt_opt = f_opt(solver, x, err, atol, rtol)
            if rec.on:
                if rec.has_a3:
                    norm = rec.norm
                else:
                    norm = float(torch.linalg.norm((err / (atol + rtol*x.abs())).flatten(), torch.inf))
                rec.control(solver, norm, dt_opt)
                if rec.redo:
                    rec.redo_try(x, err, atol, rtol, norm)
            return dt_opt

        rk.RKF45._try_step = try_step
        rk.RKF45._optimal_stepsize = optimal_stepsize

    def control(self, s, norm, dt_opt):
        """Repeat the step control of rkf45.py for one try; name its branch; check it against the tree."""
        dt, a4, x1 = self.pending
        rem = self.remaining
        c_min = c_max = c_one = False
        if norm > 1.1:
            cls = 2
            r = s._headroom / norm**(1.0/s._order)
            if r < s._minscale:
                r, c_min = s._minscale, True
        elif norm < 0.5:
            cls = 0
            if norm > 0.:
                r = s._headroom / norm**(1.0/(s._order+1.0))
            else:
                r = s._maxscale
            if r > s._maxscale:
                r, c_max = s._maxscale, True
            if r < 1.:
                r, c_one = 1., True
        else:
            cls, r = 1, 1.
        opt = s._dt * r
        c_step = opt > s._maxstep
        if c_step:
            opt = s._maxstep
        if opt != dt_opt:
            self.mismatch += 1
        past = dt > rem
        accept = not (dt > dt_opt or dt > rem)
        if accept:
            retry = 0
            self.remaining = rem - dt
            self.dt_expect = dt_opt
            self.accepted += 1
            if self.probe_on:
                self.probe_m.append(x1.detach().reshape(-1, 3)[self.probe_dev].cpu().clone())
        else:
            retry = 1 if dt_opt < rem else 2
            self.dt_expect = min(dt_opt, rem)
        self.rows.append((self.n, dt, norm, dt_opt, rem, accept, a4, past, cls, c_min, c_max, c_one, c_step, retry))

    def scaled(self, x, err, atol, rtol):
        """Norm (as the tree computes it), RMS and flat index of the largest entry of err / (atol + rtol |x|)."""
        e = (err / (atol + rtol * x.abs())).flatten()
        a = e.abs()
        norm = float(self.f_norm(err, x, atol, rtol)) if self.has_a3 else float(torch.linalg.norm(e, torch.inf))
        return norm, float(e.square().mean().sqrt()), int(a.argmax())

    def wrap_terms(self, terms):
        """--k1_redo all: wrap h of each field term. Store: copy the output into the buffer of (stage, term). Compare:
        record the largest difference to that buffer, its flat index and the largest value of the buffer."""
        self.names = [type(t).__name__ for t in terms]
        rec = self
        for j, term in enumerate(terms):
            def wrap(f, j):
                def h(state):
                    r = f(state)
                    if rec.cap is not None:
                        s = rec.calls[j]
                        rec.calls[j] += 1
                        b = rec.buf.get((s, j))
                        if rec.cap == "store":
                            if b is None or b.shape != r.shape or b.dtype != r.dtype:
                                b = rec.buf[(s, j)] = torch.empty_like(r)
                            b.copy_(r)
                        elif b is not None:
                            d = (r - b).abs()
                            rec.cmp.append((s, j, d.max(), d.argmax(), b.abs().max()))
                    return r
                return h
            term.h = wrap(term.h, j)

    def redo_try(self, x, err, atol, rtol, norm):
        """--k1_redo: add (RMS, index of the largest entry, and for a norm above 1.1 the norm, RMS and index of a
        second try with the same input) to the last row. --k1_redo all: a second try for each try; add also the
        largest difference of the error estimate in tolerance units, of x1, and the flat index of the first; record an
        event with the differences of each field term in each stage."""
        _, rms, imax = self.scaled(x, err, atol, rtol)
        n2 = r2 = float("nan")
        i2 = -1
        if self.redo == "all":
            solver, t, x0, kw = self.args
            x1 = self.pending[2]
            self.cap, self.calls, self.cmp = "compare", [0] * len(self.names), []
            x12, _, err2 = self.f_try(solver, t, x0, **kw)
            self.cap = None
            n2, r2, i2 = self.scaled(x0, err2, atol, rtol)
            de = ((err2 - err).abs() / (atol + rtol * x0.abs())).flatten()
            derr, dimax, dx1 = float(de.max()), int(de.argmax()), float((x12 - x1).abs().max())
            if derr > THR_ERR or dx1 > THR_X1:
                c = torch.stack([torch.stack([d.double(), i.double(), b.double()]) for _, _, d, i, b in self.cmp]).cpu()
                terms = [{"stage": s, "term": self.names[j], "dev": float(c[k, 0]), "index": int(c[k, 1]),
                          "scale": float(c[k, 2])} for k, (s, j, _, _, _) in enumerate(self.cmp)]
                ev = {"row": len(self.rows) - 1, "n": self.n, "dt": self.pending[0], "norm": norm, "redo_norm": n2,
                      "derr": derr, "dimax": dimax, "dx1": dx1, "imax": imax, "redo_imax": i2, "terms": terms}
                self.events.append(ev)
                print("[step_dump] redo event: n %d dt %.4g norm %.6g, second try %.6g; error estimate difference "
                      "%.3e tolerance units at %d, x1 difference %.3e" % (self.n, ev["dt"], norm, n2, derr, dimax, dx1),
                      flush=True)
                for tr in terms:
                    print("   stage %d %-26s largest difference %.3e (relative %.3e) at %d" % (
                        tr["stage"], tr["term"], tr["dev"], tr["dev"] / tr["scale"] if tr["scale"] > 0 else float("nan"),
                        tr["index"]), flush=True)
            self.rows[-1] = self.rows[-1] + (rms, imax, n2, r2, i2, derr, dx1, dimax)
            return
        if norm > 1.1:
            solver, t, x0, kw = self.args
            _, _, err2 = self.f_try(solver, t, x0, **kw)
            n2, r2, i2 = self.scaled(x0, err2, atol, rtol)
        self.rows[-1] = self.rows[-1] + (rms, imax, n2, r2, i2)


def k1_null(kind):
    """Null run: the field terms in the order of production, summed in the reverse order (kind "reverse") or in the
    order rotated by one, h1 + h2 + ... + h0 (kind "rotate"; a second null run of another order); or the Newell kernel
    times (1 + 2^-52) (kind "kernel": a coherent change at the round-off level of a fixed operator, as a set-up patch
    such as A2 makes)."""
    if kind == "kernel":
        from magnumnp.field_terms.demag import DemagField
        f_init_N = DemagField._init_N

        def init_N(self, state):
            return [[K * (1.0 + 2.0 ** -52) for K in row] for row in f_init_N(self, state)]

        DemagField._init_N = init_N
        return
    import magnumnp.solvers.llg as llg
    if hasattr(llg, "_llg_rhs_sum") or "accumulate_h(term.h(state)" not in inspect.getsource(llg.LLGSolver.dm):
        raise SystemExit("step_dump: --k1_null needs a tree without A1 (LLGSolver.dm with accumulate_h)")
    f_acc = llg.accumulate_h

    def accumulate_reordered(fields):
        hs = list(fields)
        return f_acc(hs[::-1] if kind == "reverse" else hs[1:] + hs[:1])

    llg.accumulate_h = accumulate_reordered


def k1_main(a, rl):
    import magnumnp
    import run_loops
    ref = torch.load(a.k1_arm, map_location="cpu") if a.k1_arm else None
    ck = tuple(ref["ck_list"]) if ref else tuple(a.k1_ck)
    span, full = (ref["span"], tuple(ref["full"])) if ref else (a.k1_span, tuple(a.k1_full))
    if a.k1_null:
        k1_null(a.k1_null)
    rec = TryRecord()
    rec.redo = a.k1_redo or False
    lh = {"hit": 0, "miss": 0}
    from magnumnp.field_terms.demag import DemagField
    if hasattr(DemagField, "last_h"):                     # B1: count the uses of the stored stray field
        f_last = DemagField.last_h

        def last_h(self, state):
            r = f_last(self, state)
            lh["hit" if r is not None else "miss"] += 1
            return r

        DemagField.last_h = last_h
    out = {"mode": "arm" if ref else "ref", "argv": rl, "ck_list": list(ck), "span": span, "full": list(full),
           "null": a.k1_null, "has_a3": rec.has_a3, "win": {}, "t_end": {}, "dt_end": {}, "reps": a.k1_reps}
    if ref:
        out["ref"] = str(a.k1_arm)
    else:
        out["ck"] = {}
    f_step = magnumnp.LLGSolver.step
    st = {"n": 0, "call": {}}

    def window(n):
        for k0 in ck:
            if k0 < n <= k0 + span:
                return k0
        return None

    def step(self, state, dt, **k):
        st["n"] += 1
        n = st["n"]
        st["call"][n] = (dt, k)
        solver = self._solver
        if n == 1:
            if ref:
                probe = ref["probe_idx"]
            else:
                probe = out["probe_idx"] = probe_cells(state.m)
                Ms = state.material["Ms"]
                Ms = torch.broadcast_to(Ms, state.m.shape[:3] + (1,)) if torch.is_tensor(Ms) else None
                out["particle"] = None if Ms is None else (Ms[..., 0] > Ms.max() * 1e-6).cpu()
                out["t_end"][0], out["dt_end"][0] = float(state.t), float(solver._dt)
            rec.probe_dev = probe.to(state.m.device)
            if rec.redo == "all":
                rec.wrap_terms(self._terms)
                print("[step_dump] --k1_redo all: field terms %s" % rec.names, flush=True)
        k0 = window(n)
        if ref and k0 is None:
            # outside the windows: no integration, t of the reference at the end of the interval (after the end of
            # a shortened reference, --k1_stop: t + dt)
            t1 = ref["t_end"][n] if n in ref["t_end"] else float(state.t) + float(dt)
            state.t = torch.tensor(t1, dtype=state.t.dtype, device=state.t.device)
            return None
        if k0 is not None and n == k0 + 1:
            if ref:
                c = ref["ck"][k0]
                state.m.copy_(c["m"].to(device=state.m.device, dtype=state.m.dtype))
                state.t = torch.tensor(c["t"], dtype=state.t.dtype, device=state.t.device)
                solver._dt = c["dt"]
            else:
                out["ck"][k0] = {"m": state.m.detach().cpu().clone(), "t": float(state.t), "dt": float(solver._dt)}
            out["win"][k0] = {"m": {}, "t_end": {}, "dt_end": {}}
        rec.n, rec.on = n, True
        rec.remaining, rec.dt_expect = float(dt), solver._dt
        i0, p0, acc0, step0 = len(rec.rows), len(rec.probe_m), rec.accepted, state._step
        rec.probe_on = k0 is not None
        r = f_step(self, state, dt, **k)
        rec.on = rec.probe_on = False
        if solver._dt != rec.dt_expect or state._step - step0 != rec.accepted - acc0:
            rec.mismatch += 1
        t1, dt1 = float(state.t), float(solver._dt)
        if not ref:
            out["t_end"][n], out["dt_end"][n] = t1, dt1
        if k0 is not None:
            w = out["win"][k0]
            w.setdefault("rows", []).extend(rec.rows[i0:])
            w.setdefault("probe_m", []).extend(rec.probe_m[p0:])
            del rec.probe_m[p0:]
            w["t_end"][n], w["dt_end"][n] = t1, dt1
            if n - k0 in full:
                w["m"][n - k0] = state.m.detach().cpu().clone()
            if ref and n == k0 + span and a.k1_reps > 1:
                repeat(self, state, k0, w)
        if not ref and n == a.k1_stop:
            finish()
        return r

    def repeat(self, state, k0, w):
        """--k1_reps: integrate the window again from the start state of the reference; record the rows of each
        repetition and the largest |m - m| against the first integration at the end of the window."""
        solver = self._solver
        m_end = state.m.detach().clone()
        c = ref["ck"][k0]
        w["reps"] = []
        for _ in range(1, a.k1_reps):
            state.m.copy_(c["m"].to(device=state.m.device, dtype=state.m.dtype))
            state.t = torch.tensor(c["t"], dtype=state.t.dtype, device=state.t.device)
            solver._dt = c["dt"]
            i_rep = len(rec.rows)
            for nn in range(k0 + 1, k0 + span + 1):
                dt, k = st["call"][nn]
                rec.n, rec.on = nn, True
                rec.remaining, rec.dt_expect = float(dt), solver._dt
                acc0, step0 = rec.accepted, state._step
                f_step(self, state, dt, **k)
                rec.on = False
                if solver._dt != rec.dt_expect or state._step - step0 != rec.accepted - acc0:
                    rec.mismatch += 1
            w["reps"].append({"rows": rec.rows[i_rep:], "t_end": float(state.t),
                              "dm_end": float((state.m - m_end).abs().max())})

    def finish():
        n = st["n"]
        out["n_end"] = n
        out["rows"] = rec.rows if not ref else None
        out["fields"] = ("n", "dt", "norm", "dt_opt", "remaining") + BRANCH + (REDO if a.k1_redo else ()) + (
            REDO_ALL if a.k1_redo == "all" else ())
        if a.k1_redo == "all":
            out["redo_events"], out["redo_thr"] = rec.events, (THR_ERR, THR_X1)
        out["mismatch"] = rec.mismatch
        out["last_h"] = dict(lh)
        for w in out["win"].values():
            w["probe_m"] = torch.stack(w["probe_m"]) if w.get("probe_m") else None
        torch.save(out, a.k1_out if ref else a.k1_ref)
        print("[step_dump] K1 %s: %d sample intervals, %d windows, %d tries recorded, mismatch %d, stored stray field "
              "(B1) used %d, not used %d; saved to %s"
              % (out["mode"], n, len(out["win"]), sum(len(w.get("rows", [])) for w in out["win"].values())
                 if ref else len(rec.rows), rec.mismatch, lh["hit"], lh["miss"], a.k1_out if ref else a.k1_ref),
              flush=True)
        if a.k1_redo == "all":
            print("[step_dump] --k1_redo all: %d events (error estimate above %.0e tolerance units or x1 above %.0e)"
                  % (len(rec.events), THR_ERR, THR_X1), flush=True)
        raise SystemExit(3 if rec.mismatch else 0)

    magnumnp.LLGSolver.step = step
    run_loops.main(rl)
    if st["n"] < (ref["n_end"] if ref else (a.k1_stop or 0)):
        print("[step_dump] K1: the run ended after %d sample intervals" % st["n"], flush=True)
        sys.exit(2)
    if ref or a.k1_stop is None:
        finish()                                  # after the end of the cycle: the sample rows are written


def main(argv):
    if argv[:1] == ["--compare"]:
        sys.exit(0 if compare(argv[1], argv[2]) else 1)
    if "--" not in argv:
        raise SystemExit("usage: step_dump.py --at K ... --save FILE [--rl_dir DIR] -- <run_loops arguments>")
    i = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--at", type=int, nargs="+")
    ap.add_argument("--save")
    ap.add_argument("--go_on", action="store_true", help="do not stop after the last --at (the run ends by itself)")
    ap.add_argument("--rl_dir", default=".", help="folder of run_loops.py (default: the current folder)")
    ap.add_argument("--trace", action="store_true", help="tries, accepted steps, t and probe cells after each call")
    ap.add_argument("--k1_ref", help="K1 reference run: write the checkpoints and windows to this file")
    ap.add_argument("--k1_arm", help="K1 arm: the file of the reference run")
    ap.add_argument("--k1_out", help="K1 arm: output file")
    ap.add_argument("--k1_null", nargs="?", const="reverse", choices=["reverse", "rotate", "kernel"],
                    help="K1 null run: sum the field terms in the reverse (default) or rotated order, or the kernel "
                         "times (1 + 2^-52)")
    ap.add_argument("--k1_ck", type=int, nargs="+", default=[0, 8] + list(range(16, 256, 16)),
                    help="K1 checkpoints k0 (PLAN 14a: 0, 8, 16, 32, ..., 240)")
    ap.add_argument("--k1_span", type=int, default=5, help="K1 window: sample intervals k0 + 1 ... k0 + span")
    ap.add_argument("--k1_full", type=int, nargs="+", default=[1, 2, 4], help="K1: full m after k0 + these")
    ap.add_argument("--k1_stop", type=int, default=None, help="K1 reference: stop after this sample interval")
    ap.add_argument("--k1_reps", type=int, default=1, help="K1 arm: integrate each window this many times")
    ap.add_argument("--k1_redo", nargs="?", const="rejected", choices=["rejected", "all"],
                    help="K1: error RMS per try; do each try above 1.1 again (all: each try, with the field terms)")
    a = ap.parse_args(argv[:i])
    sys.path.insert(0, str(pathlib.Path(a.rl_dir).resolve()))
    if a.k1_ref or a.k1_arm:
        if bool(a.k1_arm) != bool(a.k1_out) or (a.k1_ref and a.k1_arm):
            raise SystemExit("step_dump: --k1_ref FILE, or --k1_arm REF --k1_out FILE")
        return k1_main(a, argv[i + 1:])
    if not a.at or not a.save:
        raise SystemExit("step_dump: --at and --save are necessary")
    import magnumnp
    import run_loops
    f_step = magnumnp.LLGSolver.step
    res = {"m": {}, "t": {}, "argv": argv[i + 1:]}
    n = [0]
    if a.trace:
        from magnumnp.solvers.ode_solvers.rkf45 import RKF45
        f_try, tries = RKF45._try_step, [0]

        def try_step(self, *args, **kw):
            tries[0] += 1
            return f_try(self, *args, **kw)

        RKF45._try_step = try_step
        res["trace"] = {"tries": [], "accepted": [], "t": [], "probe_m": []}

    def step(self, state, dt, **k):
        r = f_step(self, state, dt, **k)
        n[0] += 1
        if a.trace:
            if "probe_idx" not in res:
                res["probe_idx"] = probe_cells(state.m)
            tr = res["trace"]
            tr["tries"].append(tries[0])
            tr["accepted"].append(int(state._step))
            tr["t"].append(float(state.t))
            tr["probe_m"].append(state.m.detach().reshape(-1, 3)[res["probe_idx"].to(state.m.device)].cpu().clone())
        if n[0] in a.at:
            res["m"][n[0]] = state.m.detach().cpu().clone()
            res["t"][n[0]] = float(state.t)
            if n[0] == max(a.at):
                torch.save(res, a.save)
                print("[step_dump] saved m after %s sample intervals to %s" % (sorted(res["m"]), a.save), flush=True)
                if not a.go_on:
                    raise SystemExit(0)
        return r

    magnumnp.LLGSolver.step = step
    run_loops.main(argv[i + 1:])
    print("[step_dump] the run ended after %d sample intervals" % n[0], flush=True)
    sys.exit(0 if a.go_on and n[0] >= max(a.at) else 2)


if __name__ == "__main__":
    main(sys.argv[1:])
