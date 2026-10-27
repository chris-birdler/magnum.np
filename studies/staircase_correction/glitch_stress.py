#!/usr/bin/env python
"""
Finding of the K1 pilot and of the replay (PLAN 14a, box 11): single RKF45 tries of the B code gave an error norm 4 ...
2700 times the norm of the same try in other runs; a second try with the same input gave the normal norm (transient
fault). Stress test in the production path: run run_loops.main with the given arguments; at the first call of
LLGSolver.step take the checkpoint states --k0 of the K1 reference, then evaluate each part again and again for the
given time, in turn on each state (a new tensor of m each time, as in the stages). With two or more states, a part
that reads data of the evaluation before (a stale buffer) gives a large difference; with one state it does not.
Each result is compared with the first result on the same state. An event is a result with a difference above the
threshold. For each event: the cell of the largest difference, its volume fraction (material Ms / largest Ms) and
whether the sparse matrix of the correction has a row for the cell (row_nnz: entries of the 3 rows of the
cell; row_start and row_len: position and length of the row of the index in the CSR arrays). Then stop.

Parts: try (one RKF45 try, RKF45._try_step of the tree: six stages, increment, error estimate; the try dt is the
median dt of the accepted tries in the interval k0 + 1 of the reference that are not clamped to the rest of the
interval), corr (DemagStaircaseCorrection.h), spmv (the sparse product of the correction, only with the attribute _A),
centroid (centroid_field of the correction), demag (DemagField.h), dm (LLGSolver.dm, all terms).
A part can be given more than once (for example spmv:3600 spmv:3600): each run has its own result.
Differences: field parts largest |h_r - h_0| / largest |h_0| (threshold --thr); try: the error estimate in units of
the tolerance, largest |err_r - err_0| / (atol + rtol |x|) (threshold --thr_err), and the increment, largest
|x1_r - x1_0| (threshold --thr_x1).

    CUDA_DEVICE=0 python glitch_stress.py --ref refB.pt --k0 8 144 --parts try:3600 corr:1200 --idx 1513981 \
        --json out.json -- <run_loops arguments>
"""
import argparse
import json
import pathlib
import statistics
import sys
import time

import torch


def main(argv):
    i = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--k0", type=int, nargs="+", default=[8])
    ap.add_argument("--parts", nargs="+", default=["corr:1800", "dm:900"])
    ap.add_argument("--thr", type=float, default=1e-9)
    ap.add_argument("--thr_err", type=float, default=1e-6)
    ap.add_argument("--thr_x1", type=float, default=1e-12)
    ap.add_argument("--idx", type=int, nargs="*", default=[])
    ap.add_argument("--json")
    a = ap.parse_args(argv[:i])
    sys.path.insert(0, str(pathlib.Path(".").resolve()))
    import magnumnp
    import run_loops
    ref = torch.load(a.ref, map_location="cpu", weights_only=False)
    out = {"argv": argv[i + 1:], "k0": a.k0, "thr": a.thr, "thr_err": a.thr_err, "thr_x1": a.thr_x1, "parts": {},
           "cells": {}, "try_dt": {}}

    def try_dt(k0):
        f = ref["fields"]
        i_n, i_dt, i_acc, i_rem = f.index("n"), f.index("dt"), f.index("accept"), f.index("remaining")
        v = [r[i_dt] for r in ref["win"][k0]["rows"] if r[i_n] == k0 + 1 and r[i_acc] and r[i_dt] < 0.999 * r[i_rem]]
        return statistics.median(v) if v else ref["ck"][k0]["dt"]

    def step(self, state, dt, **k):
        dev, dt_ = state.m.device, state.m.dtype
        M = [ref["ck"][k0]["m"].to(device=dev, dtype=dt_) for k0 in a.k0]
        T = [torch.tensor(ref["ck"][k0]["t"], dtype=torch.float64, device=dev) for k0 in a.k0]
        DT = [try_dt(k0) for k0 in a.k0]
        out["try_dt"] = dict(zip(a.k0, DT))
        S = len(M)
        state.t = T[0]
        names = [type(t).__name__ for t in self._terms]
        corr = self._terms[names.index("DemagStaircaseCorrection")]
        dem = self._terms[names.index("DemagField")]
        solver = self._solver
        state.m = M[0].clone()
        corr.h(state)  # set-up
        mod = sys.modules[type(corr).__module__]
        Ms = state.material["Ms"].reshape(-1).expand(state.m[..., 0].numel())
        frac = (Ms / Ms.max()).cpu()
        rows = crow = None
        if hasattr(corr, "_A"):
            crow = corr._A.crow_indices().cpu()
            rows = (crow[1:] - crow[:-1]).reshape(-1, 3).sum(-1)
        n_cut = int(((frac > 0) & (frac < 1)).sum())
        print("[stress] %s: %d cells, %d cut cells (0 < fraction < 1), %s cells with a row of the sparse matrix; "
              "states k0 %s, try dt %s" % (state.m.dtype, frac.numel(), n_cut,
                                           "-" if rows is None else int((rows > 0).sum()), a.k0,
                                           " ".join("%.4g" % x for x in DT)), flush=True)

        def where(j):
            cell = j // 3
            return {"index": j, "cell": cell, "comp": j % 3, "fraction": float(frac[cell]),
                    "row_nnz": None if rows is None else int(rows[cell]),
                "row_start": None if crow is None else int(crow[j]),
                "row_len": None if crow is None else int(crow[j + 1] - crow[j])}

        for j in a.idx:
            out["cells"][j] = where(j)
            print("[stress] index %d: %s" % (j, out["cells"][j]), flush=True)

        def f_corr(s):
            state.t, state.m = T[s], M[s].clone()
            return corr.h(state)

        def f_spmv(s):
            return corr._A @ M[s].clone().reshape(-1)

        def f_centroid(s):
            if not hasattr(dem, "_N"):
                dem._N = dem._init_N(state)
            return mod.centroid_field(dem._N, corr._v, corr._d, M[s].clone(), shape=tuple(dem._shape(state)))

        def f_demag(s):
            state.t, state.m = T[s], M[s].clone()
            return dem.h(state)

        def f_dm(s):
            state.t, state.m = T[s], M[s].clone()
            return self.dm(T[s], state.m, state=state)

        def f_try(s):
            solver._dt = DT[s]
            x1, _, err = solver._try_step(T[s], M[s].clone(), state=state)
            return x1, err

        fs = {"corr": f_corr, "spmv": f_spmv, "centroid": f_centroid, "demag": f_demag, "dm": f_dm}
        cap = 1 << 22
        for p in a.parts:
            nm, sec = p.split(":")
            if nm == "spmv" and not hasattr(corr, "_A"):
                continue
            dev_ = torch.zeros(cap, dtype=torch.float64, device=dev)
            arg = torch.zeros(cap, dtype=torch.int64, device=dev)
            dx1 = torch.zeros(cap, dtype=torch.float64, device=dev)
            if nm == "try":
                first = [tuple(y.clone() for y in f_try(s)) for s in range(S)]
                w = [1.0 / (solver._atol + solver._rtol * M[s].abs()) for s in range(S)]
            else:
                first = [fs[nm](s).clone() for s in range(S)]
                scale = [float(h.abs().max()) for h in first]
            r, t0 = 0, time.time()
            while r < cap:
                s = r % S
                if nm == "try":
                    x1, err = f_try(s)
                    d = ((err - first[s][1]).abs() * w[s]).reshape(-1)
                    dx1[r] = (x1 - first[s][0]).abs().max()
                    dev_[r], arg[r] = d.max(), d.argmax()
                else:
                    d = (fs[nm](s) - first[s]).abs().reshape(-1)
                    dev_[r], arg[r] = d.max() / scale[s], d.argmax()
                r += 1
                if r % (20 if nm == "try" else 200) == 0 and time.time() - t0 > float(sec):
                    break
            if dev.type == "cuda":
                torch.cuda.synchronize()
            wall = time.time() - t0
            dv, ag, d1 = dev_[:r].cpu(), arg[:r].cpu(), dx1[:r].cpu()
            thr = a.thr_err if nm == "try" else a.thr
            hit = (dv > thr) | (d1 > a.thr_x1) if nm == "try" else dv > thr
            ev = torch.nonzero(hit).reshape(-1)
            res = {"evaluations": r, "states": S, "wall_s": wall, "n_nonzero": int((dv > 0).sum()),
                   "n_events": int(ev.numel()), "max_dev": float(dv.max()) if r else 0.0,
                   "median_dev": float(dv.median()) if r else 0.0,
                   "events": [dict(rep=int(e), state=a.k0[int(e) % S], dev=float(dv[e]), **where(int(ag[e])))
                              for e in ev[:200]]}
            if nm == "try":
                res["max_dx1"] = float(d1.max()) if r else 0.0
                res["median_dx1"] = float(d1.median()) if r else 0.0
                for e, x in zip(ev[:200], res["events"]):
                    x["dx1"] = float(d1[e])
            out["parts"][nm if nm not in out["parts"] else "%s_%d" % (nm, len(out["parts"]))] = res
            print("[stress] %-8s %s: %d evaluations on %d states in %.0f s, nonzero %d, events %d (threshold %.0e%s), "
                  "largest difference %.3e, median %.3e%s" % (
                      nm, state.m.dtype, r, S, wall, res["n_nonzero"], res["n_events"], thr,
                      " or x1 %.0e" % a.thr_x1 if nm == "try" else "", res["max_dev"], res["median_dev"],
                      "; x1 largest %.3e, median %.3e" % (res["max_dx1"], res["median_dx1"]) if nm == "try" else ""),
                  flush=True)
            if res["n_events"] > r // 2:
                print("[stress] %s: more than half of the results are events: a first result is an outlier" % nm)
            for e in res["events"][:20]:
                print("   event %s" % e, flush=True)
            if a.json:
                pathlib.Path(a.json).write_text(json.dumps(out, indent=1))
        raise SystemExit(0)

    magnumnp.LLGSolver.step = step
    run_loops.main(argv[i + 1:])


if __name__ == "__main__":
    main(sys.argv[1:])
