#!/usr/bin/env python
"""
Finding of the K1 pilot and of the replay (PLAN 14a, box 11): single RKF45 tries of the B code gave an error norm 4 ...
2700 times the norm of the same try in other runs; a second try with the same input gave the normal norm (transient
fault). Stress test of the parts of the field on one fixed state, in the production path: run run_loops.main with the
given arguments; at the first call of LLGSolver.step set m to a checkpoint state of the K1 reference, then evaluate
each part again and again for the given time (a new tensor of m each time, as in the stages) and compare each result
with the first one: d = largest |h_r - h_0| / largest |h_0|. An event is a result with d above --thr. For each event:
the cell of the largest difference, its volume fraction (material Ms / largest Ms) and whether the sparse matrix of
the correction has a row for the cell. Then stop.

Parts: corr (DemagStaircaseCorrection.h), spmv (the sparse product of the correction, only with the attribute _A),
centroid (centroid_field of the correction), demag (DemagField.h), dm (LLGSolver.dm, all terms).

    CUDA_DEVICE=0 python glitch_stress.py --ref refB.pt --k0 8 --parts corr:1800 dm:900 --idx 1513981 876667 \
        --json out.json -- <run_loops arguments>
"""
import argparse
import json
import pathlib
import sys
import time

import torch


def main(argv):
    i = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--k0", type=int, default=8)
    ap.add_argument("--parts", nargs="+", default=["corr:1800", "dm:900"])
    ap.add_argument("--thr", type=float, default=1e-9)
    ap.add_argument("--idx", type=int, nargs="*", default=[])
    ap.add_argument("--json")
    a = ap.parse_args(argv[:i])
    sys.path.insert(0, str(pathlib.Path(".").resolve()))
    import magnumnp
    import run_loops
    ref = torch.load(a.ref, map_location="cpu")
    out = {"argv": argv[i + 1:], "k0": a.k0, "thr": a.thr, "parts": {}, "cells": {}}

    def step(self, state, dt, **k):
        c = ref["ck"][a.k0]
        m0 = c["m"].to(device=state.m.device, dtype=state.m.dtype)
        state.t = torch.tensor(c["t"], dtype=torch.float64, device=state.m.device)
        names = [type(t).__name__ for t in self._terms]
        corr = self._terms[names.index("DemagStaircaseCorrection")]
        dem = self._terms[names.index("DemagField")]
        state.m = m0.clone()
        corr.h(state)  # set-up
        mod = sys.modules[type(corr).__module__]
        Ms = state.material["Ms"].reshape(-1).expand(state.m[..., 0].numel())
        frac = (Ms / Ms.max()).cpu()
        rows = None
        if hasattr(corr, "_A"):
            rows = (corr._A.crow_indices()[1:] - corr._A.crow_indices()[:-1]).reshape(-1, 3).sum(-1).cpu()
        n_cut = int(((frac > 0) & (frac < 1)).sum())
        print("[stress] %s: %d cells, %d cut cells (0 < fraction < 1), %s cells with a row of the sparse matrix" % (
            state.m.dtype, frac.numel(), n_cut, "-" if rows is None else int((rows > 0).sum())), flush=True)

        def where(j):
            cell = j // 3
            return {"index": j, "cell": cell, "comp": j % 3, "fraction": float(frac[cell]),
                    "row_nnz": None if rows is None else int(rows[cell])}

        for j in a.idx:
            out["cells"][j] = where(j)
            print("[stress] index %d: %s" % (j, out["cells"][j]), flush=True)

        def f_corr():
            state.m = m0.clone()
            return corr.h(state)

        def f_spmv():
            return corr._A @ m0.clone().reshape(-1)

        def f_centroid():
            if not hasattr(dem, "_N"):
                dem._N = dem._init_N(state)
            return mod.centroid_field(dem._N, corr._v, corr._d, m0.clone(), shape=tuple(dem._shape(state)))

        def f_demag():
            state.m = m0.clone()
            return dem.h(state)

        def f_dm():
            state.m = m0.clone()
            return self.dm(state.t, state.m, state=state)

        fs = {"corr": f_corr, "spmv": f_spmv, "centroid": f_centroid, "demag": f_demag, "dm": f_dm}
        for p in a.parts:
            nm, sec = p.split(":")
            if nm == "spmv" and not hasattr(corr, "_A"):
                continue
            f = fs[nm]
            h0 = f().clone()
            s0 = float(h0.abs().max())
            cap = 1 << 22
            dev = torch.zeros(cap, dtype=torch.float64, device=h0.device)
            arg = torch.zeros(cap, dtype=torch.int64, device=h0.device)
            r, t0 = 0, time.time()
            while r < cap:
                d = (f() - h0).abs().reshape(-1)
                dev[r], arg[r] = d.max(), d.argmax()
                r += 1
                if r % 200 == 0 and time.time() - t0 > float(sec):
                    break
            if h0.is_cuda:
                torch.cuda.synchronize()
            wall = time.time() - t0
            dv = (dev[:r] / s0).cpu()
            ag = arg[:r].cpu()
            ev = torch.nonzero(dv > a.thr).reshape(-1)
            res = {"evaluations": r, "wall_s": wall, "scale": s0, "n_nonzero": int((dv > 0).sum()),
                   "n_events": int(ev.numel()), "max_dev": float(dv.max()) if r else 0.0,
                   "median_dev": float(dv.median()) if r else 0.0,
                   "events": [dict(rep=int(e), dev=float(dv[e]), **where(int(ag[e]))) for e in ev[:200]]}
            out["parts"][nm] = res
            print("[stress] %-8s %s: %d evaluations in %.0f s, nonzero %d, events (d > %.0e) %d, largest d %.3e, "
                  "median d %.3e" % (nm, state.m.dtype, r, wall, res["n_nonzero"], a.thr, res["n_events"],
                                     res["max_dev"], res["median_dev"]), flush=True)
            if res["n_events"] > r // 2:
                print("[stress] %s: more than half of the results are events: the first result is an outlier" % nm)
            for e in res["events"][:20]:
                print("   event %s" % e, flush=True)
            if a.json:
                pathlib.Path(a.json).write_text(json.dumps(out, indent=1))
        raise SystemExit(0)

    magnumnp.LLGSolver.step = step
    run_loops.main(argv[i + 1:])


if __name__ == "__main__":
    main(sys.argv[1:])
