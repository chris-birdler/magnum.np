#!/usr/bin/env python
"""
Finding of the K1 pilot (PLAN 14a, box 11): single RKF45 tries of the B code (with the staircase correction) gave an
error norm 10 ... 2700 times that of the other arms at the same state and dt. Test for a non-deterministic fault: run
run_loops.main with the given arguments (the production path); at the first call of LLGSolver.step set m to a
checkpoint state of the K1 reference, then evaluate the field of each term --reps times on the same m (a new tensor
each time, as in the stages) and compare each result with the first one: largest |h_i - h_0| / largest |h_0|.
Also the effective field through LLGSolver.dm. Then stop.

    CUDA_DEVICE=0 python glitch_test.py --ref refB.pt --k0 144 --reps 500 --json out.json -- <run_loops arguments>
"""
import argparse
import json
import pathlib
import sys

import torch


def main(argv):
    i = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--k0", type=int, nargs="+", default=[144])
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--json")
    a = ap.parse_args(argv[:i])
    sys.path.insert(0, str(pathlib.Path(".").resolve()))
    import magnumnp
    import run_loops
    ref = torch.load(a.ref, map_location="cpu")
    f_step = magnumnp.LLGSolver.step
    out = {"argv": argv[i + 1:], "reps": a.reps, "states": {}}

    def step(self, state, dt, **k):
        for k0 in a.k0:
            c = ref["ck"][k0]
            m0 = c["m"].to(device=state.m.device, dtype=state.m.dtype)
            state.t = torch.tensor(c["t"], dtype=torch.float64, device=state.m.device)
            res = {}
            names = ["h%d_%s" % (j, type(t).__name__) for j, t in enumerate(self._terms)] + ["dm"]
            first, dev = {}, {n: [] for n in names}
            for r in range(a.reps):
                for j, term in enumerate(self._terms):
                    state.m = m0.clone()
                    h = term.h(state)
                    if r == 0:
                        first[names[j]] = h.clone()
                    else:
                        d = float((h - first[names[j]]).abs().max()) / max(float(first[names[j]].abs().max()), 1e-300)
                        dev[names[j]].append(d)
                state.m = m0.clone()
                f = self.dm(state.t, state.m, state=state)
                if r == 0:
                    first["dm"] = f.clone()
                else:
                    dev["dm"].append(float((f - first["dm"]).abs().max()) / max(float(first["dm"].abs().max()), 1e-300))
            for n in names:
                v = dev[n]
                res[n] = {"max": max(v) if v else 0.0, "n_nonzero": sum(1 for x in v if x > 0),
                          "n_above_1e-6": sum(1 for x in v if x > 1e-6), "largest5": sorted(v)[-5:]}
                print("[glitch] k0 %d %-32s %s: max rel dev %.3e, nonzero %d, above 1e-6 %d of %d" % (
                    k0, n, state.m.dtype, res[n]["max"], res[n]["n_nonzero"], res[n]["n_above_1e-6"], len(v)),
                    flush=True)
            out["states"][k0] = res
        if a.json:
            pathlib.Path(a.json).write_text(json.dumps(out, indent=1))
        raise SystemExit(0)

    magnumnp.LLGSolver.step = step
    run_loops.main(argv[i + 1:])


if __name__ == "__main__":
    main(sys.argv[1:])
