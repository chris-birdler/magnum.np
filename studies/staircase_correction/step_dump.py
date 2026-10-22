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
"""
import argparse
import pathlib
import sys

import torch


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


def main(argv):
    if argv[:1] == ["--compare"]:
        sys.exit(0 if compare(argv[1], argv[2]) else 1)
    if "--" not in argv:
        raise SystemExit("usage: step_dump.py --at K ... --save FILE [--rl_dir DIR] -- <run_loops arguments>")
    i = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--at", type=int, nargs="+", required=True)
    ap.add_argument("--save", required=True)
    ap.add_argument("--go_on", action="store_true", help="do not stop after the last --at (the run ends by itself)")
    ap.add_argument("--rl_dir", default=".", help="folder of run_loops.py (default: the current folder)")
    ap.add_argument("--trace", action="store_true", help="tries, accepted steps, t and probe cells after each call")
    a = ap.parse_args(argv[:i])
    sys.path.insert(0, str(pathlib.Path(a.rl_dir).resolve()))
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
                mag = (state.m.detach().norm(dim=-1) > 0).reshape(-1).cpu().nonzero().squeeze(1)
                g = torch.Generator().manual_seed(0)
                res["probe_idx"] = mag[torch.randperm(len(mag), generator=g)[:2048]].sort().values
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
