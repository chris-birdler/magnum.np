"""
Read-only instrumentation of a run_loops.py run (diagnosis of the gate 1 stall of PLAN 4.21, PROTOKOLL 7.146 part 6).
It patches methods in memory only and does not change the model; it sets the thread count of torch. Every 30 s of
wall time it writes one line with the RKF45 tries and rejections, the state time, the step size (the last step, the
smallest and the largest in the interval), the smallest |s_j| of the V2 tie (and its tied cell) and the largest |h|
of the masters and of all cells.
    python diag_stall421.py <diag log> <threads> <run_loops.py arguments ...>
"""
import math
import pathlib
import runpy
import sys
import time

import torch

B = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(B))
sys.path.insert(0, str(B.parent.parent))
log = open(sys.argv[1], "a", buffering=1)
torch.set_num_threads(int(sys.argv[2]))
sys.argv = [str(B / "run_loops.py")] + sys.argv[3:]

import magnumnp.solvers.ode_solvers.rkf45 as R   # noqa: E402
import field_terms_extra as FT                    # noqa: E402

S = {"n": 0, "rej": 0, "dt": [], "t": 0.0, "ns": math.inf, "ns_idx": -1, "n_small": 0, "hm": 0.0, "ha": 0.0,
     "wall": time.time(), "t0": time.time()}


def flush():
    dt = S["dt"]
    log.write("%s wall %7.0f s  tries %7d  rejected %6d  t %.6e s  dt last %.3e min %.3e max %.3e  "
              "min|s_j| %.3e (tied cell %d, n below 0.1: %d)  max|h| masters %.3e all %.3e A/m\n" % (
                  time.strftime("%H:%M:%S"), time.time() - S["t0"], S["n"], S["rej"], S["t"],
                  dt[-1] if dt else math.nan, min(dt) if dt else math.nan, max(dt) if dt else math.nan,
                  S["ns"], S["ns_idx"], S["n_small"], S["hm"], S["ha"]))
    S.update(n=0, rej=0, dt=[], ns=math.inf, ns_idx=-1, n_small=0, hm=0.0, ha=0.0, wall=time.time())


try_step0 = R.RKF45._try_step
opt0 = R.RKF45._optimal_stepsize


def try_step(self, t, x, **kw):
    S["n"] += 1
    S["dt"].append(float(self._dt))
    S["t"] = float(t)
    if time.time() - S["wall"] > 30.0:
        flush()
    return try_step0(self, t, x, **kw)


def optimal(self, x, err, atol, rtol):
    d = opt0(self, x, err, atol, rtol)
    if self._dt > d:
        S["rej"] += 1
    return d


tie0 = FT.TiedLLGSolver.tie
field0 = FT.TiedLLGSolver.field


def tie(self, x):
    out = tie0(self, x)
    ns = out[2].reshape(-1)
    v, i = torch.min(ns, 0)
    if float(v) < S["ns"]:
        S["ns"], S["ns_idx"] = float(v), int(self._tied[int(i)])
    S["n_small"] = max(S["n_small"], int((ns < 0.1).sum()))
    return out


def field(self, state, mj, ns):
    h = field0(self, state, mj, ns)
    hf = h.reshape(-1, 3)
    S["hm"] = max(S["hm"], float(torch.linalg.norm(hf[self._masters], dim=-1).max()))
    S["ha"] = max(S["ha"], float(torch.linalg.norm(hf, dim=-1).max()))
    return h


R.RKF45._try_step = try_step
R.RKF45._optimal_stepsize = optimal
FT.TiedLLGSolver.tie = tie
FT.TiedLLGSolver.field = field
log.write("start %s, threads %d, argv %s\n" % (time.strftime("%F %T"), torch.get_num_threads(), " ".join(sys.argv)))
runpy.run_path(sys.argv[0], run_name="__main__")
