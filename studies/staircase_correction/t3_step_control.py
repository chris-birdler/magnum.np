#!/usr/bin/env python
"""
T3 of PLAN 14a: table test of the RKF45 step control, old code against A4 and A3 (local CPU, double precision).

    python t3_step_control.py OLD_TREE A4_TREE A3_TREE [--json out.json]

Each tree gives its magnumnp/solvers/ode_solvers/rkf45.py (loaded as its own module).
(a) _optimal_stepsize (and _error_norm of A3) on error arrays with set norms: 0, 1e-12, at the maxscale switch
    (0.09^5), at 0.5 and 1.1 (each with neighbours of a few ulp), at the minscale switch ((0.9/0.2)^4), 1e6; dt from
    1e-15 to maxstep. Old and A4 must give bitwise equal dt_opt; A3 bitwise equal dt_opt and norm.
(b) step() on test ODEs (rotation and damping with a pulse; a tiny start dt; a slow ODE that reaches maxstep; a jump of
    the field by 3000 that gives norms above the minscale switch), over many
    sample intervals: old and A4 must give the same accepted (t, x, dt) sequence, bitwise, except after an A4 case (a
    try past the end with dt_opt below the rest of the interval; counted); A3 against A4: the same decisions and
    branches, dt within 1e-6 relative and x within 1e-12.
The branch counts of (b) show which branches the test covers.
"""
import argparse
import importlib.util
import json
import math
import pathlib
import sys

import torch

torch.set_default_dtype(torch.float64)


def load(tree, tag):
    p = pathlib.Path(tree) / "magnumnp" / "solvers" / "ode_solvers" / "rkf45.py"
    spec = importlib.util.spec_from_file_location("rkf45_" + tag, p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class St:
    def __init__(self):
        self._dt, self._step = 0.0, 0


def branch(s, dt, norm, rem, dt_opt):
    """Branch tuple of rkf45.py (as in step_dump.TryRecord.control)."""
    if norm > 1.1:
        cls, r = 2, s._headroom / norm**(1.0/s._order)
        cmin = r < s._minscale
        cmax = cone = False
        r = s._minscale if cmin else r
    elif norm < 0.5:
        cls = 0
        r = s._headroom / norm**(1.0/(s._order+1.0)) if norm > 0. else s._maxscale
        cmax, cmin = r > s._maxscale, False
        r = s._maxscale if cmax else r
        cone = r < 1.
    else:
        cls, r, cmin, cmax, cone = 1, 1., False, False, False
    cstep = s._dt * r > s._maxstep
    past = dt > rem
    acc = not (dt > dt_opt or dt > rem)
    retry = 0 if acc else (1 if dt_opt < rem else 2)
    return (acc, past, cls, cmin, cmax, cone, cstep, retry)


def instrument(mod):
    """Record each try: (dt, norm, dt_opt, remaining, branch). The norm by the expression of the tree."""
    R = mod.RKF45
    f_opt, f_step = R._optimal_stepsize, R.step
    rec = {"rows": [], "rem": None, "k": -1}

    def opt(self, x, err, atol, rtol):
        d = f_opt(self, x, err, atol, rtol)
        if hasattr(mod, "_error_norm"):
            n = float(mod._error_norm(err, x, atol, rtol))
        else:
            n = float(torch.linalg.norm((err / (atol + rtol*x.abs())).flatten(), torch.inf))
        rem = rec["rem"][0]
        b = branch(self, self._dt, n, rem, d)
        rec["rows"].append((self._dt, n, d, rem, b, rec["k"]))
        if b[0]:
            rec["rem"][0] = rem - self._dt
        return d

    def step(self, t, x, dt, **kw):
        rec["rem"] = [float(dt)]
        rec["k"] += 1
        return f_step(self, t, x, dt, **kw)

    R._optimal_stepsize, R.step = opt, step
    return rec


def table(mods, atol=1e-7, rtol=1e-6):
    """(a): dt_opt and norm on error arrays with set norms."""
    g = torch.Generator().manual_seed(3)
    x = 2 * torch.rand((4096, 3), generator=g) - 1
    u = 2 * torch.rand((4096, 3), generator=g) - 1
    u[17, 1] = 1.0                                       # the largest |u| is exactly 1
    targets = [0.0, 1e-12, 0.09**5, 0.4999999, 0.5, 0.5000001, 1.0999999, 1.1, 1.1000001, (0.9 / 0.2)**4, 1e6]
    eps = 2.0 ** -52
    for b in (0.5, 1.1, 0.09**5, (0.9 / 0.2)**4):
        targets += [b * (1 - 2 * eps), b * (1 + 2 * eps)]
    dts = [1e-15, 1e-12, 7e-12, 1.5e-11, 9.9e-11, 1e-10]
    rows, bad = [], 0
    for tgt in targets:
        err = tgt * (atol + rtol * x.abs()) * u
        for dt in dts:
            res = []
            for name, m in mods:
                s = m.RKF45(lambda t, x, **k: x, dt=dt, atol=atol, rtol=rtol)
                d = s._optimal_stepsize(x, err, atol, rtol)
                n = (float(m._error_norm(err, x, atol, rtol)) if hasattr(m, "_error_norm")
                     else float(torch.linalg.norm((err / (atol + rtol*x.abs())).flatten(), torch.inf)))
                res.append((name, d, n))
            same_d = len({r[1] for r in res}) == 1
            same_n = len({r[2] for r in res}) == 1
            bad += not (same_d and same_n)
            rows.append({"target": tgt, "dt": dt, "dt_opt": [r[1] for r in res], "norm": [r[2] for r in res],
                         "equal": same_d and same_n})
    return rows, bad


def ode(kind):
    """Test ODEs dx/dt = f(t, x) with the scales of production (rates about 1e10 1/s)."""
    g = torch.Generator().manual_seed(5)
    n = 512
    w = 1e10 * (0.5 + torch.rand((n, 1), generator=g))
    x0 = torch.nn.functional.normalize(torch.randn((n, 3), generator=g), dim=-1)
    ez = torch.tensor([0.0, 0.0, 1.0])

    def f(t, x, **kw):
        tt = float(t)
        s = 1.0
        if kind == "pulse":
            s = 1.0 + 30.0 * math.exp(-((tt - 3e-10) / 2e-11) ** 2)
        elif kind == "slow":
            s = 1e-3
        elif kind == "kick":
            s = 1.0 if tt < 3e-10 else 3000.0
        h = ez.expand_as(x) * s
        c = torch.linalg.cross(x, h)
        return -w * c - 0.1 * w * torch.linalg.cross(x, c)

    return f, x0


def run(m, kind, n_int, dt_int, dt0):
    f, x0 = ode(kind)
    s = m.RKF45(f, dt=dt0, atol=1e-7, rtol=1e-6)
    st = St()
    t, x = torch.tensor(0.0), x0.clone()
    acc = []
    for _ in range(n_int):
        t, x = s.step(t, x, dt_int, state=st)
        acc.append((float(t), x.clone(), s._dt))
    return acc


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("trees", nargs=3, metavar=("OLD", "A4", "A3"))
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    mods = [(nm, load(tr, nm)) for nm, tr in zip(("old", "a4", "a3"), a.trees)]
    out = {}
    rows, bad = table(mods)
    out["table"] = {"cases": len(rows), "not_equal": bad, "rows": rows}
    print("T3 (a) table: %d cases (norm targets x dt), %d with a difference of dt_opt or norm between old, A4 and A3"
          % (len(rows), bad))
    recs = {nm: instrument(m) for nm, m in mods}
    ok = bad == 0
    scen = [("smooth", 64, 2.4e-11, 1e-15), ("pulse", 64, 1.2e-11, 7e-12), ("slow", 32, 3e-10, 1e-12),
            ("kick", 24, 2.5e-11, 7e-12)]
    out["step"] = {}
    for kind, n_int, dt_int, dt0 in scen:
        res, rws = {}, {}
        for nm, m in mods:
            recs[nm]["rows"].clear()
            recs[nm]["k"] = -1
            res[nm] = run(m, kind, n_int, dt_int, dt0)
            rws[nm] = list(recs[nm]["rows"])
        old_rows = rws["old"]
        cases = [r[5] for r in old_rows if r[4][1] and r[4][7] == 1]
        a4_case = len(cases)
        past = sum(1 for r in old_rows if r[4][1])
        # old against A4: the state after each interval, bitwise, up to the first interval with an A4 case
        k_end = min(cases) if cases else len(res["old"])
        eq_a4 = all(ta == tb and torch.equal(xa, xb) and da == db
                    for (ta, xa, da), (tb, xb, db) in zip(res["old"][:k_end], res["a4"][:k_end]))
        # A3 against A4: decisions and branches (without the past-end flag), dt and x close
        ka = [r[4] for r in rws["a4"]]
        k3 = [r[4] for r in rws["a3"]]
        same_branch = ka == k3
        dmax = max((abs(r3[0] / ra[0] - 1) for ra, r3 in zip(rws["a4"], rws["a3"])), default=0.0)
        xmax = max(float((xa - x3).abs().max()) for (_, xa, _), (_, x3, _) in zip(res["a4"], res["a3"]))
        cnt = {}
        for r in old_rows + rws["a4"]:
            b = r[4]
            for key, on in (("reject", not b[0]), ("past_end", b[1]), ("norm<0.5", b[2] == 0), ("norm 0.5..1.1", b[2] == 1),
                            ("norm>1.1", b[2] == 2), ("c_minscale", b[3]), ("c_maxscale", b[4]), ("c_maxstep", b[6]),
                            ("retry dt_opt", b[7] == 1), ("retry rest", b[7] == 2)):
                cnt[key] = cnt.get(key, 0) + int(on)
        pass_a4 = eq_a4 and k_end > 0
        pass_a3 = same_branch and dmax <= 1e-6 and xmax <= 1e-12
        ok &= pass_a4 and pass_a3
        out["step"][kind] = {"tries_old": len(old_rows), "tries_a4": len(rws["a4"]), "tries_a3": len(rws["a3"]),
                             "past_end_old": past, "a4_case": a4_case, "a4_case_intervals": cases,
                             "old_a4_bitwise_intervals": k_end, "old_a4_bitwise": eq_a4,
                             "a3_same_branches": same_branch, "a3_dt_rel_max": dmax, "a3_x_max": xmax,
                             "branch_counts": cnt, "pass": pass_a4 and pass_a3}
        print("T3 (b) %-6s tries old %d, A4 %d, A3 %d; past the end (old) %d, A4 cases %d; old = A4 bitwise %s; A3: same "
              "branches %s, dt %.1e, x %.1e; branches %s" % (kind, len(old_rows), len(rws["a4"]), len(rws["a3"]), past,
                                                             a4_case, eq_a4, same_branch, dmax, xmax, cnt))
    out["pass"] = ok
    print("RESULT T3: %s" % ("PASS" if ok else "FAIL"))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1, default=str))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
