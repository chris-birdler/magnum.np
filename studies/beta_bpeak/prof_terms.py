#!/usr/bin/env python
"""
Time share of each field term in a production run (Chris 2026-10-19, option 1: is the time extrapolation of the stray
field of magnum.np upstream branch feature/demag-extrapolation worth a test?). The magnum.np Timer does not wait for the
GPU, thus on a GPU it measures only the kernel launches. This wrapper synchronizes the device before and after each
field term h(), each LLGSolver.dm() (one RKF45 stage) and each LLGSolver.step(), and then runs run_loops.main with the
given arguments (same code path as production). The report (stdout and <out>/prof_terms.json) gives calls, total and mean time per method;
the h() times inside dm() (the RHS of the integrator) and outside it (samples, diagnostics, energy) are kept apart.
Counting and the wall clock start after the first --skip integrator steps (warm-up, torch.compile).

    python prof_terms.py --skip 200 -- <run_loops arguments>

Saving of the extrapolation (estimate): the stray field is evaluated in 1 of the 6 stages of a step, thus
    saved = 5/6 * T(stray field in dm) / T(wall of the run after warm-up).
"""
import json
import pathlib
import sys
import time

import torch

ACC = {}
STATE = {"in_dm": 0, "skip": 200, "t0": None, "steps": 0}


def _sync():
    if torch.cuda.is_available():
        torch.cuda.synchronize()


def wrap(cls, meth, name):
    """Time each call of cls.meth after the warm-up: counting starts with step number skip + 1 of the integrator, and
    the wall clock starts there too."""
    f = getattr(cls, meth)

    def timed(self, *a, **k):
        if meth == "step":
            STATE["steps"] += 1
            if STATE["t0"] is None and STATE["steps"] > STATE["skip"]:
                _sync()
                STATE["t0"] = time.perf_counter()
        tag = name if meth != "h" else "%s.h (%s)" % (name, "in dm" if STATE["in_dm"] else "outside dm")
        e = ACC.setdefault(tag, {"calls": 0, "counted": 0, "s": 0.0})
        e["calls"] += 1
        count = STATE["t0"] is not None
        if count:
            _sync()
            t = time.perf_counter()
        if meth == "dm":
            STATE["in_dm"] += 1
        try:
            return f(self, *a, **k)
        finally:
            if meth == "dm":
                STATE["in_dm"] -= 1
            if count:
                _sync()
                e["s"] += time.perf_counter() - t
                e["counted"] += 1

    setattr(cls, meth, timed)


def report(out):
    _sync()
    wall = time.perf_counter() - STATE["t0"] if STATE["t0"] is not None else float("nan")
    rows = sorted(ACC.items(), key=lambda kv: -kv[1]["s"])
    print("[prof] wall after the warm-up: %.1f s" % wall, flush=True)
    print("[prof] %-46s %9s %9s %10s %10s %7s" % ("method", "calls", "counted", "total s", "mean ms", "% wall"))
    for k, e in rows:
        print("[prof] %-46s %9d %9d %10.2f %10.3f %7.1f" % (k, e["calls"], e["counted"], e["s"],
                                                             1e3 * e["s"] / max(e["counted"], 1), 100 * e["s"] / wall))
    demag = sum(e["s"] for k, e in ACC.items() if k.startswith("DemagField.h (in dm)"))
    res = {"wall_s": wall, "methods": ACC, "skip": STATE["skip"], "device": torch.cuda.get_device_name(0)
           if torch.cuda.is_available() else "cpu",
           "demag_in_dm_share": demag / wall, "saved_share_estimate": 5.0 / 6.0 * demag / wall}
    print("[prof] stray field in dm: %.1f %% of the wall; saving of the extrapolation (estimate 5/6 of it): %.1f %%, "
          "speed-up %.2f" % (100 * res["demag_in_dm_share"], 100 * res["saved_share_estimate"],
                             1.0 / (1.0 - res["saved_share_estimate"])), flush=True)
    if out is not None:
        pathlib.Path(out).mkdir(parents=True, exist_ok=True)
        (pathlib.Path(out) / "prof_terms.json").write_text(json.dumps(res, indent=1))


def main(argv):
    if "--" not in argv:
        raise SystemExit("usage: prof_terms.py [--skip N] -- <run_loops arguments>")
    i = argv.index("--")
    own, rl = argv[:i], argv[i + 1:]
    if "--skip" in own:
        STATE["skip"] = int(own[own.index("--skip") + 1])
    import magnumnp
    import field_terms_extra
    import run_loops
    for cls in (magnumnp.DemagField, magnumnp.DemagFieldPBC, magnumnp.ExchangeField, magnumnp.UniaxialAnisotropyField,
                field_terms_extra.SinusoidalDrive):
        wrap(cls, "h", cls.__name__)
    wrap(magnumnp.LLGSolver, "dm", "LLGSolver.dm (one stage)")
    wrap(magnumnp.LLGSolver, "step", "LLGSolver.step")
    out = rl[rl.index("--out") + 1] if "--out" in rl else None
    try:
        run_loops.main(rl)
    finally:
        report(out)


if __name__ == "__main__":
    main(sys.argv[1:])
