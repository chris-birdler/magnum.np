#!/usr/bin/env python
"""
PLAN 15 WQ (diagnosis of the host-memory stop at q 16, N 156): the resident host memory (VmRSS, VmHWM of
/proc/self/status) and the GPU memory at the entry and the exit of each phase of the geometry and of the set-up of
DemagStaircaseCorrection as in production (pair form, PBC, centroid term, float32), and every 5 s. Phases: fcc_subcells,
centroid_offsets, NearCorrection.__init__, the loop of _targets/_coeffs in _setup (first call and the end), pair_blocks.
If the OOM killer stops the process, the last lines name the phase.

    PYTHONPATH=TREE python mem_phases.py --config RUN/config.json
"""
import argparse
import functools
import json
import pathlib
import sys
import threading
import time
import types

import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
import energy_budget as EB  # noqa: E402
import magnumnp  # noqa: E402
import magnumnp.field_terms.demag_staircase as DS  # noqa: E402
sys.path.insert(0, str(pathlib.Path(magnumnp.__file__).resolve().parents[1] / "studies" / "beta_bpeak"))
sys.modules.pop("geometry", None)
import geometry  # noqa: E402

T0 = time.perf_counter()
PHASE = ["start"]


def mem():
    s = pathlib.Path("/proc/self/status").read_text()
    kb = {line.split(":")[0]: int(line.split()[1]) for line in s.splitlines() if line.startswith(("VmRSS", "VmHWM"))}
    g = torch.cuda.memory_allocated() / 2**30 if torch.cuda.is_available() else float("nan")
    return "RSS %.1f GiB, HWM %.1f GiB, GPU %.2f GiB" % (kb["VmRSS"] / 2**20, kb["VmHWM"] / 2**20, g)


def log(msg):
    print("[mem %7.1f s] %-55s %s" % (time.perf_counter() - T0, msg, mem()), flush=True)


def wrap(owner, name, once=False):
    f = getattr(owner, name)
    seen = []

    @functools.wraps(f)
    def g(*a, **k):
        show = not (once and seen)
        seen.append(1)
        if show:
            PHASE.append(name)
            log("enter " + name)
        r = f(*a, **k)
        if show:
            log("exit " + name)
            PHASE.pop()
        return r
    setattr(owner, name, g)


def ticker():
    while True:
        time.sleep(5.0)
        log("tick (in %s)" % PHASE[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = {**EB.RL.NEW_KEYS, **json.loads(pathlib.Path(a.config).read_text())}
    args = argparse.Namespace(**cfg)
    U = EB.units(args.Js, args.A)
    a_lex, d_lex = EB.box_edge(args.d_lex, args.phi, args.a_quant)
    N = int(round(a_lex / args.dx_lex))
    d, dx = d_lex * U["l_ex"], args.dx_lex * U["l_ex"]
    q, p = args.staircase_q, args.staircase_p
    for owner, name, once in ((geometry, "fcc_subcells", False), (DS, "centroid_offsets", False),
                              (DS.NearCorrection, "__init__", False), (DS.NearCorrection, "_targets", True),
                              (DS.NearCorrection, "_coeffs", True), (DS, "pair_blocks", False)):
        wrap(owner, name, once)
    threading.Thread(target=ticker, daemon=True).start()
    log("N %d, q %d, p %d" % (N, q, p))
    n, pid, occ = geometry.fcc_subcells(N, d, dx, q)
    log("geometry done: occ %.2f GB, n %.2f GB, pid %.2f GB" % (occ.nbytes / 1e9, n.nbytes / 1e9, pid.nbytes / 1e9))
    corr = DS.DemagStaircaseCorrection(n, occ, q, p, 1.0, form="pair", periodic=True, centroid=True, demag=object())
    corr._setup(types.SimpleNamespace(device=torch.empty(0).device, dtype=torch.float32,
                                      mesh=types.SimpleNamespace(n=(N, N, N))))
    log("set-up done")


if __name__ == "__main__":
    main()
