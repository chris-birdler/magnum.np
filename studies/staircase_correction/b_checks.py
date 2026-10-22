#!/usr/bin/env python
"""
S1 of PLAN 11 for the set-up patches B5 and B4 (audit, PROTOKOLL 3.36) at the sizes of the main study:
  B5: fcc_subcells against staircase_cells(fcc_box(q N, d, dx / q)["ids"], q) on the fine mesh: bitwise, time, peak
      host memory;
  B4: NearCorrection of the code before B4 (--old_module, a copy of demag_staircase.py) against the present one, on
      the default device of magnum.np (the GPU with CUDA_DEVICE=0): B and d bitwise, time, peak GPU memory.

    CUDA_DEVICE=0 python b_checks.py --old_module /root/old_staircase.py --json out.json
"""
import argparse
import importlib.util
import json
import pathlib
import resource
import sys
import time

import numpy as np
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
import magnumnp  # noqa: E402,F401  (sets the default device)
from magnumnp.field_terms import demag_staircase as NEW  # noqa: E402
from geometry import fcc_box, fcc_subcells  # noqa: E402

LEX = 3.3421710327935706e-09


def gpu_peak():
    if not torch.cuda.is_available():
        return None
    return torch.cuda.max_memory_allocated() / 1e9


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old_module", required=True)
    ap.add_argument("--json", required=True)
    a = ap.parse_args()
    spec = importlib.util.spec_from_file_location("old_staircase", a.old_module)
    OLD = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(OLD)
    res = {"B5": [], "B4": []}
    for N, d_lex, q in [(104, 211.2, 8), (148, 300.6, 2), (148, 300.6, 4)]:
        d, dx = d_lex * LEX, 3 * LEX
        t = time.time()
        n1, p1, o1 = fcc_subcells(N, d, dx, q)
        t1 = time.time() - t
        t = time.time()
        n0, p0, o0 = NEW.staircase_cells(fcc_box(q * N, d, dx / q)["ids"], q)
        t0 = time.time() - t
        r = {"N": N, "d_lex": d_lex, "q": q, "bitwise": bool(np.array_equal(n0, n1) and np.array_equal(p0, p1)
                                                             and np.array_equal(o0, o1)),
             "cut_cells": int(((n1 > 0) & (n1 < q ** 3)).sum()), "t_fine_mesh_s": t0, "t_new_s": t1,
             "peak_host_GB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6}
        print("B5", r, flush=True)
        res["B5"].append(r)
        del n0, p0, o0, n1, p1, o1
    for q in (2, 4, 8):
        n, pid, occ = fcc_subcells(48, 92.214 * LEX, 3 * LEX, q)
        out = {}
        for name, mod in (("old", OLD), ("new", NEW)):
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
                torch.cuda.reset_peak_memory_stats()
            t = time.time()
            nc = mod.NearCorrection(n, occ, q, 3, centroid=True)
            if torch.cuda.is_available():
                torch.cuda.synchronize()
            out[name] = (nc.B.cpu(), nc.d.cpu(), time.time() - t, gpu_peak(), str(nc.B.device))
            del nc
        r = {"N": 48, "q": q, "bitwise_B": bool(torch.equal(out["old"][0], out["new"][0])),
             "bitwise_d": bool(torch.equal(out["old"][1], out["new"][1])),
             "rel_B": float((out["old"][0] - out["new"][0]).abs().max() / out["old"][0].abs().max()),
             "t_old_s": out["old"][2], "t_new_s": out["new"][2], "gpu_peak_old_GB": out["old"][3],
             "gpu_peak_new_GB": out["new"][3], "device": out["new"][4]}
        print("B4", r, flush=True)
        res["B4"].append(r)
    pathlib.Path(a.json).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
