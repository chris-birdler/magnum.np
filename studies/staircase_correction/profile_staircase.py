#!/usr/bin/env python
"""
Profile of one field evaluation with the staircase correction (PROTOKOLL 3.34): wall time of each term and of the
two parts of the correction (pair sum of the window, centroid term), set-up time, pair count, peak host memory and,
on a GPU, the peak GPU memory. The model is built as run_loops builds it (energy_budget.build); the precision is that
of the runs (single). On a CPU the times give only the parts that exist, not the GPU times of the runs.

    python profile_staircase.py --config RUN/config.json --init RUN/init.pt               # the run itself
    python profile_staircase.py --config RUN/config.json --d_lex 216 --phi 0.65 --Leff_lex 30   # other geometry
"""
import argparse
import json
import pathlib
import resource
import sys
import time

import numpy as np
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import energy_budget as EB  # noqa: E402
from magnumnp import set_precision  # noqa: E402
from magnumnp.field_terms.demag_staircase import sym_matvec, centroid_field  # noqa: E402
from anisotropy_model import box_edge  # noqa: E402


def sync():
    if torch.cuda.is_available():
        torch.cuda.synchronize()


def timed(f, reps):
    f()                                                     # warm-up (compile, caches)
    sync()
    t = time.perf_counter()
    for _ in range(reps):
        f()
    sync()
    return (time.perf_counter() - t) / reps


def gpu_gb():
    """Peak GPU memory since the last reset: (allocated, reserved) in GB, or None on a CPU."""
    if not torch.cuda.is_available():
        return None
    return torch.cuda.max_memory_allocated() / 1e9, torch.cuda.max_memory_reserved() / 1e9


def peak_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


@torch.no_grad()
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="config.json of a run with --staircase pair")
    ap.add_argument("--init", default=None, help="init.pt of that run (m); without it a random m")
    ap.add_argument("--d_lex", type=float, default=None)
    ap.add_argument("--phi", type=float, default=None)
    ap.add_argument("--Leff_lex", type=float, default=None)
    ap.add_argument("--q", type=int, default=None, help="staircase_q (default: that of the run)")
    ap.add_argument("--precision", default="single")
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    set_precision(a.precision, force=True)
    cfg = json.loads(pathlib.Path(a.config).read_text())
    for k in ("d_lex", "phi", "Leff_lex"):
        if getattr(a, k) is not None:
            cfg[k] = getattr(a, k)
    if a.q is not None:
        cfg["staircase_q"] = a.q
    a_lex, d_lex = box_edge(cfg["d_lex"], cfg["phi"], cfg.get("a_quant", 12.0))
    N = int(round(a_lex / cfg["dx_lex"]))
    if a.init and a.d_lex is None and a.phi is None:
        m = torch.load(a.init, map_location="cpu", weights_only=False)["m"]
    else:
        g = torch.Generator().manual_seed(1)
        m = torch.randn((N, N, N, 3), generator=g, dtype=torch.float64)
        m /= m.norm(dim=-1, keepdim=True)
    m = m.to("cuda" if torch.cuda.is_available() else "cpu")
    t = time.perf_counter()
    state, terms, ids, Kd = EB.build(cfg, m)
    t_build = time.perf_counter() - t
    corr, dem = terms["demag_corr"], terms["demag"]
    t = time.perf_counter()
    dem.h(state)                                            # Newell kernel (first call)
    sync()
    t_kernel = time.perf_counter() - t
    gpu_before = gpu_gb()
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()
    t = time.perf_counter()
    corr.h(state)                                           # set-up of the correction (first call)
    sync()
    t_setup = time.perf_counter() - t
    gpu_setup = gpu_gb()
    mm = state.m.reshape(-1, 3)
    parts = {
        "demag (Newell, FFT)": lambda: dem.h(state),
        "exchange": lambda: terms["exchange"].h(state),
        "anisotropy": lambda: terms["aniso_cube"].h(state),
        "correction, total": lambda: corr.h(state),
        "correction, pair sum": lambda: torch.zeros_like(mm).index_add_(0, corr._tgt, sym_matvec(corr._cf, mm[corr._src])),
        "correction, centroid term": lambda: centroid_field(dem._N, corr._v, corr._d, state.m,
                                                            shape=tuple(dem._shape(state))),
    }
    times = {k: timed(f, a.reps) for k, f in parts.items()}
    base = times["demag (Newell, FFT)"] + times["exchange"] + times["anisotropy"]
    res = {"N": N, "d_lex": d_lex, "phi": cfg["phi"], "q": cfg["staircase_q"], "precision": a.precision,
           "threads": torch.get_num_threads(), "pairs": int(corr._tgt.numel()),
           "cut_cells": int(np.count_nonzero((corr._v.cpu().numpy() > 0) & (corr._v.cpu().numpy() < 1))),
           "t_build_s": t_build, "t_kernel_s": t_kernel, "t_setup_s": t_setup, "peak_host_GB": peak_gb(),
           "device": str(state.m.device), "gpu_GB_before_setup": gpu_before, "gpu_GB_setup_peak": gpu_setup,
           "gpu_GB_after": gpu_gb(),
           "ms": {k: 1e3 * v for k, v in times.items()}, "factor_terms": (base + times["correction, total"]) / base}
    print("N %d (d/l_ex %.1f, φ %.2f), q %d, %s precision, %d threads; pairs %d, cut cells %d" % (
        N, d_lex, cfg["phi"], res["q"], a.precision, res["threads"], res["pairs"], res["cut_cells"]))
    print("device %s; build %.1f s, Newell kernel %.1f s, set-up of the correction %.1f s, peak host memory %.2f GB"
          % (res["device"], t_build, t_kernel, t_setup, res["peak_host_GB"]))
    if gpu_setup is not None:
        print("GPU memory (allocated / reserved, GB): before the set-up %.2f / %.2f, peak of the set-up %.2f / %.2f"
              % (gpu_before + gpu_setup))
    for k, v in times.items():
        print("  %-28s %8.1f ms  (%.2f × demag)" % (k, 1e3 * v, v / times["demag (Newell, FFT)"]))
    print("terms without the correction %.1f ms; with it %.1f ms: factor %.2f" % (
        1e3 * base, 1e3 * (base + times["correction, total"]), res["factor_terms"]))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
