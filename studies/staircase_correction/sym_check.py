#!/usr/bin/env python
"""
PLAN 15 WQ (information): the v-weighted symmetry of the pair operator of DemagStaircaseCorrection on a production
geometry, as test_correction_is_symmetric: with 2 random fields m1, m2 (seeds 2 and 3), a = sum v m1 . (A m2) and
b = sum v m2 . (A m1), A the PairBlocks operator of the pair form (PBC) in float64, v the volume fraction. Reports
|a - b| / (|a| + |b|) against the bound 1e-12 of the unit test, and the time and the peak GPU memory of the set-up.
The geometry is that of run_loops (fcc_subcells of the tree on PYTHONPATH), from the config.json of a run.

    PYTHONPATH=TREE python sym_check.py --config RUN/config.json [--json OUT] [--q Q --dx_lex DX (local test)]

Exit code 0 if the relative asymmetry is at most 1e-12, else 1.
"""
import argparse
import json
import pathlib
import sys
import time
import types

import numpy as np
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
import energy_budget as EB  # noqa: E402
import magnumnp  # noqa: E402
from magnumnp.field_terms.demag_staircase import DemagStaircaseCorrection  # noqa: E402
# the geometry of the tree on PYTHONPATH (fcc_subcells, audit B5), not the one of this study folder
sys.path.insert(0, str(pathlib.Path(magnumnp.__file__).resolve().parents[1] / "studies" / "beta_bpeak"))
sys.modules.pop("geometry", None)
from geometry import fcc_subcells  # noqa: E402


def random_m(shape, seed, device):
    g = torch.Generator(device="cpu").manual_seed(seed)
    m = torch.randn(shape + (3,), generator=g, dtype=torch.float64)
    return (m / m.norm(dim=-1, keepdim=True)).to(device)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="config.json of a run (geometry and staircase arguments)")
    ap.add_argument("--q", type=int, default=None, help="(local test) another staircase_q")
    ap.add_argument("--dx_lex", type=float, default=None, help="(local test) another dx_lex")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    cfg = {**EB.RL.NEW_KEYS, **json.loads(pathlib.Path(a.config).read_text())}
    if a.q is not None:
        cfg["staircase_q"] = a.q
    if a.dx_lex is not None:
        cfg["dx_lex"] = a.dx_lex
    args = argparse.Namespace(**cfg)
    U = EB.units(args.Js, args.A)
    a_lex, d_lex = EB.box_edge(args.d_lex, args.phi, args.a_quant)
    N = int(round(a_lex / args.dx_lex))
    d, dx = d_lex * U["l_ex"], args.dx_lex * U["l_ex"]
    q, p = args.staircase_q, args.staircase_p
    t = time.perf_counter()
    n, pid, occ = fcc_subcells(N, d, dx, q)
    t_geo = time.perf_counter() - t
    dev = torch.empty(0).device
    print("geometry: N %d, q %d, p %d; build %.1f s; device %s" % (N, q, p, t_geo, dev), flush=True)
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()
    t = time.perf_counter()
    corr = DemagStaircaseCorrection(n, occ, q, p, 1.0, form="pair", periodic=True, centroid=False, demag=object())
    corr._setup(types.SimpleNamespace(device=dev, dtype=torch.float64, mesh=types.SimpleNamespace(n=(N, N, N))))
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    t_setup = time.perf_counter() - t
    peak = torch.cuda.max_memory_allocated() / 2**30 if torch.cuda.is_available() else float("nan")
    A = corr._A
    v = torch.as_tensor(np.asarray(n) / q ** 3, dtype=torch.float64, device=dev).reshape(-1, 1)
    m1, m2 = random_m((N, N, N), 2, dev).reshape(-1, 3), random_m((N, N, N), 3, dev).reshape(-1, 3)
    h1 = (A @ m1.reshape(-1)).reshape(-1, 3)
    h2 = (A @ m2.reshape(-1)).reshape(-1, 3)
    ab = float((v * m1 * h2).sum())
    ba = float((v * m2 * h1).sum())
    rel = abs(ab - ba) / (abs(ab) + abs(ba))
    ok = rel <= 1e-12
    res = {"N": N, "q": q, "p": p, "blocks": int(A.V.shape[1] * A.V.shape[2]) if A.V.dim() == 3 else None,
           "a": ab, "b": ba, "rel_asymmetry": rel, "bound": 1e-12, "pass": ok, "t_geometry_s": t_geo,
           "t_setup_float64_s": t_setup, "peak_setup_float64_GiB": peak}
    print("a %.17e  b %.17e  relative asymmetry %.2e (bound 1e-12) -> %s; set-up float64 %.1f s, peak %.2f GiB"
          % (ab, ba, rel, "PASS" if ok else "FAIL", t_setup, peak), flush=True)
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(res, indent=1))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
