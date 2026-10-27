#!/usr/bin/env python
"""
PLAN 15 FQ, checks W1 and WFQ: the set-up of the staircase correction of two code trees on one production geometry,
compared bitwise. For each tree (its demag_staircase.py loaded as a separate module): the centroid offsets, the
coupling blocks B of NearCorrection and the tensors of DemagStaircaseCorrection._setup (V, S, T and rows of
PairBlocks, the centroid offsets _d and the volume fractions _v), each as the md5 of its bytes (so +0.0 and -0.0
differ), and the time and the peak GPU memory of the set-up. The new tree also computes the centroid offsets with a
small byte budget (several blocks and a partial last block). The geometry is that of run_loops (fcc_subcells, the
production path); the other modules (magnumnp, the demag kernel) come from PYTHONPATH.

    PYTHONPATH=NEW python fq_check.py --old OLD/magnumnp/field_terms/demag_staircase.py \
        --new NEW/magnumnp/field_terms/demag_staircase.py --config RUN/config.json [--q 8] [--d_lex 92 --phi 0.55]

Exit code 0 if all md5 of the two trees are equal (and the small-budget offsets equal the default ones), else 1.
"""
import argparse
import hashlib
import importlib.util
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
# the geometry of the tree on PYTHONPATH (fcc_subcells, audit B5), not the one of this study folder
sys.path.insert(0, str(pathlib.Path(magnumnp.__file__).resolve().parents[1] / "studies" / "beta_bpeak"))
sys.modules.pop("geometry", None)
from geometry import fcc_subcells  # noqa: E402


def md5(t):
    return hashlib.md5(t.detach().cpu().contiguous().numpy().tobytes()).hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def gpu_reset():
    if torch.cuda.is_available():
        torch.cuda.synchronize()
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()


def gpu_peak():
    if torch.cuda.is_available():
        torch.cuda.synchronize()
        return torch.cuda.max_memory_allocated() / 2**30, torch.cuda.max_memory_reserved() / 2**30
    return float("nan"), float("nan")


def setup(mod, n, occ, q, p, N, B_md5):
    """md5 of the tensors of one set-up of DemagStaircaseCorrection (pair form, PBC, centroid), as in production."""
    init = mod.NearCorrection.__init__

    def keep(self, *a, **k):                       # the md5 of B of the NearCorrection inside _setup
        init(self, *a, **k)
        B_md5.append(md5(self.B))
    mod.NearCorrection.__init__ = keep
    try:
        corr = mod.DemagStaircaseCorrection(n, occ, q, p, 1.0, form="pair", periodic=True, centroid=True,
                                            demag=object())
        state = types.SimpleNamespace(device=torch.empty(0).device, dtype=torch.float32,
                                      mesh=types.SimpleNamespace(n=(N, N, N)))
        corr._setup(state)
    finally:
        mod.NearCorrection.__init__ = init
    A = corr._A
    return {"V": md5(A.V), "S": md5(A.S), "T": md5(A.T), "rows": md5(A.rows), "_d": md5(corr._d),
            "_v": md5(corr._v), "B": B_md5[-1], "blocks": int(A.V.shape[0])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True, help="demag_staircase.py of the old tree")
    ap.add_argument("--new", required=True, help="demag_staircase.py of the new tree")
    ap.add_argument("--config", required=True, help="config.json of a run (geometry and staircase arguments)")
    ap.add_argument("--q", type=int, default=None)
    ap.add_argument("--d_lex", type=float, default=None)
    ap.add_argument("--phi", type=float, default=None)
    ap.add_argument("--Leff_lex", type=float, default=None)
    ap.add_argument("--small_budget_cells", type=int, default=3, help="cells per block of the small-budget run")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    cfg = json.loads(pathlib.Path(a.config).read_text())
    cfg = {**EB.RL.NEW_KEYS, **cfg}
    for k in ("q", "d_lex", "phi", "Leff_lex"):
        v = getattr(a, k)
        if v is not None:
            cfg["staircase_q" if k == "q" else k] = v
    args = argparse.Namespace(**cfg)
    U = EB.units(args.Js, args.A)
    a_lex, d_lex = EB.box_edge(args.d_lex, args.phi, args.a_quant)
    N = int(round(a_lex / args.dx_lex))
    d, dx = d_lex * U["l_ex"], args.dx_lex * U["l_ex"]
    q, p = args.staircase_q, args.staircase_p
    t = time.perf_counter()
    n, pid, occ = fcc_subcells(N, d, dx, q)
    nq = q ** 3
    cut = int(((n > 0) & (n < nq)).sum())
    print("geometry: N %d, d/l_ex %.1f, phi %.2f, q %d, p %d; cut cells %d; occ %.2f GB (host); build %.1f s; device %s"
          % (N, d_lex, args.phi, q, p, cut, occ.nbytes / 1e9, time.perf_counter() - t, torch.empty(0).device),
          flush=True)
    res = {"N": N, "q": q, "p": p, "cut_cells": cut}
    ok = True
    co = {}
    for tag, path in (("old", a.old), ("new", a.new)):
        mod = load(path, "ds_" + tag)
        gpu_reset()
        t = time.perf_counter()
        c = mod.centroid_offsets(occ, q)
        co[tag] = md5(c)
        r = {"centroid_offsets": co[tag], "t_offsets_s": time.perf_counter() - t,
             "peak_offsets_GiB": gpu_peak()[0]}
        if tag == "new":
            c2 = mod.centroid_offsets(occ, q, budget=a.small_budget_cells * 2 * 8 * nq)
            r["centroid_offsets_small_budget"] = md5(c2)
            ok &= r["centroid_offsets_small_budget"] == co[tag]
            del c2
        del c
        gpu_reset()
        t = time.perf_counter()
        r.update(setup(mod, n, occ, q, p, N, []))
        r["t_setup_s"] = time.perf_counter() - t
        r["peak_setup_GiB"], r["peak_setup_reserved_GiB"] = gpu_peak()
        res[tag] = r
        print("%s: %s" % (tag, json.dumps(r)), flush=True)
    for k in ("centroid_offsets", "V", "S", "T", "rows", "_d", "_v", "B", "blocks"):
        same = res["old"][k] == res["new"][k]
        ok &= same
        print("%-18s %s" % (k, "equal" if same else "DIFFERENT"))
    if "centroid_offsets_small_budget" in res["new"]:
        print("%-18s %s" % ("small budget", "equal" if res["new"]["centroid_offsets_small_budget"] == co["new"]
                            else "DIFFERENT"))
    res["pass"] = bool(ok)
    print("fq_check: %s" % ("PASS (all md5 equal)" if ok else "FAIL"), flush=True)
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(res, indent=1))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
