#!/usr/bin/env python
"""
Null run of PLAN 11 (kind ii): write a copy of an init.pt with a perturbed start state m + ε r. In fp64: r is a random
field with |r| 1 per cell (torch generator with the given seed), only in the cells with |m| above 0; the sum is
renormalized and cast back to the dtype of m. The void keeps m = 0. The records "init" and "config" stay, so the run
loads it with --init_from like the original. Report: RMS and maximum of |Δm| after the cast, share of changed cells.

    python perturb_init.py SRC/init.pt DST/init.pt --seed 1 --eps 1e-7
"""
import argparse
import pathlib

import torch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--eps", type=float, default=1e-7)
    a = ap.parse_args()
    c = torch.load(a.src, map_location="cpu", weights_only=False)
    m = c["m"]
    m64 = m.double()
    mag = m64.norm(dim=-1, keepdim=True) > 0
    g = torch.Generator().manual_seed(a.seed)
    r = torch.randn(m64.shape, generator=g, dtype=torch.float64)
    r = r / r.norm(dim=-1, keepdim=True)
    p = m64 + a.eps * r
    p = torch.where(mag, p / p.norm(dim=-1, keepdim=True), m64).to(m.dtype)
    d = (p.double() - m64).norm(dim=-1)
    changed = (p != m).any(dim=-1)
    print("perturbed %s -> %s: seed %d, eps %g, dtype %s; RMS |dm| %.3e, max |dm| %.3e (magnetic cells), "
          "changed cells %.4f of the magnetic cells; void unchanged %s"
          % (a.src, a.dst, a.seed, a.eps, m.dtype, float(d[mag[..., 0]].pow(2).mean().sqrt()), float(d.max()),
             float(changed[mag[..., 0]].double().mean()), bool((p[~mag[..., 0]] == m[~mag[..., 0]]).all())))
    c = dict(c, m=p)
    pathlib.Path(a.dst).parent.mkdir(parents=True, exist_ok=True)
    torch.save(c, a.dst)


if __name__ == "__main__":
    main()
