#!/usr/bin/env python
"""
PLAN 16 C1: new checkpoint states for the decision on new data. Write a copy of a K1 reference file (step_dump.py
--k1_ref; the record "ck" holds m, t and dt of each checkpoint state) with each state m perturbed to m + eps r, as
perturb_init.py: in fp64, r random with |r| 1 per cell (torch generator, seed SEED * 1000 + k0, so that each state gets
its own r), only in the cells with |m| above 0; the sum is renormalized and cast back to the dtype of m; the void keeps
m = 0. t, dt and all other records stay. Report for each state: RMS and maximum of |dm| after the cast, the share of
the magnetic cells whose float32 value changed (s1_states.py evaluates on m.float()), void unchanged.

    python perturb_ref.py SRC/refB.pt DST/refB_p3.pt --seed 3 --eps 1e-7

Exit code 0 if every state has a changed float32 share above 0.9 and an unchanged void, else 1.
"""
import argparse
import pathlib
import sys

import torch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--eps", type=float, default=1e-7)
    a = ap.parse_args()
    ref = torch.load(a.src, map_location="cpu", weights_only=False)
    ok = True
    ck = {}
    for k0, c in sorted(ref["ck"].items()):
        m = c["m"]
        m64 = m.double()
        mag = m64.norm(dim=-1, keepdim=True) > 0
        g = torch.Generator().manual_seed(a.seed * 1000 + int(k0))
        r = torch.randn(m64.shape, generator=g, dtype=torch.float64)
        r = r / r.norm(dim=-1, keepdim=True)
        p = m64 + a.eps * r
        p = torch.where(mag, p / p.norm(dim=-1, keepdim=True), m64).to(m.dtype)
        d = (p.double() - m64).norm(dim=-1)
        ch32 = (p.float() != m.float()).any(dim=-1)[mag[..., 0]].double().mean()
        void = bool((p[~mag[..., 0]] == m[~mag[..., 0]]).all())
        ok = ok and float(ch32) > 0.9 and void
        print("k0 %3d: RMS |dm| %.3e, max |dm| %.3e, changed float32 share %.4f of the magnetic cells, void unchanged %s"
              % (k0, float(d[mag[..., 0]].pow(2).mean().sqrt()), float(d.max()), float(ch32), void))
        ck[k0] = dict(c, m=p)
    out = dict(ref, ck=ck, perturbed={"src": str(a.src), "seed": a.seed, "eps": a.eps})
    pathlib.Path(a.dst).parent.mkdir(parents=True, exist_ok=True)
    torch.save(out, a.dst)
    print("perturb_ref: %d states, seed %d, eps %g -> %s: %s" % (len(ck), a.seed, a.eps, a.dst, "PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
