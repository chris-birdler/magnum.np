#!/usr/bin/env python
"""
S1 of PLAN 11 for the patch A2 (DemagField.h) on a production state and the GPU: run run_loops.main with the given
arguments; at the first call of DemagField.h after the cycles start (first LLGSolver.step), compare h() with the
reference (the loop over the components of the code before A2, the block term of true_pbc as a separate term) and
time both (20 calls each after a warm-up), then stop.

    CUDA_DEVICE=0 python s1_demag.py -- <run_loops arguments>
"""
import pathlib
import sys
import time

import torch


def reference(dem, state):
    from magnumnp.field_terms.demag import block_demag_factors
    dim = [i for i in range(3) if state.mesh.n[i] > 1]
    s = [dem._shape(state)[i] for i in dim]
    M = state.material["Ms"] * state.m
    hk = [0, 0, 0]
    for ax in range(3):
        mk = torch.fft.rfftn(M[..., ax], dim=dim, s=s)
        for i in range(3):
            hk[i] = hk[i] + dem._N[i][ax] * mk
    h = torch.stack([torch.fft.irfftn(x, s=s, dim=dim)[:state.mesh.n[0], :state.mesh.n[1], :state.mesh.n[2]]
                     for x in hk], dim=3)
    if dem._true_pbc:
        h = h + torch.tensor(block_demag_factors(state), dtype=h.dtype, device=h.device) * M.mean(dim=(0, 1, 2))
    return h


def timed(f, reps=20):
    f()
    torch.cuda.synchronize()
    t = time.perf_counter()
    for _ in range(reps):
        f()
    torch.cuda.synchronize()
    return 1e3 * (time.perf_counter() - t) / reps


def main(argv):
    i = argv.index("--")
    sys.path.insert(0, str(pathlib.Path(".").resolve()))
    import magnumnp
    import run_loops
    f_h, f_step = magnumnp.DemagField.h, magnumnp.LLGSolver.step
    st = {"steps": 0}

    def step(self, state, dt, **k):
        st["steps"] += 1
        return f_step(self, state, dt, **k)

    def h(self, state):
        r = f_h(self, state)
        if st["steps"] == 0:
            return r
        ref = reference(self, state)
        Ms = state.material["Ms"]
        dh = float((r.double() - ref.double()).abs().max() / ref.double().abs().max())
        E, E_ref = [float((Ms * state.m * x).double().sum()) for x in (r, ref)]
        t_new, t_ref = timed(lambda: f_h(self, state)), timed(lambda: reference(self, state))
        print("[s1_demag] N %s, %s: max |dh| / max |h| %.3e, relative difference of sum Ms m . h %.3e; "
              "h() %.3f ms, reference %.3f ms" % (tuple(state.mesh.n), r.dtype, dh, abs(E / E_ref - 1), t_new, t_ref),
              flush=True)
        raise SystemExit(0)

    magnumnp.DemagField.h = h
    magnumnp.LLGSolver.step = step
    run_loops.main(argv[i + 1:])


if __name__ == "__main__":
    main(sys.argv[1:])
