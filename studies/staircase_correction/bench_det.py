#!/usr/bin/env python
"""
GPU check of the fixed-order pair sum (PairBlocks, PLAN 14b) in the production path (box 11). Run run_loops.main with
the given arguments of a tree with PairBlocks; at the first call of LLGSolver.step take the checkpoint states --k0 of
the K1 reference, set up the correction and compare, in double and in single precision:

- csr: the CSR product of the same merged blocks (pair_matrix, the product of PLAN 10 with cuSPARSE);
- eager: PairBlocks @ m (block_sum without torch.compile);
- compiled: the same with torch.compile(block_sum).

For each: the largest difference to csr (relative to the largest entry of csr), the count of products that are not
bitwise equal to the first product on the same state (--reps products, on the states in turn, a new tensor of m each
time, as in the stages) and their largest difference, the time per product (median of --time products, CUDA events),
and the peak of the added GPU memory of one product.

    CUDA_DEVICE=0 python bench_det.py --ref refB.pt --k0 8 144 --json out.json -- <run_loops arguments>
"""
import argparse
import json
import pathlib
import statistics
import sys
import time

import torch


def main(argv):
    i = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--k0", type=int, nargs="+", default=[8, 144])
    ap.add_argument("--reps", type=int, default=20000)
    ap.add_argument("--time", type=int, default=300)
    ap.add_argument("--modes", nargs="+", default=["csr", "eager", "compiled"])
    ap.add_argument("--dtypes", nargs="+", default=["float64", "float32"])
    ap.add_argument("--json")
    a = ap.parse_args(argv[:i])
    sys.path.insert(0, str(pathlib.Path(".").resolve()))
    import magnumnp
    import run_loops
    ref = torch.load(a.ref, map_location="cpu", weights_only=False)
    out = {"argv": argv[i + 1:], "k0": a.k0, "reps": a.reps, "torch": torch.__version__, "res": {}}

    def step(self, state, dt, **k):
        dev = state.m.device
        names = [type(t).__name__ for t in self._terms]
        corr = self._terms[names.index("DemagStaircaseCorrection")]
        mod = sys.modules[type(corr).__module__]
        state.m = ref["ck"][a.k0[0]]["m"].to(device=dev, dtype=state.m.dtype).clone()
        corr.h(state)  # set-up
        B = corr._A
        if not isinstance(B, getattr(mod, "PairBlocks", ())):
            raise SystemExit("bench_det: the correction of this tree has no PairBlocks (%s)" % type(B).__name__)
        # the merged blocks back from the layout (slots of real chunks with a nonzero coefficient), for pair_matrix;
        # on the CPU (magnum.np sets the default device to the GPU)
        n_chunk, C = B.S.shape[0] - 1, B.S.shape[1]
        with torch.device("cpu"):
            T, rows = B.T.cpu().long(), B.rows.cpu().long()
            row_of_chunk = torch.full((n_chunk + 1,), -1, dtype=torch.long)
            ok = T < n_chunk
            row_of_chunk[T[ok]] = rows[:, None].expand_as(T)[ok]
            Vc, Sc = B.V.cpu().reshape(6, -1), B.S.cpu().reshape(-1).long()
            slot = torch.nonzero((row_of_chunk.repeat_interleave(C) >= 0) & (Vc != 0).any(0)).reshape(-1)
            tgt, src, cf = row_of_chunk.repeat_interleave(C)[slot], Sc[slot], Vc[:, slot].T.contiguous()
        out["blocks"] = {"n_rows": len(rows), "n_chunk": n_chunk, "chunk": C, "L": int(T.shape[1]),
                         "n_blocks": len(slot), "slots": (n_chunk + 1) * C}
        print("[bench] %d rows, %d blocks, %d chunks of %d (%d slots), at most %d chunks per row" % (
            len(rows), len(slot), n_chunk, C, (n_chunk + 1) * C, T.shape[1]), flush=True)
        cb = torch.compile(mod.block_sum, dynamic=False)
        for dn in a.dtypes:
            dt_ = getattr(torch, dn)
            A = mod.pair_matrix(tgt, src, cf, B.n_cells, dtype=dt_).to(dev)
            P = mod.PairBlocks(B.V.to(dt_), B.S, B.T, B.rows, B.n_cells)
            M = [ref["ck"][k0]["m"].to(device=dev, dtype=dt_) for k0 in a.k0]

            def f_compiled(x, P=P):
                m = x.reshape(-1, 3)
                return torch.zeros_like(m).index_copy_(0, P.rows, cb(P.V, P.S, P.T, m)).reshape(-1)

            fs = {"csr": lambda x, A=A: A @ x, "eager": lambda x, P=P: P @ x, "compiled": f_compiled}
            csr = [(A @ m.clone().reshape(-1)) for m in M]
            scale = [float(h.abs().max()) for h in csr]
            res = out["res"][dn] = {}
            for nm in a.modes:
                f = fs[nm]
                torch.cuda.synchronize()
                base = torch.cuda.memory_allocated()
                torch.cuda.reset_peak_memory_stats()
                first = [f(m.clone().reshape(-1)) for m in M]
                torch.cuda.synchronize()
                peak = torch.cuda.max_memory_allocated() - base
                d_csr = max(float((h - c).abs().max()) / s for h, c, s in zip(first, csr, scale))
                n_neq = torch.zeros((), dtype=torch.long, device=dev)
                d_max = torch.zeros((), dtype=torch.float64, device=dev)
                t0 = time.time()
                for r in range(a.reps):
                    s = r % len(M)
                    d = (f(M[s].clone().reshape(-1)) - first[s]).abs()
                    n_neq += (d.max() > 0).long()
                    d_max = torch.maximum(d_max, d.max().double())
                torch.cuda.synchronize()
                wall = time.time() - t0
                ev = [(torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)) for _ in range(a.time)]
                x = M[0].clone().reshape(-1)
                for e0, e1 in ev:
                    e0.record()
                    f(x)
                    e1.record()
                torch.cuda.synchronize()
                ms = statistics.median(e0.elapsed_time(e1) for e0, e1 in ev)
                res[nm] = {"rel_diff_to_csr": d_csr, "reps": a.reps, "not_bitwise_equal": int(n_neq),
                           "largest_repeat_diff": float(d_max), "ms_per_product": ms, "peak_extra_MiB": peak / 2 ** 20,
                           "wall_s_reps": wall}
                print("[bench] %s %-8s: relative difference to csr %.3e; %d of %d repeats not bitwise equal (largest "
                      "difference %.3e); %.3f ms per product; added peak memory %.0f MiB" % (
                          dn, nm, d_csr, int(n_neq), a.reps, float(d_max), ms, peak / 2 ** 20), flush=True)
                if a.json:
                    pathlib.Path(a.json).write_text(json.dumps(out, indent=1))
            if "eager" in res and "compiled" in res:
                pe = [P @ m.clone().reshape(-1) for m in M]
                pc = [f_compiled(m.clone().reshape(-1)) for m in M]
                res["compiled_vs_eager"] = max(float((x - y).abs().max()) / s for x, y, s in zip(pc, pe, scale))
                print("[bench] %s compiled against eager: relative difference %.3e" % (dn, res["compiled_vs_eager"]),
                      flush=True)
            del A, P, M, csr
            torch.cuda.empty_cache()
        if a.json:
            pathlib.Path(a.json).write_text(json.dumps(out, indent=1))
        raise SystemExit(0)

    magnumnp.LLGSolver.step = step
    run_loops.main(argv[i + 1:])


if __name__ == "__main__":
    main(sys.argv[1:])
