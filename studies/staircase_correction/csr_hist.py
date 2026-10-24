#!/usr/bin/env python
"""
Row lengths of the sparse matrix of the correction (PLAN 10) in the production path, for the outlier finding of
PLAN 14a: run run_loops.main with the given arguments; at the first call of LLGSolver.step set up the correction,
then print the histogram of the row lengths of the CSR matrix (the most frequent lengths, with their share of the
rows and of the entries), and for each given flat index its row length, the share of the rows with this length, the
volume fraction of its cell and the share of the cells with a row and with this fraction. Then stop.

    CUDA_DEVICE=-1 python csr_hist.py 3060786 3615150 -- <run_loops arguments>
"""
import pathlib
import sys

import torch


def main(argv):
    i = argv.index("--")
    idx = [int(x) for x in argv[:i]]
    sys.path.insert(0, str(pathlib.Path(".").resolve()))
    import magnumnp
    import run_loops

    def step(self, state, dt, **k):
        names = [type(t).__name__ for t in self._terms]
        corr = self._terms[names.index("DemagStaircaseCorrection")]
        corr.h(state)  # set-up
        crow = corr._A.crow_indices().cpu()
        L = crow[1:] - crow[:-1]
        nz = L[L > 0]
        n, nnz = nz.numel(), int(nz.sum())
        print("[hist] %d rows, %d rows with entries, %d entries; row length min %d, median %d, max %d" % (
            L.numel(), n, nnz, int(nz.min()), int(nz.median()), int(nz.max())), flush=True)
        vals, cnt = torch.unique(nz, return_counts=True)
        for o in cnt.argsort(descending=True)[:20]:
            v, c = int(vals[o]), int(cnt[o])
            print("[hist] length %4d: %8d rows (%.4f of the rows, %.4f of the entries)" % (v, c, c / n, v * c / nnz))
        for b in (256, 512, 1024):
            print("[hist] rows longer than %d: %.4f of the rows; rows of length %d: %.4f; of length %d: %.4f" % (
                b, float((nz > b).double().mean()), b, float((nz == b).double().mean()),
                b + 1, float((nz == b + 1).double().mean())))
        Ms = state.material["Ms"].reshape(-1).expand(state.m[..., 0].numel())
        frac = (Ms / Ms.max()).cpu()
        has = (L.reshape(-1, 3).sum(-1) > 0)
        for j in idx:
            cell, l = j // 3, int(L[j])
            f = float(frac[cell])
            same_f = int(((frac - f).abs() < 1e-12)[has].sum())
            print("[hist] index %d: cell %d comp %d, row length %d (rows of this length: %.4f of the rows with "
                  "entries), fraction %.6f (cells with a row and this fraction: %.4f)" % (
                      j, cell, j % 3, l, float((nz == l).double().mean()), f, same_f / int(has.sum())), flush=True)
        raise SystemExit(0)

    magnumnp.LLGSolver.step = step
    run_loops.main(argv[i + 1:])


if __name__ == "__main__":
    main(sys.argv[1:])
