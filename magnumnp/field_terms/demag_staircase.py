#
# This file is part of the magnum.np distribution
# (https://gitlab.com/magnum.np/magnum.np).
# Copyright (c) 2023 magnum.np team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
"""
Near-field correction of the stray field for cut surface cells, after Donahue & McMichael, IEEE Trans. Magn. 43, 2878
(2007).

A fine mesh with q sub-cells per axis and cell carries the body (the sub-cells whose centre is inside it). A cell with
n_i > 0 body sub-cells is magnetic, with the volume fraction v_i = n_i / q^3. DemagField gives the far field with the
moments M_s v_j m_j (material "Ms" = v M_s). Inside the window W(i) (max |offset| <= p, minimum image under PBC) the
correction replaces the coarse coupling v_j K^c(i - j) by the coupling L_ij of the body sub-cells:

    h_i = sum_{j in W(i)} (L_ij - v_j K^c(i - j)) m_j                 (pairwise form, their eq. 5)
    h_i = [sum_{j in W(i)} (L_ij - v_j K^c(i - j))] m_i               (cell-local form, their eqs. 8 and 9)

L_ij = (1 / n_i) sum_{a in i} sum_{b in j} K^f(a - b) over the body sub-cells. Sign convention of magnum.np: H = K M
(self term about -1/3). The kernels come from demag_f and demag_g with the Newell radius of DemagField. v_i L_ij =
v_j L_ji, so both forms give a symmetric operator and a consistent energy.
"""
import numpy as np
import torch

from magnumnp.common import constants
from magnumnp.field_terms.field_terms import LinearFieldTerm
from magnumnp.field_terms.demag import demag_f, demag_g, _kernel_product

__all__ = ["DemagStaircaseCorrection", "DemagStaircaseFine", "NearCorrection", "staircase_cells", "centroid_offsets",
           "centroid_field", "centroid_field_fourier", "pair_matrix", "pair_blocks", "PairBlocks"]

NEWELL_P = 20          # Newell radius of DemagField (in cell diagonals), the same on both meshes
SYM = ((0, 1, 2), (1, 3, 4), (2, 4, 5))   # 3x3 index -> component of (xx, xy, xz, yy, yz, zz)


def kernel_table(R, p=NEWELL_P):
    """Free-space kernel K(offset) [2R+1]^3 x 6 (xx, xy, xz, yy, yz, zz) for integer offsets in [-R, R]^3 of unit
    cubes, from demag_f / demag_g with the same component permutations as DemagField._init_N."""
    r = torch.arange(-R, R + 1, dtype=torch.float64)
    X, Y, Z = torch.meshgrid(r, r, r, indexing="ij")
    one = (1.0, 1.0, 1.0)
    comps = [demag_f(X, Y, Z, *one, *one, p), demag_g(X, Y, Z, *one, *one, p), demag_g(X, Z, Y, *one, *one, p),
             demag_f(Y, Z, X, *one, *one, p), demag_g(Y, Z, X, *one, *one, p), demag_f(Z, X, Y, *one, *one, p)]
    return torch.stack(comps, dim=-1)


def sym_matvec(c, m):
    """c [..., 6] (xx, xy, xz, yy, yz, zz), m [..., 3] -> [..., 3]."""
    return torch.stack([c[..., 0] * m[..., 0] + c[..., 1] * m[..., 1] + c[..., 2] * m[..., 2],
                        c[..., 1] * m[..., 0] + c[..., 3] * m[..., 1] + c[..., 4] * m[..., 2],
                        c[..., 2] * m[..., 0] + c[..., 4] * m[..., 1] + c[..., 5] * m[..., 2]], dim=-1)


def pair_matrix(tgt, src, cf, n_cells, dtype=None):
    """Sparse matrix (CSR, int32 indices) of the pair sum h_t += C m_s over the pairs (tgt, src, cf), C the symmetric
    3x3 matrix of the 6 components cf (SYM): rows 3 t + a, columns 3 s + b, size 3 n_cells. Built on the CPU; pairs
    with the same (t, s) are summed. One sparse product replaces the gather, sym_matvec and the atomic index_add of
    each evaluation (PROTOKOLL 3.36). All work is on the CPU, also when magnum.np set the default device to a GPU."""
    with torch.device("cpu"):
        return _pair_matrix_cpu(tgt, src, cf, n_cells, dtype)


def _merge_pairs(tgt, src, cf, n_cells, dtype):
    """The pairs sorted by (t, s); the coefficients of pairs with the same (t, s) summed (CPU)."""
    tgt, src = tgt.cpu().long(), src.cpu().long()
    cf = cf.cpu()
    key, order = torch.sort(tgt * n_cells + src)
    cf = cf[order]
    key, inv = torch.unique_consecutive(key, return_inverse=True)
    if len(key) < len(inv):
        cf = torch.zeros(len(key), 6, dtype=cf.dtype).index_add_(0, inv, cf)
    cf = cf.to(dtype if dtype is not None else cf.dtype)
    return key // n_cells, key % n_cells, cf


def _pair_matrix_cpu(tgt, src, cf, n_cells, dtype):
    t, s, cf = _merge_pairs(tgt, src, cf, n_cells, dtype)
    nb = len(t)
    if 9 * nb >= 2 ** 31:
        raise ValueError("pair_matrix: %d entries do not fit int32 indices" % (9 * nb))
    K = torch.bincount(t, minlength=n_cells)                         # blocks per block row
    j = torch.arange(nb) - (torch.cumsum(K, 0) - K)[t]               # place of each block in its block row
    crow = torch.zeros(3 * n_cells + 1, dtype=torch.int64)
    crow[1:] = torch.cumsum(torch.repeat_interleave(3 * K, 3), 0)    # rows 3 t + a hold 3 K_t entries each
    col = torch.empty(9 * nb, dtype=torch.int32)
    val = torch.empty(9 * nb, dtype=cf.dtype)
    c3 = torch.arange(3)
    for a in range(3):
        pos = ((crow[3 * t + a] + 3 * j)[:, None] + c3).reshape(-1)  # row 3 t + a, block j, column b
        col[pos] = (3 * s[:, None] + c3).reshape(-1).to(torch.int32)
        val[pos] = cf[:, list(SYM[a])].reshape(-1)
    return torch.sparse_csr_tensor(crow.to(torch.int32), col, val, size=(3 * n_cells, 3 * n_cells))


def pair_blocks(tgt, src, cf, n_cells, dtype=None, chunk=16):
    """The pair sum of pair_matrix as PairBlocks: the same merged pairs, summed in a fixed order. Built on the CPU,
    also when magnum.np set the default device to a GPU."""
    with torch.device("cpu"):
        return _pair_blocks_cpu(tgt, src, cf, n_cells, dtype, chunk)


def _pair_blocks_cpu(tgt, src, cf, n_cells, dtype, C):
    t, s, cf = _merge_pairs(tgt, src, cf, n_cells, dtype)
    rows, K = torch.unique_consecutive(t, return_counts=True)       # block rows with blocks, blocks per row
    if len(torch.unique(rows)) != len(rows):
        raise ValueError("pair_blocks: a block row occurs twice")  # the write h[rows] must not accumulate
    nch = (K + C - 1) // C                                           # chunks per row
    n_chunk, n_rows = int(nch.sum()), len(rows)
    if (n_chunk + 1) * C >= 2 ** 31:
        raise ValueError("pair_blocks: %d blocks do not fit int32 indices" % ((n_chunk + 1) * C))
    r = torch.repeat_interleave(torch.arange(n_rows), K)             # row of each block
    j = torch.arange(len(t)) - (torch.cumsum(K, 0) - K)[r]           # place of each block in its row
    c0 = torch.cumsum(nch, 0) - nch                                  # first chunk of each row
    slot = (c0[r] + j // C) * C + j % C
    S = torch.zeros((n_chunk + 1) * C, dtype=torch.int32)            # padding: source cell 0 with zero coefficients
    V = torch.zeros(6, (n_chunk + 1) * C, dtype=cf.dtype)            # chunk n_chunk: zeros only (padding of T)
    S[slot] = s.to(torch.int32)
    V[:, slot] = cf.T
    k = torch.arange(n_chunk) - torch.repeat_interleave(c0, nch)     # place of each chunk in its row
    T = torch.full((n_rows, int(nch.max()) if n_rows else 1), n_chunk, dtype=torch.int32)
    T[torch.repeat_interleave(torch.arange(n_rows), nch), k] = torch.arange(n_chunk, dtype=torch.int32)
    return PairBlocks(V.reshape(6, n_chunk + 1, C), S.reshape(n_chunk + 1, C), T, rows, n_cells)


@torch.compile
def block_sum(V, S, T, m):
    """The two stages of PairBlocks (compiled: in eager mode the products of sym_matvec make about 15 temporary arrays
    of the size of the blocks; 60 against 9 ms per product, float64, geometry B, PLAN 14b D1). V [6, n_chunk + 1, C]
    coefficients, S [n_chunk + 1, C] source cells, T [n_rows, L] chunks of each row, m [n_cells, 3] -> h [n_rows, 3]."""
    g = m.index_select(0, S.reshape(-1)).reshape(S.shape + (3,))
    p = sym_matvec(V.movedim(0, -1), g).sum(-2)                     # stage 1: sum over the blocks of each chunk
    return p.index_select(0, T.reshape(-1)).reshape(T.shape + (3,)).sum(-2)   # stage 2: over the chunks of each row


class PairBlocks:
    """The pair sum h_t = sum_s C_ts m_s of pair_matrix, summed in a fixed order, without cuSPARSE and without atomic
    additions (PLAN 14b; the cuSPARSE product had transient faults, PROTOKOLL 3.40). The merged blocks of each block
    row t (target cell) are cut into chunks of C blocks; the last chunk of a row is padded with blocks of zero
    coefficients. Stage 1 sums the blocks of each chunk, stage 2 the chunks of each row (T is padded with one added
    chunk of zeros). Only gather, element-wise operations and sums over a fixed axis, and a write to unique rows: the
    result is bitwise the same on repetition. A @ x as for the matrix (x and the result flat, size 3 n_cells), so the
    tools that wrap the product of the correction (step_dump.SparsePart) see this product."""

    def __init__(self, V, S, T, rows, n_cells):
        self.V, self.S, self.T, self.rows, self.n_cells = V, S, T, rows, n_cells

    def to(self, device):
        return PairBlocks(self.V.to(device), self.S.to(device), self.T.to(device), self.rows.to(device), self.n_cells)

    def __matmul__(self, x):
        m = x.reshape(-1, 3)
        return torch.zeros_like(m).index_copy_(0, self.rows, block_sum(self.V, self.S, self.T, m)).reshape(-1)


def blocks(fine, q):
    """fine [qN]^3 (any dtype) -> [N, N, N, q^3]: the q^3 sub-cells of each coarse cell, sub-cell order (a, b, c)."""
    Nf = fine.shape[0]
    N = Nf // q
    return fine.reshape(N, q, N, q, N, q).transpose(0, 2, 4, 1, 3, 5).reshape(N, N, N, q ** 3)


def staircase_cells(ids_fine, q):
    """Cut cells from the fine particle map ids_fine ([qN]^3, -1 = void, 0..3 = particle).
    Returns (n [N,N,N] int, pid [N,N,N] int, occ [N,N,N,q^3] bool). Reading R0: raises if a coarse cell holds fine
    cells of two particles."""
    b = blocks(np.asarray(ids_fine), q).astype(np.int16)
    occ = b >= 0
    n = occ.sum(-1).astype(np.int64)
    pid = b.max(-1).astype(np.int64)
    lo = np.where(occ, b, 99).min(-1)
    if np.any((n > 0) & (lo != pid)):
        raise ValueError("R0: a coarse cell holds fine cells of two particles")
    pid[n == 0] = -1
    return n, pid, occ


def centroid_offsets(occ, q, budget=1 << 30):
    """Offset of the centroid of the body sub-cells of each cell from its centre, in cells: [N, N, N, 3] (0 where a
    cell has no body sub-cell or all of them). Only the cut cells (0 < n < q^3) are computed, in blocks whose two live
    float64 blocks take at most budget bytes; the map occ stays where it is (at q 16 the full map has 4096 bytes per
    cell: 15.5 GB at N 156; PLAN 15 FQ). Each sum is exact in float64 (the positions are multiples of 1/(2q), at most
    q^3 terms), so the result does not depend on the blocks or on the device."""
    occ_t = occ if torch.is_tensor(occ) else torch.as_tensor(occ, device="cpu")
    dev = occ.device if torch.is_tensor(occ) else torch.empty(0).device     # as before: the default device for numpy
    nq = occ_t.shape[-1]
    shape = occ_t.shape[:-1]
    occ_t = occ_t.reshape(-1, nq)
    # the count of the body sub-cells in blocks of cells: torch.sum of a bool map makes an int64 copy of its input
    # (124 GB at q 16, N 156); integer sums are exact, so the blocks do not change the result (PLAN 16 Q)
    cb = max(1, int(budget) // (8 * nq))
    n = torch.cat([occ_t[i:i + cb].sum(-1) for i in range(0, occ_t.shape[0], cb)])
    rows = torch.nonzero((n > 0) & (n < nq)).reshape(-1)
    aloc = torch.stack(torch.meshgrid(*(torch.arange(q, dtype=torch.float64),) * 3, indexing="ij"), -1).reshape(-1, 3)
    pos = ((aloc + 0.5) / q - 0.5).to(dev)                           # sub-cell centres relative to the cell centre
    d = torch.zeros(occ_t.shape[0], 3, dtype=torch.float64, device=dev)
    chunk = max(1, int(budget) // (2 * 8 * nq))
    for i in range(0, len(rows), chunk):
        r = rows[i:i + chunk]
        o = occ_t[r].to(device=dev, dtype=torch.float64)
        d[r.to(dev)] = (o @ pos) / o.sum(-1, keepdim=True).clamp(min=1.0)
        del o
    return d.reshape(*shape, 3)


def centroid_field(Nhat, v, d, m, shape=None, h1=None):
    """First-moment part of the far field of the cut cells, with the rfft kernel Nhat of DemagField (list of lists
    [3][3]; shape: its FFT shape, (N, N, N) with PBC, (2N, 2N, 2N) for open boundaries):
        h_i = d_i . (D K * (v m))_i - sum_k (D_k K * (v d_k m))_i,
    D_k K = central difference of the kernel along k (i sin w_k in Fourier space). Arrays [N, N, N(, 3)].

    With PBC the FFT is circular: i sin w_k is then exactly the periodic central difference (f[x+1] - f[x-1]) / 2,
    and the convolution commutes with it. Thus h = sum_k d_k CD_k(K * v m) - K * sum_k CD_k(v d_k m): two
    convolutions (12 FFTs) instead of 24 (PROTOKOLL 3.36). Open boundaries use centroid_field_fourier.

    h1: K * (v m), if the caller has it (the stray field of DemagField divided by M_s; a uniform part does not
    matter, the differences remove it): then one convolution (6 FFTs; audit B1)."""
    N = v.shape[0]
    s = tuple(shape) if shape is not None else (N, N, N)
    if s != (N, N, N):
        return centroid_field_fourier(Nhat, v, d, m, shape=s)

    def conv(src):
        """K * src for the 3 components [N, N, N, 3]: one batched rfftn and one batched irfftn."""
        f = torch.fft.rfftn(src, s=s, dim=(0, 1, 2))
        return torch.fft.irfftn(torch.stack([sum(Nhat[b][a] * f[..., a] for a in range(3)) for b in range(3)], -1),
                                s=s, dim=(0, 1, 2))

    def cd(f, k):
        """Periodic central difference along axis k: (f[x+1] - f[x-1]) / 2."""
        return 0.5 * (torch.roll(f, -1, k) - torch.roll(f, 1, k))

    vm = v[..., None] * m
    if h1 is None:
        h1 = conv(vm)
    h = sum(d[..., k:k + 1] * cd(h1, k) for k in range(3))
    return h - conv(sum(cd(vm * d[..., k:k + 1], k) for k in range(3)))


@torch.compile
def _centroid_source(v, d, m):
    """sum_k CD_k(v d_k m), the source of the second convolution of centroid_field, component first [3, N, N, N]:
    one fused kernel (PLAN 15 F1; the periodic central differences become index shifts)."""
    vm = v[..., None] * m
    s = sum(0.5 * (torch.roll(vm * d[..., k:k + 1], -1, k) - torch.roll(vm * d[..., k:k + 1], 1, k)) for k in range(3))
    return s.permute(3, 0, 1, 2).contiguous()


@torch.compile
def _centroid_total(hp, Ms, hd, d, c0, c1, c2):
    """Ms h_pair + Ms (sum_k d_k CD_k(h1) - K * source) with h1 = hd / Ms, as h() of the pair form with the centroid
    term and B1 computes it: one fused kernel (PLAN 15 F1). c0, c1, c2: the components of K * source."""
    h1 = hd / Ms
    h = sum(d[..., k:k + 1] * (0.5 * (torch.roll(h1, -1, k) - torch.roll(h1, 1, k))) for k in range(3))
    return Ms * hp + Ms * (h - torch.stack([c0, c1, c2], -1))


def centroid_field_fourier(Nhat, v, d, m, shape=None):
    """centroid_field with the differences in Fourier space (24 FFTs); any FFT shape, also open boundaries."""
    N = v.shape[0]
    s = tuple(shape) if shape is not None else (N, N, N)
    dev, dt = Nhat[0][0].device, torch.real(Nhat[0][0]).dtype
    w = [2 * np.pi * torch.fft.fftfreq(s[0], dtype=dt, device=dev),
         2 * np.pi * torch.fft.fftfreq(s[1], dtype=dt, device=dev),
         2 * np.pi * torch.fft.rfftfreq(s[2], dtype=dt, device=dev)]
    S = [torch.sin(w[0]).reshape(-1, 1, 1), torch.sin(w[1]).reshape(1, -1, 1), torch.sin(w[2]).reshape(1, 1, -1)]
    vm = (v[..., None] * m)

    def conv(src):
        """K * src in Fourier space: list of the 3 field components (complex)."""
        f = [torch.fft.rfftn(src[..., a], s=s, dim=(0, 1, 2)) for a in range(3)]
        return [sum(Nhat[b][a] * f[a] for a in range(3)) for b in range(3)]

    def back(x):
        return torch.fft.irfftn(x, s=s, dim=(0, 1, 2))[:N, :N, :N]

    hv = conv(vm)
    h = torch.zeros_like(m)
    src = [0, 0, 0]
    for k in range(3):
        h += d[..., k:k + 1] * torch.stack([back(1j * S[k] * hv[b]) for b in range(3)], -1)
        hk = conv(vm * d[..., k:k + 1])
        src = [src[b] + 1j * S[k] * hk[b] for b in range(3)]   # source terms summed in Fourier space
    h -= torch.stack([back(src[b]) for b in range(3)], -1)
    return h


class NearCorrection:
    """The correction terms of the pairwise and the cell-local form for the V1 cells (n, occ) of a box of N^3 cells.

    periodic=True: the window wraps across the box edge (minimum image under PBC). periodic=False: offsets that leave
    the box are dropped (open boundaries; used by the unit tests).

    centroid=True: the far field also carries the first moment of the cut cells (each moment at the centroid d_i of
    the body part of its cell; centroid_field). The window then removes the coarse coupling
    v_j [K^c(r) + ((d_i - d_j) . D) K^c(r)] with the central difference D of the coarse kernel."""

    def __init__(self, n, occ, q, p, periodic=True, batch=2048, centroid=False, table_bytes=1e9):
        self.N = N = n.shape[0]
        self.q, self.p, self.periodic, self.centroid = q, p, periodic, centroid
        nq = q ** 3
        self.n = torch.as_tensor(n, dtype=torch.float64).reshape(-1)
        self.v = self.n / nq
        self.mag = self.n > 0
        self.full = self.n == nq
        mod = self.mag & ~self.full                                  # modified cells (0 < n < q^3)
        self.J = torch.nonzero(mod).reshape(-1)                      # flat indices of the modified cells
        self.jxyz = torch.stack([self.J // (N * N), (self.J // N) % N, self.J % N], dim=-1)
        r = torch.arange(-p, p + 1)
        D = torch.stack(torch.meshgrid(r, r, r, indexing="ij"), dim=-1).reshape(-1, 3)
        self.offsets = D                                             # [n_off, 3]
        Kc = kernel_table(p + 1)                                     # coarse kernel for |offset| <= p + 1
        i0, i1, i2 = D[:, 0] + p + 1, D[:, 1] + p + 1, D[:, 2] + p + 1
        self.Kc = Kc[i0, i1, i2]                                     # [n_off, 6]
        # central differences of the coarse kernel along x, y, z: [n_off, 3, 6]
        self.DKc = torch.stack([(Kc[i0 + 1, i1, i2] - Kc[i0 - 1, i1, i2]) / 2,
                                (Kc[i0, i1 + 1, i2] - Kc[i0, i1 - 1, i2]) / 2,
                                (Kc[i0, i1, i2 + 1] - Kc[i0, i1, i2 - 1]) / 2], dim=1)
        self.d = centroid_offsets(occ, q).reshape(-1, 3) if centroid else torch.zeros(N ** 3, 3, dtype=torch.float64)
        Rf = (p + 1) * q - 1                                         # largest fine offset between cells in the window
        Kf = kernel_table(Rf)
        aloc = torch.stack(torch.meshgrid(*(torch.arange(q),) * 3, indexing="ij"), dim=-1).reshape(-1, 3)  # [q^3, 3]
        # the rows of occ only for the modified cells (PLAN 15 FQ: no full map on the GPU); rows(idx) gives the float64
        # rows of any cells, with ones for a full cell and zeros for an empty cell (the values of the full map)
        occ_c = occ if torch.is_tensor(occ) else torch.as_tensor(occ, device="cpu")
        occ_J = occ_c.reshape(-1, nq)[self.J.to(occ_c.device)].to(self.J.device)   # [M, q^3] bool
        slot = torch.full((N ** 3,), -1, dtype=torch.int64)
        slot[self.J] = torch.arange(len(self.J))

        def rows(idx):
            s = slot[idx]
            o = occ_J[s.clamp(min=0)].to(torch.float64)
            return torch.where((s >= 0)[..., None], o, self.full[idx][..., None].to(torch.float64))
        # B[j, offset, 6] = L_{t, j}, target t = j + offset (zero if the target is not magnetic, or if the offset
        # leaves the box in the open case); computed in blocks of offsets (the table K^f(q offset + a - b) of all
        # offsets has 343 q^6 6 values: 4.3 GB at q 8) and in blocks of source cells (bounded memory for any q)
        M = len(self.J)
        self.B = torch.zeros(M, len(D), 6, dtype=torch.float64)
        per_off = nq * nq * 6 * 8
        n_chunk = max(1, min(len(D), int(table_bytes) // per_off))
        t_all, ok_all = self._targets(self.jxyz)                     # [M, n_off]
        for d0 in range(0, len(D), n_chunk):
            dc = slice(d0, min(d0 + n_chunk, len(D)))
            base = q * D[dc, None, None, :] + aloc[None, :, None, :] - aloc[None, None, :, :] + Rf   # [n_d, a, b, 3]
            Kab = Kf[base[..., 0], base[..., 1], base[..., 2]]       # [n_d, a, b, 6]: K^f(q offset + a - b)
            nb = max(1, min(batch, int(2e8) // ((dc.stop - dc.start) * nq * 6 * 8)))
            for s in range(0, M, nb):
                sl = slice(s, min(s + nb, M))
                S = torch.einsum("jb,dabc->jdac", occ_J[sl].to(torch.float64), Kab)  # sum over source b
                t, ok = t_all[sl, dc], ok_all[sl, dc]                # [b, n_d]
                L = torch.einsum("bda,bdac->bdc", rows(t), S)
                nt = self.n[t]
                L = torch.where((ok & (nt > 0))[..., None], L / nt.clamp(min=1.0)[..., None], torch.zeros_like(L))
                self.B[sl, dc] = L
            del Kab, base
        del occ_J, slot
        self._kloc = None

    def _targets(self, xyz):
        """Flat target index and validity for each source cell (rows of xyz) and each offset."""
        N = self.N
        T = xyz[:, None, :] + self.offsets[None, :, :]
        if self.periodic:
            ok = torch.ones(T.shape[:2], dtype=torch.bool)
            T = T % N
        else:
            ok = ((T >= 0) & (T < N)).all(-1)
            T = T.clamp(0, N - 1)
        return (T[..., 0] * N + T[..., 1]) * N + T[..., 2], ok

    def _coeffs(self, o, t, ok):
        """Coefficients of offset o: CA (source = modified cell j, target t) and CB (target = modified cell i,
        source t; only for unmodified magnetic sources)."""
        vJ = self.v[self.J]
        dJ = self.d[self.J]
        # coarse coupling of the far-field model for target t, source J (offset r = o):
        #     K(r) + ((d_t - d_J) . D) K(r)
        KA = self.Kc[o][None, :] + torch.einsum("bk,kc->bc", self.d[t] - dJ, self.DKc[o])
        CA = self.B[:, o] - vJ[:, None] * KA
        CA = torch.where((ok & self.mag[t])[:, None], CA, torch.zeros_like(CA))
        # target J, full source t (d_t = 0), offset -o: K(-r) + (d_J . D) K(-r) = K(r) - (d_J . D) K(r)
        KB = self.Kc[o][None, :] - torch.einsum("bk,kc->bc", dJ, self.DKc[o])
        CB = self.B[:, o] / vJ[:, None] - KB
        CB = torch.where((ok & self.full[t])[:, None], CB, torch.zeros_like(CB))
        return CA, CB

    def h_pair(self, m):
        """Correction of the pairwise form; m [N, N, N, 3] -> h [N, N, N, 3]."""
        N = self.N
        mf = m.reshape(-1, 3).to(torch.float64)
        h = torch.zeros_like(mf)
        t_all, ok_all = self._targets(self.jxyz)
        mJ = mf[self.J]
        for o in range(len(self.offsets)):
            t, ok = t_all[:, o], ok_all[:, o]
            CA, CB = self._coeffs(o, t, ok)
            h.index_add_(0, t, sym_matvec(CA, mJ))
            h.index_add_(0, self.J, sym_matvec(CB, mf[t]))
        return h.reshape(N, N, N, 3)

    def k_local(self):
        """Cell-local tensors [N^3, 6]: sum over the window of (L_ij - v_j K^c(i - j))."""
        if self._kloc is None:
            k = torch.zeros(self.N ** 3, 6, dtype=torch.float64)
            t_all, ok_all = self._targets(self.jxyz)
            for o in range(len(self.offsets)):
                t, ok = t_all[:, o], ok_all[:, o]
                CA, CB = self._coeffs(o, t, ok)
                k.index_add_(0, t, CA)
                k.index_add_(0, self.J, CB)
            self._kloc = k
        return self._kloc

    def h_local(self, m):
        """Correction of the cell-local form."""
        N = self.N
        return sym_matvec(self.k_local(), m.reshape(-1, 3).to(torch.float64)).reshape(N, N, N, 3)

    def energy(self, m, h):
        """Energy per box volume in units of K_d: -(1 / N^3) sum_i v_i m_i . h_i."""
        v = self.v.reshape(self.N, self.N, self.N, 1)
        return float(-(v * m.to(torch.float64) * h).sum() / self.N ** 3)


class DemagStaircaseCorrection(LinearFieldTerm):
    r"""
    Near-field correction of the stray field for cut surface cells (Donahue & McMichael 2007). Use it together with
    DemagField, with the material "Ms" set to v M_s (v: volume fraction of each cell from the same q sub-cells).

    :param n: number of body sub-cells of each cell, array [nx, nx, nx] (cubic mesh, 0 ... q^3)
    :param occ: body sub-cells of each cell, bool array [nx, nx, nx, q^3] (order of staircase_cells)
    :param q: sub-cells per axis and cell
    :param p: half width of the window in cells
    :param Ms: saturation magnetisation of the full material [A/m]
    :param form: "pair" (eq. 5) or "local" (eqs. 8 and 9)
    :param periodic: window across the box edge (PBC)
    :param centroid: far field also with the first moment of the cut cells (centroid_field); needs demag
    :param demag: the DemagField term of the far field (its kernel is used for the centroid term)
    """
    def __init__(self, n, occ, q, p, Ms, form="pair", periodic=True, centroid=False, demag=None):
        if form not in ("pair", "local"):
            raise ValueError("form must be 'pair' or 'local'")
        if centroid and demag is None:
            raise ValueError("centroid needs the DemagField term (demag)")
        self._n, self._occ, self._q, self._p = n, occ, q, p
        self._Ms, self._form, self._periodic = float(Ms), form, periodic
        self._centroid, self._demag = centroid, demag
        self._ready = False
        self._Nk = None

    def _setup(self, state):
        nc = NearCorrection(self._n, self._occ, self._q, self._p, periodic=self._periodic, centroid=self._centroid)
        dev, dt = state.device, state.dtype
        if self._centroid:
            self._v = torch.as_tensor(np.asarray(self._n) / self._q ** 3, device=dev, dtype=dt)
            self._d = nc.d.reshape(tuple(np.asarray(self._n).shape) + (3,)).to(device=dev, dtype=dt)   # computed once
        if self._form == "local":
            k = nc.k_local()
            idx = torch.nonzero((k != 0).any(-1)).reshape(-1)
            self._idx, self._k = idx.to(dev), k[idx].to(device=dev, dtype=dt)
        else:
            t_all, ok_all = nc._targets(nc.jxyz)
            tg, sr, cf = [], [], []
            for o in range(len(nc.offsets)):
                t, ok = t_all[:, o], ok_all[:, o]
                CA, CB = nc._coeffs(o, t, ok)
                a = (CA != 0).any(-1)
                b = (CB != 0).any(-1)
                tg += [t[a], nc.J[b]]
                sr += [nc.J[a], t[b]]
                cf += [CA[a], CB[b]]
            self._A = pair_blocks(torch.cat(tg), torch.cat(sr), torch.cat(cf), nc.N ** 3, dtype=dt).to(dev)
            del tg, sr, cf
        self._shape = tuple(state.mesh.n)
        self._ready = True
        del nc

    def h(self, state):
        if not self._ready:
            self._setup(state)
        m = state.m.reshape(-1, 3)
        if self._form == "local":
            h = torch.zeros_like(m)
            h[self._idx] = sym_matvec(self._k, m[self._idx])
        else:
            h = (self._A @ m.reshape(-1)).reshape(-1, 3)
        if self._centroid:
            dem = self._demag
            if not hasattr(dem, "_N"):
                dem._N = dem._init_N(state)
            # K * (v m) from the stray field of this m, if DemagField has it (audit B1): material Ms = v M_s, so
            # h_demag / M_s = K * (v m) + the field of the void moments (ms_void, about 1e-8) + the uniform block
            # term of true_pbc, which the central differences remove
            hd = dem.last_h(state) if hasattr(dem, "last_h") else None
            s = tuple(dem._shape(state))
            if hd is not None and s == tuple(self._shape):           # PBC with B1 (production): fused (PLAN 15 F1)
                return self._centroid_fused(h.reshape(self._shape + (3,)), hd, state.m, s)
            h = (self._Ms * h).reshape(self._shape + (3,))
            return h + self._Ms * centroid_field(dem._N, self._v, self._d, state.m, shape=s,
                                                 h1=None if hd is None else hd / self._Ms)
        return (self._Ms * h).reshape(self._shape + (3,))

    def _centroid_fused(self, hp, hd, m, s):
        """h() of the pair form with the centroid term and B1, as centroid_field with h1 = hd / Ms: the work outside
        the 6 FFTs in 3 fused kernels (_centroid_source, _kernel_product of DemagField, _centroid_total)."""
        if self._Nk is None:
            N = self._demag._N
            self._Nk = [K.unsqueeze(-1) for K in (N[0][0], N[0][1], N[0][2], N[1][1], N[1][2], N[2][2])]
        src = _centroid_source(self._v, self._d, m)
        f = [torch.view_as_real(torch.fft.rfftn(src[a], s=s, dim=(0, 1, 2))) for a in range(3)]
        c = [torch.fft.irfftn(torch.view_as_complex(k), s=s, dim=(0, 1, 2)) for k in _kernel_product(self._Nk, f)]
        return _centroid_total(hp, self._Ms, hd, self._d, *c)

    def E(self, state, domain = Ellipsis): # without torch.compile, as DemagField.E
        E = -0.5 * constants.mu_0 * state.material["Ms"] * state.m * self.h(state) * state.mesh.cell_volumes
        return E[domain].sum(dtype=torch.float64)


class DemagStaircaseFine(LinearFieldTerm):
    r"""
    Reference model of the staircase correction (studies/staircase_correction PROTOKOLL 3.18): the exact stray field
    of the fine staircase for the coarse magnetisation, without a window. Each body sub-cell takes m of its cell; the
    field of a cell is the mean of the fine field over its body sub-cells. This is the gradient of the stray-field
    energy of the fine body, so the term is variational. One FFT on the fine mesh (q^3 times the cells) per call.
    Use it instead of DemagField, with the material "Ms" set to v M_s (v: volume fraction from the same sub-cells).

    :param occ_fine: body of the fine mesh, bool array [q nx, q ny, q nz]
    :param q: sub-cells per axis and cell
    :param Ms: saturation magnetisation of the full material [A/m]
    :param demag_kwargs: arguments of the DemagField of the fine mesh (e.g. true_pbc, extrapolate)
    """
    def __init__(self, occ_fine, q, Ms, **demag_kwargs):
        from magnumnp.field_terms.demag import DemagField
        self._occ, self._q, self._Ms = np.asarray(occ_fine, dtype=bool), int(q), float(Ms)
        self._demag = DemagField(**demag_kwargs)
        self._fine = None

    def _setup(self, state):
        from magnumnp.common.mesh import Mesh
        from magnumnp.common.state import State
        q, n = self._q, tuple(state.mesh.n)
        if self._occ.shape != tuple(q * k for k in n):
            raise ValueError("occ_fine has the shape %s, not q times the mesh %s" % (self._occ.shape, n))
        fine = State(Mesh(tuple(q * k for k in n), tuple(d / q for d in state.mesh.dx), pbc=state.mesh.pbc))
        occ = torch.as_tensor(self._occ, device=state.device, dtype=state.dtype).reshape(n[0], q, n[1], q, n[2], q)
        fine.material["Ms"] = (self._Ms * occ.reshape(fine.mesh.n))[..., None]
        self._w = occ[..., None]
        self._cnt = occ.sum(dim=(1, 3, 5)).clamp(min=1.0)[..., None]
        self._fine = fine

    def h(self, state):
        if self._fine is None:
            self._setup(state)
        q, (nx, ny, nz) = self._q, state.mesh.n
        self._fine.m = state.m.repeat_interleave(q, 0).repeat_interleave(q, 1).repeat_interleave(q, 2)
        hf = self._demag.h(self._fine).reshape(nx, q, ny, q, nz, q, 3)
        return (hf * self._w).sum(dim=(1, 3, 5)) / self._cnt

    def E(self, state, domain = Ellipsis): # without torch.compile, as DemagField.E
        E = -0.5 * constants.mu_0 * state.material["Ms"] * state.m * self.h(state) * state.mesh.cell_volumes
        return E[domain].sum(dtype=torch.float64)
