"""
Near-field correction of the stray field for the staircase surface of voxel spheres, after Donahue & McMichael,
IEEE Trans. Magn. 43, 2878 (2007) (PROTOKOLL 3.1 ... 3.4, PLAN.md section 2).

V1 cells (cut cells): a fine mesh with q sub-cells per axis and coarse cell carries the fine body (fine cells whose
centre is inside a sphere). A coarse cell with at least one fine cell of the body is magnetic, with the volume
fraction v_i = n_i / q^3. The coarse FFT field uses the moments v_j m_j; the correction replaces, inside the window
W(i) (max |offset| <= p, minimum image), the coarse coupling v_j K^c(i - j) by the coupling L_ij of the fine cells:

    h_i = h_i(FFT, v_j m_j) + sum_{j in W(i)} (L_ij - v_j K^c(i - j)) m_j            (pairwise form, eq. 5)
    h_i = h_i(FFT, v_j m_j) + [sum_{j in W(i)} (L_ij - v_j K^c(i - j))] m_i          (cell-local form, eqs. 8, 9)

L_ij = (1 / n_i) sum_{a in i} sum_{b in j} K^f(a - b). Sign convention of magnum.np: H = K M (self term about -1/3).
Units: lengths in cells of the respective mesh (the kernel is scale-free), fields h = H / M_s.
"""
import numpy as np
import torch

from magnumnp.field_terms.demag import demag_f, demag_g

NEWELL_P = 20          # Newell radius of DemagField (in cell diagonals), the same on both meshes
SYM = ((0, 1, 2), (1, 3, 4), (2, 4, 5))   # 3x3 index -> component of (xx, xy, xz, yy, yz, zz)


def kernel_table(R, p=NEWELL_P, periodic_n=None):
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


def blocks(fine, q):
    """fine [qN]^3 (any dtype) -> [N, N, N, q^3]: the q^3 sub-cells of each coarse cell, sub-cell order (a, b, c)."""
    Nf = fine.shape[0]
    N = Nf // q
    return fine.reshape(N, q, N, q, N, q).transpose(0, 2, 4, 1, 3, 5).reshape(N, N, N, q ** 3)


def v1_cells(ids_fine, q):
    """V1 cells from the fine particle map ids_fine ([qN]^3, -1 = void, 0..3 = particle).
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


class NearCorrection:
    """The correction terms of the pairwise and the cell-local form for the V1 cells (n, occ) of a box of N^3 cells.

    periodic=True: the window wraps across the box edge (minimum image under PBC). periodic=False: offsets that leave
    the box are dropped (open boundaries; used by the unit tests)."""

    def __init__(self, n, occ, q, p, periodic=True, batch=2048):
        self.N = N = n.shape[0]
        self.q, self.p, self.periodic = q, p, periodic
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
        Kc = kernel_table(p)                                         # coarse kernel for |offset| <= p
        self.Kc = Kc[D[:, 0] + p, D[:, 1] + p, D[:, 2] + p]          # [n_off, 6]
        Rf = (p + 1) * q - 1                                         # largest fine offset between cells in the window
        Kf = kernel_table(Rf)
        aloc = torch.stack(torch.meshgrid(*(torch.arange(q),) * 3, indexing="ij"), dim=-1).reshape(-1, 3)  # [q^3, 3]
        occ_t = torch.as_tensor(occ.reshape(-1, nq))
        # unique occupancy patterns of the modified source cells
        pats, pat_idx = torch.unique(occ_t[self.J], dim=0, return_inverse=True)
        # S[pattern, offset, a, 6] = sum_{b in pattern} K^f(q * offset + a - b)
        S = torch.zeros(len(pats), len(D), nq, 6, dtype=torch.float64)
        base = q * D[:, None, None, :] + aloc[None, :, None, :] - aloc[None, None, :, :] + Rf   # [n_off, a, b, 3]
        Kab = Kf[base[..., 0], base[..., 1], base[..., 2]]           # [n_off, a, b, 6]
        for k in range(len(pats)):
            S[k] = Kab[:, :, pats[k], :].sum(dim=2)
        # B[j, offset, 6] = L_{t, j}, target t = j + offset (zero if the target is not magnetic, or if the offset
        # leaves the box in the open case)
        M = len(self.J)
        self.B = torch.zeros(M, len(D), 6, dtype=torch.float64)
        occ_f = occ_t.to(torch.float64)
        for s in range(0, M, batch):
            sl = slice(s, min(s + batch, M))
            t, ok = self._targets(self.jxyz[sl])                     # [b, n_off]
            o_t = occ_f[t]                                           # [b, n_off, q^3]
            L = torch.einsum("bda,bdac->bdc", o_t, S[pat_idx[sl]])
            nt = self.n[t]
            L = torch.where((ok & (nt > 0))[..., None], L / nt.clamp(min=1.0)[..., None], torch.zeros_like(L))
            self.B[sl] = L
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
        CA = self.B[:, o] - vJ[:, None] * self.Kc[o][None, :]
        CA = torch.where((ok & self.mag[t])[:, None], CA, torch.zeros_like(CA))
        CB = self.B[:, o] / vJ[:, None] - self.Kc[o][None, :]
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


class FFTField:
    """Stray field h = H / M_s of DemagField on one mesh, with the kernel set up once.
    setting: "open" (no PBC), or ("newell", k) / ("newell_ext", k) with k image shells and true_pbc (PROTOKOLL 7.105)."""

    def __init__(self, N, setting, p=NEWELL_P):
        from magnumnp import Mesh, State, DemagField
        if setting == "open":
            pbc, term = (0, 0, 0), DemagField(p=p)
        else:
            kind, k = setting
            pbc, term = (k, k, k), DemagField(p=p, true_pbc=True, extrapolate=(kind == "newell_ext"))
        self.N = N
        self.state = State(Mesh((N, N, N), (1e-9, 1e-9, 1e-9), pbc=pbc))
        self.term = term

    def h(self, ms, m):
        """ms [N, N, N] (moment weight per cell), m [N, N, N, 3] -> h [N, N, N, 3] (double)."""
        st = self.state
        st.material["Ms"] = torch.as_tensor(ms, dtype=st.dtype).reshape(self.N, self.N, self.N, 1)
        st.m = torch.as_tensor(m, dtype=st.dtype)
        with torch.no_grad():
            return self.term.h(st).to(torch.float64)


def fine_owner_m(m_coarse, q):
    """Fine magnetisation of the V1 cells: each fine cell takes m of its coarse cell (the caller masks the body)."""
    m = torch.as_tensor(m_coarse)
    return m.repeat_interleave(q, 0).repeat_interleave(q, 1).repeat_interleave(q, 2)


def coarse_mean(h_fine, occ_fine, q):
    """Mean of the fine field over the body fine cells of each coarse cell (0 where a cell has none)."""
    Nf = h_fine.shape[0]
    N = Nf // q
    w = torch.as_tensor(occ_fine, dtype=torch.float64).reshape(N, q, N, q, N, q)
    hf = h_fine.reshape(N, q, N, q, N, q, 3)
    s = (hf * w[..., None]).sum(dim=(1, 3, 5))
    n = w.sum(dim=(1, 3, 5))
    return s / n.clamp(min=1.0)[..., None]
