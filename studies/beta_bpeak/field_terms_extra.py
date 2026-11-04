"""
Additional field terms for the beta(B_peak) study.

SinusoidalDrive
    Homogeneous AC field  H(t) = H_amp * sin(2 pi f (t - t0)) * e.
    H_amp, f and t0 are plain Python floats that the driver changes between
    llg.step() calls (at zero crossings of the field). The time is taken
    from state.t (a device tensor), so there is no device synchronisation.

CutPairExchange
    Diagnostic exchange term of PLAN 4.18 (E1, --diag_k_cut): the exchange
    of magnumnp with the factor K on the coupling of each neighbour pair
    with at least one cut cell (0 < f < 1). Not a production model.

tie_thin_face, tie_ghost_mirror, TiedLLGSolver
    PLAN 4.20 E6 (--tie_cells): tied cut cells. A tied cell j has no degree of
    freedom: m_j = s_j / |s_j|, s_j = sum_k w_jk m_k over its masters k (free
    cells of the same particle). V1 thin_face ties each thin cell to a face
    neighbour, V2 ghost_mirror ties each cut cell with the centre outside its
    sphere to the mirror point of its centre. TiedLLGSolver integrates the LLG
    of the free cells with the field of the tied energy.
"""
import math
import numpy as np
import torch
from magnumnp import constants, timedmethod, normalize, LLGSolver
from magnumnp.field_terms.field_terms import LinearFieldTerm
from magnumnp.solvers.llg import _llg_rhs_sum

__all__ = ["SinusoidalDrive", "CutPairExchange", "TieMap", "tie_thin_face", "tie_ghost_mirror", "TiedLLGSolver"]


class SinusoidalDrive(object):
    def __init__(self, state, direction=(1.0, 0.0, 0.0)):
        e = torch.tensor(direction, dtype=state.dtype)
        self._e = (e / torch.linalg.norm(e)).reshape(1, 1, 1, 3)
        self._n = state.mesh.n
        self.H_amp = 0.0   # [A/m]
        self.freq = 1.0    # [Hz]
        self.t0 = 0.0      # [s] start of the current segment (phase 0)
        self.hold = None   # if set: constant field H_amp * hold (ring-down test), no time dependence

    def value(self, t):
        """Scalar field value [A/m] at time t (Python float)."""
        if self.hold is not None:
            return self.H_amp * self.hold
        return self.H_amp * math.sin(2.0 * math.pi * self.freq * (t - self.t0))

    @timedmethod
    def h(self, state):
        if self.hold is not None:
            return (self.H_amp * self.hold * self._e).expand(self._n + (3,))
        s = torch.sin(2.0 * math.pi * self.freq * (state.t - self.t0))
        return (self.H_amp * s * self._e).expand(self._n + (3,))

    def E(self, state, domain=Ellipsis):
        E = -constants.mu_0 * state.material["Ms"] * state.m * self.h(state) * state.mesh.cell_volumes
        return E[domain].sum(dtype=torch.float64)


class CutPairExchange(LinearFieldTerm):
    r"""ExchangeField of magnumnp (periodic branch, the same formula) with the coupling of a pair multiplied by
    1 + (K - 1) s_ij, s_ij = 1 if cell i or cell j is a cut cell (PLAN 4.18 E1). The factor is set on the pair
    coefficient before it is used for both cells, so the field stays the gradient of a sum over pairs.

    :param cut: bool array [nx, ny, nz], True at the cut cells (magnetic, 0 < f < 1)
    :param K: factor on the coupling of the pairs with a cut cell
    """
    parameters = ["A"]

    def __init__(self, cut, K, **kwargs):
        self._cut = cut
        self._K = float(K)
        self._fac = None
        super().__init__(**kwargs)

    def factors(self, device, dtype):
        """Per axis: the factor of the pair (i - 1, i) along the axis, stored at cell i (as A_avg of ExchangeField
        before its roll)."""
        if self._fac is None:
            c = torch.as_tensor(self._cut, device=device)
            self._fac = [(1.0 + (self._K - 1.0) * (c | torch.roll(c, 1, dims=ax)).to(dtype))[..., None]
                         for ax in range(3)]
        return self._fac

    def n_pairs_changed(self, A):
        """Number of pairs with a cut cell and a nonzero coupling (both cells carry A)."""
        c = torch.as_tensor(self._cut).cpu()            # the default device of magnumnp can be the GPU
        a = torch.as_tensor(A).detach().cpu().reshape(c.shape) > 0
        return int(sum(int(((c | torch.roll(c, 1, dims=ax)) & a & torch.roll(a, 1, dims=ax)).sum()) for ax in range(3)))

    @timedmethod
    def h(self, state):
        if any(int(p) == 0 for p in state.mesh.pbc):
            raise ValueError("CutPairExchange: only for a periodic mesh in all three axes")
        A = state.material[self.A]
        Ms = state.material["Ms"]
        m = state.m
        fac = self.factors(m.device, m.dtype)
        h = torch.zeros_like(m)
        for ax in range(3):
            d = state.mesh.dx_tensor[ax].reshape([-1 if k == ax else 1 for k in range(3)] + [1])
            A_next = torch.roll(A, +1, dims=ax)
            d_next = torch.roll(d, +1, dims=ax)
            A_avg = fac[ax] * 2. * A_next * A / (A_next * d + A * d_next)
            h += A_avg * (torch.roll(m, +1, dims=ax) - m) / d            # m_i-1 - m_i
            A_avg = torch.roll(A_avg, -1, dims=ax)
            h += A_avg * (torch.roll(m, -1, dims=ax) - m) / d            # m_i+1 - m_i
        h *= 2. / (constants.mu_0 * Ms)
        return h.nan_to_num_(posinf=0, neginf=0)

    def E(self, state, domain=Ellipsis):      # without torch.compile: the factors are made at the first call of h
        E = -0.5 * constants.mu_0 * state.material["Ms"] * state.m * self.h(state) * state.mesh.cell_volumes
        return E[domain].sum(dtype=torch.float64)


# --- PLAN 4.20 E6: tied cut cells ------------------------------------------------------------------------------------

_FACES = ((0, +1), (0, -1), (1, +1), (1, -1), (2, +1), (2, -1))     # (axis, side); the first largest face wins


class TieMap(object):
    """The tie of PLAN 4.20 E6 on a mesh with N^3 cells (flat indices i N^2 + j N + k):
    tied [n_t] the tied cells; src [n_t, K] their masters, w [n_t, K] the weights (padded with the first master and the
    weight 0); info: the counts for the record of the run."""

    def __init__(self, kind, N, tied, src, w, info):
        self.kind, self.N = kind, N
        self.tied = np.asarray(tied, np.int64)
        self.src = np.asarray(src, np.int64).reshape(len(self.tied), -1)
        self.w = np.asarray(w, np.float64).reshape(len(self.tied), -1)
        self.info = dict(info, kind=kind, n_tied=int(len(self.tied)))
        if np.isin(self.src[self.w > 0], self.tied).any():
            raise ValueError("tie %s: a master is a tied cell" % kind)

    def master_table(self):
        """(masters [n_m], rows [n_m, L], weights [n_m, L]): for each master the rows of the tied cells that use it,
        padded with the row n_t (a zero row) and the weight 0. Sorted by master and row, so that the sums over a
        master are a fixed gather (no scatter with atomics: the same result on each run)."""
        j, k = np.nonzero(self.w > 0)
        mk = self.src[j, k]
        order = np.lexsort((j, mk))
        j, k, mk = j[order], k[order], mk[order]
        masters, start, cnt = np.unique(mk, return_index=True, return_counts=True)
        L = int(cnt.max()) if len(cnt) else 1
        rows = np.full((len(masters), L), len(self.tied), np.int64)
        wts = np.zeros((len(masters), L), np.float64)
        pos = np.arange(len(mk)) - np.repeat(start, cnt)
        g = np.repeat(np.arange(len(masters)), cnt)
        rows[g, pos] = j
        wts[g, pos] = self.w[j, k]
        return masters, rows, wts


def _cut_and_thin(ids, frac, dx_lex):
    cut = (ids >= 0) & (frac < 1.0)
    return cut, cut & (frac < 1.0 / dx_lex)


def _nearest_cell(point, ok, rmax=6):
    """The cell with ok True whose centre is nearest to point (centre-grid units: the centre of cell (i, j, k) is at
    (i, j, k); periodic), the first in the order of the window offsets on equal distance; None within rmax cells."""
    N = ok.shape[0]
    pt = np.asarray(point, np.float64)
    c0 = np.round(pt).astype(np.int64)
    rng = np.arange(-rmax, rmax + 1)
    off = np.stack(np.meshgrid(rng, rng, rng, indexing="ij"), -1).reshape(-1, 3)
    d2 = ((c0[None, :] + off - pt[None, :]) ** 2).sum(1)
    c = (c0[None, :] + off) % N
    sel = ok[c[:, 0], c[:, 1], c[:, 2]]
    if not sel.any():
        return None
    order = np.lexsort((np.arange(len(off)), d2))
    best = order[np.argmax(sel[order])]
    return tuple(int(v) for v in c[best])


def tie_thin_face(ids, frac, dx_lex, occ, q, max_chain=3):
    """V1 thin_face (PLAN 4.20): each thin cell (cut cell with f below l_ex/dx = 1/dx_lex) is tied to the face neighbour
    of the same particle with the largest shared occupied face: the number of pairs of occupied sub-cells across the
    face (occ [N, N, N, q^3], sub-cell order (a, b, c) of geometry.fcc_subcells). If that neighbour is thin, the tie
    follows its master. A thin cell without an occupied shared face, a cycle, or a chain longer than max_chain cells
    takes the nearest free cell of the same particle. Weight 1."""
    N = ids.shape[0]
    cut, thin = _cut_and_thin(ids, frac, dx_lex)
    free = (ids >= 0) & ~thin
    th = np.argwhere(thin)
    o = occ[th[:, 0], th[:, 1], th[:, 2]].reshape(-1, q, q, q)
    best_n = np.zeros(len(th), np.int64)
    best_c = -np.ones((len(th), 3), np.int64)
    for ax, side in _FACES:
        nb = th.copy()
        nb[:, ax] = (nb[:, ax] + side) % N
        on = occ[nb[:, 0], nb[:, 1], nb[:, 2]].reshape(-1, q, q, q)
        same = ids[nb[:, 0], nb[:, 1], nb[:, 2]] == ids[th[:, 0], th[:, 1], th[:, 2]]
        own = np.take(o, q - 1 if side > 0 else 0, axis=ax + 1)
        oth = np.take(on, 0 if side > 0 else q - 1, axis=ax + 1)
        n = (own & oth).sum(axis=(1, 2)) * same
        better = n > best_n
        best_n[better] = n[better]
        best_c[better] = nb[better]
    flat = lambda c: (c[..., 0] * N + c[..., 1]) * N + c[..., 2]
    free_p = {p: free & (ids == p) for p in np.unique(ids[thin])}
    nxt = dict(zip(flat(th).tolist(), [int(v) if n > 0 else -1 for v, n in zip(flat(best_c).tolist(), best_n)]))
    thin_set = set(nxt)
    src, n_fall_face, n_fall_chain, longest = [], 0, 0, 0
    for j in flat(th).tolist():
        k, chain, seen = nxt[j], 1, {j}
        while k >= 0 and k in thin_set and chain <= max_chain and k not in seen:
            seen.add(k)
            k, chain = nxt[k], chain + 1
        if nxt[j] < 0 or k < 0:
            fb, n_fall_face = True, n_fall_face + 1
        elif k in thin_set or chain > max_chain:
            fb, n_fall_chain = True, n_fall_chain + 1
        else:
            fb = False
            longest = max(longest, chain)
        if fb:
            ci = (j // (N * N), (j // N) % N, j % N)
            c = _nearest_cell(ci, free_p[ids[ci]])
            if c is None:
                raise ValueError("tie thin_face: no free cell of the particle near cell %d" % j)
            k = (c[0] * N + c[1]) * N + c[2]
        src.append(k)
    tied = flat(th)
    info = {"n_cut": int(cut.sum()), "n_thin": int(len(th)), "sum_f_tied": float(frac[thin].astype(np.float64).sum()),
            "longest_chain": int(longest), "n_fallback_no_face": int(n_fall_face),
            "n_fallback_chain": int(n_fall_chain), "q": int(q), "thin_below_f": 1.0 / dx_lex}
    return TieMap("thin_face", N, tied, np.array(src, np.int64)[:, None], np.ones((len(tied), 1)), info)


def tie_ghost_mirror(ids, frac, N, d, dx):
    """V2 ghost_mirror (PLAN 4.20; the idea from Garcia-Cervera, Gimbutas & E, J. Comput. Phys. 184, 37 (2003), Sect.
    3.1): each cut cell whose centre lies outside its sphere (minimum image, the arithmetic of geometry.fcc_box) is tied
    to the mirror point of its centre across the sphere surface along the radial normal (r -> 2R - r): trilinear
    weights of the 8 cell centres around the mirror point, only cells of the same particle with the centre inside;
    the weights renormalised; without such a cell the nearest one (centre distance to the mirror point)."""
    from geometry import FCC_SITES
    a = N * dx
    R = 0.5 * d
    cut = (ids >= 0) & (frac < 1.0)
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rr = np.full((N, N, N), np.inf)
    for p, s in enumerate(FCC_SITES * a):
        sel = ids == p
        ddx = X[sel] - s[0]; ddx -= a * np.round(ddx / a)
        ddy = Y[sel] - s[1]; ddy -= a * np.round(ddy / a)
        ddz = Z[sel] - s[2]; ddz -= a * np.round(ddz / a)
        rr[sel] = ddx**2 + ddy**2 + ddz**2
    inside = (ids >= 0) & (rr <= R**2)
    tc = np.argwhere(cut & ~inside)
    sites = FCC_SITES * a
    src = np.zeros((len(tc), 8), np.int64)
    w = np.zeros((len(tc), 8))
    corners = np.stack(np.meshgrid((0, 1), (0, 1), (0, 1), indexing="ij"), -1).reshape(-1, 3)
    n_near, w_min_sum = 0, 1.0
    inside_p = {p: inside & (ids == p) for p in range(len(FCC_SITES))}
    for n, (i, j, k) in enumerate(tc):
        p = ids[i, j, k]
        c = np.array([X[i, j, k], Y[i, j, k], Z[i, j, k]]) - sites[p]
        c -= a * np.round(c / a)
        r = math.sqrt(float((c**2).sum()))
        mp = (sites[p] + c * (2.0 * R - r) / r) / dx - 0.5            # mirror point in units of the centre grid
        i0 = np.floor(mp).astype(np.int64)
        t = mp - i0
        cc = (i0[None, :] + corners) % N
        wt = np.prod(np.where(corners == 1, t[None, :], 1.0 - t[None, :]), axis=1)
        ok = inside[cc[:, 0], cc[:, 1], cc[:, 2]] & (ids[cc[:, 0], cc[:, 1], cc[:, 2]] == p) & (wt > 0)
        flat = (cc[:, 0] * N + cc[:, 1]) * N + cc[:, 2]
        if ok.any():
            w_min_sum = min(w_min_sum, float(wt[ok].sum()))
            m = int(ok.sum())
            src[n, :m] = flat[ok]
            w[n, :m] = wt[ok] / wt[ok].sum()
            src[n, m:] = flat[ok][0]
        else:
            nc = _nearest_cell(mp, inside_p[p])
            if nc is None:
                raise ValueError("tie ghost_mirror: no inside cell of particle %d near the mirror point" % p)
            n_near += 1
            src[n, :] = (nc[0] * N + nc[1]) * N + nc[2]
            w[n, 0] = 1.0
    info = {"n_cut": int(cut.sum()), "n_centre_outside": int(len(tc)),
            "sum_f_tied": float(frac[cut & ~inside].astype(np.float64).sum()), "n_fallback_nearest": int(n_near),
            "smallest_weight_sum_before_norm": w_min_sum}
    return TieMap("ghost_mirror", N, (tc[:, 0] * N + tc[:, 1]) * N + tc[:, 2], src, w, info)


class TiedLLGSolver(LLGSolver):
    r"""LLGSolver of magnumnp with tied cells (PLAN 4.20 E6). dm ties its input x (m_j = s_j / |s_j|,
    s_j = \sum_k w_jk \hat x_k with the normalised masters), evaluates all field terms on the tied state and replaces
    the field of each master k by

        h_k^* = (\mu_k h_k + \sum_j w_jk \mu_j P_j h_j / |s_j|) / D_k,   P_j = I - m_j m_j^T,   D_k = \mu_k + \sum_j w_jk \mu_j

    with \mu = M_s f (state.material["Ms"]). The numerator is the chain rule, -(1/\mu_0 V) dE/dm_k of the tied energy;
    D_k is the moment of the group. dm of a tied cell is 0. step: the step of the ODE solver, normalise, tie (the state
    is then a fixed point of the tie). The relaxation loop of run_loops.py calls the ODE solver directly and ties with
    tie_state at the end of each chunk.

    :param terms: the field terms (as LLGSolver)
    :param tie: TieMap
    :param Ms: state.material["Ms"] (mesh shape [N, N, N, 1]); its device and dtype are used
    """

    def __init__(self, terms, tie, Ms, **kwargs):
        super().__init__(terms, **kwargs)
        dev, dt = Ms.device, Ms.dtype
        mu = Ms.detach().reshape(-1).to(torch.float64).cpu().numpy()
        masters, rows, wts = tie.master_table()
        mu_t = mu[tie.tied]
        mu_pad = np.append(mu_t, 0.0)
        D = mu[masters] + (wts * mu_pad[rows]).sum(1)
        weight = mu.copy()
        weight[masters] = D
        weight[tie.tied] = 0.0
        L = lambda a: torch.as_tensor(a, dtype=torch.int64, device=dev)
        F = lambda a: torch.as_tensor(a, dtype=dt, device=dev)
        self.tie_map = tie
        self._tied, self._src, self._w = L(tie.tied), L(tie.src), F(tie.w)[..., None]
        self._masters, self._rows, self._wts = L(masters), L(rows), F(wts)[..., None]
        self._mu_t, self._mu_m, self._D = F(mu_t)[:, None], F(mu[masters])[:, None], F(D)[:, None]
        self.p_weight = F(weight).reshape(Ms.shape)          # M_s of the free cells, D_k of the masters, 0 if tied

    def tie(self, x):
        """(x with the tied cells set, out of place; m_j [n_t, 3]; |s_j| [n_t, 1])"""
        xf = x.reshape(-1, 3)
        mk = xf[self._src]
        mk = mk / torch.linalg.norm(mk, dim=-1, keepdim=True)
        sj = (self._w * mk).sum(1)
        ns = torch.linalg.norm(sj, dim=-1, keepdim=True)
        mj = sj / ns
        return xf.index_copy(0, self._tied, mj).reshape(x.shape), mj, ns

    def tie_state(self, state):
        state.m = self.tie(state.m)[0]

    def field(self, state, mj, ns):
        """h^* of the masters, h of the other free cells, 0 at the tied cells (state.m tied)."""
        hs = [term.h(state) for term in self._terms]
        h = hs[0]
        for f in hs[1:]:
            h = h + f
        hf = h.reshape(-1, 3)
        ht = hf[self._tied]
        g = self._mu_t * (ht - mj * (mj * ht).sum(-1, keepdim=True)) / ns
        g = torch.cat([g, g.new_zeros(1, 3)])
        G = (self._wts * g[self._rows]).sum(1)
        hm = (self._mu_m * hf[self._masters] + G) / self._D
        return hf.index_copy(0, self._masters, hm).index_fill(0, self._tied, 0.0).reshape(h.shape)

    def field_p(self, state):
        """(h^*, the moment weight) for the dissipation p = p_coef weight |m x h^*|^2 of run_loops.diagnostics; ties
        state.m first (a tied state does not change)."""
        state.m, mj, ns = self.tie(state.m)
        return self.field(state, mj, ns), self.p_weight

    def dm(self, t, x, state, alpha=None):
        state.t = t
        state.m, mj, ns = self.tie(x)
        alpha = alpha or state.material["alpha"]
        if isinstance(alpha, torch.Tensor) and type(alpha) is not torch.Tensor:
            alpha = alpha.as_subclass(torch.Tensor)
        d = _llg_rhs_sum(state.m, [self.field(state, mj, ns)], alpha, self._no_precession)
        return d.reshape(-1, 3).index_fill(0, self._tied, 0.0).reshape(x.shape)

    @timedmethod
    def step(self, state, dt, **kwargs):
        state.t, state.m = self._solver.step(state.t, state.m, dt, state=state, **kwargs)
        normalize(state.m)
        self.tie_state(state)

    def relax(self, *args, **kwargs):
        raise NotImplementedError("TiedLLGSolver: use the converge relaxation of run_loops.py")
