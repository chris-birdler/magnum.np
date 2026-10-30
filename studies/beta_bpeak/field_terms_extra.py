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
"""
import math
import torch
from magnumnp import constants, timedmethod
from magnumnp.field_terms.field_terms import LinearFieldTerm

__all__ = ["SinusoidalDrive", "CutPairExchange"]


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
        c = torch.as_tensor(self._cut)
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
