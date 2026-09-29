"""
Additional field terms for the beta(B_peak) study.

StressAnisotropyField
    Magnetoelastic energy of an isotropic magnetostrictive material in a
    prescribed (frozen) stress field:

        e_me = -(3/2) * lambda_s * sum_ij sigma_ij m_i m_j        [J/m^3]
        h    =  3 * lambda_s / (mu0 * Ms) * (sigma . m)            [A/m]

    Isotropic lambda_s is an approximation (single crystals have
    lambda_100 != lambda_111). The hydrostatic part of sigma adds a constant
    (|m| = 1) and has no effect. The term is linear in m, so the energy of
    LinearFieldTerm (-1/2 mu0 Ms m.h) is correct.

SinusoidalDrive
    Homogeneous AC field  H(t) = H_amp * sin(2 pi f (t - t0)) * e.
    H_amp, f and t0 are plain Python floats that the driver changes between
    llg.step() calls (at zero crossings of the field). The time is taken
    from state.t (a device tensor), so there is no device synchronisation.
"""
import math
import torch
from magnumnp import constants, LinearFieldTerm, timedmethod

__all__ = ["StressAnisotropyField", "SinusoidalDrive"]


class StressAnisotropyField(LinearFieldTerm):
    def __init__(self, sigma, lambda_s, **kwargs):
        """
        :param sigma:    stress field, tensor [nx,ny,nz,6] in Pa (Voigt order
                         xx, yy, zz, yz, xz, xy); set it to 0 in the void
        :param lambda_s: saturation magnetostriction (float or [nx,ny,nz,1])
        """
        super().__init__(**kwargs)
        self._sigma = sigma
        self._lambda_s = lambda_s
        self._c = None

    def _init(self, state):
        Ms = state.material["Ms"]
        c = 3.0 * self._lambda_s * self._sigma.to(dtype=state.dtype) / (constants.mu_0 * Ms)
        self._c = c.nan_to_num(posinf=0., neginf=0.).contiguous()

    @timedmethod
    def h(self, state):
        if self._c is None:
            self._init(state)
        c = self._c
        m = state.m
        mx, my, mz = m[..., 0], m[..., 1], m[..., 2]
        # c: [..., 0]=xx 1=yy 2=zz 3=yz 4=xz 5=xy
        return torch.stack([c[..., 0]*mx + c[..., 5]*my + c[..., 4]*mz,
                            c[..., 5]*mx + c[..., 1]*my + c[..., 3]*mz,
                            c[..., 4]*mx + c[..., 3]*my + c[..., 2]*mz], dim=-1)


class SinusoidalDrive(object):
    def __init__(self, state, direction=(1.0, 0.0, 0.0)):
        e = torch.tensor(direction, dtype=state.dtype)
        self._e = (e / torch.linalg.norm(e)).reshape(1, 1, 1, 3)
        self._n = state.mesh.n
        self.H_amp = 0.0   # [A/m]
        self.freq = 1.0    # [Hz]
        self.t0 = 0.0      # [s] start of the current segment (phase 0)

    def value(self, t):
        """Scalar field value [A/m] at time t (Python float)."""
        return self.H_amp * math.sin(2.0 * math.pi * self.freq * (t - self.t0))

    @timedmethod
    def h(self, state):
        s = torch.sin(2.0 * math.pi * self.freq * (state.t - self.t0))
        return (self.H_amp * s * self._e).expand(self._n + (3,))

    def E(self, state, domain=Ellipsis):
        E = -constants.mu_0 * state.material["Ms"] * state.m * self.h(state) * state.mesh.cell_volumes
        return E[domain].sum()
