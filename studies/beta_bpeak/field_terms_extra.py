"""
Additional field terms for the beta(B_peak) study.

SinusoidalDrive
    Homogeneous AC field  H(t) = H_amp * sin(2 pi f (t - t0)) * e.
    H_amp, f and t0 are plain Python floats that the driver changes between
    llg.step() calls (at zero crossings of the field). The time is taken
    from state.t (a device tensor), so there is no device synchronisation.
"""
import math
import torch
from magnumnp import constants, timedmethod

__all__ = ["SinusoidalDrive"]


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
        return E[domain].sum()
