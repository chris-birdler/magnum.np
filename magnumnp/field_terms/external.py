from magnumnp.common import constants
import torch

__all__ = ["ExternalField"]

class ExternalField(object):
    def __init__(self, state, h):
        self._state = state
        self._mesh = state._mesh
        self._Ms = state._material["Ms"]
        self._h = state.zeros(self._mesh.n + (3,))
        self._h[:,:,:,:] = torch.DoubleTensor(h)

    def h(self, t, m):
        return self._h

    def E(self, t, m):
        return - constants.mu_0 * self._mesh.cell_volume \
               * torch.sum(self._Ms * m * self.h(t, m))
