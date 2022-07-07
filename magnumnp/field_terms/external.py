from magnumnp.common import constants
import torch

__all__ = ["ExternalField"]

class ExternalField(object):
    def __init__(self, state, h):
        self._state = state
        self._mesh = state._mesh
        self._Ms = state._material["Ms"]
        self._h = state.zeros(self._mesh.n + (3,))
        if isinstance(h, list):
            self._h[:,:,:,:] = state.tensor(h)
        elif isinstance(h, torch.Tensor):
            self._h[:,:,:,:] = h
        elif callable(h):
            self._lambda_h = h
        else:
            raise TypeError("h needs to be 'list', 'torch.Tensor', or 'function'!")

    def h(self, t, m):
        if hasattr(self, "_lambda_h"):
            self._h[:,:,:,:] = self._state.tensor(self._lambda_h(t))
        return self._h

    def E(self, t, m):
        return - constants.mu_0 * self._mesh.cell_volume \
               * torch.sum(self._Ms * m * self.h(t, m))
