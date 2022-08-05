from magnumnp.common import constants
import torch

__all__ = ["ExternalField"]

class ExternalField(object):
    def __init__(self, state, h):
        self._state = state
        self._mesh = state._mesh
        self._h = state.Tensor(h)

    def h(self, t, m):
        return self._h(t)

    def __setattr__(self, name, value):
        if name == "h":
            self._h = self._state.Tensor(value)
        else:
            super().__setattr__(name, value)

    def E(self, t, m):
        return - constants.mu_0 * self._mesh.cell_volume \
               * torch.sum(self._state._material["Ms"] * m * self.h(t, m))
