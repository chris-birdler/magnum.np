from magnumnp.common import constants
import torch

__all__ = ["ExternalField"]

class ExternalField(object):
    def __init__(self, state, h): # TODO: state should be removed
        self._h = state.Tensor(h)

    def h(self, state):
        return self._h(state.t)

    def __setattr__(self, name, value):
        if name == "h":
            self._h = self._state.Tensor(value)
        else:
            super().__setattr__(name, value)

    def E(self, state):
        return - constants.mu_0 * state.mesh.cell_volume \
               * torch.sum(state.material["Ms"] * state.m * self.h(state))
