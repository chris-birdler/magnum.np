from magnumnp.common import timedmethod, constants
import torch

__all__ = ["ExternalField"]

class ExternalField(object):
    def __init__(self, h):
        self._h = h

    @timedmethod
    def h(self, state):
        return state.Tensor(self._h)(state.t)

    def __setattr__(self, name, value):
        if name == "h":
            self._h = value
        else:
            super().__setattr__(name, value)

    def E(self, state):
        return - constants.mu_0 * state.mesh.cell_volume \
               * torch.sum(state.material["Ms"] * state.m * self.h(state))
