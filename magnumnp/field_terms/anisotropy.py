from magnumnp.common import timedmethod, constants
import torch

__all__ = ["AnisotropyField"]

class AnisotropyField(object):
    def __init__(self, state):
        self._state = state

        self._mesh = state._mesh
        self._Ms = state._material["Ms"]
        self._K = state._material["K"]
        self._K_axis = state._zeros(state._mesh.n + (3,)) #TODO: fix this and add unit tests for anisotropy field
        self._K_axis[:,:,:,:] = state.tensor(state._material["K_axis"])

    @timedmethod
    def h(self, t, m):
        h = 2. * self._K * self._K_axis / (constants.mu_0 * self._Ms) * torch.sum(self._K_axis * m, dim=3, keepdim=True)
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, t, m):
        return -0.5 * constants.mu_0 * self._mesh.cell_volume * torch.sum(self._Ms * m * self.h(t, m))
