from magnumnp.common import timedmethod, constants
import torch

__all__ = ["AnisotropyField"]

class AnisotropyField(object):
    def __init__(self, state):
        self._state = state

        self._mesh = state._mesh
        self._Ms = state._material["Ms"]
        self._K = state._material["K"]
        self._K_axis = state.function(state._mesh.n + (3,))
        self._K_axis[:,:,:,:] = torch.DoubleTensor(material["K_axis"])

        # initialize scratch space
        self._h = torch.zeros(mesh.n + (3,), dtype=torch.float64, device = cuda)

    @timedmethod
    def h(self, t, m):
        self._h[:,:,:,:] = 2. * self._K * self._K_axis / (constants.mu_0 * self._Ms) * torch.sum(self._K_axis * m, dim=3).unsqueeze(-1)
        self._h = torch.nan_to_num(self._h, posinf=0, neginf=0)
        return self._h

    def E(self, t, m):
        return - 0.5 * constants.mu_0 * self._mesh.cell_volume \
               * torch.sum(self._Ms * m * self.h(t, m))
