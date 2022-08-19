from magnumnp.common import timedmethod, constants
import torch
from torch import sin, cos

__all__ = ["UniaxialAnisotropyField", "CubicAnisotropyField"]

class UniaxialAnisotropyField(object):
    def __init__(self, state):
        self._state = state

        self._mesh = state._mesh
        self._Ms = state._material["Ms"]
        self._Ku = state._material["Ku"]
        self._Ku_axis = state.Tensor(state._material["Ku_axis"])

    @timedmethod
    def h(self, t, m):
        h = 2. * self._Ku * self._Ku_axis / (constants.mu_0 * self._Ms) * torch.sum(self._Ku_axis * m, dim=3, keepdim=True)
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, t, m):
        return -0.5 * constants.mu_0 * self._mesh.cell_volume * torch.sum(self._Ms * m * self.h(t, m))


class CubicAnisotropyField(object):
    def __init__(self, state):
        self._state = state
        self._mesh = state._mesh
        self._Ms = state._material["Ms"]
        self._Kc1 = state._material["Kc1"]
        self._Kc2 = state._material["Kc2"]
        self._K_alpha = state._material["Kc_alpha"]
        self._K_beta  = state._material["Kc_beta"]
        self._K_gamma = state._material["Kc_gamma"]

    def _R(self, state):
        a = self._K_alpha
        b = self._K_beta
        g = self._K_gamma

        return torch.stack([torch.concat([-sin(a)*sin(g) + cos(a)*cos(b)*cos(g), cos(a)*sin(g) + sin(a)*cos(b)*cos(g), -sin(b)*cos(g)], dim = -1),
                            torch.concat([-sin(a)*cos(g) - cos(a)*cos(b)*sin(g), cos(a)*cos(g) - sin(a)*cos(b)*sin(g),  sin(b)*sin(g)], dim = -1),
                            torch.concat([ cos(a)*sin(b)                       , sin(a)*sin(b)                       ,  cos(b)], dim = -1)], dim = -1)

    @timedmethod
    def h(self, t, m):
        R = self._R(self._state)
        mx, my, mz = torch.einsum('...a, ...ab-> ...b', m, R).unbind(dim=-1) # matmult

        h =  2. * self._Kc1 * torch.stack([mx * (my**2 + mz**2), my*(mz**2 + mx**2), mz*(mx**2 + my**2)], dim = -1) + \
             2. * self._Kc2 * torch.stack([mx * my**2. * mz**2., mx**2. * my * mz**2., mx**2. * my**2. * mz], dim = -1)
        h = torch.einsum('...a, ...ba-> ...b', h, R) # matmult transpose
        return torch.nan_to_num(-1./constants.mu_0/self._Ms * h, posinf=0, neginf=0)

    def E(self, t, m):
        R = self._R(self._state)
        mx, my, mz = torch.einsum('...a, ...ab-> ...b', m, R).unbind(dim=-1) # matmult

        return (self._Kc1 * (mx**2 * my**2 + mx**2 * mz**2 + my**2 * mz**2).sum() +
                self._Kc2 * (mx**2 * my**2 * mz**2).sum()) * self._mesh.cell_volume
