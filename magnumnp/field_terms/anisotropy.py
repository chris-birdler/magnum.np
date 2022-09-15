from magnumnp.common import timedmethod, constants
import torch
from torch import sin, cos

__all__ = ["UniaxialAnisotropyField", "CubicAnisotropyField"]

class UniaxialAnisotropyField(object):
    @timedmethod
    def h(self, state):
        h = 2. * state.material["Ku"] * state.material["Ku_axis"] / (constants.mu_0 * state.material["Ms"]) * torch.sum(state.material["Ku_axis"] * state.m, dim=3, keepdim=True)
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))


class CubicAnisotropyField(object):
    def _R(self, state):
        a = state.material["Kc_alpha"]
        b = state.material["Kc_beta"]
        g = state.material["Kc_gamma"]

        R = torch.stack([torch.concat([-sin(a)*sin(g) + cos(a)*cos(b)*cos(g), cos(a)*sin(g) + sin(a)*cos(b)*cos(g), -sin(b)*cos(g)], dim = -1),
                         torch.concat([-sin(a)*cos(g) - cos(a)*cos(b)*sin(g), cos(a)*cos(g) - sin(a)*cos(b)*sin(g),  sin(b)*sin(g)], dim = -1),
                         torch.concat([ cos(a)*sin(b)                       , sin(a)*sin(b)                       ,  cos(b)], dim = -1)], dim = -1)
        return R

    @timedmethod
    def h(self, state):
        R = self._R(state)
        mx, my, mz = torch.einsum('...a, ...ab-> ...b', state.m, R).unbind(dim=-1) # matmult

        h =  2. * state.material["Kc1"] * torch.stack([mx * (my**2 + mz**2), my*(mz**2 + mx**2), mz*(mx**2 + my**2)], dim = -1) + \
             2. * state.material["Kc2"] * torch.stack([mx * my**2. * mz**2., mx**2. * my * mz**2., mx**2. * my**2. * mz], dim = -1)
        h = torch.einsum('...a, ...ba-> ...b', h, R) # matmult transpose
        return torch.nan_to_num(-1./constants.mu_0/state.material["Ms"] * h, posinf=0, neginf=0)

    def E(self, state):
        R = self._R(state)
        mx, my, mz = torch.einsum('...a, ...ab-> ...b', state.m, R).unbind(dim=-1) # matmult

        return (state.material["Kc1"] * (mx**2 * my**2 + mx**2 * mz**2 + my**2 * mz**2).sum() +
                state.material["Kc2"] * (mx**2 * my**2 * mz**2).sum()) * state.mesh.cell_volume
