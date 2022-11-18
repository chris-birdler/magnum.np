from magnumnp.common import timedmethod, constants
import torch
from .field_terms import FieldTerm, LinearFieldTerm
from torch import sin, cos

__all__ = ["UniaxialAnisotropyField", "CubicAnisotropyField"]

class UniaxialAnisotropyField(LinearFieldTerm):
    parameters = ["Ku", "Ku_axis"]

    @timedmethod
    def h(self, state):
        Ku = state.material[self.Ku]
        Ku_axis = state.material[self.Ku_axis]

        h = 2. * Ku * Ku_axis / (constants.mu_0 * state.material["Ms"]) * torch.sum(Ku_axis * state.m, dim=3, keepdim=True)
        return torch.nan_to_num(h, posinf=0, neginf=0)


class CubicAnisotropyField(FieldTerm):
    parameters = ["Kc1", "Kc2", "Kc_alpha", "Kc_beta", "Kc_gamma"]

    def _R(self, state):
        a = state.material[self.Kc_alpha]
        b = state.material[self.Kc_beta]
        g = state.material[self.Kc_gamma]

        R = torch.stack([torch.concat([-sin(a)*sin(g) + cos(a)*cos(b)*cos(g), cos(a)*sin(g) + sin(a)*cos(b)*cos(g), -sin(b)*cos(g)], dim = -1),
                         torch.concat([-sin(a)*cos(g) - cos(a)*cos(b)*sin(g), cos(a)*cos(g) - sin(a)*cos(b)*sin(g),  sin(b)*sin(g)], dim = -1),
                         torch.concat([ cos(a)*sin(b)                       , sin(a)*sin(b)                       ,  cos(b)], dim = -1)], dim = -1)
        return R

    @timedmethod
    def h(self, state):
        Kc1 = state.material[self.Kc1]
        Kc2 = state.material[self.Kc2]

        R = self._R(state)
        mx, my, mz = torch.einsum('...a, ...ab-> ...b', state.m, R).unbind(dim=-1) # matmult

        h =  2. * Kc1 * torch.stack([mx * (my**2 + mz**2), my*(mz**2 + mx**2), mz*(mx**2 + my**2)], dim = -1) + \
             2. * Kc2 * torch.stack([mx * my**2. * mz**2., mx**2. * my * mz**2., mx**2. * my**2. * mz], dim = -1)
        h = torch.einsum('...a, ...ba-> ...b', h, R) # matmult transpose
        return torch.nan_to_num(-1./constants.mu_0/state.material["Ms"] * h, posinf=0, neginf=0)

    def E(self, state):
        R = self._R(state)
        mx, my, mz = torch.einsum('...a, ...ab-> ...b', state.m, R).unbind(dim=-1) # matmult

        return (state.material[self.Kc1] * (mx**2 * my**2 + mx**2 * mz**2 + my**2 * mz**2).sum() +
                state.material[self.Kc2] * (mx**2 * my**2 * mz**2).sum()) * state.mesh.cell_volume
