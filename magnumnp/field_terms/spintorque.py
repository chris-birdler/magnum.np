from magnumnp.common import timedmethod, constants
import torch

__all__ = ["SpinOrbitTorque", "SpinTorqueZhangLi"]

class SpinOrbitTorque(object):
    @timedmethod
    def h(self, state):
        p = state.material["p"]
        if len(p.shape) == 1: # fixes API change in torch.linalg.cross! Maybe better to expand all parameters in state.material
            p = p.reshape((1,1,1,-1))
        h = state.material["eta_damp"] * torch.cross(state.m, p) + state.material["eta_field"] * p
        h *= -state.material["je"] * constants.hbar / (2. * constants.e * state.material["Ms"] * constants.mu_0 * state.material["d"])
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, state):
        raise NotImplemented()


class SpinTorqueZhangLi(object):
    @timedmethod
    def h(self, state):
        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        dx = [state.mesh.dx[i] for i in range(3) if state.mesh.n[i] > 1]

        j = state.j(state.t)
        jgradm = torch.einsum('...a,...ba-> ...b', j[...,dim], torch.stack(torch.gradient(state.m, spacing=dx, dim=dim), dim=-1)) # matmult

        return state.material["b"] / constants.gamma * (torch.cross(state.m, jgradm) + state.material["xi"] * jgradm)

    def E(self, state):
        raise NotImplemented()
