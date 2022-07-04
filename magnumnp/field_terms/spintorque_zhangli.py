from magnumnp.common import timedmethod, constants
import torch

__all__ = ["SpinTorqueZhangLi"]

class SpinTorqueZhangLi(object):
    def __init__(self, state):
        self._state = state
        self._mesh = state._mesh
        self._b = self._state._material["b"]
        self._xi = self._state._material["xi"]
        self.j = self._state.j

    @timedmethod
    def h(self, t, m):
        dim = [i for i in range(3) if self._mesh.n[i] > 1]
        dx = [self._mesh.dx[i] for i in range(3) if self._mesh.n[i] > 1]
        jgradm = (torch.stack(torch.gradient(m, spacing=dx, dim=dim), dim=-1)*self.j[dim]).sum(axis=-1)

        return self._b / constants.gamma * (torch.cross(m, jgradm) + self._xi * jgradm)

    def E(self, t, m):
        raise NotImplemented()
