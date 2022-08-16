from magnumnp.common import timedmethod, constants
import torch

__all__ = ["SpinOrbitTorque", "SpinTorqueZhangLi"]

class SpinOrbitTorque(object):
    def __init__(self, state):
        self._state = state
        self._mesh = self._state._mesh
        self.eta_damp = self._state._material["eta_damp"]
        self.eta_field = self._state._material["eta_field"]
        self.Ms = self._state._material["Ms"]
        self.je = self._state._material["je"]
        self.d = self._state._material["d"] #thickness of FM
        self._p = state.Tensor(state._material["p"])
         
    @timedmethod
    def h(self, t, m):
        p = self._p(t)
        h = self.eta_damp * torch.cross(m, p) + self.eta_field * p
        h *= -self.je *constants.hbar / (2. * constants.e * self.Ms * constants.mu_0 * self.d) 
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, t, m):
        raise NotImplemented()


class SpinTorqueZhangLi(object):
    def __init__(self, state):
        self._state = state
        self._mesh = state._mesh
        self._b = self._state._material["b"]
        self._xi = self._state._material["xi"]
        if isinstance(self._state.j, list):
            self._j = state.tensor(self._state.j)
        elif isinstance(self._state.j, torch.Tensor):
            if len(self._state.j.shape) == 4:
                self._j = self._state.j[:,:,:,None,:]
            else:
                self._j = self._state.j
        elif callable(self._state.j):
            self._lambda_j = self._state.j
        else:
            raise TypeError("state.j needs to be 'list', 'torch.Tensor', or 'function'!")

    @timedmethod
    def h(self, t, m):
        if hasattr(self, "_lambda_j"):
            j = self._state.tensor(self._lambda_j(t))
        else:
            j = self._j

        dim = [i for i in range(3) if self._mesh.n[i] > 1]
        dx = [self._mesh.dx[i] for i in range(3) if self._mesh.n[i] > 1]
        jgradm = (torch.stack(torch.gradient(m, spacing=dx, dim=dim), dim=-1)*j[...,dim]).sum(axis=-1)

        return self._b / constants.gamma * (torch.cross(m, jgradm) + self._xi * jgradm)

    def E(self, t, m):
        raise NotImplemented()
