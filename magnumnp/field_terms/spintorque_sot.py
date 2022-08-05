from magnumnp.common import timedmethod, constants
import torch

__all__ = ["SpinOrbitTorque"]

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

