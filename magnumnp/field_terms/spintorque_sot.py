from magnumnp.common import timedmethod, constants
import torch

__all__ = ["SpinOrbitTorque"]

class SpinOrbitTorque(object):
    def __init__(self, state):
        self._state = state
        self._mesh = self._state._mesh
        self.gamma = constants.gamma
        self.alpha = self._state._material["alpha"]
        self.eta_damp = self._state._material["eta_damp"]
        self.eta_field = self._state._material["eta_field"]
        self.Ms = self._state._material["Ms"]
        self.je = self._state._material["je"]
        self.d = self._state._material["d"] #thickness of FM
        p = self._state._material["p"] #thickness of FM
         
        self._p = self._state.zeros(self._mesh.n + (3,))
        if isinstance(p, list):
            self._p[:,:,:,:] = self._state.tensor(p)
        elif isinstance(p, torch.Tensor):
            self._p[:,:,:,:] = p
        elif callable(p):
            self._lambda_p = p
        else:
            raise TypeError("p needs to be 'list', 'torch.Tensor', or 'function'!")

    @timedmethod
    def h(self, t, m):
        if hasattr(self, "_lambda_h"):
            p = self._lambda_p(t)
        else:
            p = self._p
        h = self.eta_damp * torch.cross(m, p) + self.eta_field * p
        h *= -self.je *constants.hbar / (2. * constants.e * self.Ms * constants.mu_0 * self.d) 
        return h

    def E(self, t, m):
        raise NotImplemented()

