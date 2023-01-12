from magnumnp.common import timedmethod, constants
import torch
from .field_terms import LinearFieldTerm

__all__ = ["ExchangeFieldPBC"]

class ExchangeFieldPBC(LinearFieldTerm):
    parameters = ["A"]

    def __init__(self, domain=None, **kwargs):
        self._domain = domain
        super().__init__(**kwargs)

    @timedmethod
    def h(self, state):
        h = state._zeros(state.mesh.n + (3,))

        A = state.material[self.A]
        if self._domain != None:
            A = A * self._domain[:,:,:,None]

        for dim in range(3):
            if isinstance(A, torch.Tensor) and A.dim() == 4: # TODO: A could be a 1D tensor instead of a 4D tensor field
                A_next = torch.roll(A, +1, dim) # N
                A_avg = 2.*A_next*A/(A_next+A)
                h += A_avg * (torch.roll(state.m, +1, dim) - state.m) / state.mesh.dx[dim]**2 # m_i+1 - m_i
          
                A_avg = torch.roll(A_avg, -1, dim)
                h += A_avg * (torch.roll(state.m, -1, dim) - state.m) / state.mesh.dx[dim]**2 # m_i-1 - m_i
            else:
                h += A * (torch.roll(state.m, +1, dim) - state.m) / state.mesh.dx[dim]**2 # m_i+1 - m_i
                h += A * (torch.roll(state.m, -1, dim) - state.m) / state.mesh.dx[dim]**2 # m_i-1 - m_i

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)
