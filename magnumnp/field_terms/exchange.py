from magnumnp.common import timedmethod, constants
import torch

__all__ = ["ExchangeField"]

class ExchangeField(object):
    def __init__(self, state, domain=None):
        self._state = state
        self._mesh = state._mesh
        self._Ms = state._material["Ms"]
        self._A = state._material["A"]
        if domain != None:
            self._A *= domain[:,:,:,None]

    @timedmethod
    def h(self, t, m):
        h = self._state._zeros(self._mesh.n + (3,))

        A = self._A
        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        for dim in range(3):
            if isinstance(self._A, torch.Tensor):
                h[current] += (2.*A[next]*A[current]) / (A[next]+A[current]) * (m[next] - m[current]) / self._mesh.dx[dim]**2 # m_i+1 - m_i
                h[next]    += (2.*A[next]*A[current]) / (A[next]+A[current]) * (m[current] - m[next]) / self._mesh.dx[dim]**2 # m_i-1 - m_i
            else:
                h[current] += A * (m[next] - m[current]) / self._mesh.dx[dim]**2 # m_i+1 - m_i
                h[next]    += A * (m[current] - m[next]) / self._mesh.dx[dim]**2 # m_i-1 - m_i

            # rotate dimension
            current = current[-1:] + current[:-1]
            next = next[-1:] + next[:-1]

        h *= 2. / (constants.mu_0 * self._Ms)
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return h

    def E(self, t, m):
        return -0.5 * constants.mu_0 * self._mesh.cell_volume * torch.sum(self._Ms * m * self.h(t, m))
