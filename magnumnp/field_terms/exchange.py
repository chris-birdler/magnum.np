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
        self._h = state.zeros(self._mesh.n + (3,))

    @timedmethod
    def h(self, t, m):
        A = self._A
        self._h[:,:,:,:] = 0.

        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        for dim in range(3):
            if isinstance(self._A, torch.Tensor):
                self._h[current] += (2.*A[next]*A[current]) / (A[next]+A[current]) * (m[next] - m[current]) / self._mesh.dx[dim]**2 # m_i+1 - m_i
                self._h[next]    += (2.*A[next]*A[current]) / (A[next]+A[current]) * (m[current] - m[next]) / self._mesh.dx[dim]**2 # m_i-1 - m_i
            else:
                self._h[current] += A * (m[next] - m[current]) / self._mesh.dx[dim]**2 # m_i+1 - m_i
                self._h[next]    += A * (m[current] - m[next]) / self._mesh.dx[dim]**2 # m_i-1 - m_i

            # rotate dimension
            current = current[-1:] + current[:-1]
            next = next[-1:] + next[:-1]

        self._h *= 2. / (constants.mu_0 * self._Ms)
        self._h = torch.nan_to_num(self._h, posinf=0, neginf=0)
        return self._h

    def E(self, t, m):
        return -0.5 * constants.mu_0 * self._mesh.cell_volume * torch.sum(self._Ms * m * self.h(t, m))
