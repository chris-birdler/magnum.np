import torch
import numpy as np
from scipy import ndimage, constants
import os
CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
cuda = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")

class ExchangeField(object):
    def __init__(self, mesh, material):
        self._mesh = mesh
        self._A = material["A"]
        self._Ms = material["Ms"]

        # initialize scratch space
        self._h = torch.zeros(mesh.n + (3,), dtype=torch.float64, device = cuda)

    def h(self, t, m):
        self._h[:,:,:,:] = 0.

        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

    #    for dim in range(3):
    #        self._h += torch.diff(m, dim=dim, append=m[bnd2[-dim:]+bnd2[:-dim]]) / self._mesh.dx[dim]**2
    #        self._h -= torch.diff(m, dim=dim, prepend=m[bnd1[-dim:]+bnd1[:-dim]]) / self._mesh.dx[dim]**2

        for dim in range(3):
            #self._h[:-1,:,:,:] += (m[1:,:,:,:] - m[:-1,:,:,:]) / self._mesh.dx[dim]**2 # m_i+1 - m_i
            #self._h[1:,:,:,:] += (m[:-1,:,:,:] - m[1:,:,:,:]) / self._mesh.dx[dim]**2 # m_i-1 - m_i

            self._h[current] += (m[next] - m[current]) / self._mesh.dx[dim]**2 # m_i+1 - m_i
            self._h[next]    += (m[current] - m[next]) / self._mesh.dx[dim]**2 # m_i-1 - m_i

            # rotate dimension
            current = current[-1:] + current[:-1]
            next = next[-1:] + next[:-1]

        self._h *= 2. * self._A / (constants.mu_0 * self._Ms)
        self._h = torch.nan_to_num(self._h, posinf=0, neginf=0)
        return self._h

    def E(self, t, m):
        return - 0.5 * constants.mu_0 * self._mesh.cell_volume \
               * np.sum(self._Ms * m * self.h(t, m))
