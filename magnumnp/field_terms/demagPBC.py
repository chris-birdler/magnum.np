#
# This file is part of the magnum.np distribution
# (https://gitlab.com/magnum.np/magnum.np).
# Copyright (c) 2023 magnum.np team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

from magnumnp.common import timedmethod, constants
from .field_terms import LinearFieldTerm
from math import pi
import torch

__all__ = ["DemagFieldPBC"]

class DemagFieldPBC(LinearFieldTerm):
    r"""
    True-periodic demagnetization field.

    Solves the discrete Poisson equation for the magnetic scalar potential in
    Fourier space, treating the mesh as periodically continued in all
    dimensions. The k = 0 mode of the potential is set to zero, which
    corresponds to the infinite periodic continuation with exactly canceled
    surface charges (zero average field).
    """
    def _init_factors(self, state):
        # all spectral factors depend only on the mesh -> compute them once.
        # rfftn halves the last transformed dimension.
        n = state.mesh.n
        dims = [i for i in range(3) if n[i] > 1]
        n_fft = list(n)
        if len(dims) > 0:
            n_fft[dims[-1]] = n[dims[-1]] // 2 + 1

        dx = state.mesh.dx
        k = []
        for i in range(3):
            ki = 2. * pi * torch.arange(n_fft[i]) / n[i]
            k.append(ki.reshape([-1 if j == i else 1 for j in range(3)]))

        # discrete divergence/gradient stencil factors (1-D per axis)
        self._div = [(1. - torch.exp(-1j*k[i])) / dx[i] for i in range(3)]
        self._grad = [(1. - torch.exp(1j*k[i])) / dx[i] for i in range(3)]

        # negative inverse of the discrete Laplacian; the k = 0 entry is set to
        # zero (canceled surface charges -> zero average field)
        L = 4./dx[0]**2 * torch.sin(k[0]/2.)**2 \
          + 4./dx[1]**2 * torch.sin(k[1]/2.)**2 \
          + 4./dx[2]**2 * torch.sin(k[2]/2.)**2
        L[0,0,0] = 1.
        self._inv_L = -1. / L
        self._inv_L[0,0,0] = 0.

        self._dims = dims
        self._s = [n[i] for i in dims]

    @timedmethod
    def h(self, state):
        if not hasattr(self, "_dims"):
            self._init_factors(state)

        m_fft = torch.fft.rfftn(state.material["Ms"] * state.m, dim = self._dims)

        div_fft = self._div[0] * m_fft[:,:,:,0] \
                + self._div[1] * m_fft[:,:,:,1] \
                + self._div[2] * m_fft[:,:,:,2]

        u_fft = div_fft * self._inv_L

        h_fft = torch.stack([self._grad[0] * u_fft,
                             self._grad[1] * u_fft,
                             self._grad[2] * u_fft], dim=-1)

        return torch.fft.irfftn(h_fft, dim = self._dims, s = self._s)

    def E(self, state, domain = Ellipsis): # TODO: remove as soon as @compile works for DemagField
        E = -0.5 * constants.mu_0 * state.material["Ms"] * state.m * self.h(state) * state.mesh.cell_volumes
        return E[domain].sum()
