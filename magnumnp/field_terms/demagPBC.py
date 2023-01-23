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
from math import pi
import torch

__all__ = ["DemagFieldPBC"]

class DemagFieldPBC(object):
    @timedmethod
    def h(self, state):
        m_fft = torch.fft.fftn(state.material["Ms"] * state.m, dim = [i for i in range(3) if state.mesh.n[i] > 1]).squeeze(-1) #TODO: use rfftn -> kz should be size N//2+1
        dx, dy, dz = state.mesh.dx

        kx = (2. * pi * state._arange(state.mesh.n[0]) / state.mesh.n[0]).reshape(-1,1,1)
        ky = (2. * pi * state._arange(state.mesh.n[1]) / state.mesh.n[1]).reshape(1,-1,1)
        kz = (2. * pi * state._arange(state.mesh.n[2]) / state.mesh.n[2]).reshape(1,1,-1)

        div_fft = (1.-torch.exp(-1j*kx)) * m_fft[:,:,:,0] / dx \
                + (1.-torch.exp(-1j*ky)) * m_fft[:,:,:,1] / dy \
                + (1.-torch.exp(-1j*kz)) * m_fft[:,:,:,2] / dz

        u_fft = -div_fft / (4./dx**2*torch.sin(kx/2.)**2 + \
                            4./dy**2*torch.sin(ky/2.)**2 + \
                            4./dz**2*torch.sin(kz/2.)**2)
        u_fft[0,0,0] = 0

        h_fft = torch.empty_like(m_fft)
        h_fft[:,:,:,0] = (1.-torch.exp(1j*kx)) * u_fft / dx
        h_fft[:,:,:,1] = (1.-torch.exp(1j*ky)) * u_fft / dy
        h_fft[:,:,:,2] = (1.-torch.exp(1j*kz)) * u_fft / dz

        h = torch.fft.ifftn(h_fft, dim = [i for i in range(3) if state.mesh.n[i] > 1])
        return h.real

    def E(self, state):
        return - 0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))
