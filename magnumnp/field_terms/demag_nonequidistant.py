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

from magnumnp.common import logging, timedmethod, constants, Timer
from .field_terms import LinearFieldTerm
from . import newell_f, newell_g, dipole_f, dipole_g, newell_N
import numpy as np
import torch
from time import time
import os

__all__ = ["DemagFieldNonEquidistant"]

class DemagFieldNonEquidistant(LinearFieldTerm):
    r"""
    Demagneti_dstation Field:

    The dipole-dipole interaction gives rise to a long-range interaction.
    The integral formulation of the corresponding Maxwell equations can
    be represented as convolution of the magneti_dstation :math:`\vec{M} = M_s \; \vec{m}` with a proper
    demagneti_dstation kernel :math:`\vec{N}`

    .. math::
        \vec{h}^\text{dem}_{\vec{i}} = \sum\limits_{\vec{j}} \vec{N}_{\vec{i} - \vec{j}} \, \vec{M}_{\vec{j}},

    The convolution can be evaluated efficiently using an FFT method.

    :param p: number of next neighbors for near field via Newell's equation (default = 20)
    :type p: int, optional
    """
    def __init__(self, p = 20):
        self._p = p

    def _init_N_component(self, state, i_dst, i_src, perm, func_near, func_far):
        # rescale dx to avoid NaNs when using single precision
        # TODO: add scale to state and rescale like in DemagField
        dx = state.mesh.dx
        z = (torch.cumsum(state.dx[2], dim=0) - state.dx[2][0])
        dx_dst = state._tensor([[dx[0], dx[1], state.dx[2][i_dst]][ind] for ind in perm])
        dx_src = state._tensor([[dx[0], dx[1], state.dx[2][i_src]][ind] for ind in perm])

        # dipole far-field
        shape = [1 if n==1 else 2*n for n in state.mesh.n[:2] + (1,)]

        ij = [torch.fft.fftshift(state._arange(n)) - n//2 for n in shape]
        ij = torch.meshgrid(*ij,indexing='ij')

        xyz = [[ij[0]*dx[0], ij[1]*dx[1], ij[2]*0. + z[i_dst]+state.dx[2][i_dst]/2. - (z[i_src]+state.dx[2][i_src]/2.)][ind] for ind in perm] # diff of cell centers
        Nc = func_far(*xyz) * torch.prod(dx_src) / (4.*torch.pi)

        # newell near-field
        n_near = np.minimum(state.mesh.n, self._p)
        n_near[2] = 1
        shape = [1 if n==1 else 2*n for n in n_near]
        N_near = state._zeros(shape)
        ij = [torch.fft.fftshift(state._arange(n)) - n//2 for n in shape]
        ij = torch.meshgrid(*ij,indexing='ij')

        xyz = [[ij[0]*dx[0], ij[1]*dx[1], z[i_dst] - z[i_src]][ind] for ind in perm] # diff of cell centers origins

        N_near = -newell_N(func_near, *xyz, *dx_dst, *dx_src) / (4.*torch.pi*torch.prod(dx_dst))

        Nc[:n_near[0]   ,:n_near[1]   ,:n_near[2]   ] = N_near[:n_near[0]   ,:n_near[1]   ,:n_near[2]   ]
        Nc[:n_near[0]   ,:n_near[1]   ,-n_near[2]+1:] = N_near[:n_near[0]   ,:n_near[1]   ,-n_near[2]+1:]
        Nc[:n_near[0]   ,-n_near[1]+1:,:n_near[2]   ] = N_near[:n_near[0]   ,-n_near[1]+1:,:n_near[2]   ]
        Nc[:n_near[0]   ,-n_near[1]+1:,-n_near[2]+1:] = N_near[:n_near[0]   ,-n_near[1]+1:,-n_near[2]+1:]
        Nc[-n_near[0]+1:,:n_near[1]   ,:n_near[2]   ] = N_near[-n_near[0]+1:,:n_near[1]   ,:n_near[2]   ]
        Nc[-n_near[0]+1:,:n_near[1]   ,-n_near[2]+1:] = N_near[-n_near[0]+1:,:n_near[1]   ,-n_near[2]+1:]
        Nc[-n_near[0]+1:,-n_near[1]+1:,:n_near[2]   ] = N_near[-n_near[0]+1:,-n_near[1]+1:,:n_near[2]   ]
        Nc[-n_near[0]+1:,-n_near[1]+1:,-n_near[2]+1:] = N_near[-n_near[0]+1:,-n_near[1]+1:,-n_near[2]+1:]

        Nc = torch.fft.rfftn(Nc, dim = [i for i in range(2) if state.mesh.n[i] > 1])#.real#.clone()
        return Nc

    def _init_N(self, state):
        if isinstance(state.mesh.dx[2], float):
            logging.warning("mesh.dx[2] should not be constant when using DemagFieldNonEquidistant! Use the equidistant DemagField otherwise!")
        if not all([isinstance(dx, float) for dx in state.mesh.dx[:2]]):
            raise ValueError("Demag field only implemented for non-equidistant z-spacings. mesh.dx[0] and mesh.dx[1] need to be constant!")

        dtype = state._dtype
        state._dtype = torch.float64 # always use double precision

        time_kernel = time()
        self._N = [None]*state.mesh.n[2]
        for i_dst in range(state.mesh.n[2]):
            self._N[i_dst] = [None]*state.mesh.n[2]
            for i_src in range(i_dst+1):
                Nxx = self._init_N_component(state, i_dst, i_src, [0,1,2], newell_f, dipole_f)
                Nxy = self._init_N_component(state, i_dst, i_src, [0,1,2], newell_g, dipole_g)
                Nxz = self._init_N_component(state, i_dst, i_src, [0,2,1], newell_g, dipole_g)
                Nyy = self._init_N_component(state, i_dst, i_src, [1,2,0], newell_f, dipole_f)
                Nyz = self._init_N_component(state, i_dst, i_src, [1,2,0], newell_g, dipole_g)
                Nzz = self._init_N_component(state, i_dst, i_src, [2,0,1], newell_f, dipole_f)

                self._N[i_dst][i_src] = [[ Nxx,  Nxy,  Nxz],
                                         [ Nxy,  Nyy,  Nyz],
                                         [ Nxz,  Nyz,  Nzz]]
                self._N[i_src][i_dst] = [[ Nxx,  Nxy, -Nxz],
                                         [ Nxy,  Nyy, -Nyz],
                                         [-Nxz, -Nyz,  Nzz]]
        logging.info(f"[DEMAG]: Time calculation of demag kernel = {time() - time_kernel} s")
        state._dtype = dtype # restore dtype

    @timedmethod
    def h(self, state):
        if not hasattr(self, "_N"):
            self._init_N(state)

        m_pad_fft = torch.fft.rfftn(state.material["Ms"] * state.m, dim = [i for i in range(2) if state.mesh.n[i] > 1], s = [2*state.mesh.n[i] for i in range(2) if state.mesh.n[i] > 1])
        hx = state._zeros(m_pad_fft.shape[:-1], dtype=state.complex_dtype)
        hy = state._zeros(m_pad_fft.shape[:-1], dtype=state.complex_dtype)
        hz = state._zeros(m_pad_fft.shape[:-1], dtype=state.complex_dtype)

        for i_dst in range(state.mesh.n[2]):
            for i_src in range(state.mesh.n[2]):
                for ax in range(3):
                    hx[:,:,i_src] += self._N[i_src][i_dst][0][ax][:,:,0]*m_pad_fft[:,:,i_dst,ax]
                    hy[:,:,i_src] += self._N[i_src][i_dst][1][ax][:,:,0]*m_pad_fft[:,:,i_dst,ax]
                    hz[:,:,i_src] += self._N[i_src][i_dst][2][ax][:,:,0]*m_pad_fft[:,:,i_dst,ax]

        hx = torch.fft.irfftn(hx, dim = [i for i in range(2) if state.mesh.n[i] > 1])
        hy = torch.fft.irfftn(hy, dim = [i for i in range(2) if state.mesh.n[i] > 1])
        hz = torch.fft.irfftn(hz, dim = [i for i in range(2) if state.mesh.n[i] > 1])

        return torch.stack([hx[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hy[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hz[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]]], dim=3)
