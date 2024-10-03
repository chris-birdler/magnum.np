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

from magnumnp.common import logging, timedmethod, constants, Timer, complex_dtype
from .field_terms import LinearFieldTerm
from .demag import demag_f, demag_g
import numpy as np
import torch
from time import time
import os

__all__ = ["DemagFieldNonEquidistant"]

class DemagFieldNonEquidistant(LinearFieldTerm):
    r"""
    Demagnetization Field:

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

    def _shape(self, state):
        s = [1,1,1]
        for i in range(2):
            if state.mesh.n[i] == 1:
                continue
            if state.mesh.pbc[i] == 0:
                s[i] = 2*state.mesh.n[i]
            else:
                s[i] = state.mesh.n[i] # no need to pad if nonzero pbc
        return s

    def _init_N_component(self, state, i_dst, i_src, perm, func):
        # rescale dx to avoid NaNs when using single precision
        # TODO: add scale to state and rescale like in DemagField
        dx = state.mesh.dx
        z = (torch.cumsum(state.mesh.dx_tensor[2], dim=0) - state.mesh.dx_tensor[2][0])
        dx_dst = torch.tensor([[dx[0], dx[1], state.mesh.dx_tensor[2][i_dst]][ind] for ind in perm])
        dx_src = torch.tensor([[dx[0], dx[1], state.mesh.dx_tensor[2][i_src]][ind] for ind in perm])

        shape = self._shape(state)
        ij = [torch.fft.fftfreq(n,1/n) for n in shape] # local indices
        ij[2] = ij[2]*0. + z[i_dst] - z[i_src] # use fixed distance for z-direction
        ij = torch.meshgrid(*ij,indexing='ij')
        x, y, z = [[ij[0]*dx[0], ij[1]*dx[1], ij[2].clone()][ind] for ind in perm]

        Lx = [state.mesh.n[0]*dx[0], state.mesh.n[1]*dx[1], torch.cumsum(state.mesh.dx_tensor[2], dim=0)[-1]]
        Lx = [Lx[ind] for ind in perm]

        offsets = [torch.arange(-state.mesh.pbc[ind], state.mesh.pbc[ind]+1) for ind in perm] # offset of pseudo PBC images
        offsets = torch.stack(torch.meshgrid(*offsets, indexing="ij"), dim=-1).flatten(end_dim=-2)

        Nc = torch.zeros(shape)
        for offset in offsets:
            Nc += func(x + offset[0]*Lx[0], y + offset[1]*Lx[1], z + offset[2]*Lx[2], *dx_dst, *dx_src, self._p)

        dim = [i for i in range(2) if state.mesh.n[i] > 1]
        if len(dim) > 0:
            Nc = torch.fft.rfftn(Nc, dim = dim) # seems to be complex!!
        return Nc # .real.clone()

    def _init_N(self, state):
        if isinstance(state.mesh.dx[2], float):
            logging.warning("mesh.dx[2] should not be constant when using DemagFieldNonEquidistant! Use the equidistant DemagField otherwise!")
        if not all([isinstance(dx, float) for dx in state.mesh.dx[:2]]):
            raise ValueError("Demag field only implemented for non-equidistant z-spacings. mesh.dx[0] and mesh.dx[1] need to be constant!")

        dtype = torch.get_default_dtype()
        torch.set_default_dtype(torch.float64) # always use double precision

        time_kernel = time()
        self._N = [None]*state.mesh.n[2]
        for i_dst in range(state.mesh.n[2]):
            self._N[i_dst] = [None]*state.mesh.n[2]
            for i_src in range(i_dst+1):
                Nxx = self._init_N_component(state, i_dst, i_src, [0,1,2], demag_f).to(dtype=complex_dtype[dtype])
                Nxy = self._init_N_component(state, i_dst, i_src, [0,1,2], demag_g).to(dtype=complex_dtype[dtype])
                Nxz = self._init_N_component(state, i_dst, i_src, [0,2,1], demag_g).to(dtype=complex_dtype[dtype])
                Nyy = self._init_N_component(state, i_dst, i_src, [1,2,0], demag_f).to(dtype=complex_dtype[dtype])
                Nyz = self._init_N_component(state, i_dst, i_src, [1,2,0], demag_g).to(dtype=complex_dtype[dtype])
                Nzz = self._init_N_component(state, i_dst, i_src, [2,0,1], demag_f).to(dtype=complex_dtype[dtype])

                self._N[i_dst][i_src] = [[ Nxx,  Nxy,  Nxz],
                                         [ Nxy,  Nyy,  Nyz],
                                         [ Nxz,  Nyz,  Nzz]]
                self._N[i_src][i_dst] = [[ Nxx,  Nxy, -Nxz],
                                         [ Nxy,  Nyy, -Nyz],
                                         [-Nxz, -Nyz,  Nzz]]
        logging.info(f"[DEMAG]: Time calculation of demag kernel = {time() - time_kernel} s")
        torch.set_default_dtype(dtype) # restore dtype

    @timedmethod
    def h(self, state):
        if not hasattr(self, "_N"):
            self._init_N(state)
        dim = [i for i in range(2) if state.mesh.n[i] > 1]
        shape = self._shape(state)
        s = [shape[i] for i in dim]

        m_pad_fft = torch.fft.rfftn(state.material["Ms"] * state.m, dim = dim, s = s)
        hx = torch.zeros_like(m_pad_fft[...,0])
        hy = torch.zeros_like(m_pad_fft[...,0])
        hz = torch.zeros_like(m_pad_fft[...,0])

        for i_dst in range(state.mesh.n[2]):
            for i_src in range(state.mesh.n[2]):
                for ax in range(3):
                    hx[:,:,i_src] += self._N[i_src][i_dst][0][ax][:,:,0]*m_pad_fft[:,:,i_dst,ax]
                    hy[:,:,i_src] += self._N[i_src][i_dst][1][ax][:,:,0]*m_pad_fft[:,:,i_dst,ax]
                    hz[:,:,i_src] += self._N[i_src][i_dst][2][ax][:,:,0]*m_pad_fft[:,:,i_dst,ax]

        hx = torch.fft.irfftn(hx, dim = dim)
        hy = torch.fft.irfftn(hy, dim = dim)
        hz = torch.fft.irfftn(hz, dim = dim)

        return torch.stack([hx[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hy[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hz[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]]], dim=3)
