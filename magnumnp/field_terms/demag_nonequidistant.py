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
        # the kernel is always evaluated on the CPU in double precision and
        # only the final stacked tensor is moved to the target device
        dx = state.mesh.dx
        dz = state.mesh.dx_tensor[2].detach().to(device="cpu", dtype=torch.float64)
        z = torch.cumsum(dz, dim=0) - dz[0]
        dx_dst = [[dx[0], dx[1], float(dz[i_dst])][ind] for ind in perm]
        dx_src = [[dx[0], dx[1], float(dz[i_src])][ind] for ind in perm]

        shape = self._shape(state)
        ij = [torch.fft.fftfreq(n, 1/n, dtype=torch.float64, device="cpu").reshape([n if j == i else 1 for j in range(3)])
              for i, n in enumerate(shape)] # local indices (1-D, broadcast inside func)
        ij[2] = ij[2]*0. + z[i_dst] - z[i_src] # use fixed distance for z-direction
        x, y, z = [[ij[0]*dx[0], ij[1]*dx[1], ij[2].clone()][ind] for ind in perm]

        Lx = [state.mesh.n[0]*dx[0], state.mesh.n[1]*dx[1], float(torch.cumsum(dz, dim=0)[-1])]
        Lx = [Lx[ind] for ind in perm]

        offsets = [torch.arange(-state.mesh.pbc[ind], state.mesh.pbc[ind]+1, device="cpu") for ind in perm] # offset of pseudo PBC images
        offsets = torch.stack(torch.meshgrid(*offsets, indexing="ij"), dim=-1).flatten(end_dim=-2)

        Nc = torch.zeros(shape, dtype=torch.float64, device="cpu")
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
        time_kernel = time()
        # kernel is stored as one stacked tensor N[i_dst, i_src, comp, ax, x, y],
        # so h() can contract all layer pairs in a single einsum instead of a
        # O(nz^2 * 3) python loop
        n_z = state.mesh.n[2]
        N = None
        for i_dst in range(n_z):
            for i_src in range(i_dst+1):
                Nxx = self._init_N_component(state, i_dst, i_src, [0,1,2], demag_f).squeeze(-1).to(dtype=complex_dtype[dtype])
                Nxy = self._init_N_component(state, i_dst, i_src, [0,1,2], demag_g).squeeze(-1).to(dtype=complex_dtype[dtype])
                Nxz = self._init_N_component(state, i_dst, i_src, [0,2,1], demag_g).squeeze(-1).to(dtype=complex_dtype[dtype])
                Nyy = self._init_N_component(state, i_dst, i_src, [1,2,0], demag_f).squeeze(-1).to(dtype=complex_dtype[dtype])
                Nyz = self._init_N_component(state, i_dst, i_src, [1,2,0], demag_g).squeeze(-1).to(dtype=complex_dtype[dtype])
                Nzz = self._init_N_component(state, i_dst, i_src, [2,0,1], demag_f).squeeze(-1).to(dtype=complex_dtype[dtype])

                if N is None:
                    N = torch.zeros((n_z, n_z, 3, 3) + Nxx.shape, dtype=complex_dtype[dtype], device="cpu")
                N[i_dst,i_src] = torch.stack([torch.stack([ Nxx,  Nxy,  Nxz]),
                                              torch.stack([ Nxy,  Nyy,  Nyz]),
                                              torch.stack([ Nxz,  Nyz,  Nzz])])
                if i_src != i_dst:
                    N[i_src,i_dst] = torch.stack([torch.stack([ Nxx,  Nxy, -Nxz]),
                                                  torch.stack([ Nxy,  Nyy, -Nyz]),
                                                  torch.stack([-Nxz, -Nyz,  Nzz])])
        self._N = N.to(state.device)
        logging.info(f"[DEMAG]: Time calculation of demag kernel = {time() - time_kernel} s")

    @timedmethod
    def h(self, state):
        if not hasattr(self, "_N"):
            self._init_N(state)
        dim = [i for i in range(2) if state.mesh.n[i] > 1]
        shape = self._shape(state)
        s = [shape[i] for i in dim]

        if len(dim) == 0: # single spin   TODO: remove this when torch issue #96518 has been solved
            m_pad_fft = state.material["Ms"] * state.m
        else:
            m_pad_fft = torch.fft.rfftn(state.material["Ms"] * state.m, dim = dim, s = s)

        # contract kernel N[a, b, comp, ax, x, y] with m_fft[x, y, b, ax]
        # over all layer pairs and components in one batched operation
        if not m_pad_fft.is_complex():
            m_pad_fft = m_pad_fft.to(self._N.dtype) # single-spin branch passes a real tensor
        h_fft = torch.einsum("abcdxy,xybd->xyac", self._N, m_pad_fft)

        if len(dim) == 0: # single spin   TODO: remove this when torch issue #96518 has been solved
            h = h_fft.real.clone()
        else:
            h = torch.fft.irfftn(h_fft, dim = dim)

        return h[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]].contiguous()
