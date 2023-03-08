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
import numpy as np
import torch
import torch.fft
from torch import asinh, atan, sqrt, log, abs
from time import time
import os

__all__ = ["DemagField", "newell_f", "newell_g", "dipole_f", "dipole_g", "newell_N"]

def newell_f(x, y, z):
    x, y, z = abs(x), abs(y), abs(z)

    result = 1.0 / 6.0 * (2*x**2 - y**2 - z**2) * sqrt(x**2 + y**2 + z**2)

    mask = (x**2 + z**2).gt(0)
    result[mask] += (y / 2.0 * (z**2 - x**2) * asinh(y / sqrt(x**2 + z**2)))[mask]

    mask = (x**2 + y**2).gt(0)
    result[mask] += (z / 2.0 * (y**2 - x**2) * asinh(z / sqrt(x**2 + y**2)))[mask]

    mask = (x * (x**2 + y**2 + z**2)).gt(0)
    result[mask] -= (x * y * z * atan(y*z / (x * sqrt(x**2 + y**2 + z**2))))[mask]

    return result

def newell_g(x, y, z):
    z = abs(z)

    result = - x*y * sqrt(x**2 + y**2 + z**2) / 3.0

    mask = (x**2 + y**2).gt(0) # x**2 + y**2 > 0
    result[mask] += (x*y*z * asinh(z / sqrt(x**2 + y**2)))[mask]

    mask = (y**2 + z**2).gt(0)
    result[mask] += (y / 6.0 * (3.0 * z**2 - y**2) * asinh(x / sqrt(y**2 + z**2)))[mask]

    mask = (x**2 + z**2).gt(0)
    result[mask] += (x / 6.0 * (3.0 * z**2 - x**2) * asinh(y / sqrt(x**2 + z**2)))[mask]

    mask = (z * (x**2 + y**2 + z**2)).ne(0)
    result[mask] -= ( z**3 / 6.0 * atan(x*y / (z * sqrt(x**2 + y**2 + z**2))))[mask]

    mask = (y * (x**2 + y**2 + z**2)).ne(0)
    result[mask] -= (z * y**2 / 2.0 * atan(x*z / (y * sqrt(x**2 + y**2 + z**2))))[mask]

    mask = (x * (x**2 + y**2 + z**2)).ne(0)
    result[mask] -= (z * x**2 / 2.0 * atan(y*z / (x * sqrt(x**2 + y**2 + z**2))))[mask]

    return result

def F1(func, x, y, z, dz, dZ):
    return func(x, y, z      + dZ) \
         - func(x, y, z          ) \
         - func(x, y, z - dz + dZ) \
         + func(x, y, z - dz     )

def F0(func, x, y, z, dy, dY, dz, dZ):
    return F1(func, x, y      + dY, z, dz, dZ) \
         - F1(func, x, y,           z, dz, dZ) \
         - F1(func, x, y - dy + dY, z, dz, dZ) \
         + F1(func, x, y - dy,      z, dz, dZ)

def newell_N(func, x, y, z, dx, dy, dz, dX, dY, dZ):
    return F0(func, x,           y, z, dy, dY, dz, dZ) \
         - F0(func, x - dx,      y, z, dy, dY, dz, dZ) \
         - F0(func, x + dX,      y, z, dy, dY, dz, dZ) \
         + F0(func, x - dx + dX, y, z, dy, dY, dz, dZ)


def dipole_f(x, y, z):
    result = (2.*x**2 - y**2 - z**2) * pow(x**2 + y**2 + z**2, -5./2.)
    result[0,0,0] = 0.
    return result

def dipole_g(x, y, z):
    result = 3.*x*y * pow(x**2 + y**2 + z**2, -5./2.)
    result[0,0,0] = 0.
    return result


class DemagField(LinearFieldTerm):
    r"""
    Demagnetization Field:

    The dipole-dipole interaction gives rise to a long-range interaction.
    The integral formulation of the corresponding Maxwell equations can
    be represented as convolution of the magnetization :math:`\vec{M} = M_s \; \vec{m}` with a proper
    demagnetization kernel :math:`\vec{N}`

    .. math::
        \vec{h}^\text{dem}_{\vec{i}} = \sum\limits_{\vec{j}} \vec{N}_{\vec{i} - \vec{j}} \, \vec{M}_{\vec{j}},

    The convolution can be evaluated efficiently using an FFT method.

    :param p: number of next neighbors for near field via Newell's equation (default = 20)
    :type p: int, optional
    """
    def __init__(self, p = 20):
        self._p = p

    def _init_N_component(self, state, perm, func_near, func_far):
        # rescale dx to avoid NaNs when using single precision
        dx = np.array(state.mesh.dx)
        dx /= dx.min()

        # dipole far-field
        shape = [1 if n==1 else 2*n for n in state.mesh.n]
        ij = [torch.fft.fftshift(state._arange(n)) - n//2 for n in shape]
        ij = torch.meshgrid(*ij,indexing='ij')

        xyz = [ij[ind]*dx[ind] for ind in perm]
        Nc = func_far(*xyz) * np.prod(dx) / (4.*np.pi)

        # newell near-field
        n_near = np.minimum(state.mesh.n, self._p)
        N_near = state._zeros([1 if n==1 else 2*n for n in n_near])
        ij = [torch.fft.fftshift(state._arange(n)) - n//2 for n in N_near.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        xyz = [ij[ind]*dx[ind] for ind in perm]
        dx = [dx[ind] for ind in perm]
        N_near = -newell_N(func_near, *xyz, *dx, *dx) / (4.*np.pi*np.prod(dx))

        Nc[:n_near[0]   ,:n_near[1]   ,:n_near[2]   ] = N_near[:n_near[0]   ,:n_near[1]   ,:n_near[2]   ]
        Nc[:n_near[0]   ,:n_near[1]   ,-n_near[2]+1:] = N_near[:n_near[0]   ,:n_near[1]   ,-n_near[2]+1:]
        Nc[:n_near[0]   ,-n_near[1]+1:,:n_near[2]   ] = N_near[:n_near[0]   ,-n_near[1]+1:,:n_near[2]   ]
        Nc[:n_near[0]   ,-n_near[1]+1:,-n_near[2]+1:] = N_near[:n_near[0]   ,-n_near[1]+1:,-n_near[2]+1:]
        Nc[-n_near[0]+1:,:n_near[1]   ,:n_near[2]   ] = N_near[-n_near[0]+1:,:n_near[1]   ,:n_near[2]   ]
        Nc[-n_near[0]+1:,:n_near[1]   ,-n_near[2]+1:] = N_near[-n_near[0]+1:,:n_near[1]   ,-n_near[2]+1:]
        Nc[-n_near[0]+1:,-n_near[1]+1:,:n_near[2]   ] = N_near[-n_near[0]+1:,-n_near[1]+1:,:n_near[2]   ]
        Nc[-n_near[0]+1:,-n_near[1]+1:,-n_near[2]+1:] = N_near[-n_near[0]+1:,-n_near[1]+1:,-n_near[2]+1:]

        Nc = torch.fft.rfftn(Nc, dim = [i for i in range(3) if state.mesh.n[i] > 1])
        return Nc.real.clone()

    def _init_N(self, state):
        dtype = state._dtype
        state._dtype = torch.float64 # always use double precision
        time_kernel = time()

        Nxx = self._init_N_component(state, [0,1,2], newell_f, dipole_f).to(dtype=dtype)
        Nxy = self._init_N_component(state, [0,1,2], newell_g, dipole_g).to(dtype=dtype)
        Nxz = self._init_N_component(state, [0,2,1], newell_g, dipole_g).to(dtype=dtype)
        Nyy = self._init_N_component(state, [1,2,0], newell_f, dipole_f).to(dtype=dtype)
        Nyz = self._init_N_component(state, [1,2,0], newell_g, dipole_g).to(dtype=dtype)
        Nzz = self._init_N_component(state, [2,0,1], newell_f, dipole_f).to(dtype=dtype)

        self._N = [[Nxx, Nxy, Nxz],
                   [Nxy, Nyy, Nyz],
                   [Nxz, Nyz, Nzz]]
        logging.info(f"[DEMAG]: Time calculation of demag kernel = {time() - time_kernel} s")
        state._dtype = dtype # restore dtype

    @timedmethod
    def h(self, state):
        if not hasattr(self, "_N"):
            self._init_N(state)

        hx = state._zeros(self._N[0][0].shape, dtype=state.complex_dtype)
        hy = state._zeros(self._N[0][0].shape, dtype=state.complex_dtype)
        hz = state._zeros(self._N[0][0].shape, dtype=state.complex_dtype)
        for ax in range(3):
            m_pad_fft1D = torch.fft.rfftn(state.material["Ms"] * state.m[:,:,:,(ax,)], dim = [i for i in range(3) if state.mesh.n[i] > 1], s = [2*state.mesh.n[i] for i in range(3) if state.mesh.n[i] > 1]).squeeze(-1)

            hx += self._N[0][ax] * m_pad_fft1D
            hy += self._N[1][ax] * m_pad_fft1D
            hz += self._N[2][ax] * m_pad_fft1D

        hx = torch.fft.irfftn(hx, dim = [i for i in range(3) if state.mesh.n[i] > 1])
        hy = torch.fft.irfftn(hy, dim = [i for i in range(3) if state.mesh.n[i] > 1])
        hz = torch.fft.irfftn(hz, dim = [i for i in range(3) if state.mesh.n[i] > 1])

        return torch.stack([hx[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hy[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hz[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]]], dim=3)

    def E(self, state):
        return - 0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))
