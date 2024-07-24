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
import numpy as np
import torch
import torch.fft
from torch import asinh, atan, sqrt, log, abs, pi
from time import time
import os

__all__ = ["DemagField"]

def f(x, y, z):
    x, y, z = abs(x), abs(y), abs(z)
    x2, y2, z2 = x**2, y**2, z**2
    r = sqrt(x2 + y2 + z2)
    res = 1.0 / 6.0 * (2*x2 - y2 - z2) * r
    res += (y / 2.0 * (z2 - x2) * asinh(y / sqrt(x2 + z2))).nan_to_num(posinf=0, neginf=0)
    res += (z / 2.0 * (y2 - x2) * asinh(z / sqrt(x2 + y2))).nan_to_num(posinf=0, neginf=0)
    res -= (x * y * z * atan(y*z / (x * r))).nan_to_num(posinf=0, neginf=0)
    return res

def g(x, y, z):
    z = abs(z)
    x2, y2, z2 = x**2, y**2, z**2
    r = sqrt(x2 + y2 + z2)
    res = -x * y * r / 3.0
    res += (x * y * z * asinh(z / sqrt(x2 + y2))).nan_to_num(posinf=0, neginf=0)
    res += (y / 6.0 * (3.0 * z2 - y2) * asinh(x / sqrt(y2 + z2))).nan_to_num(posinf=0, neginf=0)
    res += (x / 6.0 * (3.0 * z2 - x2) * asinh(y / sqrt(x2 + z2))).nan_to_num(posinf=0, neginf=0)
    res -= (z**3 / 6.0 * atan(x * y / (z * r))).nan_to_num(posinf=0, neginf=0)
    res -= (z * y2 / 2.0 * atan(x * z / (y * r))).nan_to_num(posinf=0, neginf=0)
    res -= (z * x2 / 2.0 * atan(y * z / (x * r))).nan_to_num(posinf=0, neginf=0)
    return res

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

def newell(func, x, y, z, dx, dy, dz, dX, dY, dZ):
    res = F0(func, x,           y, z, dy, dY, dz, dZ) \
        - F0(func, x - dx,      y, z, dy, dY, dz, dZ) \
        - F0(func, x + dX,      y, z, dy, dY, dz, dZ) \
        + F0(func, x - dx + dX, y, z, dy, dY, dz, dZ)
    return -res / (4.*pi*dx*dy*dz)

def dipole_f(x, y, z, dx, dy, dz, dX, dY, dZ):
    z = z + dZ/2. - dz/2. # diff of cell centers for non-equidistant demag
    res = (2.*x**2 - y**2 - z**2) * pow(x**2 + y**2 + z**2, -5./2.)
    res[0,0,0] = 0.
    return res * dx*dy*dz / (4.*pi)

def dipole_g(x, y, z, dx, dy, dz, dX, dY, dZ):
    z = z + dZ/2. - dz/2. # diff of cell centers for non-equidistant demag
    res = 3.*x*y * pow(x**2 + y**2 + z**2, -5./2.)
    res[0,0,0] = 0.
    return res * dx*dy*dz / (4.*pi)

def demag_f(x, y, z, dx, dy, dz, dX, dY, dZ, p):
    res = dipole_f(x, y, z, dx, dy, dz, dX, dY, dZ)
    near = (x**2 + y**2 + z**2) / max(dx**2 + dy**2 + dz**2, dX**2 + dY**2 + dZ**2) < p**2
    res[near] = newell(f, x[near], y[near], z[near], dx, dy, dz, dX, dY, dZ)
    return res

def demag_g(x, y, z, dx, dy, dz, dX, dY, dZ, p):
    res = dipole_g(x, y, z, dx, dy, dz, dX, dY, dZ)
    near = (x**2 + y**2 + z**2) / max(dx**2 + dy**2 + dz**2, dX**2 + dY**2 + dZ**2) < p**2
    res[near] = newell(g, x[near], y[near], z[near], dx, dy, dz, dX, dY, dZ)
    return res


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
    def __init__(self, p = 20, cache_dir = None):
        self._p = p
        self._cache_dir = cache_dir

    def _shape(self, state): # TODO: try padding to 2N-1 for small N like mumax does
        s = [1,1,1]
        for i in range(3):
            if state.mesh.n[i] == 1:
                continue
            if state.mesh.pbc[i] == 0:
                s[i] = 2*state.mesh.n[i]
            else:
                s[i] = state.mesh.n[i] # no need to pad if nonzero pbc
        return s

    def _init_N_component(self, state, perm, func):
        dx = np.array(state.mesh.dx)
        dx /= dx.min() # rescale dx to avoid NaNs when using single precision

        shape = self._shape(state)
        ij = [torch.fft.fftfreq(n,1/n) for n in shape] # local indices
        ij = torch.meshgrid(*ij,indexing='ij')
        x, y, z = [ij[ind]*dx[ind] for ind in perm]
        Lx = [state.mesh.n[ind]*dx[ind] for ind in perm]
        dx = [dx[ind] for ind in perm]

        offsets = [torch.arange(-state.mesh.pbc[ind], state.mesh.pbc[ind]+1) for ind in perm] # offset of pseudo PBC images
        offsets = torch.stack(torch.meshgrid(*offsets, indexing="ij"), dim=-1).flatten(end_dim=-2)

        Nc = torch.zeros(shape)
        for offset in offsets:
            Nc += func(x + offset[0]*Lx[0], y + offset[1]*Lx[1], z + offset[2]*Lx[2], *dx, *dx, self._p)

        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        if len(dim) > 0:
            Nc = torch.fft.rfftn(Nc, dim = dim)
        return Nc.real.clone()


    def _init_N(self, state):
        name = "/N_%s.pt" % str(state.mesh).replace(" ","")
        if self._cache_dir != None and os.path.isfile(self._cache_dir + name):
            [Nxx,Nxy,Nxz,Nyy,Nyz,Nzz] = torch.load(self._cache_dir + name, map_location=state.device)
            logging.info("[DEMAG]: Use cached demag kernel from '%s'" % (self._cache_dir + name))
        else:
            dtype = torch.get_default_dtype()
            torch.set_default_dtype(torch.float64) # always use double precision
            time_kernel = time()

            Nxx = self._init_N_component(state, [0,1,2], demag_f).to(dtype=dtype)
            Nxy = self._init_N_component(state, [0,1,2], demag_g).to(dtype=dtype)
            Nxz = self._init_N_component(state, [0,2,1], demag_g).to(dtype=dtype)
            Nyy = self._init_N_component(state, [1,2,0], demag_f).to(dtype=dtype)
            Nyz = self._init_N_component(state, [1,2,0], demag_g).to(dtype=dtype)
            Nzz = self._init_N_component(state, [2,0,1], demag_f).to(dtype=dtype)

            logging.info(f"[DEMAG]: Time calculation of demag kernel = {time() - time_kernel} s")
            torch.set_default_dtype(dtype) # restore dtype

            # cache demag tensor
            if self._cache_dir != None:
                if not os.path.isdir(self._cache_dir):
                    os.makedirs(self._cache_dir)
                torch.save([Nxx,Nxy,Nxz,Nyy,Nyz,Nzz], self._cache_dir + name)
                logging.info("[DEMAG]: Save demag kernel to '%s'" % (self._cache_dir + name))

        return [[Nxx, Nxy, Nxz],
                [Nxy, Nyy, Nyz],
                [Nxz, Nyz, Nzz]]


    @timedmethod
    def h(self, state):
        if not hasattr(self, "_N"):
            self._N = self._init_N(state)

        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        shape = self._shape(state)
        s = [shape[i] for i in dim]

        if len(dim) == 0: # single spin   TODO: remove this when torch issue #96518 has been solved
            N = torch.stack([torch.stack(self._N[0], dim=-1),
                             torch.stack(self._N[1], dim=-1),
                             torch.stack(self._N[2], dim=-1)], dim=-1)
            return (N * state.m).sum(dim=-1)

        hx = torch.zeros(self._N[0][0].shape, dtype=complex_dtype[self._N[0][0].dtype], device=state.device)
        hy = torch.zeros(self._N[0][0].shape, dtype=complex_dtype[self._N[0][0].dtype], device=state.device)
        hz = torch.zeros(self._N[0][0].shape, dtype=complex_dtype[self._N[0][0].dtype], device=state.device)

        for ax in range(3):
            m_pad_fft1D = torch.fft.rfftn(state.material["Ms"] * state.m[:,:,:,(ax,)], dim = dim, s = s).squeeze(-1)

            hx += self._N[0][ax] * m_pad_fft1D
            hy += self._N[1][ax] * m_pad_fft1D
            hz += self._N[2][ax] * m_pad_fft1D

        hx = torch.fft.irfftn(hx, dim = dim)
        hy = torch.fft.irfftn(hy, dim = dim)
        hz = torch.fft.irfftn(hz, dim = dim)

        return torch.stack([hx[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hy[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hz[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]]], dim=3)
