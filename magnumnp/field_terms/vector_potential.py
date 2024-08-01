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

__all__ = ["VectorPotential"]

def v(x, y, z):
    x2, y2, z2 = x**2, y**2, z**2
    r = sqrt(x2 + y2 + z2)
    return y * z * log(x + r) - x2 / 2. * atan((y*z) / (x*r))

def V(x, y, z):
    return v(x,y,z) + v(y,z,x) + v(z,x,y)

def V0(x, y, z, dz):
    return V(x, y,  dz/2.-z) \
         - V(x, y, -dz/2.-z)

def V1(x, y, z, dy, dz):
    return V0(x,  dy/2.-y, z, dz) \
         - V0(x, -dy/2.-y, z, dz)

def seidov(x, y, z, dx, dy, dz):
    res = V1( dx/2.-x, y, z, dy, dz) \
        - V1(-dx/2.-x, y, z, dy, dz)
    return res / (4*pi)

def vector_far(x, y, z, dx, dy, dz):
    res = (1./sqrt(x**2 + y**2 + z**2)).nan_to_num(posinf=0, neginf=0)
    return res * dx*dy*dz / (4.*np.pi)

def vector_func(x, y, z, dx, dy, dz, p):
    res = vector_far(x, y, z, dx, dy, dz)
    near = (x**2 + y**2 + z**2) / (dx**2 + dy**2 + dz**2) < p**2
    res[near] = seidov(x[near], y[near], z[near], dx, dy, dz)
    return res


class VectorPotential(object):
    r"""
    Vector Potential originating from current density state.j

    .. math::

        \vec{A}(\vec{x}) = \frac{1}{4 \pi} \int \vec{j}(\vec{x}') \; \frac{1}{\vert \vec{x}-\vec{x}'\vert} \, d\vec{x}'.

    Intergrals over a cuboid source region are given in [seidov], and the occuring convolution can be efficiently calculated by means of an FFT method. The integrals are evaluated for target points in the centers of each cell.

    :param p: number of next neighbors for near field via Seidov's equations (default = 20)
    :type p: int, optional
    """
    def __init__(self, p = 20, cache_dir = None):
        self._p = p
        self._cache_dir = cache_dir

    def _init_A_component(self, state):
        dx = np.array(state.mesh.dx)

        shape = [1 if n==1 else 2*n for n in state.mesh.n]
        ij = [torch.fft.fftfreq(n,1/n) for n in shape] # local indices
        ij = torch.meshgrid(*ij,indexing='ij')
        x, y, z = [ij[ind]*dx[ind] for ind in range(3)]

        Ac = vector_func(x, y, z, *dx, self._p) # TODO: handle PBCs and non-equidistant grids

        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        if len(dim) > 0:
            Ac = torch.fft.rfftn(Ac, dim = dim)

        return Ac #.real.clone() ## TODO: check symmetry


    def _init_A(self, state):
        if state.mesh.pbc[0] != 0 or state.mesh.pbc[1] != 0 or state.mesh.pbc[2] != 0:
            logging.warning(f"[VectorPotential]: PBCs are not used by VectorPotential! (mesh.pbc = %s)" % str(state.mesh.pbc))

        name = "/A_%s.pt" % str(state.mesh).replace(" ","")
        if self._cache_dir != None and os.path.isfile(self._cache_dir + name):
            Axx = torch.load(self._cache_dir + name, map_location=state.device)
            logging.info("[DEMAG]: Use cached VectorPotential kernel from '%s'" % (self._cache_dir + name))
        else:
            dtype = torch.get_default_dtype()
            torch.set_default_dtype(torch.float64) # always use double precision
            time_kernel = time()

            Axx = self._init_A_component(state).to(dtype=complex_dtype[dtype])

            logging.info(f"[DEMAG]: Time calculation of VectorPotential kernel = {time() - time_kernel} s")
            torch.set_default_dtype(dtype) # restore dtype

            # cache VectorPotential tensor
            if self._cache_dir != None:
                if not os.path.isdir(self._cache_dir):
                    os.makedirs(self._cache_dir)
                torch.save(Axx, self._cache_dir + name)
                logging.info("[DEMAG]: Save VectorPotential kernel to '%s'" % (self._cache_dir + name))

        return Axx

    @timedmethod
    def A(self, state):
        if not hasattr(self, "_Axx"):
            self._Axx = self._init_A(state)

        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        s = [2*state.mesh.n[i] for i in range(3) if state.mesh.n[i] > 1]

        Ax = torch.zeros_like(self._Axx)
        Ay = torch.zeros_like(self._Axx)
        Az = torch.zeros_like(self._Axx)

        j = state.j # state calls j(state) if j is a function
        Ai = [torch.fft.irfftn(self._Axx * torch.fft.rfftn(j[:,:,:,ax], dim=dim, s=s), dim=dim) for ax in range(3)]

        return torch.stack([Ai[0][:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            Ai[1][:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            Ai[2][:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]]], dim=3)
