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

from magnumnp.common import logging, timedmethod, constants, complex_dtype
from .field_terms import FieldTerm
import numpy as np
import torch
import torch.fft
from torch import asinh, atan, sqrt, log, abs
import os
from time import time

__all__ = ["OerstedField"]

def g(x, y, z):
    R = sqrt(x**2 + y**2 + z**2)

    res = (3.*x**2 + 3.*y**2 - 2.*z**2)*z*R/24.
    res += np.pi*z/4.*abs(x*y*z)
    res += ((x**4 - 6.*x**2*y**2 + y**4)/24. * log(z+R)).nan_to_num(posinf=0, neginf=0)
    res += (x*y/6. * (y**2 - 3.*z**2) * atan(x*z/(y*R))).nan_to_num(posinf=0, neginf=0)
    res += (x*y/6. * (x**2 - 3.*z**2) * atan(y*z/(x*R))).nan_to_num(posinf=0, neginf=0)
    res += (z/6. * x * (z**2 - 3.*y**2) * log(x+R)).nan_to_num(posinf=0, neginf=0)
    res += (z/6. * y * (z**2 - 3.*x**2) * log(y+R)).nan_to_num(posinf=0, neginf=0)

    return res

def G1(x, y, z, dz):
    return - 1.*g(x, y, z+dz) \
           + 2.*g(x, y, z   ) \
           - 1.*g(x, y, z-dz)

def G0(x, y, z, dy, dz):
    return - 1.*G1(x, y+dy, z, dz) \
           + 2.*G1(x, y   , z, dz) \
           - 1.*G1(x, y-dy, z, dz)

def krueger_g(x, y, z, dx, dy, dz):
    ret = - 1.*G0(x+dx, y, z, dy, dz) \
          + 2.*G0(x   , y, z, dy, dz) \
          - 1.*G0(x-dx, y, z, dy, dz)
    return ret / (4.*np.pi*dx*dy*dz)


def dipole_g(x, y, z, dx, dy, dz, zero_origin = True):
    R = sqrt(x**2 + y**2 + z**2)
    res = -z/R**3
    if zero_origin:
        res[0,0,0] = 0.
    return res * dx*dy*dz / (4.*np.pi)

def oersted_g(x, y, z, dx, dy, dz, p, zero_origin = True):
    x, y, z = torch.broadcast_tensors(x, y, z) # no-op for already dense inputs
    res = dipole_g(x, y, z, dx, dy, dz, zero_origin)
    near = (x**2 + y**2 + z**2) / (dx**2 + dy**2 + dz**2) < p**2
    res[near] = krueger_g(x[near], y[near], z[near], dx, dy, dz)
    return res


class OerstedField(FieldTerm):
    r"""
    The Oersted field created by some current density :math:`\vec{j}` can be calculated by means of the Biot-Savart law

    .. math::

        \vec{h}^\text{oersted}(\vec{x}) = \frac{1}{4 \pi} \int \vec{j}(\vec{x}') \times \frac{\vec{x}-\vec{x}'}{\vert \vec{x}-\vec{x}'\vert^3} \, d\vec{x}'.

    The occuring equations [krueger] look very similar to those of the demag field, and the occuring convolution can be efficiently calculated by means of an FFT method.

    :param p: number of next neighbors for near field via Krueger's equations (default = 20)
    :type p: int, optional
    :param chunk_cells: chunk size (in cells) for the kernel setup; bounds the size of
        the temporaries during initialization (default = 4*1024**2)
    :type chunk_cells: int, optional
    """
    def __init__(self, p = 20, cache_dir = None, chunk_cells = 4<<20):
        self._p = p
        self._cache_dir = cache_dir
        self._chunk_cells = chunk_cells

    def _init_K_component(self, state, perm, func):
        # dipole far-field
        dx = np.array(state.mesh.dx)

        # the kernel is always evaluated on the CPU in double precision and
        # only the final rfft components are moved to the target device by
        # the caller; 1-D coordinate vectors + broadcasting replace the full
        # padded meshgrids, chunking bounds the elementwise temporaries
        shape = [1 if n==1 else 2*n for n in state.mesh.n]
        coords = []
        for i, n in enumerate(shape):
            v = torch.fft.fftfreq(n, 1/n, dtype=torch.float64, device="cpu") * dx[i] # local indices
            coords.append(v.reshape([n if j == i else 1 for j in range(3)]))
        x, y, z = [coords[ind] for ind in perm]
        dx = [dx[ind] for ind in perm]

        Kc = torch.zeros(shape, dtype=torch.float64, device="cpu")
        nc = max(1, int(self._chunk_cells) // max(1, shape[1]*shape[2]))
        for i0 in range(0, shape[0], nc):
            xs, ys, zs = [c[i0:i0+nc] if ind == 0 else c for ind, c in zip(perm, (x, y, z))]
            Kc[i0:i0+nc] = func(xs, ys, zs, *dx, self._p, zero_origin = (i0 == 0)) # TODO: handle PBCs and non-equidistant grids

        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        if len(dim) > 0:
            Kc = torch.fft.rfftn(Kc, dim = dim)

        return Kc #.real.clone() ## Oersted Kernel is not symmetric

    def _init_K(self, state):
        if state.mesh.pbc[0] != 0 or state.mesh.pbc[1] != 0 or state.mesh.pbc[2] != 0:
            logging.warning(f"[OERSTED]: PBCs are not used by OerstedField! (mesh.pbc = %s)" % str(state.mesh.pbc))

        name = "/K_%s.pt" % str(state.mesh).replace(" ","")
        if self._cache_dir != None and os.path.isfile(self._cache_dir + name):
            [Kxy, Kyz, Kxz] = torch.load(self._cache_dir + name, map_location=state.device)
            logging.info("[OERSTED]: Use cached Oersted kernel from '%s'" % (self._cache_dir + name))
        else:
            dtype = torch.get_default_dtype()
            time_kernel = time()

            Kxy = self._init_K_component(state, [0,1,2], oersted_g).to(dtype=complex_dtype[dtype], device=state.device)
            Kyz = self._init_K_component(state, [1,2,0], oersted_g).to(dtype=complex_dtype[dtype], device=state.device)
            Kxz = self._init_K_component(state, [2,0,1], oersted_g).to(dtype=complex_dtype[dtype], device=state.device)

            logging.info(f"[OERSTED]: Time calculation of oersted kernel = {time() - time_kernel} s")

            # cache Oersted tensor
            if self._cache_dir != None:
                if not os.path.isdir(self._cache_dir):
                    os.makedirs(self._cache_dir)
                torch.save([Kxy, Kyz, Kxz], self._cache_dir + name)
                logging.info("[OERSTED]: Save Oersted kernel to '%s'" % (self._cache_dir + name))

        return [[  0., -Kxy, +Kxz],
                [+Kxy,   0., -Kyz],
                [-Kxz, +Kyz,   0.]]

    @timedmethod
    def h(self, state):
        if not hasattr(self, "_K"):
            self._K = self._init_K(state)

        hx = torch.zeros_like(self._K[0][1])
        hy = torch.zeros_like(self._K[0][1])
        hz = torch.zeros_like(self._K[0][1])

        j = state.j # state calls j(state) if j is a function

        for ax in range(3):
            j_pad_fft1D = torch.fft.rfftn(j[:,:,:,ax], dim = [i for i in range(3) if state.mesh.n[i] > 1], s = [2*state.mesh.n[i] for i in range(3) if state.mesh.n[i] > 1])

            hx += self._K[0][ax] * j_pad_fft1D
            hy += self._K[1][ax] * j_pad_fft1D
            hz += self._K[2][ax] * j_pad_fft1D

        hx = torch.fft.irfftn(hx, dim = [i for i in range(3) if state.mesh.n[i] > 1])
        hy = torch.fft.irfftn(hy, dim = [i for i in range(3) if state.mesh.n[i] > 1])
        hz = torch.fft.irfftn(hz, dim = [i for i in range(3) if state.mesh.n[i] > 1])

        return torch.stack([hx[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hy[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hz[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]]], dim=3)
