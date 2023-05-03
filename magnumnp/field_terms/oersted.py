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

from magnumnp.common import logging, timedmethod, constants
import numpy as np
import torch
import torch.fft
from torch import asinh, atan, sqrt, log, abs
import os
from time import time

__all__ = ["OerstedField"]

def krueger_g(points):
    x = points[:,:,:,0]
    y = points[:,:,:,1]
    z = points[:,:,:,2]

    R = sqrt(x**2 + y**2 + z**2)

    res = (3.*x**2 + 3.*y**2 - 2.*z**2)*z*R/24.
    res += np.pi*z/4.*abs(x*y*z)

    mask = (x**2 + y**2).gt(0)
    res[mask] += ((x**4 - 6.*x**2*y**2 + y**4)/24. * log(z+R))[mask]

    mask = (y**2).gt(0)
    res[mask] += (x*y/6. * (y**2 - 3.*z**2) * atan(x*z/(y*R)))[mask]

    mask = (x**2).gt(0)
    res[mask] += (x*y/6. * (x**2 - 3.*z**2) * atan(y*z/(x*R)))[mask]

    mask = (y**2+z**2).gt(0)
    res[mask] += (z/6. * x * (z**2 - 3.*y**2) * log(x+R))[mask]

    mask = (x**2+z**2).gt(0)
    res[mask] += (z/6. * y * (z**2 - 3.*x**2) * log(y+R))[mask]

    return res


def dipole_g(points):
    x = points[:,:,:,0]
    y = points[:,:,:,1]
    z = points[:,:,:,2]

    R = sqrt(x**2 + y**2 + z**2)
    res = -z/R**3
    res[0,0,0] = 0.
    return res


class OerstedField(object):
    r"""
    The Oersted field created by some current density :math:`\vec{j}` can be calculated by means of the Biot-Savart law

    .. math::

        \vec{h}^\text{oersted}(\vec{x}) = \frac{1}{4 \pi} \int \vec{j}(\vec{x}') \times \frac{\vec{x}-\vec{x}'}{\vert \vec{x}-\vec{x}'\vert^3} \, d\vec{x}'.

    The occuring equations look very similar to those of the demag field [krueger], and the occuring convolution can be efficiently calculated by means of an FFT method.

    :param p: number of next neighbors for near field via Krueger's equations (default = 20)
    :type p: int, optional
    """
    def __init__(self, p = 20):
        self._p = p

    def _init_K_component(self, state, perm, func_near, func_far):
        # dipole far-field
        shape = [1 if n==1 else 2*n for n in state.mesh.n]
        ij = [torch.fft.fftshift(state._arange(n)) - n//2 for n in shape]
        ij = torch.meshgrid(*ij,indexing='ij')

        r = torch.stack([ij[ind]*state.mesh.dx[ind] for ind in perm], dim=-1)
        Kc = func_far(r) * np.prod(state.mesh.dx) / (4.*np.pi)

        # newell near-field
        n_near = np.minimum(state.mesh.n, self._p)
        K_near = state._zeros([1 if i==1 else 2*i for i in n_near])
        ij = [torch.fft.fftshift(state._arange(n)) - n//2 for n in K_near.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        for k in np.rollaxis(np.indices((3,)*3), 0, 4).reshape(27, -1) - 1:
            r = torch.stack([(ij[ind] + k[ind])*state.mesh.dx[ind] for ind in perm], dim=-1)
            K_near[:,:,:] += np.prod(2.-3*np.abs(k)) * func_near(r) / (4.*np.pi*np.prod(state.mesh.dx))

        Kc[:n_near[0]   ,:n_near[1]   ,:n_near[2]   ] = K_near[:n_near[0]   ,:n_near[1]   ,:n_near[2]   ]
        Kc[:n_near[0]   ,:n_near[1]   ,-n_near[2]+1:] = K_near[:n_near[0]   ,:n_near[1]   ,-n_near[2]+1:]
        Kc[:n_near[0]   ,-n_near[1]+1:,:n_near[2]   ] = K_near[:n_near[0]   ,-n_near[1]+1:,:n_near[2]   ]
        Kc[:n_near[0]   ,-n_near[1]+1:,-n_near[2]+1:] = K_near[:n_near[0]   ,-n_near[1]+1:,-n_near[2]+1:]
        Kc[-n_near[0]+1:,:n_near[1]   ,:n_near[2]   ] = K_near[-n_near[0]+1:,:n_near[1]   ,:n_near[2]   ]
        Kc[-n_near[0]+1:,:n_near[1]   ,-n_near[2]+1:] = K_near[-n_near[0]+1:,:n_near[1]   ,-n_near[2]+1:]
        Kc[-n_near[0]+1:,-n_near[1]+1:,:n_near[2]   ] = K_near[-n_near[0]+1:,-n_near[1]+1:,:n_near[2]   ]
        Kc[-n_near[0]+1:,-n_near[1]+1:,-n_near[2]+1:] = K_near[-n_near[0]+1:,-n_near[1]+1:,-n_near[2]+1:]

        return torch.fft.rfftn(Kc, dim = [i for i in range(3) if state.mesh.n[i] > 1])#.real.clone()

    def _init_K(self, state):
        dtype = state._dtype
        state._dtype = torch.float64 # always use double precision
        time_kernel = time()
        Kxy = self._init_K_component(state, [0,1,2], krueger_g, dipole_g).to(dtype=dtype)
        Kyz = self._init_K_component(state, [1,2,0], krueger_g, dipole_g).to(dtype=dtype)
        Kxz = self._init_K_component(state, [2,0,1], krueger_g, dipole_g).to(dtype=dtype)

        self._K = [[  0., -Kxy, +Kxz],
                   [+Kxy,   0., -Kyz],
                   [-Kxz, +Kyz,   0.]]

        logging.info(f"[OERSTED]: Time calculation of oersted kernel = {time() - time_kernel} s")
        state._dtype = dtype # restore dtype

    @timedmethod
    def h(self, state):
        if not hasattr(self, "_K"):
            self._init_K(state)

        hx = state._zeros(list(self._K[0][1].shape), dtype=state.complex_dtype)
        hy = state._zeros(list(self._K[0][1].shape), dtype=state.complex_dtype)
        hz = state._zeros(list(self._K[0][1].shape), dtype=state.complex_dtype)

        for ax in range(3):
            j_pad_fft1D = torch.fft.rfftn(state.j[:,:,:,ax], dim = [i for i in range(3) if state.mesh.n[i] > 1], s = [2*state.mesh.n[i] for i in range(3) if state.mesh.n[i] > 1])

            hx += self._K[0][ax] * j_pad_fft1D
            hy += self._K[1][ax] * j_pad_fft1D
            hz += self._K[2][ax] * j_pad_fft1D

        hx = torch.fft.irfftn(hx, dim = [i for i in range(3) if state.mesh.n[i] > 1])
        hy = torch.fft.irfftn(hy, dim = [i for i in range(3) if state.mesh.n[i] > 1])
        hz = torch.fft.irfftn(hz, dim = [i for i in range(3) if state.mesh.n[i] > 1])

        return torch.stack([hx[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hy[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]],
                            hz[:state.mesh.n[0],:state.mesh.n[1],:state.mesh.n[2]]], dim=3)
