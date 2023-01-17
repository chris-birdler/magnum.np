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

import torch
import os
import numpy as np
from magnumnp.common import logging, DecoratedTensor, Material

__all__ = ["State"]

class State(object):
    def __init__(self, mesh, t0=0., device=None, dtype=None):
        self.mesh = mesh
        if device == None:
            CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
            self._device = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")
        else:
            self._device = device
        self._material = Material(self)
        self._dtype = dtype
        self.t = t0
        logging.info_green("[State] running on device:%s" % self._device)
        logging.info_green("[Mesh] %dx%dx%d (size= %g x %g x %g)" % (mesh.n + mesh.dx))

    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, value):
        self._t = self._tensor(value)

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, values):
        if isinstance(values, dict):
            self._material = Material(self)
            for key, value in values.items():
                self._material[key] = value
        else:
            raise ValueError("Dictionary needs to be provided to set material")

    def _zeros(self, size, dtype=None, **kwargs):
        dtype = dtype or self._dtype or torch.get_default_dtype()
        return torch.zeros(size, dtype=dtype, device=self._device, **kwargs)

    def _arange(self, start, end=None, step=1, dtype=None, **kwargs):
        dtype = dtype or self._dtype or torch.get_default_dtype()
        if end == None:
           end = start
           start = 0
        return torch.arange(start, end, step, dtype=dtype, device=self._device, **kwargs)

    def _linspace(self, start, end, steps, dtype=None, **kwargs):
        dtype = dtype or self._dtype or torch.get_default_dtype()
        return torch.linspace(start, end, steps, dtype=dtype, device=self._device, **kwargs)

    # _tensor for internal use only
    def _tensor(self, data, dtype=None):
        dtype = dtype or self._dtype or torch.get_default_dtype()
        if isinstance(data, torch.Tensor):
            return data
        else:
            return torch.tensor(data, dtype=dtype, device=self._device)

    def Tensor(self, data, dtype=None, requires_grad=False):
        if isinstance(data, list) or isinstance(data, tuple) or isinstance(data, float) or isinstance(data, int) or isinstance(data, np.ndarray):
            dtype = dtype or self._dtype or torch.get_default_dtype()
            t = torch.tensor(data, dtype=dtype, device=self._device).as_subclass(DecoratedTensor)
            t.requires_grad = requires_grad
            return t
        elif isinstance(data, torch.Tensor):
            requires_grad = requires_grad or data.requires_grad
            return data.requires_grad_(requires_grad).as_subclass(DecoratedTensor)
        elif callable(data):
            return lambda t: self.Tensor(data(t))
        else:
            raise TypeError("Unknown data of type '%s' (needs to be 'list', 'tuple', 'torch.Tensor', or 'function')!" % type(data))

    def Constant(self, c, dtype=None):
        dtype = dtype or self._dtype or torch.get_default_dtype()
        c = self.Tensor(c, dtype=dtype)
        x = self._zeros(self.mesh.n + c.shape, dtype=dtype).as_subclass(DecoratedTensor)
        x[...] = c
        return x

    def SpatialCoordinate(self):
        x = self._arange(self.mesh.n[0]) * self.mesh.dx[0] + self.mesh.dx[0]/2. + self.mesh.origin[0]
        y = self._arange(self.mesh.n[1]) * self.mesh.dx[1] + self.mesh.dx[1]/2. + self.mesh.origin[1]
        z = self._arange(self.mesh.n[2]) * self.mesh.dx[2] + self.mesh.dx[2]/2. + self.mesh.origin[2]

        XX, YY, ZZ = torch.meshgrid(x, y, z, indexing = "ij")
        return DecoratedTensor(XX), DecoratedTensor(YY), DecoratedTensor(ZZ)
