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
import subprocess
import numpy as np
from magnumnp.common import logging, DecoratedTensor, Material
from magnumnp.common.io import write_vti, write_vtr

__all__ = ["State", "complex_dtype"]

complex_dtype = {
    torch.float: torch.complex,
    torch.float32: torch.complex64,
    torch.float64: torch.complex128
    }

class State(object):
    def __init__(self, mesh, t0 = 0., device = None, dtype = None):
        if device == None:
            device_id = os.environ.get('CUDA_DEVICE')
            if device_id == None:
                device_id = get_gpu_with_least_memory()
            self._device = torch.device(f"cuda:{device_id}" if int(device_id) >= 0 else "cpu")
        else:
            self._device = device

        #TODO: add scale parameter to fix paraview issue, and use characteristic length scales
        self._dtype = dtype or torch.get_default_dtype()
        self.mesh = mesh

        self._is_equidistant = all([isinstance(dx, (float, int)) for dx in mesh.dx])
        self.dx = [self._tensor(dx).expand(n) for n, dx in zip(mesh.n, mesh.dx)] # use state.dx when a torch.tensor is needed

        # compute cell_volumes (use expand for equidistant dimentions)
        dx, dy, dz = torch.meshgrid([self._tensor(dx) for dx in mesh.dx], indexing = "ij")
        self._cell_volumes = (dx*dy*dz).expand(mesh.n).unsqueeze(-1)

        self._material = Material(self)
        self.t = t0
        self._step = 0
        self._dt = 0.

        dtype_str = str(self._dtype).split('.')[1]
        logging.info_green("[State] running on device: %s (dtype = %s)" % (self._device, dtype_str))
        logging.info_green("[Mesh] %s" % mesh)

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

    def zeros(self, size, dtype = None, **kwargs):
        dtype = dtype or self._dtype
        return torch.zeros(size, dtype=dtype, device=self._device, **kwargs)

    def arange(self, start, end = None, step=1, dtype = None, **kwargs):
        dtype = dtype or self._dtype
        if end == None:
           end = start
           start = 0
        return torch.arange(start, end, step, dtype=dtype, device=self._device, **kwargs)

    def linspace(self, start, end, steps, dtype = None, **kwargs):
        dtype = dtype or self._dtype
        return torch.linspace(start, end, steps, dtype=dtype, device=self._device, **kwargs)

    def _normal(self, mean, std, dtype = None, **kwargs):
        if not hasattr(self, "_rng"):
            self._rng = torch.Generator(device=self._device)
            self._rng.manual_seed(2147483647) # fixed seed
        dtype = dtype or self._dtype or torch.get_default_dtype()
        return torch.normal(mean, std, dtype=dtype, device=self._device, **kwargs)

    # _tensor for internal use only
    def _tensor(self, data, dtype = None):
        dtype = dtype or self._dtype
        if isinstance(data, torch.Tensor):
            return data.to(dtype=dtype, device=self._device)
        else:
            return torch.tensor(data, dtype=dtype, device=self._device)

    # TODO: avoid unneeded DecoratedTensors (e.g. state.Tensor(0.))
    def Tensor(self, data, dtype = None, requires_grad = False):
        dtype = dtype or self._dtype
        if isinstance(data, DecoratedTensor) and data.dtype == dtype:
            return data

        if isinstance(data, list) or isinstance(data, tuple) or isinstance(data, float) or isinstance(data, int) or isinstance(data, np.ndarray):
            t = DecoratedTensor(torch.tensor(data, dtype=dtype, device=self._device), self.cell_volumes)
            t.requires_grad = requires_grad
            return t
        elif isinstance(data, torch.Tensor):
            requires_grad = requires_grad or data.requires_grad
            return DecoratedTensor(data.requires_grad_(requires_grad), self.cell_volumes)
        elif callable(data):
            return lambda t: self.Tensor(data(t))
        else:
            raise TypeError("Unknown data of type '%s' (needs to be 'list', 'tuple', 'torch.Tensor', or 'function')!" % type(data))

    def Constant(self, c, dtype = None, requires_grad = False):
        dtype = dtype or self._dtype
        c = self.Tensor(c, dtype=dtype)
        x = DecoratedTensor(self.zeros(self.mesh.n + c.shape, dtype=dtype), self.cell_volumes)
        x[...] = c
        x.requires_grad = requires_grad
        return x

    def SpatialCoordinate(self):
        x = self.dx[0].cumsum(0) - self.dx[0]/2. + self.mesh.origin[0]
        y = self.dx[1].cumsum(0) - self.dx[1]/2. + self.mesh.origin[1]
        z = self.dx[2].cumsum(0) - self.dx[2]/2. + self.mesh.origin[2]

        XX, YY, ZZ = torch.meshgrid(x, y, z, indexing = "ij")
        return DecoratedTensor(XX, self.cell_volumes), DecoratedTensor(YY, self.cell_volumes), DecoratedTensor(ZZ, self.cell_volumes)

    def convert_tensorfield(self, value):
        ''' convert arbitrary input to tensor-fields '''
        value = self.Tensor(value)
        if len(value.shape) == 0: # convert dim=0 tensor into dim=1 tensor
            value = value.reshape(1)
        if len(value.shape) < 3: # expand homogeneous material to [nx,ny,nz,...] tensor-field
            shape = value.shape
            value = value.reshape((1,1,1) + tuple(shape))
            value = value.expand(self.mesh.n + tuple(shape))
            value._expanded = True # annotate expanded tensor (clone will be before individual items are modified)
        elif len(value.shape) == 3: # scalar-field should have dimension [nx,ny,nz,1]
            value = value.unsqueeze(-1)
        else: # otherwise assume the dimention is correct!
            pass
        return value

    def write_vtk(self, fields, filename):
        if self._is_equidistant:
            write_vti(fields, filename + ".vti", self)
        else:
            write_vtr(fields, filename + ".vtr", self)

    @property
    def dtype(self):
        return self._dtype

    @property
    def cell_volumes(self):
        return self._cell_volumes



def get_gpu_with_least_memory():
    if not torch.cuda.is_available():
        return -1

    import pynvml
    pynvml.nvmlInit()
    num_gpus = pynvml.nvmlDeviceGetCount()

    gpu_memory = []
    for i in range(num_gpus):
        handle = pynvml.nvmlDeviceGetHandleByIndex(i)
        mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
        gpu_memory.append(mem_info.used)

    pynvml.nvmlShutdown()
    return gpu_memory.index(min(gpu_memory))
