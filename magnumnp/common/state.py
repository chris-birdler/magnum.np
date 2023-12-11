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
from magnumnp.common import logging, Material
from magnumnp.common.io import write_vti, write_vtr

__all__ = ["State", "complex_dtype", "avg", "normalize"]

complex_dtype = {
    torch.float: torch.complex,
    torch.float32: torch.complex64,
    torch.float64: torch.complex128
    }


def avg(data, cell_volumes = None, dim=(0,1,2)):
    if cell_volumes == None:
        if data.dim() <= 1: # e.g. [0,0,1]
            return data
        elif data.dim() == 2: # state.m[domain]
            return data.mean(dim=0)
        else:                 # [nx,ny,nz,...]
            return data.mean(dim=dim)
    else: # non-equidistant
        if data.dim() <= 1: # e.g. [0,0,1]
            return data
        elif data.dim() == 2: # state.m[domain]
            return (data * cell_volumes).sum(dim=0) / cell_volumes.sum(dim=0)
        elif data.dim() == 3: # [nx,ny,nz]
            return (data * cell_volumes.squeeze(-1)).sum(dim=dim) / cell_volumes.sum()
        else:                 # [nx,ny,nz,...]
            return (data * cell_volumes).sum(dim=dim) / cell_volumes.sum()


def normalize(data):
    data /= torch.linalg.norm(data, dim = -1, keepdim = True)
    data[...] = torch.nan_to_num(data, posinf=0, neginf=0)
    return data

class State(object):
    def __init__(self, mesh, t0 = 0.):
        #TODO: add scale parameter to fix paraview issue, and use characteristic length scales
        self.mesh = mesh

        self._is_equidistant = all([isinstance(dx, (float, int)) for dx in mesh.dx])
        self.dx = [torch.tensor(dx).expand(n) for n, dx in zip(mesh.n, mesh.dx)] # use state.dx when a torch.tensor is needed

        # compute cell_volumes (use expand for equidistant dimentions)
        dx, dy, dz = torch.meshgrid([torch.tensor(dx) for dx in mesh.dx], indexing = "ij")
        self._cell_volumes = (dx*dy*dz).expand(mesh.n).unsqueeze(-1)

        self._material = Material(self)
        self.t = t0
        self._step = 0
        self._dt = 0.

        x = torch.tensor(1.)
        dtype_str = str(x.dtype).split('.')[1]
        device = x.device
        logging.info_green("[State] running on device: %s (dtype = %s)" % (device, dtype_str))
        logging.info_green("[Mesh] %s" % mesh)

    @property
    def t(self): # TODO: should t be a float? would this break the inverse code?
        return self._t

    @t.setter
    def t(self, value):
        self._t = torch.tensor(value)

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

    def _normal(self, mean, std, **kwargs):
        if not hasattr(self, "_rng"):
            self._rng = torch.Generator()
            self._rng.manual_seed(2147483647) # fixed seed
        return torch.normal(mean, std, **kwargs)


    def Constant(self, c, dtype = None):
        c = torch.tensor(c, dtype = dtype)
        x = torch.zeros(self.mesh.n + c.shape, dtype = dtype)
        x[...] = c
        return x

    def SpatialCoordinate(self):
        x = self.dx[0].cumsum(0) - self.dx[0]/2. + self.mesh.origin[0]
        y = self.dx[1].cumsum(0) - self.dx[1]/2. + self.mesh.origin[1]
        z = self.dx[2].cumsum(0) - self.dx[2]/2. + self.mesh.origin[2]

        XX, YY, ZZ = torch.meshgrid(x, y, z, indexing = "ij")
        return XX, YY, ZZ

    def write_vtk(self, fields, filename):
        if self._is_equidistant:
            write_vti(fields, filename + ".vti", self)
        else:
            write_vtr(fields, filename + ".vtr", self)

    @property
    def cell_volumes(self):
        return self._cell_volumes
