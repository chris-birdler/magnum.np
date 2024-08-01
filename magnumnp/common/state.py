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
from magnumnp.common.utils import randM

__all__ = ["State"]


class State(object):
    def __init__(self, mesh):
        self.mesh = mesh

        self._material = Material(self)
        self._t = torch.tensor(0.)
        self._step = 0
        self._dt = 0.

        dtype_str = str(self.dtype).split('.')[1]
        logging.info_green("[State] running on device: %s (dtype = %s)" % (self.device, dtype_str))
        logging.info_green("[Mesh] %s" % mesh)

    # Time
    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, value):
        if isinstance(value, (int, float)):
            self._t = torch.tensor(float(value))
        else:
            self._t = value

    # Material
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

    # Current density
    @property
    def j(self):
        return self._j(self)

    @j.setter
    def j(self, value):
        if callable(value):
            self._j = value
        else:
            self._j = lambda state: value

    # Temperature
    @property
    def T(self):
        return self._T(self)

    @T.setter
    def T(self, value):
        if callable(value):
            self._T = value
        else:
            self._T = lambda state: value

    @property
    def dtype(self):
        return self.mesh.dx_tensor[0].dtype

    @property
    def device(self):
        return self.mesh.dx_tensor[0].device

    def Constant(self, c, dtype = None, requires_grad = False):
        if not isinstance(c, torch.Tensor):
            c = torch.tensor(c, dtype = dtype, device = self.device)
        if c.dim() == 0 and c.dtype != torch.bool:
            c = c.reshape(1)
        x = torch.zeros(self.mesh.n + c.shape, dtype = dtype, device = self.device)
        x[...] = c
        if requires_grad == True:
            x.requires_grad = requires_grad
        return x

    def RandM(self):
        x = self.Constant([0.,0.,0.])
        randM(x)
        return x

    def SpatialCoordinate(self):
        logging.warning("State.SpatialCoordinate() is deprecated! Use mesh.SpatialCoordinate() instead!")
        return self.mesh.SpatialCoordinate()

    def convert_tensorfield(self, value):
        ''' convert arbitrary input to tensor-fields '''
        if not isinstance(value, torch.Tensor):
            value = torch.tensor(value, dtype=self.dtype)

        if len(value.shape) == 0: # convert dim=0 tensor into dim=1 tensor
            value = value.reshape(1)
        if len(value.shape) < 3: # expand homogeneous material to [nx,ny,nz,...] tensor-field
            shape = value.shape
            value = value.reshape((1,1,1) + tuple(shape))
            value = value.expand(self.mesh.n + tuple(shape))
            #value._expanded = True # annotate expanded tensor (clone will be before individual items are modified)
            value = value.clone()
        elif len(value.shape) == 3: # scalar-field should have dimension [nx,ny,nz,1]
            value = value.unsqueeze(-1)
        else: # otherwise assume the dimention is correct!
            pass
        return value

    def write_vtk(self, fields, filename, scale = 1.):
        if self.mesh.is_equidistant:
            if not filename.endswith(".vti"):
                filename += ".vti"
            write_vti(fields, filename, self, scale)
        else:
            if not filename.endswith(".vtr"):
                filename += ".vtr"
            write_vtr(fields, filename, self, scale)


    def avg(self, data, cell_volumes = None, dim=(0,1,2)):
        r"""
        Average over spatial dimensions of tensor fields.

        :param data: tensor field to average
        :type A: :class:`Tensor`
        :param dim: dimensions to average over
        :type dim: tuple, optional
        :param cell_volumes: volume of each cell (required only in case of non-equidistant meshes)
        :type cell_volumes: :class:`Tensor`, optional

        :Examples:

        .. code::
            Ms_avg = avg(state.material["Ms"])
            m_avg = avg(state.m)
        """
        if self.mesh.is_equidistant:
            if data.dim() <= 1: # e.g. [0,0,1]
                return data
            elif data.dim() == 2: # state.m[domain]
                return data.mean(dim=0)
            else:                 # [nx,ny,nz,...]
                return data.mean(dim=dim)
        else: # non-equidistant
            if cell_volumes == None:
                cell_volumes = self.mesh.cell_volumes
            if data.dim() <= 1: # e.g. [0,0,1]
                return data
            if data.shape[:3] != cell_volumes.shape[:3]:
                raise ValueError("Data shape (%s) does not match cell_volumes shape (%s). When averaging over slices of non-equidistant tensors you have to provide a sliced version of state.mesh.cell_volumes!" % (str(data.shape), str(cell_volumes.shape)))
            if data.dim() == 2: # state.m[domain]
                return (data * cell_volumes).sum(dim=0) / cell_volumes.sum(dim=0)
            if data.dim() == 3: # [nx,ny,nz]
                return (data * cell_volumes.squeeze(-1)).sum(dim=dim) / cell_volumes.sum()
            # [nx,ny,nz,...]
            return (data * cell_volumes).sum(dim=dim) / cell_volumes.sum()
