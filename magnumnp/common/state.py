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

__all__ = ["State"]


class State(object):
    def __init__(self, mesh):
        #TODO: add scale parameter to fix paraview issue, and use characteristic length scales
        self.mesh = mesh

        self._material = Material(self)
        self.t = 0.
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
        if c.dim() == 0 and c.dtype != torch.bool:
            c = c.reshape(1)
        x = torch.zeros(self.mesh.n + c.shape, dtype = dtype)
        x[...] = c
        return x

    def SpatialCoordinate(self):
        logging.warning("State.SpatialCoordinate() is deprecated! Use mesh.SpatialCoordinate() instead!")
        return self.mesh.SpatialCoordinate()

    def write_vtk(self, fields, filename):
        if self.mesh.is_equidistant:
            write_vti(fields, filename + ".vti", self)
        else:
            write_vtr(fields, filename + ".vtr", self)
