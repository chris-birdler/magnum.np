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

from magnumnp.common import logging
import torch

__all__ = ["Mesh"]

class Mesh(object):
    def __init__(self, n, dx, origin=(0,0,0), pbc=(0,0,0)):
        self.n = tuple(n)
        self.dx = tuple(dx)
        self.origin = tuple(origin)
        self.pbc = tuple(pbc)

        self.is_equidistant = all([isinstance(dx, (float, int)) for dx in dx])
        self.dx_tensor = [torch.tensor(dx, dtype=torch.get_default_dtype()).expand(n) for n, dx in zip(n, dx)]

        # compute cell_volumes
        if self.is_equidistant:
            self.cell_volumes = dx[0] * dx[1] * dx[2]
        else:
            dx, dy, dz = torch.meshgrid([torch.tensor(dx) for dx in dx], indexing = "ij") # use expand for equidistant dimentions
            self.cell_volumes = (dx*dy*dz).expand(self.n).unsqueeze(-1)

    def __str__(self):
        str_dx = ["%g" % dx if isinstance(dx, (int,float)) else "XX" for dx in self.dx]
        if self.pbc[0] != 0 or self.pbc[1] != 0 or self.pbc[2] != 0:
            str_pbc = ", pbc=[%d,%d,%d])" % self.pbc
        else:
            str_pbc = ")"
        return "%dx%dx%d (dx= %s x %s x %s%s" % (*self.n, *str_dx, str_pbc)

    def SpatialCoordinate(self):
        x = self.dx_tensor[0].cumsum(0) - self.dx_tensor[0]/2. + self.origin[0]
        y = self.dx_tensor[1].cumsum(0) - self.dx_tensor[1]/2. + self.origin[1]
        z = self.dx_tensor[2].cumsum(0) - self.dx_tensor[2]/2. + self.origin[2]

        XX, YY, ZZ = torch.meshgrid(x, y, z, indexing = "ij")
        return XX, YY, ZZ
