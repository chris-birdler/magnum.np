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

__all__ = ["Voronoi"]

class Voronoi(object):
    def __init__(self, mesh, num_points=None, seed_points=None):
        if seed_points == None:
            L = torch.Tensor(mesh.dx_tuple)*torch.Tensor(mesh.n)
            offset = torch.Tensor(mesh.origin)
            self._points = L * torch.rand((num_points, 3)) + offset
        self._grid = torch.stack(mesh.SpatialCoordinate(),dim=-1).reshape(-1, 3)
        self._mesh = mesh

    def _update_points(self):
        domains = self._update_domains().reshape(-1)
        for i in range(self._points.shape[0]):
            print("i:", i)
            mask = domains == i
            if mask.sum() > 0:
                centroid = self._grid[mask].mean(dim=0)
                self._points[i] = centroid
        return self._points

    def _update_domains(self):
        distances = torch.cdist(self._grid, self._points)
        self._domains = distances.argmin(dim=1).reshape(self._mesh.n)
        return self._domains

    def relax(self, it = 10):
        for i in range(it):
            self._update_points()
        return self._update_domains()
