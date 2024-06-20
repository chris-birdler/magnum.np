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
from magnumnp.common import logging

__all__ = ["Voronoi"]

class Voronoi(object):
    def __init__(self, mesh, num_points, seed_points=None):
        if seed_points == None:
            L = torch.Tensor(mesh.dx_tuple)*torch.Tensor(mesh.n)
            offset = torch.Tensor(mesh.origin)
            self._points = L * torch.rand((num_points, 3)) + offset
        self._grid = torch.stack(mesh.SpatialCoordinate(),dim=-1).reshape(-1, 3)
        self._mesh = mesh
        self._update_domains()

        logging.info_green("[Voronoi] Setup initial Tesselation (num_points = %d)" % (num_points))

    @property
    def domains(self):
        return self._domains.reshape(self._mesh.n)

    @property
    def points(self):
        return self._points

    def _update_points(self):
        #for i in range(self._points.shape[0]):
        #    mask = self._domains == i
        #    if mask.sum() > 0:
        #        centroid = self._grid[mask].mean(dim=0)
        #        self._points[i] = centroid

        new_points = torch.zeros_like(self._points)
        counts = torch.zeros(self.points.shape[0])
        
        # Scatter add the grid points to their corresponding centroid accumulators
        domains = self._domains.unsqueeze(-1).expand(-1, 3)
        new_points = new_points.scatter_add_(0, domains, self._grid)
        
        # Count the number of points in each Voronoi cell
        counts.scatter_add_(0, domains[:,0], torch.ones(self._grid.size(0)))
        
        # Avoid division by zero
        valid_mask = counts > 0
        new_points[valid_mask] /= counts[valid_mask].unsqueeze(1)
    
        self._points = new_points

    def _update_domains(self):
        distances = torch.cdist(self._grid, self._points)
        self._domains = distances.argmin(dim=1)

    def relax(self, it = 10):
        for i in range(it):
            self._update_points()
            self._update_domains()
            logging.info_blue("[Voronoi] Relax Tesselation")
        return self.domains
