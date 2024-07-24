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
from scipy.spatial import KDTree

__all__ = ["Voronoi"]

class Voronoi(object):
    def __init__(self, mesh, num_grains, seed_points=None):
        """
        Creates a Voronoi Tesselation with a given number of seed points.
        The seed points will be randomly selected with a uniform distribution.
        Optional it is possible to manually provide seed_points.

        *Example*
            .. code:: python

                # create simple Voronoi tesselation
                voi = Voronoi(mesh, 10)
                state.write_vtk(voi.domains, "data/domains.vti")

                # perform Llloyd's iteration to improve tesselation
                domains = voi.relax()
                state.write_vtk(domains, "data/domains.vti")

                # add an intergrain phase
                voi.add_intergrain_phase(2)
                state.write_vtk(voi.domains, "domains.vti")

                # set parameters material 
                Ms = 8e5
                Ms_values = torch.normal(Ms, 0.1*Ms, (11,))
                Ms_values[-1] = 0. # set Ms=0 for intergrain phase
                state.material['Ms'] = Ms_values.take(voi.domains)

        *Arguments*
            mesh ([:class:`Mesh`])
                Mesh object which should be tesselated
            num_grains (int)
                Number of grains to be created
            seed_points ([:class:`torch.Tensor`])
                User provided seed points of size (num_grains, 3)
        """
        if seed_points == None:
            L = torch.tensor(mesh.dx)*torch.tensor(mesh.n)
            offset = torch.tensor(mesh.origin)
            self._points = L * torch.rand((num_grains, 3)) + offset
        self._grid = torch.stack(mesh.SpatialCoordinate(),dim=-1).reshape(-1, 3)
        self._mesh = mesh

        self._points = self._points.to(dtype=torch.float32)
        self._grid = self._grid.to(dtype=torch.float32)
        self._update_domains()

        logging.info_green("[Voronoi] Setup initial Tesselation (num_grains = %d)" % (num_grains))

    @property
    def domains(self):
        return self._domains.reshape(self._mesh.n)

    @property
    def points(self):
        return self._points

    def _update_points(self):
        ## orginal non-vectorized code (factor 2 slower)
        #for i in range(self._points.shape[0]):
        #    mask = self._domains == i
        #    if mask.sum() > 0:
        #        centroid = self._grid[mask].mean(dim=0)
        #        self._points[i] = centroid

        new_points = torch.zeros_like(self._points)
        counts = torch.zeros(self.points.shape[0])

        # Scatter add the grid points to their corresponding centroid accumulators
        new_points.index_add_(0, self._domains, self._grid)

        # Count the number of points in each Voronoi cell
        ones = torch.ones(self._grid.size(0))
        counts.index_add_(0, self._domains, ones)

        # Avoid division by zero
        valid_mask = counts > 0
        new_points[valid_mask] /= counts[valid_mask].unsqueeze(1)

        self._points = new_points

    def _update_domains(self):
        ## original vectorized code (requries N_mesh * N_points memory)
        #distances = torch.cdist(self._grid, self._points)
        #self._domains = distances.argmin(dim=1)

        tree = KDTree(self._points.cpu().numpy())
        dd, ii = tree.query(self._grid.cpu().numpy())
        self._domains = torch.tensor(ii)

    def relax(self, it = 5):
        for i in range(it):
            self._update_points()
            self._update_domains()
            logging.info_blue("[Voronoi] Relax Tesselation (Lloyd's Iteration %d)" % i)
        return self.domains

    def add_intergrain_phase(self, thickness):
        domains = self.domains
        intergrain_id = domains.max() + 1
        grad = domains
        dim = [i for i in range(3) if domains.shape[i] > 1]

        for i in range(thickness):
            grad = torch.gradient(grad, dim=dim)
            grad = sum([g.abs() for g in grad])
            domains[grad > 1e-15] = intergrain_id
