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
import torch

__all__ = ["OhmSolver"]

class OhmSolver(object):
    def __init__(self, dirichlet_bc_nodes):
        self._dirichlet_bc_nodes = dirichlet_bc_nodes

    def u(self, state, **kwargs):
        dx = state.mesh.dx
        sigma = state.material["sigma"]
        rhs = state.Constant(0.)
        u0 = state.u.clone()

        def _M(u):
            res = torch.zeros_like(u)

            # homogeneous Neumann conditions for now
            if u.shape[0] > 2:
                u[ 0,:,:] = u[ 1,:,:]
                u[-1,:,:] = u[-2,:,:]
            if u.shape[1] > 2:
                u[:, 0,:] = u[:, 1,:]
                u[:,-1,:] = u[:,-2,:]
            if u.shape[2] > 2:
                u[:,:, 0] = u[:,:, 1]
                u[:,:,-1] = u[:,:,-2]

            # set Dirichlet conditions
            u[self._dirichlet_bc_nodes] = u0[self._dirichlet_bc_nodes]

            # assemble laplace
            # x
            sigma_avg = 2. * sigma[1:,:,:] * sigma[:-1,:,:] / (sigma[1:,:,:] + sigma[:-1,:,:])
            sigma_avg = sigma_avg.nan_to_num(posinf=0, neginf=0)
            du = u[1:,:,:] - u[:-1,:,:]
            du *= sigma_avg / dx[0]
            res[1:-1,:,:] += (du[1:,:,:] - du[:-1,:,:]) / dx[0]

            # y
            sigma_avg = 2. * sigma[:,1:,:] * sigma[:,:-1,:] / (sigma[:,1:,:] + sigma[:,:-1,:])
            sigma_avg = sigma_avg.nan_to_num(posinf=0, neginf=0)
            du = u[:,1:,:] - u[:,:-1,:]
            du *= sigma_avg / dx[1]
            res[:,1:-1,:] += (du[:,1:,:] - du[:,:-1,:]) / dx[1]

            # z
            sigma_avg = 2. * sigma[:,:,1:] * sigma[:,:,:-1] / (sigma[:,:,1:] + sigma[:,:,:-1])
            sigma_avg = sigma_avg.nan_to_num(posinf=0, neginf=0)
            du = u[:,:,1:] - u[:,:,:-1]
            du *= sigma_avg / dx[2]
            res[:,:,1:-1] += (du[:,:,1:] - du[:,:,:-1]) / dx[2]

            return res

        state.u = conjugate_gradient(_M, state.u, rhs, **kwargs)
        return state.u

    def j(self, state, **kwargs):
        sigma = state.material["sigma"]
        u = self.u(state, **kwargs).squeeze(-1)
        grad = [torch.gradient(u, dim=i)[0] if state.mesh.n[i] > 2 else torch.zeros_like(u) for i in range(3)]
        return -sigma * torch.stack(grad, dim=-1)

def conjugate_gradient(A, x, b, tol=1e-6, max_iter=1000):
    r = b - A(x)
    p = r.clone()
    rsold = (r*r).sum()

    for i in range(max_iter):
        Ap = A(p)
        pAp = (p*Ap).sum()
        alpha = rsold / pAp
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = (r*r).sum()
        norm = torch.norm(r)
        if norm < tol:
            break
        beta = rsnew / rsold
        p = r + beta * p
        rsold = rsnew
        if i % 100 ==0:
            logging.info_blue("[OhmSolver] it= %d, norm= %.5e" % (i, norm))
    return x
