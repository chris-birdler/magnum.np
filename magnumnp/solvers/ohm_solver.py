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

from magnumnp.common import logging, timedmethod, constants, DecoratedTensor
import torch

__all__ = ["OhmSolver"]

class OhmSolver(object):
    def __init__(self, dirichlet_bc_nodes, dirichlet_bc_u0):
        self._dirichlet_bc_u0 = dirichlet_bc_u0
        self._dirichlet_bc_nodes = dirichlet_bc_nodes

    def u(self, state):
        sigma = state.material["sigma"]
        u0 = state._zeros(state.mesh.n)
        u_pad = state._zeros([n+2 for n in state.mesh.n])
        sigma_pad = state._zeros([n+2 for n in state.mesh.n])
        sigma_pad[1:-1,1:-1,1:-1] = sigma[:,:,:,0]
        rhs = state._zeros(state.mesh.n)

        def _M(u):
            lap = torch.zeros_like(u)

            # setup padding magnetization
            u_pad[1:-1,1:-1,1:-1] = u
            #u_pad[0,1:-1,1:-1] = u[0,:,:,:] + g[0,:,:,:] * dx^? * sigma^?

            # homogeneous Neumann conditions for now
            sigma_pad[ 0,1:-1,1:-1] = sigma[ 0,:,:,0]
            sigma_pad[-1,1:-1,1:-1] = sigma[-1,:,:,0]
            sigma_pad[1:-1, 0,1:-1] = sigma[:, 0,:,0]
            sigma_pad[1:-1,-1,1:-1] = sigma[:,-1,:,0]
            sigma_pad[1:-1,1:-1, 0] = sigma[:,:, 0,0]
            sigma_pad[1:-1,1:-1,-1] = sigma[:,:,-1,0]

            u_pad[ 0,1:-1,1:-1] = u[ 0,:,:]
            u_pad[-1,1:-1,1:-1] = u[-1,:,:]
            u_pad[1:-1, 0,1:-1] = u[:, 0,:]
            u_pad[1:-1,-1,1:-1] = u[:,-1,:]
            u_pad[1:-1,1:-1, 0] = u[:,:, 0]
            u_pad[1:-1,1:-1,-1] = u[:,:,-1]

            # set Dirichlet conditions
            u_pad[-1,1:-1,1:-1] = u[-1,:,:]
            u_pad[self._dirichlet_bc_nodes] = self._dirichlet_bc_u0[self._dirichlet_bc_nodes]

            # assemble laplace
            sigma_avg = 2. * sigma_pad[1:,:,:] * sigma_pad[:-1,:,:] / (sigma_pad[1:,:,:] + sigma_pad[:-1,:,:])
            sigma_avg = sigma_avg.nan_to_num(posinf=0, neginf=0)
            du = u_pad[1:,:,:] - u_pad[:-1,:,:]
            du *= sigma_avg / state.mesh.dx[0]
            lap += (du[1:,1:-1,1:-1] - du[:-1,1:-1,1:-1]) / state.mesh.dx[0]

            sigma_avg = 2. * sigma_pad[:,1:,:] * sigma_pad[:,:-1,:] / (sigma_pad[:,1:,:] + sigma_pad[:,:-1,:])
            sigma_avg = sigma_avg.nan_to_num(posinf=0, neginf=0)
            du = u_pad[:,1:,:] - u_pad[:,:-1,:]
            du *= sigma_avg / state.mesh.dx[1]
            lap += (du[1:-1,1:,1:-1] - du[1:-1,:-1,1:-1]) / state.mesh.dx[1]

            sigma_avg = 2. * sigma_pad[:,:,1:] * sigma_pad[:,:,:-1] / (sigma_pad[:,:,1:] + sigma_pad[:,:,:-1])
            sigma_avg = sigma_avg.nan_to_num(posinf=0, neginf=0)
            du = u_pad[:,:,1:] - u_pad[:,:,:-1]
            du *= sigma_avg / state.mesh.dx[2]
            lap += (du[1:-1,1:-1,1:] - du[1:-1,1:-1,:-1]) / state.mesh.dx[2]

            return lap

        return conjugate_gradient(_M, u0, rhs)

#    def j(self, state):
#        sigma = state.material["sigma"] # TODO: add string parameter
#        u = self.u(state)
#
#        # x
#        du = u[1:,:,:] - u[:-1,:,:] # nx, ny+1, nz+1, 1]
#        jx = du[:, 1:, 1:]
#           + du[:, 1:,:-1]
#           + du[:,:-1, 1:]
#           + du[:,:-1,:-1]
#        jx /= 4. * state.mesh.dx[0]
#
#        # y
#        du = u[:,1:,:] - u[:,:-1,:] # nx+1, ny, nz+1, 1]
#        jy = du[ 1:,:, 1:]
#           + du[ 1:,:,:-1]
#           + du[:-1,:, 1:]
#           + du[:-1,:,:-1]
#        jy /= 4. * state.mesh.dx[1]
#
#        # z
#        du = u[:,:,1:] - u[:,:,:-1] # nx+1, ny+1, nz, 1]
#        jz = du[ 1:, 1:,:]
#           + du[ 1:,:-1,:]
#           + du[:-1, 1:,:]
#           + du[:-1,:-1,:]
#        jz /= 4. * state.mesh.dx[2]
#
#        return -state.material.sigma * torch.stack([jx, jy, jz], dim = -1)


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
        print("i:", i, "norm:", norm)
    return x
