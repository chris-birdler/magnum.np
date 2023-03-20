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
    def u(self, state):
        sigma = state.matrial["sigma"]
        u0 = state._zeros([n+1 for n in state.mesh.n])
        rhs = torch.zeros_like(u0)

        def _M(u):
            # u[self.bnd_dirichlet] = self.u_dirichlet[self.bnd_dirichlet]
            rhs = torch.zeros_like(u)

            # x
            du = u[1:,:,:] - u[:-1,:,:]
            du /= state.mesh.dx[0]
            

            res[1:,:,:] += sigma * du[ #  m_i-1 - m_i
            
            res += simga_{i+1/2} 
            # res
            
            return res

        return self.conjugate_gradient(_M, u0, rhs)

    def j(self, state):
        sigma = state.material["sigma"] # TODO: add string parameter
        u = self.u(state) 

        # x
        du = u[1:,:,:] - u[:-1,:,:] # nx, ny+1, nz+1, 1]
        jx = du[:, 1:, 1:]
           + du[:, 1:,:-1]
           + du[:,:-1, 1:]
           + du[:,:-1,:-1]
        jx /= 4. * state.mesh.dx[0]
        
        # y
        du = u[:,1:,:] - u[:,:-1,:] # nx+1, ny, nz+1, 1]
        jy = du[ 1:,:, 1:]
           + du[ 1:,:,:-1]
           + du[:-1,:, 1:]
           + du[:-1,:,:-1]
        jy /= 4. * state.mesh.dx[1]
        
        # z
        du = u[:,:,1:] - u[:,:,:-1] # nx+1, ny+1, nz, 1]
        jz = du[ 1:, 1:,:]
           + du[ 1:,:-1,:]
           + du[:-1, 1:,:]
           + du[:-1,:-1,:]
        jz /= 4. * state.mesh.dx[2]

        return -state.material.sigma * torch.stack([jx, jy, jz], dim = -1)

j = OhmSolver.j(state) 


def conjugate_gradient(A, x, b, tol=1e-6, max_iter=1000):
    r = b - A(x)
    p = r.clone()
    rsold = torch.dot(r, r)

    for i in range(max_iter):
        Ap = A(p)
        pAp = torch.dot(p, Ap)
        alpha = rsold / pAp
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = torch.dot(r, r)
        if torch.norm(r) < tol:
            break
        beta = rsnew / rsold
        p = r + beta * p
        rsold = rsnew
    return x
