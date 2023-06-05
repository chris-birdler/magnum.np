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

from magnumnp.common import logging, constants, write_vti, complex_dtype
import torch
import numpy as np
from scipy.sparse.linalg import LinearOperator, aslinearoperator, eigs
from scipy.linalg import eig

__all__ = ["EigenSolver", "EigenResult"]

class EigenSolver(object):
    def __init__(self, state, linear_terms, constant_terms):
        self._linear_terms = linear_terms 
        self._state = state
        self._m0 = state.m
        self._h0 = torch.sum(sum([term.h(state) for term in self._linear_terms + constant_terms])*self._m0, dim=-1, keepdim=True)

        ez = state.Constant([1e-15,0.,1.])
        self._e1 = torch.linalg.cross(ez, self._m0)
        self._e1 = self._e1 / torch.linalg.norm(self._e1, axis=3, keepdim=True)
        self._e0 = -torch.linalg.cross(self._e1, self._m0)
        self._e0 = self._e0 / torch.linalg.norm(self._e0, axis=3, keepdim=True)

        self._it = 0

    def _C(self, m):
        self._state.m = m
        return sum([term.h(self._state) for term in self._linear_terms])
        
    def _D0(self, vv):
        self._it += 1
        if self._it % 500 == 0:
            logging.info_blue("[Eigensolver] it= %d" % self._it)

        vv = torch.from_numpy(vv).to(dtype=complex_dtype[self._state._dtype], device=self._state._device)
        vv = vv.reshape(self._m0.shape[:3] + (2,))

        # apply R
        vvv = vv[:,:,:,(0,)]*self._e0 + vv[:,:,:,(1,)]*self._e1

        # calculate A0 v = (C+I H0)*v
        hr = self._C(vvv.real)
        hi = self._C(vvv.imag)
        h = hr + 1j*hi
        h -= self._h0*vvv

        # apply B0
        rrr = -constants.gamma * torch.linalg.cross(self._m0.to(dtype = h.dtype), h)

        # apply R^T (reuse vv tensor)
        vv[:,:,:,0] = (rrr*self._e0).sum(dim=-1)
        vv[:,:,:,1] = (rrr*self._e1).sum(dim=-1)

        return vv.reshape(-1).detach().cpu().numpy()

    def solve(self, k=10, tol=0):
        N = np.prod(self._m0.shape[:3])
        D0 = LinearOperator((2*N,2*N), self._D0, dtype=np.complex128)

        evals, evecs2D = eigs(D0, k = 2*k, which = 'SM', tol = tol)
        #evals, evecs2D = eigs(D0, k = 2*k, sigma = 0, which = 'LM', tol = tol)

        evalvecs_sorted = sorted(zip(evals,evecs2D.T), key=lambda x: np.abs(x[0].imag))
        evals = np.array([x[0] for x in evalvecs_sorted if x[0].imag > 1000.])
        evecs2D = np.array([x[1] for x in evalvecs_sorted if x[0].imag > 1000.]).transpose()

        omega = self._state.Tensor(evals.imag)
        evecs2D = self._state.Tensor(torch.from_numpy(evecs2D)).reshape(self._m0.shape[:3] + (2,-1))
        return EigenResult(omega, evecs2D, self._state, m0 = self._m0, e0 = self._e0, e1 = self._e1, D0 = D0)


class EigenResult(object):
    def __init__(self, omega, evecs2D, state, **kwargs):
        self.omega = omega
        self._evecs2D = evecs2D
        self._state = state
        self.__dict__.update(kwargs)

    def store(self, filename):
        torch.save({"m0":self.m0, "omega":self.omega, "evecs2D":self._evecs2D}, filename)
        logging.info_green("[Eigensolver] Stored %d eigenvalues to '%s'it= %d" % (len(self.omega), filename))

#    @staticmethod
#    def load(state, filename):
#        stored = np.load(filename)
#        omega, evecs2D = stored['omega'], stored['evecs2D']
#
#        m0 = fd.Function(state.m.function_space())
#        m0.vector()[:] = stored['m0']
#        state.m = m0
#
#        solver = EigenSolverBase(state)
#        logging.info_green("%s: Loaded %d eigenvalues from '%s'" % (__class__.__name__, len(omega), filename))
#        return EigenResult(omega, evecs2D, state, m0 = m0, R = PETSc2CSR(solver.R), A0 = PETSc2Scipy(solver.A0), B0 = 1j*PETSc2Scipy(solver.B0))

    def evecs(self, N = slice(None)):
        vvv = self._evecs2D[:,:,:,(0,),:]*self.e0[:,:,:,:,None] + self._evecs2D[:,:,:,(1,),:]*self.e1[:,:,:,:,None]
        return vvv[...,N]

    def save_evecs3D(self, filename, which = "abs", N = slice(None)):
        if which == "abs":
            op = torch.abs
        elif which == "real":
            op = torch.real
        elif which == "imag":
            op = torch.imag
        else:
            op = which

        evecs_list = [op(v.squeeze(-1)) for v in self.evecs(N = N).split(1,dim=-1)]
        write_vti(evecs_list, filename, self._state) # TODO: write pvd instead of vti
