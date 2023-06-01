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

from magnumnp.common import constants
import torch
import numpy as np
from scipy.sparse.linalg import LinearOperator, aslinearoperator, eigs
from scipy.linalg import eig

__all__ = ["EigenSolver"]

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
            print("D0(2D):", self._it)

        vv = self._state.Tensor(torch.from_numpy(vv))
        vv = vv.reshape(self._m0.shape[:3] + (2,))

        # apply R
        vvv = vv[:,:,:,(0,)]*self._e0 + vv[:,:,:,(1,)]*self._e1

        # calculate A0 v = (C+I H0)*v
        hr = self._C(vvv.real)
        hi = self._C(vvv.imag)
        h = hr + 1j*hi
        h -= self._h0*vvv

        # apply B0
        rrr = -constants.gamma * torch.linalg.cross(self._m0.to(dtype = torch.complex128), h)

        # apply R^T
        rr = self._state.zeros(vv.shape[:3] + (2,), dtype=torch.complex128)
        rr[:,:,:,0] = (rrr*self._e0).sum(dim=3)
        rr[:,:,:,1] = (rrr*self._e1).sum(dim=3)

        return rr.reshape(-1).detach().cpu().numpy()

    def solve(self, k=10, tol=0, method="eigs"):
        N = np.prod(self._m0.shape[:3])
        D0 = LinearOperator((2*N,2*N), self._D0, dtype=np.complex128)

        if method == "eigs":
            evals, evecs2D = eigs(D0, k = 2*k, which = 'SM', tol = tol)
            #evals, evecs2D = eigs(D0, k = 2*k, sigma = 0, which = 'LM', tol = tol)

            evalvecs_sorted = sorted(zip(evals,evecs2D.T), key=lambda x: np.abs(x[0].imag))
            evals = np.array([x[0] for x in evalvecs_sorted if x[0].imag > 1000.])
            evecs2D = np.array([x[1] for x in evalvecs_sorted if x[0].imag > 1000.]).transpose()
        elif method == "eig":
            D0_op = matrix_from_operator(D0)
            evals, evecs2D = eig(D0_op)

            evalvecs_sorted = sorted(zip(evals,evecs2D.T), key=lambda x: np.abs(x[0].imag))
            evals = np.array([x[0] for x in evalvecs_sorted if x[0].imag > 1000.])[:k]
            evecs2D = np.array([x[1] for x in evalvecs_sorted if x[0].imag > 1000.]).transpose()[:,:k]

        return evals, evecs2D

    def save_evecs3D(self, evecs2D, filename):
        vv = evecs2D.reshape(self.m0.shape[:3] + (2,-1))
        for i in list(range(vv.shape[-1])):
            vvv = vv[:,:,:,(0,),i]*self._e0 + vv[:,:,:,(1,),i]*self._e1
            write_vtr.write_vtr(np.abs(vvv), "%s_%02d" % (filename, i))


##############################
### some utility functions ###
##############################
def matrix_from_operator(op, N=None): # NOTE: this may change m0
    if N is None:
        N = aslinearoperator(op).shape[1]

    mat = np.zeros((op.shape[0],N), dtype=np.complex128)
    for i in range(N):
        m = np.zeros(op.shape[1])
        m[i] = 1.0
        mat[:,i] = aslinearoperator(op)(m)
    return mat[:N,:N]

def plot_operator(op, filename, N=None, log=False):
    from matplotlib import pyplot as plt
    from matplotlib import cm as cm

    if N == None:
        N = op.shape[0]
    X = matrix_from_operator(op, N=N)
    #X = X.imag.copy()
    X = np.absolute(X)

    if log is True:
        X = np.log(np.abs(X))

    fig, ax = plt.subplots()
    img = ax.imshow(X, cmap=cm.jet, interpolation='nearest')#, vmin=-50, vmax=30)
    numrows, numcols = X.shape

    def format_coord(x, y):
        col = int(x + 0.5)
        row = int(y + 0.5)
        if col >= 0 and col < numcols and row >= 0 and row < numrows:
            z = X[row, col]
            return 'x=%1.4g, y=%1.4g, z=%1.4g' % (x, y, z)
        else:
            return 'x=%1.4g, y=%1.4g' % (x, y)

    ax.format_coord = format_coord
    fig.colorbar(img)
    plt.savefig(filename)
    return X

