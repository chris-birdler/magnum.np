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

from magnumnp.common import timedmethod, constants
import torch
from .field_terms import LinearFieldTerm

__all__ = ["AtomisticRuXExchangeField"]

class AtomisticRuXExchangeField(LinearFieldTerm):
    r"""
    Noncollinear coupling caused by Ru-X Spacer Layers:
    see: https://www.science.org/doi/full/10.1126/sciadv.abd8861
         https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.054401

    :param Jij: List of Coupling Constants
    :type Jij: list, optional
    """
    parameters = ["J"]

    def __init__(self, Jij, hom = None, **kwargs):
        self._hom = hom
        self._Jij = Jij
        super().__init__(**kwargs)

    def _init_coefficients(self, state):
        # everything except state.m is time-invariant, so the interaction
        # coefficients (material prefactor * J-lookup * spacer masks) are
        # computed once and cached (trades ~6*len(J) stored fields for the
        # per-step recomputation; this term targets small atomistic meshes)
        J = self._Jij
        mu = state.material["Ms"]*state.mesh.cell_volumes
        f = 2. / (constants.mu_0 * mu)
        mat = state.material["RuxDistribution"].repeat_interleave(3).reshape(state.mesh.n + (3,))
        coeffs = []
        for i in range(J.size(1)):
            smx, smy, smz = 1., 1., 1.
            for k in range(i):
                smx *= (1 - torch.sign(mat[1+k:-i+k,:,:]))
                smy *= (1 - torch.sign(mat[:,1+k:-i+k,:]))
                smz *= (1 - torch.sign(mat[:,:,1+k:-i+k]))

            Jx = J[mat[+i+1:,:,:,:] + mat[:-i-1,:,:,:], i] * smx
            Jy = J[mat[:,+i+1:,:,:] + mat[:,:-i-1,:,:], i] * smy
            Jz = J[mat[:,:,+i+1:,:] + mat[:,:,:-i-1,:], i] * smz

            coeffs.append((f[+i+1:,:,:,:] * Jx, f[:-i-1,:,:,:] * Jx,
                           f[:,+i+1:,:,:] * Jy, f[:,:-i-1,:,:] * Jy,
                           f[:,:,+i+1:,:] * Jz, f[:,:,:-i-1,:] * Jz))
        return coeffs

    @timedmethod
    def h(self, state):
        # identity catches parameter reassignment, the version counters catch
        # in-place/domain writes (which keep the object identity)
        Ms = state.material["Ms"]
        dist = state.material["RuxDistribution"]
        key = (id(Ms), Ms._version, id(dist), dist._version)
        if getattr(self, "_coeffs_key", None) != key:
            self._coeffs = self._init_coefficients(state)
            self._coeffs_key = key

        h = torch.zeros_like(state.m)
        m = state.m
        for i, (cxa, cxb, cya, cyb, cza, czb) in enumerate(self._coeffs):
            h[+i+1:,:,:,:] += cxa * m[:-i-1,:,:,:]
            h[:-i-1,:,:,:] += cxb * m[+i+1:,:,:,:]
            h[:,+i+1:,:,:] += cya * m[:,:-i-1,:,:]
            h[:,:-i-1,:,:] += cyb * m[:,+i+1:,:,:]
            h[:,:,+i+1:,:] += cza * m[:,:,:-i-1,:]
            h[:,:,:-i-1,:] += czb * m[:,:,+i+1:,:]

        if self._hom is not None:
            for i in self._hom:
                h[i[0]:i[1],i[2]:i[3],i[4]:i[5],0] = torch.mean(h[i[0]:i[1],i[2]:i[3],i[4]:i[5], 0])
                h[i[0]:i[1],i[2]:i[3],i[4]:i[5],1] = torch.mean(h[i[0]:i[1],i[2]:i[3],i[4]:i[5], 1])
                h[i[0]:i[1],i[2]:i[3],i[4]:i[5],2] = torch.mean(h[i[0]:i[1],i[2]:i[3],i[4]:i[5], 2])
        return h
