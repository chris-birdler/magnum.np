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

    @timedmethod
    def h(self, state):
        h = torch.zeros_like(state.m)
        J = self._Jij
        mu = state.material["Ms"]*state.mesh.cell_volumes
        mat = torch.tensor(state.material["RuxDistribution"].repeat_interleave(3).reshape(state.mesh.n + (3,)))
        Jsize = list(J.size())
        for i in range(Jsize[1]):
            f = 2. / (constants.mu_0 * mu)
            smx, smy, smz = 1., 1., 1.
            for k in range(i):
                smx *= (1 - torch.sign(mat[1+k:-i+k,:,:]))
                smy *= (1 - torch.sign(mat[:,1+k:-i+k,:]))
                smz *= (1 - torch.sign(mat[:,:,1+k:-i+k]))

            h[+i+1:,:,:,:] += f[+i+1:,:,:,:] * J[mat[+i+1:,:,:,:] + mat[:-i-1,:,:,:], i] * smx * state.m[:-i-1,:,:,:]
            h[:-i-1,:,:,:] += f[:-i-1,:,:,:] * J[mat[+i+1:,:,:,:] + mat[:-i-1,:,:,:], i] * smx * state.m[+i+1:,:,:,:]
            h[:,+i+1:,:,:] += f[:,+i+1:,:,:] * J[mat[:,+i+1:,:,:] + mat[:,:-i-1,:,:], i] * smy * state.m[:,:-i-1,:,:]
            h[:,:-i-1,:,:] += f[:,:-i-1,:,:] * J[mat[:,+i+1:,:,:] + mat[:,:-i-1,:,:], i] * smy * state.m[:,+i+1:,:,:]
            h[:,:,+i+1:,:] += f[:,:,+i+1:,:] * J[mat[:,:,+i+1:,:] + mat[:,:,:-i-1,:], i] * smz * state.m[:,:,:-i-1,:]
            h[:,:,:-i-1,:] += f[:,:,:-i-1,:] * J[mat[:,:,+i+1:,:] + mat[:,:,:-i-1,:], i] * smz * state.m[:,:,+i+1:,:]

        if self._hom is not None:
            for i in self._hom:
                h[i[0]:i[1],i[2]:i[3],i[4]:i[5],0] = torch.mean(h[i[0]:i[1],i[2]:i[3],i[4]:i[5], 0])
                h[i[0]:i[1],i[2]:i[3],i[4]:i[5],1] = torch.mean(h[i[0]:i[1],i[2]:i[3],i[4]:i[5], 1])
                h[i[0]:i[1],i[2]:i[3],i[4]:i[5],2] = torch.mean(h[i[0]:i[1],i[2]:i[3],i[4]:i[5], 2])
        return h
