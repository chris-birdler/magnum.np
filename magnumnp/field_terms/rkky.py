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

__all__ = ["RKKYField"]

class RKKYField(object):
    r"""
    Effective field contribution corresponding to the RKKY interaction between two interface layers.

    .. math::

        E_\text{RKKY} = \int_\Gamma J_\text{RKKY} \, \vec{m}_1 \codt \vec{m}_2 d\vec{A}

    :param state: current state
    :type state: State

    :Example:

      .. code::

        # create state with named domains from mesh
        state = State(mesh)

        # create domains as bool arrays, e.g:
        domain1 = state.zeros(n, dtype=torch.bool)
        domain1[n[0]//2:,:,:] = True

        domain2 = state.zeros(n, dtype=torch.bool)
        domain2[:-n[0]//2:,:,:] = True

        # rotate magnetization within one subdomain
        state.m[domain1] = state.Tensor([np.cos(phi), np.sin(phi), 0])

        # without interface layer, two seperate exchange fields need to be defined
        exchange1 = ExchangeField(Aex1, domain1)
        exchange2 = ExchangeField(Aex2, domain2)
        rkky = RKKYField(J_rkky, "z", id1, id2)
    """
    def __init__(self, J_rkky, dir, id1, id2, order = 0):
        self._J_rkky = J_rkky
        self._dir = dir #TODO: dir is ignored
        self._id1 = min(id1,id2)
        self._id2 = max(id1,id2)
        self._order = order

    @timedmethod
    def h(self, state):
        h = state._zeros(state.mesh.n + (3,))
        m1 = state.m[:,:,(self._id1,),:]
        m2 = state.m[:,:,(self._id2,),:]
        if self._order == 1:
            m1 += 0.5 * (m[:,:,(self._id1,),:] - m[:,:,(self._id1-1,),:])
            m2 += 0.5 * (m[:,:,(self._id2,),:] - m[:,:,(self._id2+1,),:])
        elif self._order == 2:
            m1 += 0.25 * (3*m[:,:,(self._id1,),:] - 4*m[:,:,(self._id1-1,),:] + m[:,:,(self._id1-2,),:])
            m2 += 0.25 * (3*m[:,:,(self._id2,),:] - 4*m[:,:,(self._id2+1,),:] + m[:,:,(self._id2+2,),:])

        h[:,:,(self._id1,),:] = self._J_rkky * (m2 - (m1*m2).sum(axis = 3, keepdim=True) * m1)
        h[:,:,(self._id2,),:] = self._J_rkky * (m1 - (m1*m2).sum(axis = 3, keepdim=True) * m2)

        h /= constants.mu_0 * state.material["Ms"] * state.mesh.dx[2]
        return torch.nan_to_num(h)

    def E(self, state):
        m1 = state.m[:,:,(self._id1,),:]
        m2 = state.m[:,:,(self._id2,),:]
        if self._order == 1:
            m1 += 0.5 * ( state.m[:,:,(self._id1,),:] - state.m[:,:,(self._id1-1,),:])
            m2 += 0.5 * (-state.m[:,:,(self._id2,),:] + state.m[:,:,(self._id2+1,),:])
        elif self._order == 2:
            m1 += 0.25 * (3*state.m[:,:,(self._id1,),:] - 4*state.m[:,:,(self._id1-1,),:] + state.m[:,:,(self._id1-2,),:])
            m2 += 0.25 * (3*state.m[:,:,(self._id2,),:] - 4*state.m[:,:,(self._id2+1,),:] + state.m[:,:,(self._id2+2,),:])

        E = (m2*m1).sum()
        E *= -state.mesh.dx[0] * state.mesh.dx[1] * self._J_rkky
        return E


