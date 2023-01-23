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

__all__ = ["ExchangeField", "ExchangeDMIField"]

class ExchangeField(LinearFieldTerm):
    r"""
    Exchange Field

    .. math::
        \vec{h}^\text{ex}_i = \frac{2}{\mu_0 \, M_{s,i}} \; \sum_{k=\pm x, \pm y,\pm z} \frac{2}{\Delta_k} \frac{A_{i+\vec{e}_k} \; A_i}{A_{i+\vec{e}_k} + A_i} \; \left( \vec{m}_{i+\vec{e}_k} - \vec{m}_i \right),

    with the vacuum permeability :math:`\mu_0`, the saturation magnetization :math:`M_s`, and the exchange constant :math:`A`. :math:`\Delta_k` and :math:`\vec{e}_k` represent the grid spacing and the unit vector in direction :math:`k`, respectively.

    :param A: Name of the material parameter for the exchange constant :math:`A`, defaults to "A"
    :type A: str, optional
    """
    parameters = ["A"]

    def __init__(self, domain=None, **kwargs):
        self._domain = domain
        super().__init__(**kwargs)

    @timedmethod
    def h(self, state):
        h = state._zeros(state.mesh.n + (3,))

        A = state.material[self.A]
        if self._domain != None:
            A = A * self._domain[:,:,:,None]
        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        for dim in range(3):
            A_avg = 2.*A[next]*A[current]/(A[next]+A[current])
            h[current] += A_avg * (state.m[next] - state.m[current]) / state.mesh.dx[dim]**2 # m_i+1 - m_i
            h[next]    += A_avg * (state.m[current] - state.m[next]) / state.mesh.dx[dim]**2 # m_i-1 - m_i

            # rotate dimension
            current = current[-1:] + current[:-1]
            next = next[-1:] + next[:-1]

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)


# deprecated version of ExchangeDMI field
class ExchangeDMIField(object):
    def __init__(self, domain=None):
        self._domain = domain

    @timedmethod
    def h(self, state):
        h = state._zeros(state.mesh.n + (3,))

        A = state.material["A"]
        Di = state.material["Di"]
        if self._domain != None:
            A = A * self._domain[:,:,:,None]
        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        # assue 1. dimension
        # calculate m_1.5:
        for dim in [0]: #range(3): # TODO: omit singelton dimentions # TODO: material needs to be a TensorField (use expand)
            A0, A1 = A.pad(dim, -1), A.pad(dim, +1)
            Di0, Di1 = Di.pad(dim, -1), Di.pad(dim, +1)
            m0, m1 = state.m.pad(dim, -1), state.m.pad(dim, +1)

            rhs = A0*m0 + A1*m1

            m_15 = torch.concat([+rhs[...,(0,)] *                       (A0 + A1) / ((A0 + A1)**2 + (state.mesh.dx[0]/4.*(Di0-Di1))**2)
                                 -rhs[...,(2,)] * state.mesh.dx[dim]/4.*(Di0-Di1) / ((A0 + A1)**2 + (state.mesh.dx[0]/4.*(Di0-Di1))**2),

                                 rhs[...,(1,)] / (A0 + A1),

                                 +rhs[...,(0,)] * state.mesh.dx[dim]/4.*(Di0-Di1) / ((A0 + A1)**2 + (state.mesh.dx[0]/4.*(Di0-Di1))**2)
                                 +rhs[...,(2,)] *                       (A0 + A1) / ((A0 + A1)**2 + (state.mesh.dx[0]/4.*(Di0-Di1))**2)], dim=-1)

            # exchange field
            h += 2.*A * (m_15[current] - 2.*state.m + m_15[next]) / state.mesh.dx[dim]**2

            # rotate dimension
            #current = current[-1:] + current[:-1]
            #next = next[-1:] + next[:-1]

        # DMI field
        dmdx = (m_15[next] - m_15[current]) / state.mesh.dx[dim]
        h += 2.*Di * torch.stack((dmdx[...,2], 0.*dmdx[...,1], -dmdx[...,0]), dim=-1)

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))
