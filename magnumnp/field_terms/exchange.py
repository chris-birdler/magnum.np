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

__all__ = ["ExchangeField"]

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
        A = state.material[self.A].torch_tensor
        Ms = state.material["Ms"].torch_tensor
        m = state.m.torch_tensor
        if self._domain != None:
            A = A * self._domain[:,:,:,None]
        dx = state.dx[0].reshape(-1,1,1,1)
        dy = state.dx[1].reshape(1,-1,1,1)
        dz = state.dx[2].reshape(1,1,-1,1)
        h = self._h(m, A, Ms, dx, dy, dz, state)
        return state.Tensor(h)

    @torch.compile
    def _h(self, m, A, Ms, dx, dy, dz, state):
        h = state._zeros(state.mesh.n + (3,))

        A_avg = 2.*A[1:,:,:]*A[:-1,:,:] / (A[1:,:,:]*dx[:-1,:,:,:] + A[:-1,:,:]*dx[1:,:,:,:])
        h[:-1,:,:,:] += A_avg * (m[ 1:,:,:,:]-m[:-1,:,:,:]) / dx[:-1,:,:,:] # m_i-1 - m_i
        h[ 1:,:,:,:] += A_avg * (m[:-1,:,:,:]-m[ 1:,:,:,:]) / dx[ 1:,:,:,:] # m_i+1 - m_i

        A_avg = 2.*A[:,1:,:]*A[:,:-1,:] / (A[:,1:,:]*dy[:,:-1,:,:] + A[:,:-1,:]*dy[:,1:,:,:])
        h[:,:-1,:,:] += A_avg * (m[:, 1:,:,:]-m[:,:-1,:,:]) / dy[:,:-1,:,:] # m_i-1 - m_i
        h[:, 1:,:,:] += A_avg * (m[:,:-1,:,:]-m[:, 1:,:,:]) / dy[:, 1:,:,:] # m_i+1 - m_i

        A_avg = 2.*A[:,:,1:]*A[:,:,:-1] / (A[:,:,1:]*dz[:,:,:-1,:] + A[:,:,:-1]*dz[:,:,1:,:])
        h[:,:,:-1,:] += A_avg * (m[:,:, 1:,:]-m[:,:,:-1,:]) / dz[:,:,:-1,:] # m_i-1 - m_i
        h[:,:, 1:,:] += A_avg * (m[:,:,:-1,:]-m[:,:, 1:,:]) / dz[:,:, 1:,:] # m_i+1 - m_i

        h *= 2. / (constants.mu_0 * Ms)
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return h
