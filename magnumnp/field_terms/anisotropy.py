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
from .field_terms import FieldTerm, LinearFieldTerm
from torch import sin, cos

__all__ = ["UniaxialAnisotropyField", "CubicAnisotropyField"]

class UniaxialAnisotropyField(LinearFieldTerm):
    r"""
    Uniaxial Anisotropy Field:

    .. math::
        \vec{h}^\text{u} = \frac{2 K_\text{u}}{\mu_0 \, M_s} \; \vec{e}_\text{u} \; (\vec{e}_\text{u} \cdot \vec{m}),

    with the anisotropy constant :math:`K_\text{u}` given in units of :math:`\text{J/m}^3`.

    :param Ku: Name of the material parameter for the anisotropy constant :math:`K_\text{u}`, defaults to "Ku"
    :type Ku: str, optional
    :param Ku_axis: Name of the material parameter for the anisotropy axis :math:`\vec{e}_\text{u}`, defaults to "Ku_axis"
    :type Ku_axis: str, optional
    """
    parameters = ["Ku", "Ku_axis"]

    @timedmethod
    def h(self, state):
        Ku = state.material[self.Ku]
        Ku_axis = state.material[self.Ku_axis]

        h = 2. * Ku * Ku_axis / (constants.mu_0 * state.material["Ms"]) * torch.sum(Ku_axis * state.m, dim=3, keepdim=True)
        return torch.nan_to_num(h, posinf=0, neginf=0)


class CubicAnisotropyField(FieldTerm):
    r"""
    Cubic Anisotropy Field:

    .. math::
        \vec{h}^\text{c} = -\frac{2 K_\text{c1}}{\mu_0 \, M_s} \; \begin{pmatrix} m_1 \, m_2^2 + m_1 \, m_3^2 \\ m_2 \, m_3^2 + m_2 \, m_1^2 \\ m_3 \, m_1^2 + m_3 \, m_2^2\end{pmatrix}
                           -\frac{2 K_\text{c2}}{\mu_0 \, M_s} \; \begin{pmatrix} m_1 \, m_2^2 \, m_3^2 \\ m_1^2 \, m_2 \, m_3^2 \\ m_1^2 \, m_2^2 \, m_3\end{pmatrix},

    with the anisotropy constants :math:`K_{c1}` and :math:`K_{c2}` given in units of :math:`\text{J/m}^3`.
    If Euler angles :math:`\alpha`, :math:`\beta` and :math:`\gamma` are provided, the cubic anisotropy axes are rotated such that the effective magnetization components :math:`m_i` with :math:`i \in \{x,y,z\}` read

    .. math::
      m_i = (\mat{A} \vec{e}_i) \cdot \vec{m}

    with

    .. math::
      \small
      \mat{A} = \begin{pmatrix}
            \cos(\alpha) \cos(\gamma) - \cos(\beta) \sin(\alpha) \sin(\gamma) & -\cos(\beta) \cos(\gamma) \sin(\alpha) - \cos(\alpha) \sin(\gamma) &  \sin(\alpha) \sin(\beta)\\
            \cos(\gamma) \sin(\alpha) + \cos(\alpha) \cos(\beta) \sin(\gamma) &  \cos(\alpha) \cos(\beta) \cos(\gamma) - \sin(\alpha) \sin(\gamma) & -\cos(\alpha) \sin(\beta)\\
            \sin(\beta) \sin(\gamma) & \cos(\gamma) \sin(\beta) & \cos(\beta)
            \end{pmatrix}

    :param Kc1: Name of the material parameter for the anisotropy constant :math:`K_\text{c1}`, defaults to "Kc1"
    :type Kc1: str, optional
    :param Kc2: Name of the material parameter for the anisotropy constant :math:`K_\text{c2}`, defaults to "Kc2"
    :type Kc2: str, optional
    :param Kc_alpha: Euler angle :math:`\alpha` for the rotation of the anisotropy axes
    :type Kc_alpha: str, optional
    :param Kc_beta: Euler angle :math:`\beta` for the rotation of the anisotropy axes
    :type Kc_beta: str, optional
    :param Kc_gamma: Euler angle :math:`\gamma` for the rotation of the anisotropy axes
    :type Kc_gamma: str, optional
    """
    parameters = ["Kc1", "Kc2", "Kc_alpha", "Kc_beta", "Kc_gamma"]

    def _R(self, state):
        a = state.material[self.Kc_alpha]
        b = state.material[self.Kc_beta]
        g = state.material[self.Kc_gamma]

        R = torch.stack([torch.concat([cos(a)*cos(g) - cos(b)*sin(a)*sin(g), -cos(b)*cos(g)*sin(a) - cos(a)*sin(g),  sin(a)*sin(b)], dim = -1),
                         torch.concat([cos(g)*sin(a) + cos(a)*cos(b)*sin(g),  cos(a)*cos(b)*cos(g) - sin(a)*sin(g), -cos(a)*sin(b)], dim = -1),
                         torch.concat([sin(b)*sin(g),                         cos(g)*sin(b),                         cos(b)       ], dim = -1)], dim = -1)
        return R

    @timedmethod
    def h(self, state):
        Kc1 = state.material[self.Kc1]
        Kc2 = state.material[self.Kc2]

        R = self._R(state)
        mx, my, mz = torch.einsum('...a, ...ab-> ...b', state.m, R).unbind(dim=-1) # matmult

        h =  2. * Kc1 * torch.stack([mx * (my**2 + mz**2), my*(mz**2 + mx**2), mz*(mx**2 + my**2)], dim = -1) + \
             2. * Kc2 * torch.stack([mx * my**2. * mz**2., mx**2. * my * mz**2., mx**2. * my**2. * mz], dim = -1)
        h = torch.einsum('...a, ...ba-> ...b', h, R) # matmult transpose
        return torch.nan_to_num(-1./constants.mu_0/state.material["Ms"] * h, posinf=0, neginf=0)

    def E(self, state):
        R = self._R(state)
        mx, my, mz = torch.einsum('...a, ...ab-> ...b', state.m, R).unbind(dim=-1) # matmult

        return (state.material[self.Kc1] * (mx**2 * my**2 + mx**2 * mz**2 + my**2 * mz**2).sum() +
                state.material[self.Kc2] * (mx**2 * my**2 * mz**2).sum()) * state.mesh.cell_volume
