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

__all__ = ["SpinOrbitTorque", "SpinTorqueZhangLi"]

#TODO: generalize interface (ZhangLi uses state.j, SOT uses material["je"]
class SpinOrbitTorque(object):
    r"""
    General spin torque contributions can be described by the following field

    .. math::
        \vec{h}^\text{sot} = -\frac{j_e \hbar}{2 e \mu_0 M_s} \left[\eta_\text{damp} \, \vec{m} \times \vec{p} + \eta_\text{field} \, \vec{p} \right],

    with the current density :math:`j_e`, the reduced Planck constant :math:`\hbar`,
    the elementary charge :math:`e`, and the polarization of the electrons :math:`\vec{p}`.
    :math:`\eta_\text{damp}` and :math:`\eta_\text{field}` are material parameters which
    describe the amplitude of damping- and field-like torque.

    In case of Spin-Orbit-Torqe (SOT) :math:`\eta_\text{field}` and :math:`\eta_\text{damp}` are constant material parameters.
    """
    @timedmethod
    def h(self, state):
        p = state.material["p"].expand_as(state.m)
        h = state.material["eta_damp"] * torch.linalg.cross(state.m, p) + state.material["eta_field"] * p
        h *= -state.material["je"] * constants.hbar / (2. * constants.e * state.material["Ms"] * constants.mu_0 * state.material["d"])
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, state):
        raise NotImplemented()


class SpinTorqueZhangLi(object):
    r"""
    Zhang Lie spin torque contributions can be described by the following field

    .. math::
        \vec{h}^\text{stt,zl} = \frac{b}{\gamma} \left[\vec{m} \times (\vec{j}_e \cdot \nabla) \vec{m} + \xi \; (\vec{j}_e \cdot \nabla) \vec{m} \right],

    with the reduced gyromagnetic ratio :math:`\gamma`, the degree of nonadiabacity :math:`\xi`. :math:`b` is the polarization rate of the conducting electrons and can be written as

    .. math::
        b = \frac{\beta \mu_B}{e M_s (1+\xi^2)},

    with the Bohr magneton :math:`\mu_B`, and the dimensionless polarization rate :math:`\beta`.
    """
    @timedmethod
    def h(self, state):
        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        dx = [state.mesh.dx[i] for i in range(3) if state.mesh.n[i] > 1]

        j = state.j(state.t)
        jgradm = torch.einsum('...a,...ba-> ...b', j[...,dim], torch.stack(torch.gradient(state.m, spacing=dx, dim=dim), dim=-1)) # matmult

        return state.material["b"] / constants.gamma * (torch.linalg.cross(state.m, jgradm) + state.material["xi"] * jgradm)

    def E(self, state):
        raise NotImplemented()
