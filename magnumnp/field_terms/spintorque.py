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

class SpinOrbitTorque(object):
    @timedmethod
    def h(self, state):
        p = state.material["p"].expand_as(state.m)
        h = state.material["eta_damp"] * torch.cross(state.m, p) + state.material["eta_field"] * p
        h *= -state.material["je"] * constants.hbar / (2. * constants.e * state.material["Ms"] * constants.mu_0 * state.material["d"])
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, state):
        raise NotImplemented()


class SpinTorqueZhangLi(object):
    @timedmethod
    def h(self, state):
        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        dx = [state.mesh.dx[i] for i in range(3) if state.mesh.n[i] > 1]

        j = state.j(state.t)
        jgradm = torch.einsum('...a,...ba-> ...b', j[...,dim], torch.stack(torch.gradient(state.m, spacing=dx, dim=dim), dim=-1)) # matmult

        return state.material["b"] / constants.gamma * (torch.cross(state.m, jgradm) + state.material["xi"] * jgradm)

    def E(self, state):
        raise NotImplemented()
