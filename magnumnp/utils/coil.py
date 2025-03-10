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

from magnumnp.utils.logging_helpers import LogDt
from magnumnp.field_terms import OerstedField, VectorPotential

__all__ = ["Coil"]

class Coil(object):
    def __init__(self, state, I0, h_oersted0 = None, A0 = None):
        self.I0 = I0
        self.h_oersted0 = h_oersted0 or OerstedField().h(state)
        self.A0 = A0 or VectorPotential().A(state)

        # calculate induced voltage by numerical time-derivative
        self.Ui_m = LogDt(self.psi_m)
        self.Ui_j = LogDt(self.psi_j)
        self.Ui = LogDt(self.psi)

    def psi_m(self, state):
        return (state.material["Ms"] * state.m * self.h_oersted0).sum() / self.I0 * constants.mu_0 * state.mesh.cell_volumes

    def psi_j(self, state):
        return (self.A0*state.j).sum() * constants.mu_0 * state.mesh.cell_volumes

    def psi(self, state):
        return self.psi_m(state) + self.psi_j(state)

    def L(self, state):
        return self.psi_j(state) / self.I0

# state.j = ...
# h_oersted = OerstedField().h(state)
# coil = Coil(h_oersted, I=j0*A)
# 
# logger = Logger("data", ["t", "m", coil.U_m, coil.U_j, coil.U])

