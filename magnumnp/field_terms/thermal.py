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
from .field_terms import FieldTerm

__all__ = ["ThermalField"]

class ThermalField(FieldTerm):
    r"""
    """
    parameters = ["T"]
    def __init__(self, domain=None, **kwargs):
        self._step = None
        super().__init__(**kwargs)

    @timedmethod
    def h(self, state):
        if state._step != self._step: # update random field
            self._sigma = state._normal(0., 1., size = state.m.shape)
            self._step = state._step

        alpha = state.material["alpha"].torch_tensor
        Ms = state.material["Ms"].torch_tensor
        return self._h(self._sigma, alpha, Ms, state)

    @torch.compile
    def _h(self, sigma, alpha, Ms, state):
        return sigma * torch.sqrt(2. * alpha * constants.kb * state.T / (constants.mu_0 * Ms * constants.gamma * state.cell_volumes * state._dt))
