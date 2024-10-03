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
    #@torch.compile
    def h(self, state):
        if state._step != self._step: # update random field
            #self._sigma = torch.normal(0., 1., size = state.m.shape) # creates tensor on CPU by default
            self._sigma = torch.normal(0., 1., size = state.m.shape).to(state.device)
            self._step = state._step

        h = self._sigma * torch.sqrt(2. * state.material["alpha"]  * constants.kb * state.T / (constants.mu_0 * state.material["Ms"] * constants.gamma * state.mesh.cell_volumes * state._dt))
        return h.nan_to_num(posinf=0, neginf=0)
