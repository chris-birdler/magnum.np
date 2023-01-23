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

from magnumnp.common import constants
import torch

__all__ = ["FieldTerm", "LinearFieldTerm"]

class FieldTerm(object):
    parameters = []

    def __init__(self, **kwargs):
        unknown_params = set(kwargs.keys()) - set(self.parameters)
        if unknown_params:
            raise Warning("Got unknown parameters '%s'. Ignoring!" % unknown_params)
        params = {key:key for key in self.parameters}
        params.update(kwargs)
        for key, value in params.items():
            setattr(self, key, value)

class LinearFieldTerm(FieldTerm):
    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))
