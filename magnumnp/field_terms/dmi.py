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
import numpy as np

__all__ = ["InterfaceDMIField", "BulkDMIField", "D2dDMIField"]

class DMIField(LinearFieldTerm):
    parameters = ["D"]
    def __init__(self, dmi_vector, **kwargs):
        self._dmi_vector = dmi_vector
        super().__init__(**kwargs)

    @timedmethod
    def h(self, state):
        h = state._zeros(state.mesh.n + (3,))
        D = state.material[self.D]

        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        for dim in range(3):
            v = state.Tensor(self._dmi_vector[dim]).expand(state.m[next].shape)
            D_avg = torch.where(D[next]*D[current] < 0,
                                torch.sqrt(torch.sqrt(-D[next]*D[current])*torch.abs(D[next]+D[current]) / 2.),
                                2.*D[next]*D[current]/(D[next]+D[current]))
            h[current] += D_avg * torch.linalg.cross(v, state.m[next]   ) / (2.*state.mesh.dx[dim])
            h[next]    -= D_avg * torch.linalg.cross(v, state.m[current]) / (2.*state.mesh.dx[dim])

            # rotate dimension
            current = current[-1:] + current[:-1]
            next = next[-1:] + next[:-1]

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)


class InterfaceDMIField(DMIField):
    def __init__(self, Di = "Di"):
        dmi_vector = [[ 0, 1, 0], # x
                      [-1, 0, 0], # y
                      [ 0, 0, 0]] # z
        super().__init__(dmi_vector, D = Di)


class BulkDMIField(DMIField):
    def __init__(self, Db = "Db"):
        dmi_vector = [[1, 0, 0], # x
                      [0, 1, 0], # y
                      [0, 0, 1]] # z
        super().__init__(dmi_vector, D = Db)


class D2dDMIField(DMIField):
    def __init__(self, DD2d = "DD2d"):
        dmi_vector = [[-1, 0, 0], # x
                      [ 0, 1, 0], # y
                      [ 0, 0, 0]] # z
        super().__init__(dmi_vector, D = DD2d)
