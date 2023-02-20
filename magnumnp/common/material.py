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

import torch
from . import DecoratedTensor

__all__ = ["Material"]

class Material(dict):
    def __init__(self, state):
        self._state = state

    def _convert(self, value):
        ''' convert arbitrary input to tensor-fields '''
        value = self._state.Tensor(value)
        if len(value.shape) == 0: # convert dim=0 tensor into dim=1 tensor
            value = value.reshape(1)
        if len(value.shape) < 3: # expand homogeneous material to [nx,ny,nz,...] tensor-field
            shape = value.shape
            value = value.reshape((1,1,1) + tuple(shape))
            value = value.expand(self._state.mesh.n + tuple(shape))
        elif len(value.shape) == 3: # scalar-field should have dimension [nx,ny,nz,1]
            value = value.unsqueeze(-1)
        else: # otherwise assume the dimention is correct!
            pass
        return value

    def __getitem__(self, key):
        return super().__getitem__(key)(self._state.t)

    def __setitem__(self, key, value):
        if callable(value) and not isinstance(value, DecoratedTensor):
            super().__setitem__(key, lambda t: self._convert(value(t)))
        else:
            super().__setitem__(key, self._convert(value))
