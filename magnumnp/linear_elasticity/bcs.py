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

__all__ = ["BC", "PlaneBC", "Plane"]

class BC:
    def __init__(self, mask, condition):
        self._mask = mask
        
        if callable(condition):
            self._condition = condition
        elif isinstance(condition, torch.Tensor):
            self._condition = lambda state : condition
        else:
            raise Exception("BC condition has invalid type " + type(condition) + "!/n condition has to be either torch.Tensor or a function of the state that returns torch.Tensor.")

    @property
    def mask(self):
        return self._mask
    
    def condition(self, state):
        return self._condition(state)

class PlaneBC:
    def __init__(self, plane, condition):
        self.plane = plane

        if callable(condition):
            self._condition = condition
        elif isinstance(condition, torch.Tensor):
            self._condition = lambda state : condition
        else:
            raise Exception("BC condition has invalid type " + type(condition) + "!/n condition has to be either torch.Tensor or a function of the state that returns torch.Tensor.")

    def condition(self, state):
        return self._condition(state)
    
class Plane:
    def __init__(self, 
                 normal_dim, # (int) direction of the normal vector (0=x, 1=y, 2=z)
                 normal_position_index, # (int) positional index in the normal dimension
                 normal_sign, # (int), 1 or -1, direction in which the normal vector points, points away from the interiour
                 trans_lim1=[0,None], # bounds of the transversal direction with lowest index
                 trans_lim2=[0,None]): # bounds of the transversal direction with the higgest index
        
        assert isinstance(normal_dim,int)
        self.dim = normal_dim

        assert isinstance(normal_position_index,int)
        self.pos = normal_position_index

        assert isinstance(normal_sign, int)
        if (self.pos == 0):
            assert normal_sign == -1
        if (self.pos == -1):
            assert normal_sign == 1
            
        self.sign = normal_sign

        assert isinstance(trans_lim1, list) or isinstance(trans_lim1, tuple) 
        assert isinstance(trans_lim1[0], int) or isinstance(trans_lim1[0], None)
        self.trans_slice1 = slice(trans_lim1[0], trans_lim1[1])

        assert isinstance(trans_lim2, list) or isinstance(trans_lim2, tuple) 
        assert isinstance(trans_lim2[0], int) or isinstance(trans_lim2[0], None)
        self.trans_slice2 = slice(trans_lim2[0], trans_lim2[1])

        if self.dim == 0:
            self.trans_dim1 = 1
            self.trans_dim2 = 2
        elif self.dim == 1:
            self.trans_dim1 = 0
            self.trans_dim2 = 2
        elif self.dim == 2:
            self.trans_dim1 = 0
            self.trans_dim2 = 1

        self.slices = [None,None,None]
        self.slices[self.dim] = self.pos
        self.slices[self.trans_dim1] = self.trans_slice1
        self.slices[self.trans_dim2] = self.trans_slice2
        self.slices = tuple(self.slices)

    def get_mask(self, state):
        mask = torch.zeros(state.mesh.n, dtype=torch.bool)
        mask[self.slices] = True
        return mask