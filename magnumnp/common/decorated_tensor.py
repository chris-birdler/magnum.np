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

__all__ = ["DecoratedTensor"]

class DecoratedTensor(torch.Tensor):
    @staticmethod
    def __new__(cls, x, *args, **kwargs): # TODO: is this needed?
        return super().__new__(cls, x, *args, **kwargs)

    def avg(self, dim=(0,1,2)):
        if self.dim() <= 1: # e.g. [0,0,1]
            return self
        elif self.dim() == 2: # state.m[domain]
            return self.mean(dim=0)
        else:
            return self.mean(dim=dim)

    def average(self, dim=(0,1,2)):
        return self.avg(dim)

    def normalize(self):
        self /= torch.linalg.norm(self, dim = -1, keepdim = True)
        self[...] = torch.nan_to_num(self, posinf=0, neginf=0)
        return self

    def __call__(self, t):
        return self

    def __getitem__(self, key):
        if self.dim() == 0 or (self.dim() == 1 and self.shape[0] == 1 and key != 0):
            return self
        else:
            return super().__getitem__(key)
