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
    def __new__(self, x, cell_volumes, *args, **kwargs): # TODO: is this needed?
        self._cell_volumes = cell_volumes
        return super().__new__(self, x, *args, **kwargs)

    def avg(self, dim=(0,1,2)):
        if self.dim() <= 1: # e.g. [0,0,1]
            return self
        elif self.dim() == 2: # state.m[domain]
            return (self * self._cell_volumes.unsqueeze(-1)).sum(dim=0) / self._cell_volumes.sum(dim=0)
        elif self.dim() == 3: # [nx,ny,nz]
            return (self * self._cell_volumes).sum(dim=dim) / self._cell_volumes.sum()
        else:                 # [nx,ny,nz,...]
            return (self * self._cell_volumes.unsqueeze(-1)).sum(dim=dim) / self._cell_volumes.sum()

    def average(self, dim=(0,1,2)):
        return self.avg(dim)

    def normalize(self):
        self /= torch.linalg.norm(self, dim = -1, keepdim = True)
        self[...] = torch.nan_to_num(self, posinf=0, neginf=0)
        return self

    def __call__(self, t):
        return self

    def __getitem__(self, idx):
        item = super().__getitem__(idx)

        # update cell_volumes
        if isinstance(idx, tuple): # apply slicing
            if idx[0] == Ellipsis:
                item._cell_volumes = self._cell_volumes
            else:
                item._cell_volumes = self._cell_volumes[idx[:3]]
        elif isinstance(idx, int):
            return item
        else: # apply fancy indexing
            item._cell_volumes = self._cell_volumes.__getitem__(idx)

        return item

    @property
    def torch_tensor(self):
        ''' return original tensor (in order to apply pytorch.compile) '''
        return self.as_subclass(torch.Tensor)
