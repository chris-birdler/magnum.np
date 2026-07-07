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

__all__ = ["CoWTensor"]


class CoWTensor(torch.Tensor):
    r"""
    Copy-on-write tensor used for homogeneous material parameters.

    Homogeneous parameters are stored as broadcasted (stride-0) views of a
    single value, so they occupy O(1) instead of O(nx*ny*nz) memory. In all
    operations a :class:`CoWTensor` behaves exactly like a plain tensor
    (torch-function dispatch is disabled, like for :class:`torch.nn.Parameter`),
    but an indexed write such as

    .. code::
        state.material["Ms"][domain] = 0.

    first materializes a dense private copy. The storage swap happens in
    place (via :code:`Tensor.set_`), so every existing reference to the
    parameter observes the write.

    Note that in-place operations other than :code:`__setitem__` (e.g.
    :code:`mul_` or :code:`+=`) on a still-expanded tensor raise the usual
    "multiple elements refer to a single memory location" error; call
    :code:`materialize()` first in that case.
    """
    __torch_function__ = torch._C._disabled_torch_function_impl

    @staticmethod
    def wrap(tensor):
        t = tensor.as_subclass(CoWTensor)
        t._materialized = False
        return t

    def materialize(self):
        if not getattr(self, "_materialized", True):
            self.set_(self.clone()) # clone of the expanded view -> dense tensor; same python object
            self._materialized = True
        return self

    def __setitem__(self, idx, value):
        self.materialize()
        super().__setitem__(idx, value)
