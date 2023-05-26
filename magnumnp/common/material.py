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

    def __getitem__(self, key):
        return super().__getitem__(key)(self._state.t)

    def __setitem__(self, key, value):
        if callable(value) and not isinstance(value, DecoratedTensor):
            super().__setitem__(key, lambda t: self._state.convert_tensorfield(value(t)))
        else:
            super().__setitem__(key, self._state.convert_tensorfield(value))

    def set(self, material, domain=None):
        r"""
        Setting several material parameters at once

        :param materials: material that should be set
        :type materials:  :class:`Material`

        :param domains: domain where the materials should be set
        :type domains:  :class:`torch.Tensor`, optional

        :Example:

            .. code::

            # set material everywhere
            state.material.set(material0)

            # set material in certain domain
            state.material.set(material1, domain1)
        """
        for key, value in material.items():
            if domain == None:
                self[key] = value
            else:
                if key not in self.keys():
                    self[key] = 0.
                self[key][domain] = value
