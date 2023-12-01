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

__all__ = ["ExternalField"]

class ExternalField(object):
    r"""
    External Field

    :param h: External Field
    :type h: list or tuple or :class:`Tensor` or function

    :Examples:

    .. code::

        # homogenious, constant field
        external = ExternalField([Hx, 0, 0])

        # homogenious, time-dependent field
        external = ExternalField(lambda t: [Hx*t, 0, 0])

        # inhomogenious, constant field
        x, y, z = SpatialCoordinate(state)
        h = torch.stack([x,y,z], dim=-1)
        external = ExternalField(h)
    """
    def __init__(self, h):
        self._h = h

    @timedmethod
    @torch.compile
    def h(self, state):
        if callable(self._h):
            return state.Tensor(self._h(state.t))
        else:
            return state.Tensor(self._h)

    def __setattr__(self, name, value):
        if name == "h":
            self._h = value
        else:
            super().__setattr__(name, value)

    @torch.compile
    def E(self, state, domain = Ellipsis):
        E = - constants.mu_0 * state.material["Ms"] * state.m * self.h(state) * state.cell_volumes
        return E[domain].sum()
