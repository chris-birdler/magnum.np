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
        external = ExternalField(lambda state: [Hx*state.t, 0, 0])

        # inhomogenious, constant field
        x, y, z = SpatialCoordinate(state)
        external = ExternalField(Expression([x,y,z]))
    """
    def __init__(self, h = None):
        self._h_cache = None
        self._h_is_const = False
        if h != None:
            self.__setattr__("h", h)

    @timedmethod
    def h(self, state):
        if self._h_is_const:
            # cache the converted field for constant h; converting on every call
            # would expand+copy a full [nx,ny,nz,3] tensor per evaluation.
            # Callers must not modify the returned tensor (accumulate_h complies).
            if self._h_cache is None:
                self._h_cache = state.convert_tensorfield(self._h(state))
            return self._h_cache
        return state.convert_tensorfield(self._h(state))

    def __setattr__(self, name, value):
        if name == "h":
            super().__setattr__("_h_is_const", not callable(value))
            super().__setattr__("_h_cache", None) # invalidate cache when h is reassigned
            if callable(value):
                self._h = value
            else:
                self._h = lambda state: value
        else:
            super().__setattr__(name, value)

    def E(self, state, domain = Ellipsis):
        E = - constants.mu_0 * state.material["Ms"] * state.m * self.h(state) * state.mesh.cell_volumes
        return E[domain].sum()
