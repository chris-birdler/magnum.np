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

from magnumnp.common import logging

__all__ = ["Mesh"]

class Mesh(object):
    def __init__(self, n, dx, origin=(0,0,0)):
        self.n = tuple(n)
        self.dx = tuple(dx)
        self.origin = tuple(origin)

    @property
    def cell_volume(self):
        return self.dx[0] * self.dx[1] * self.dx[2]

    @property
    def volume(self):
        return self.n[0] * self.n[1] * self.n[2] * self.cell_volume

    def __str__(self):
        str_dx = ["%g" % dx if isinstance(dx, (float, int)) else "XX" for dx in self.dx]
        return "%dx%dx%d (size= %s x %s x %s)" % (*self.n, *str_dx)
