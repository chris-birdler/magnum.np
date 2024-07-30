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

import numpy as np
import torch
from scipy import interpolate

__all__ = ["TimeInterpolator"]

class TimeInterpolator(object):
    def __init__(self, state, points, order = 1):
        self._tp = torch.tensor(list(points.keys()))
        self._fp = [state.convert_tensorfield(f) for f in points.values()]
        self._state = state
        self._order = order

    def __call__(self, state):
        i = torch.searchsorted(self._tp, state.t) # upper index
        tp = self._tp
        fp = self._fp
        if self._order == 0:
            i = torch.clamp(i, min=1, max=len(self._tp)) # extrapolate only on lower bound
            return fp[i-1]
        elif self._order == 1:
            i = torch.clamp(i, min=1, max=len(self._tp)-1) # extrapolate on bounds
            return fp[i-1] + (state.t-tp[i-1]) / (tp[i]-tp[i-1]) * (fp[i] - fp[i-1])
        else:
            raise ValueError("Order must be 0 or 1(default)")

    @property
    def final_time(self):
        return self._tp[-1]
