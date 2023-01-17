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
    def __init__(self, state, points):
        self._tp = state.Tensor(list(points.keys()))
        self._fp = state._tensor(list(points.values()))
        self._state = state

    def __call__(self, t):
        i = torch.searchsorted(self._tp, t) # upper index
        i = torch.clamp(i, min=1, max=len(self._tp)-1) # extrapolate on bounds
        tp = self._tp
        fp = self._fp
        return fp[i-1] + (t-tp[i-1]) / (tp[i]-tp[i-1]) * (fp[i] - fp[i-1])

    @property
    def final_time(self):
        return self._tp[-1]
