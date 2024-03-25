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
import numpy as np

__all__ = ["LogDt", "LogMovingAverage", "LogCumSum"]

class LogDt():
    ''' Discrete difference '''
    def __init__(self, func):
        self._func = func

    def __call__(self, state):
        try:
            f0 = self._f
            t0 = self._t
        except:
            f0 = 0
            t0 = 0
        self._f = self._func(state)
        self._t = state.t
        return (self._f - f0) / (self._t - t0)

class LogMovingAverage():
    def __init__(self, func, N = 10):
        self._func = func
        self._vals = np.zeros(N)
        self._fill = np.zeros(N)

    def __call__(self, state):
        self._vals = np.roll(self._vals, -1)
        self._fill = np.roll(self._fill, -1)

        self._vals[-1] = self._func(state)
        self._fill[-1] = 1
        return self._vals.sum() / self._fill.sum()

class LogCumSum():
    def __init__(self, func):
        self._func = func
        self._val = 0.

    def __call__(self, state):
        self._val += self._func(state)
        return self._val
