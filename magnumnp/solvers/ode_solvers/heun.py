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
from magnumnp.common import logging

__all__ = ["Heun"]

class Heun(object):
    def __init__(self, f, update_thermal_field, get_thermal_field, dt = 1e-15):
        self._f = f
        self._dt = dt
        self._update_thermal_field = update_thermal_field
        self._get_thermal_field = get_thermal_field

        logging.info_green("[LLGSolver] using Heun solver (dt = %g)" % dt)

    def _f_wrapper(self, state, t, m, **llg_args):
        t0 = state.t
        m0 = state.m.detach()
        state.t = t
        state.m = m
        f = self._f(state, **llg_args)
        state.t = t0
        state.m = m0
        return f

    def step(self, state, dt, **llg_args):
        f, m, t = self._f_wrapper, state.m, state.t

        self._update_thermal_field(state)
        self._get_thermal_field(state, dt)
        k1 = dt * f(state, t,      m,           **llg_args)
        self._update_thermal_field(state)
        self._get_thermal_field(state, dt)
        k2 = dt * f(state, t + dt, m + k1 + k1, **llg_args)
    
        state.m = (k1 + k2) / 2.
        state.t = t + dt
