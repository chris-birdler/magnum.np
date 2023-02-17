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
from magnumnp.common import logging, DecoratedTensor
from scipy.integrate import odeint

__all__ = ["ScipyOdeint"]

class ScipyOdeint(object):
    def __init__(self, f, hmax = 0.0, hmin = 0.0, mxordn = 12, mxords = 5, rtol = 1e-4, atol = 1e-4):
        self._f = f
        self._hmax = hmax
        self._hmin = hmin
        self._mxordn = mxordn
        self._mxords = mxords
        self._rtol = rtol
        self._atol = atol

        logging.info_green("[LLGSolver] using Scipy odeint Solver (rtol = %g, atol = %g)" % (rtol, atol))

    def _f_wrapper(self, t, m, state, llg_args):
        state.t = t * 1e-9 # scale time by 1e9 to prevent underflow error
        state.m = state.Tensor(m.reshape(state.mesh.n + (3,), order = "F"))
        f = self._f(state, **llg_args) * 1e-9
        return f.detach().cpu().numpy().flatten(order = "F")

    def step(self, state, dt, rtol = None, atol = None, **llg_args):
        m0 = state.m.numpy().reshape(-1, order = 'F')

        t1 = state.t + dt
        m1 = odeint(self._f_wrapper,
                    m0,
                    [state.t*1e9, t1*1e9],
                    args = (state, llg_args),
                    rtol = rtol or self._rtol,
                    atol = atol or self._atol,
                    tfirst = True)[1]

        state.m = state.Tensor(m1.reshape(state.mesh.n + (3,), order = "F"))
        state.t = t1
