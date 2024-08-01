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

    def _f_wrapper(self, t, x, kwargs):
        x_torch = torch.tensor(x.reshape(kwargs["state"].mesh.n + (-1,), order = "F"))
        f = self._f(t * 1e-9, x_torch, **kwargs) * 1e-9 # scale time by 1e9 to prevent underflow error
        return f.detach().cpu().numpy().flatten(order = "F")

    def step(self, t, x_torch, dt, rtol = None, atol = None, **kwargs):
        x0 = x_torch.detach().cpu().numpy().reshape(-1, order = 'F')
        t1 = t + dt
        x1 = odeint(self._f_wrapper,
                    x0,
                    [(t*1e9).detach().cpu(), (t1*1e9).detach().cpu()],
                    args = (kwargs,),
                    rtol = rtol or self._rtol,
                    atol = atol or self._atol,
                    tfirst = True)[1]

        return t1, torch.tensor(x1.reshape(x_torch.shape, order = "F"))

    def solve(self, tt, x_torch, rtol = None, atol = None, **kwargs):
        x0 = x_torch.detach().cpu().numpy().reshape(-1, order = 'F')
        res = odeint(self._f_wrapper,
                     x0,
                     tt.detach().cpu().numpy()*1e9,
                     args = (kwargs,),
                     rtol = rtol or self._rtol,
                     atol = atol or self._atol,
                     tfirst = True)

        return torch.tensor(res.reshape(tt.shape + x_torch.shape, order = "F"))
