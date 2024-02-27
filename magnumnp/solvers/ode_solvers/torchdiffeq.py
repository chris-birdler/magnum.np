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
from torchdiffeq import odeint, odeint_adjoint

__all__ = ["TorchDiffEq", "TorchDiffEqAdjoint"]

class TorchDiffEq(object):
    def __init__(self, f, method = "dopri5", rtol = 1e-5, atol = 1e-5, options = {}):
        self._f = f
        self._method = method
        self._rtol = rtol
        self._atol = atol
        self._options = options
        logging.info_green("[LLGSolver] using TorchDiffEq solver (method = '%s', rtol = %g, atol = %g)" % (method, rtol, atol))

    def _f_wrapper(self, t, x, **kwargs): # TODO: move scaling to dm?
        return self._f(t * 1e-9, x, **kwargs) * 1e-9 # scale time by 1e9 to prevent underflow error

    def step(self, t, x, dt, rtol = None, atol = None, **kwargs):
        t1 = t + dt
        res = odeint(lambda t, x: self._f_wrapper(t, x, **kwargs),
                     x,
                     torch.tensor([t*1e9, t1*1e9]),
                     method = self._method,
                     rtol = rtol or self._rtol,
                     atol = atol or self._atol,
                     options = self._options) # TODO: reuse solver object?
        return t1, res[1]

class TorchDiffEqAdjoint(object):
    def __init__(self, f, adjoint_parameters, method = "dopri5", rtol = 1e-5, atol = 1e-5, options = {}):
        self._f = f
        self._adjoint_parameters = adjoint_parameters
        self._method = method
        self._rtol = rtol
        self._atol = atol
        self._options = options
        logging.info_green("[LLGSolver] using TorchDiffEq adjoint solver (method = '%s', rtol = %g, atol = %g)" % (method, rtol, atol))

    def _f_wrapper(self, t, x, **kwargs): # TODO: move scaling to dm?
        return self._f(t * 1e-9, x, **kwargs) * 1e-9 # scale time by 1e9 to prevent underflow error

    def step(self, t, x, dt, rtol = None, atol = None, **kwargs):
        t1 = t + dt
        res = odeint_adjoint(lambda t, m: self._f_wrapper(t, x, **kwargs),
                     x,
                     torch.tensor([t*1e9, t1*1e9]),
                     method = self._method,
                     rtol = rtol or self._rtol,
                     atol = atol or self._atol,
                     adjoint_params = self._adjoint_parameters,
                     options = self._options) # TODO: reuse solver object?
        return t1, res[1]
