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

    def _f_wrapper(self, t, m, state, **llg_args):
        state.t = t * 1e-9 # scale time by 1e9 to prevent underflow error
        state.m = m
        return self._f(state, **llg_args) * 1e-9

    def step(self, state, dt, rtol = None, atol = None, **llg_args):
        t1 = state.t + dt
        res = odeint(lambda t, m: self._f_wrapper(t, m, state, **llg_args),
                     state.m,
                     state.Tensor([state.t*1e9, t1*1e9]),
                     method = self._method,
                     rtol = rtol or self._rtol,
                     atol = atol or self._atol,
                     options = self._options) # TODO: reuse solver object?
        state.m = state.Tensor(res[1])
        state.t = t1

class TorchDiffEqAdjoint(object):
    def __init__(self, f, adjoint_parameters, method = "dopri5", rtol = 1e-5, atol = 1e-5, options = {}):
        self._f = f
        self._adjoint_parameters = adjoint_parameters
        self._method = method
        self._rtol = rtol
        self._atol = atol
        self._options = options
        logging.info_green("[LLGSolver] using TorchDiffEq adjoint solver (method = '%s', rtol = %g, atol = %g)" % (method, rtol, atol))

    def _f_wrapper(self, t, m, state, **llg_args):
        state.t = t * 1e-9 # scale time by 1e9 to prevent underflow error
        state.m = state.Tensor(m)
        return self._f(state, **llg_args) * 1e-9

    def step(self, state, dt, rtol = None, atol = None, **llg_args):
        t1 = state.t + dt
        res = odeint_adjoint(lambda t, m: self._f_wrapper(t, m, state, **llg_args),
                     state.m,
                     state.Tensor([state.t*1e9, t1*1e9]),
                     method = self._method,
                     rtol = rtol or self._rtol,
                     atol = atol or self._atol,
                     adjoint_params = self._adjoint_parameters,
                     options = self._options) # TODO: reuse solver object?
        state.m = state.Tensor(res[1])
        state.t = t1
