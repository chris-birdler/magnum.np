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
from scipy.integrate import ode

__all__ = ["ScipyODE"]

class ScipyODE(object):
    def __init__(self, f, name = "vode", method = "BDF", rtol = 1e-5, atol = 1e-5):
        self._f = f
        self._solver = ode(self._f_wrapper)
        self._solver.set_integrator(name = name,
                                    method = method,
                                    rtol = rtol,
                                    atol = atol)
        self._initialized = False
        logging.info_green("[LLGSolver] using Scipy ODE solver '%s' (method = '%s', rtol = %g, atol = %g)" % (name, method, rtol, atol))

    def _f_wrapper(self, t, m, state, llg_args):
        state.t = t
        state.m = state.Tensor(m.reshape(state.mesh.n + (3,), order = "F"))
        f = self._f(state, **llg_args)
        return f.detach().cpu().numpy().flatten(order = "F")

    def step(self, state, dt, rtol = None, atol = None, **llg_args):
        if not self._initialized:
            m = state.m.numpy().reshape(-1, order = 'F')
            self._solver.set_initial_value(m, state.t.item())
            self._initialized = True

        self._solver.set_f_params(state, llg_args)

        t1 = self._solver.t + dt
        m1 = self._solver.integrate(t1)
        if not self._solver.successful():
            logging.warning("[LLGSolver] Scipy ODE solver: integration not successful!")

        state.m = state.Tensor(m1.reshape(state.mesh.n + (3,), order = "F"))
        state.t = self._solver.t
