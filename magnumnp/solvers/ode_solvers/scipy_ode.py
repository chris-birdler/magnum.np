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

    def _f_wrapper(self, t, x, kwargs): # TODO: should x be a 1D array?
        x_torch = torch.tensor(x.reshape(kwargs["state"].mesh.n + (-1,), order = "F"))
        f = self._f(t, x_torch, **kwargs)
        return f.detach().cpu().numpy().flatten(order = "F")

    def step(self, t, x_torch, dt, rtol = None, atol = None, **kwargs):
        if not self._initialized:
            x = x_torch.detach().cpu().numpy().reshape(-1, order = 'F')
            self._solver.set_initial_value(x, t)
            self._initialized = True

        self._solver.set_f_params(kwargs)

        t1 = self._solver.t + dt
        x1 = self._solver.integrate(t1)
        if not self._solver.successful():
            logging.warning("[LLGSolver] Scipy ODE solver: integration not successful!")

        return self._solver.t, torch.tensor(x1.reshape(x_torch.shape, order = "F"))
