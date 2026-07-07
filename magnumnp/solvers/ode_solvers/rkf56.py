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

import math
import torch
from magnumnp.common import logging, normalize

__all__ = ["RKF56"]

#Runge-Kutta-Fehlberg method with stepsize control
class RKF56(object):
    def __init__(self, f, dt = 1e-15, atol = 1e-5, rtol = 1e-5):
        self._f = f
        self._dt = dt
        self._order = 5

        # Numerical Recipies 3rd Edition suggests these values:
        self._headroom = 0.9
        self._maxstep = 1e-10
        self._minscale = 0.2
        self._maxscale = 10.
        self._atol = atol
        self._rtol = rtol # currently rtol is ignored since m is normalized!
        logging.info_green("[LLGSolver] using RKF56 solver (atol = %g)" % atol)

    def _f_wrapper(self, state, t, m, **llg_args):
        t0 = state.t
        m0 = state.m.detach()
        state.t = t
        state.m = m
        f = self._f(state, **llg_args)
        state.t = t0
        # m0 shares storage with the stage base held by _try_step, so the
        # normalization must act on a copy - normalizing in place would
        # silently change the base point of the remaining RK stages
        state.m = normalize(m0.clone())
        return f

    def _try_step(self, state, **llg_args):
        f, m, t, dt = self._f_wrapper, state.m, state.t, self._dt
        state._dt = dt  # update current dt in state used by thermal field class
        k1 = dt * f(state, t,               m, **llg_args)
        k2 = dt * f(state, t +  1./ 6.*dt,  m +    1. / 6.*k1, **llg_args)
        k3 = dt * f(state, t +  4./ 15.*dt, m +    4. / 75.*k1    + 16./75.*k2, **llg_args)
        k4 = dt * f(state, t +  2./ 3.*dt,  m +    5. / 6.*k1     - 8./3.*k2    + 5./2.*k3, **llg_args)
        k5 = dt * f(state, t +  4./ 5.*dt,  m -    8. / 5.*k1     + 144./25.*k2 - 4.*k3        + 16./25.*k4, **llg_args)
        k6 = dt * f(state, t +  1.*dt,      m +    361. / 320.*k1 - 18/5.*k2    + 407./128.*k3 - 11./80.*k4  + 55./128.*k5, **llg_args)
        k7 = dt * f(state, t,               m -    11. / 640.*k1                + 11./256.*k3  - 11./160.*k4 + 11./256.*k5, **llg_args)
        k8 = dt * f(state, t +  1.*dt,      m +    93. / 640.*k1  - 18./5.*k2   + 803./256.*k3 - 11./160.*k4 + 99./256.*k5 + 1.*k7, **llg_args)

        dm =             31./384.*k1 + 1125./2816.*k3 + 9./32.*k4 + 125./768.*k5 + 5./66.*k6
        rk_error = dm - (7./1408.*k1 + 1125./2816.*k3 + 9./32.*k4 + 125./768.*k5             + 5./66.*k7 + 5./66.*k8)
        return m+dm, t+dt, rk_error

    def _optimal_stepsize(self, rk_error, atol):
        # single device sync per attempted step; all step-size control
        # afterwards is plain python float arithmetic
        norm = float(torch.linalg.norm(rk_error.flatten() / atol, torch.inf))
        if math.isnan(norm):
            raise RuntimeError("Unexpected error norm= %.5g!" % norm)

        if norm > 1.1:
            # decrease step, no more than factor of 5, but a fraction S more
            # than scaling suggests (for better accuracy)
            r = self._headroom / norm**(1.0/self._order)
            if (r < self._minscale):
                r = self._minscale
        elif norm < 0.5:
            # increase step, but no more than by a factor of 5
            if norm > 0.:
                r = self._headroom / norm**(1.0/(self._order+1.0))
            else:
                r = self._maxscale
            if r > self._maxscale: # increase no more than factor of 5
                r = self._maxscale
            if r < 1.: # don't allow any decrease caused by S<1
                r = 1.
        else: # no change
            r = 1.

        dt_opt = self._dt * r
        if dt_opt > self._maxstep:
            dt_opt = self._maxstep
        return dt_opt

    def step(self, state, dt, rtol = None, atol = None, **llg_args):
        remaining = float(dt) # track remaining time in python to avoid device syncs
        while remaining > 0.:
            _m1, _t1, err = self._try_step(state, **llg_args)
            dt_opt = self._optimal_stepsize(err, atol or self._atol)
            if self._dt > dt_opt or self._dt > remaining:
                # step size was too large, retry with optimal stepsize
                # also rescale the thermal field accordingly
                self._dt = min(dt_opt, remaining)
                logging.debug("REVERT step: %g, new step size: %g", self._dt, dt_opt)
            else:
                # accept step, adapt stepsize for next step
                remaining -= self._dt
                state.m = _m1
                state.t = _t1
                logging.debug("ACCEPT step: %g, new step size: %g", self._dt, dt_opt)
                self._dt = dt_opt
                state._step += 1
