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

from magnumnp.common import logging, timedmethod
import torch
from torchdiffeq import odeint

__all__ = ["Minimizer"]

class Minimizer(object):
    def __init__(self, state, terms):
        self._state = state
        self._terms = terms

    def _dm(self, t, m):
        h = sum([term.h(t, m) for term in self._terms])
        return -torch.cross(m, torch.cross(m, h))

    @timedmethod
    def minimize(self, rtol = 1e-4, dt = 1e-4, maxiter = 1000, method = 'dopri5', options = {'first_step':1e-7}):
        m = self._state.m
        m_old = torch.zeros_like(m)
        t     = 0.
        dmdt  = torch.inf

        for i in range(maxiter):
            m_old[:,:,:,:] = m

            m = odeint(lambda t, m: self._dm(t, m), m, self._state.tensor([t, t + dt]), method=method, options=options)[1] # TODO: reuse Solver object?
            t += dt
            dmdt = torch.linalg.norm((m - m_old).reshape(-1), ord = float("Inf")) / dt
            E = sum([term.E(t, m) for term in self._terms])
            logging.info("[MIN]: i=%g dmdt=%g E=%g" % (i, dmdt.item(), E.item()))
            if dmdt < rtol:
                break

    @timedmethod
    def minimize_event(self, tol = 1e-2, T = 10e-4, dt = 1e-4, maxiter = 1000, method = 'dopri5', options = {'first_step':1e-7}):
        m = self._state.m
        m_old = torch.zeros_like(m)
        t     = 0.
        t_old = 0.
        t_old_log = -1.
        dmdt  = torch.inf

        def event_fn(t, m):
            nonlocal m_old, t_old, t_old_log
            dmdt = torch.linalg.norm((m - m_old).reshape(-1), ord = float("Inf")) / (t-t_old)
            if t - t_old_log > 1e-4:
                logging.info("[MIN]: t= %.04g   dmdt=%.04g" % (t, dmdt.item()))
                t_old_log = t
            m_old[:,:,:,:] = m
            t_old = t
            return t < T

        t, m_opt = odeint(lambda t, m: self._dm(t, m), m, self._state.tensor([t, t+T]), method=method, options=options, event_fn = event_fn)
        self._state.m[:,:,:,:] = m_opt[-1,:,:,:,:]
