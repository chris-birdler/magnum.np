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

from magnumnp.common import logging, timedmethod, constants, normalize
from .ode_solvers import RKF45
import torch

__all__ = ["LLGSolver"]

class LLGSolver(object):
    def __init__(self, terms, solver = RKF45, no_precession = False, **kwargs):
        """
        This class implements the LLG dm term as well as the corresponding energy. 
        It also provides the interface for time-integration and allows switching between
        different ODE solvers.

        *Example*
            .. code:: python

            llg = LLGSolver([demag, exchange, external])
            logger = Logger("data", ['t', 'm'])
            while state.t < 1e-9-eps:
                llg.step(state, 1e-11)
                logger << state

        *Arguments*
            terms ([:class:`LLGTerm`])
                List of LLG contributions to be considered for time integration
            solver ([:class:`Solver`])
                ODE solver to be used (chose one of RKF45 (default), RKF56, ScipyODE, ScipyOdeint, TorchDiffEq, TorchDiffEqAdjoint)
            no_precession (bool)
                integrate without precession term (default: False)
        """
        self._terms = terms
        self._solver = solver(self.dm, **kwargs)
        self._no_precession = no_precession

    def dm(self, t, x, state, alpha = None):
        state.t = t
        state.m = x
        alpha = alpha or state.material["alpha"]

        gamma_prime = constants.gamma / (1. + alpha**2)
        alpha_prime = alpha * gamma_prime

        h = sum([term.h(state) for term in self._terms])

        cross_m_h = torch.linalg.cross(state.m, h)

        dm = -alpha_prime * torch.linalg.cross(state.m, cross_m_h)
        if not self._no_precession:
            dm -= gamma_prime * cross_m_h

        return dm

    def E(self, state):
        return sum([term.E(state) for term in self._terms])

    @timedmethod
    def step(self, state, dt, **kwargs):
        state.t, state.m = self._solver.step(state.t, state.m, dt, state=state, **kwargs)
        normalize(state.m)
        logging.info_blue("[LLG] step: dt= %g  t=%g" % (dt, state.t))

    @timedmethod
    def solve(self, state, tt, **kwargs):
        logging.info_blue("[LLG] solve: t0=%g  t1=%g Integrating ..." % (tt[0].cpu().numpy(), tt[-1].cpu().numpy()))
        res = self._solver.solve(tt, state.m, state=state, **kwargs)
        logging.info_blue("[LLG] solve: t0=%g  t1=%g Finished" % (tt[0].cpu().numpy(), tt[-1].cpu().numpy()))

        state.t = tt[-1]
        state.m = res[-1]
        return res

    @timedmethod
    def relax(self, state, maxiter = 500, dm_tol = 1e2, dt = 1e-11):
        t0 = state.t

        for i in range(maxiter):
            state.t, state.m = self._solver.step(state.t, state.m, dt, state=state, alpha = 1.0) #, no_precession = True) # no_precession requires more iterations for SP4 demo!?

            dm = self.dm(state.t, state.m, state=state, alpha = 1.0).abs().max() / constants.gamma # use same scaling as within minimizer
            logging.info_blue("[LLG] relax: t=%g |dm|=%g" % (state.t-t0, dm))
            if dm < dm_tol:
                break

        state.t = t0
