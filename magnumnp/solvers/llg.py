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

from magnumnp.common import logging, timedmethod, constants
from .rkf45 import RKF45
import torch

__all__ = ["LLGSolver"]

class LLGSolver(object):
    def __init__(self, terms, atol = 1e-5):
        self._terms = terms
        def dm(state, t, m):
            t0 = state.t
            m0 = state.m.detach()
            state.t = t
            state.m = m
            dm = self._dm(state)
            state.t = t0
            state.m = m0
            return dm
        self._solver = RKF45(dm, atol=atol)

    def _dm(self, state):
        gamma_prime = constants.gamma / (1. + state.material["alpha"]**2)
        alpha_prime = state.material["alpha"] * gamma_prime

        h = sum([term.h(state) for term in self._terms])
        return - gamma_prime * torch.cross(state.m, h) \
               - alpha_prime * torch.cross(state.m, torch.cross(state.m, h))

    @timedmethod
    def step(self, state, dt):
        self._solver.step(state, dt)
        logging.info_blue("[LLG] step: dt= %g  t=%g" % (dt, state.t))

    @timedmethod
    def relax(self, state, maxiter = 500, rtol = 1e-5, dt = 1e-11):
        alpha0 = state.material["alpha"]
        t0 = state.t
        state.material["alpha"] = 1.0
        E0 = sum([term.E(state) for term in self._terms])
        for i in range(maxiter) :
            self._solver.step(state, dt)
            E = sum([term.E(state) for term in self._terms])
            dE = torch.linalg.norm(((E - E0)/E).reshape(-1), ord = float("Inf"))
            logging.info_blue("[LLG] relax: t=%g dE=%g E=%g" % (state.t-t0, dE, E))
            if dE < rtol:
                break
            E0 = E
        state.t = t0
        state.material["alpha"] = alpha0
