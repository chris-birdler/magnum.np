from magnumnp.common import logging, timedmethod, constants
from .rkf45 import RKF45
import torch

__all__ = ["LLGSolver"]

class LLGSolver(object):
    def __init__(self, terms, atol = 1e-5):
        self._terms = terms
        def dm(state, t, m):
            t0 = state.tensor(state.t)
            m0 = state.m.detach()
            state.t = t
            state.m = m # TODO: m needs to be reset after changing
            dm = self._dm(state)
            state.t = t0
            state.m = m0
            return dm
        #self._solver = RKF45(lambda state: self._dm(state), atol=atol)
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
