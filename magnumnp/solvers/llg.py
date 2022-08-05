from magnumnp.common import logging, timedmethod, constants
from .rkf45 import RKF45
import torch
from torchdiffeq import odeint

__all__ = ["LLGSolver"]

class LLGSolver(object):
    def __init__(self, state, terms, atol = 1e-5):
        self._state = state
        self._terms = terms
        self._gamma_prime = constants.gamma / (1. + state._material["alpha"]**2)
        self._alpha_prime = state._material["alpha"] * self._gamma_prime
        self._solver = RKF45(lambda t, m: self._dm(t, m), state, atol=atol)

    def _dm(self, t, m):
        h = sum([term.h(t, m) for term in self._terms])
        return - self._gamma_prime * torch.cross(m, h) \
               - self._alpha_prime * torch.cross(m, torch.cross(m, h))

    @timedmethod
    def step(self, dt):
        self._solver.step(dt)
        logging.info_blue("[LLG] step: dt= %g  t=%g" % (dt, self._state.t))
