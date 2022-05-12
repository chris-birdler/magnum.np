from magnumnp.common import logging, timedmethod
import torch
from torchdiffeq import odeint

__all__ = ["LLGSolver"]

class LLGSolver(object):
    def __init__(self, state, terms):
        self._state = state
        self._terms = terms
        self._gamma_prime = state._material["gamma"] / (1. + state._material["alpha"]**2)
        self._alpha_prime = state._material["alpha"] * self._gamma_prime

    def _dm(self, t, m):
        h = sum([term.h(t, m) for term in self._terms])
        return - self._gamma_prime * torch.cross(m, h) \
               - self._alpha_prime * torch.cross(m, torch.cross(m, h))

    def step(self, dt, method = 'dopri5', options = {}):
        self.m = odeint(lambda t, m: self._dm(t, m), self._state.m, torch.DoubleTensor([self.t, self.t + dt]), method=method, options=options)[1] # TODO: reuse Solver object?
        self._state.t += dt
        logging.info("[LLG]: t=%g" % self._state.t)

    @timedmethod
    def solve(self, t_final, dt, method = 'dopri5', options = {}):
        tt = self._state.arange(self._state.t, t_final, dt)
        res = odeint(lambda t, m: self._dm(t, m), self._state.m, tt, method=method, options=options)
        self._state.t = t_final
        self._state.m[:,:,:,:] = res[-1,:,:,:,:]
        logging.info("[LLG]: t=%g" % self._state.t)
        return tt, res
