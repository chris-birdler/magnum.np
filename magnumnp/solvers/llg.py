from magnumnp.common import logging, timedmethod, constants
from .rkf45 import RKF45
import torch
from torchdiffeq import odeint

__all__ = ["LLGSolver"]

class LLGSolver(object):
    def __init__(self, state, terms, rtol = 1e-4, atol = 1e-4):
        self._state = state
        self._terms = terms
        self._gamma_prime = constants.gamma / (1. + state._material["alpha"]**2)
        self._alpha_prime = state._material["alpha"] * self._gamma_prime
        self._rtol = rtol
        self._atol = atol
        self._solver = RKF45(lambda t, m: self._dm(t, m), state)

    def _dm(self, t, m):
        h = sum([term.h(t, m) for term in self._terms])
        return - self._gamma_prime * torch.cross(m, h) \
               - self._alpha_prime * torch.cross(m, torch.cross(m, h))

    @timedmethod
    def step(self, dt):
        self._solver.step(dt)
        logging.info("[LLG] step: dt= %g  t=%g" % (dt, self._state.t))



    # Deprecated Functions (will be removed)
    @timedmethod
    def step_torchdiffeq(self, dt, method = 'dopri5', options = {}):
        self._state.m = odeint(lambda t, m: self._dm(t, m), self._state.m, state.DoubleTensor([self._state.t, self._state.t + dt]), method=method, options=options)[1] # TODO: reuse Solver object?
        self._state.t += dt
        logging.info("[LLG]: t=%g" % self._state.t)

    @timedmethod
    def solve_torchdiffeq(self, t_final, dt=None, method = 'dopri5', options = {}):
        if dt == None:
            tt = self._state.linspace(self._state.t, t_final, steps=2)
            print("tt:", tt)
        else:
            tt = self._state.arange(self._state.t, t_final, dt)
        res = odeint(lambda t, m: self._dm(t, m), self._state.m, tt, rtol=self._rtol, atol=self._atol, method=method, options=options)
        self._state.t = t_final
        self._state.m[:,:,:,:] = res[-1,:,:,:,:]

        logging.info("[LLG]: t=%g" % self._state.t)
        return tt, res

