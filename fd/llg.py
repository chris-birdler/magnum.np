import torch
from torchdiffeq import odeint
import os
CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
cuda = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")

import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

class LLGSolver(object):
    def __init__(self, terms, material, m0, t0 = 0.0):
        self._terms = terms
        self.t = t0
        self.m = m0
        self._gamma_prime = material["gamma"] / (1. + material["alpha"]**2)
        self._alpha_prime = material["alpha"] * self._gamma_prime

    def _dm(self, t, m):
        h = sum([term.h(t, m) for term in self._terms])
        return - self._gamma_prime * torch.cross(m, h) \
               - self._alpha_prime * torch.cross(m, torch.cross(m, h))

    def step(self, dt, method = 'dopri5', options = {}):
        self.m = odeint(lambda t, m: self._dm(t, m), self.m, torch.DoubleTensor([self.t, self.t + dt], device=cuda), method=method, options=options)[1] # TODO: reuse Solver object?
        self.t += dt
        logging.info("[LLG]: t=%g" % self.t)
        return self.m
