import torch
from torchdiffeq import odeint
import os
CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
cuda = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")

import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

class Minimizer(object):
    def __init__(self, terms):
        self._terms = terms

    def _dm(self, t, m):
        h = sum([term.h(t, m) for term in self._terms])
        return -torch.cross(m, torch.cross(m, h))

    def minimize(self, m, rtol = 1e-4, dt = 1e-4, maxiter = 1000, method = 'dopri5', options = {'first_step':1e-7}):
        m_old = torch.zeros_like(m)
        t     = 0.
        dmdt  = torch.inf

        for i in range(maxiter):
            m_old[:,:,:,:] = m

            m = odeint(lambda t, m: self._dm(t, m), m, torch.DoubleTensor([t, t + dt], device=cuda), method=method, options=options)[1] # TODO: reuse Solver object?
            t += dt
            dmdt = torch.linalg.norm((m - m_old).reshape(-1), ord = float("Inf")) / dt
            E = sum([term.E(t, m) for term in self._terms])
            logging.info("[MIN]: i=%g dmdt=%g E=%g" % (i, dmdt.item(), E.item()))
            if dmdt < rtol:
                break
        return m
