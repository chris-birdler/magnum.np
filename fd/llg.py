import numpy as np
from scipy import integrate

class LLG(object):
    def __init__(self, terms, material, m0, t0 = 0.0):
        self._terms = terms
        self.t = t0
        self._gamma_prime = material["gamma"] / (1. + material["alpha"]**2)
        self._alpha_prime = material["alpha"] * self._gamma_prime

        # initialize ode integrator
        self._shape = m0.shape
        self._ode = integrate.ode(lambda t, m: self._dmdt(t, m.reshape(self._shape)).reshape(-1))
        self._ode.set_integrator("dopri5")
        self._ode.set_initial_value(m0.reshape(-1), t0)

    def _dmdt(self, t, m):
        h = np.sum([term.h(t, m) for term in self._terms], axis=0)

        return - self._gamma_prime * np.cross(m, h) \
               - self._alpha_prime * np.cross(m, np.cross(m, h))

    def step(self, dt = 1e-12):
        self.t += dt
        return self._ode.integrate(self.t).reshape(self._shape)
