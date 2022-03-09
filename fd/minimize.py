import numpy as np
from scipy import integrate

class Minimizer(object):
    def __init__(self, terms):
        self._terms = terms

    def _dm(self, t, m):
        h = np.sum([term.h(t, m) for term in self._terms], axis=0)
        result =- np.cross(m, np.cross(m, h))
        return result

    def minimize(self, m, rtol = 1e-4, dt = 1e-4, maxiter = np.inf):
        shape = m.shape
        ode = integrate.ode(lambda t, m: self._dm(t, m.reshape(shape)).reshape(-1))
        ode.set_initial_value(m.reshape(-1), 0.)
        ode.set_integrator("dopri5")

        m_current = np.zeros(m.shape)
        m_next    = m.copy()
        t         = 0.
        dmdt      = np.inf
        i         = 0

        with open('data/m0.log', 'w') as f:
            while dmdt > rtol and i < maxiter:
                t += dt
                m_current[:] = m_next
                m_next = ode.integrate(t).reshape(shape)
                f.write("%g %g %g %g\n" % ((t,) + tuple(np.mean(m_next, axis=(0,1,2)))))
                f.flush()
                # TODO renormalize

                # compute max norm of dm/dt
                dmdt = np.linalg.norm((m_current - m_next).reshape(-1), ord = np.inf) / dt
                # compute total energy
                E = np.sum([term.E(t, m_next) for term in self._terms])

                print("i=%d dmdt=%g E=%g" % (i, dmdt, E))
                i += 1

        return m_next
