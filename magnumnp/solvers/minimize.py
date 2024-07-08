from magnumnp.common import logging, timedmethod, constants
import torch

__all__ = ["MinimizerBB"]

class MinimizerBB(object):
    def __init__(self, terms):
        """ This class implements the direct energy minimizing algorithm introduced in [Exl2014]_.

        .. note:: This feature is experimental.

        *Example*
          .. code:: python

            state = State(mesh)
            minimizer = MinimizerBB([ExchangeField()])
            minimizer.minimize(state)

        *Arguments*
            terms ([:class:`LLGTerm`])
                List of LLG contributions to be considered for energy minimization
        """
        self._terms = terms

    def E(self, state):
        return sum([term.E(state) for term in self._terms])

    def h(self, state):
        return sum([term.h(state) for term in self._terms])

    def dm(self, state, h):
        return torch.linalg.cross(state.m, torch.linalg.cross(state.m, h))

    def _midpoint(self, m, h, tau):
        """
        Solving the the semi-implicit midpoint scheme:
            m_i+1 = m_i + tau * (m_i + m_i+1)/2 x (m_i x Heff[m_i])

        see "Abert, 'Efficient Energyminimization in Finite-Difference Micromagnetics', 2014"
        see "Goldfarb, 'A Curvilinear Search Method for p-Harmonic Flows on Spheres', 2009"
        """
        mxh = torch.linalg.cross(m, h)
        mx, my, mz = m.unbind(-1)
        mxh_x, mxh_y, mxh_z = mxh.unbind(-1)

        N = 4 + tau*tau * (mxh_x*mxh_x + mxh_y*mxh_y + mxh_z*mxh_z)
        return torch.stack([(4*mx + 4*tau * (mxh_y*mz - mxh_z*my) + tau*tau*mx * (+ mxh_x*mxh_x - mxh_y*mxh_y - mxh_z*mxh_z) + 2*tau*tau*mxh_x * (mxh_y*my + mxh_z*mz)) / N,
                            (4*my + 4*tau * (mxh_z*mx - mxh_x*mz) + tau*tau*my * (- mxh_x*mxh_x + mxh_y*mxh_y - mxh_z*mxh_z) + 2*tau*tau*mxh_y * (mxh_z*mz + mxh_x*mx)) / N,
                            (4*mz + 4*tau * (mxh_x*my - mxh_y*mx) + tau*tau*mz * (- mxh_x*mxh_x - mxh_y*mxh_y + mxh_z*mxh_z) + 2*tau*tau*mxh_z * (mxh_x*mx + mxh_y*my)) / N], dim=-1)
        # mumax code
        #// m = 1 / (4 + τ²(m x H)²) [{4 - τ²(m x H)²} m - 4τ(m x m x H)]
        #// note: torque from LLNoPrecess has negative sign

    def _linesearch(self, state, m0, h0, dm0, tau):
        r = 0.5  # Reduction factor
        c = 0.5  # Sufficient decrease parameter
        m = -(constants.mu_0*state.material["Ms"]*state.mesh.cell_volumes*dm0*dm0).sum()
        t = -c*m
        E0 = self.E(state)

        for j in range(10):
            state.m = self._midpoint(m0, h0, tau)
            E = self.E(state)
            if E0 - E >= tau * t:
                break
            tau *= r
            logging.info_blue("[MinimizerBB] Linesearch: %d, E=%g" % (j, E))

    @timedmethod
    def minimize(self, state, maxiter = 2000, dm_tol = 1., tau_min = 1e-13, tau_max = 1e-5):
        tau = tau_min
        steps = 0
        dm_max = 1e18
        m0 = state.m.clone()
        h0 = self.h(state)
        dm0 = self.dm(state, h0)

        for i in range(maxiter):
            #self._linesearch(state, m0, h0, dm0, tau)
            state.m = self._midpoint(m0, h0, tau)
            h = self.h(state)

            # compute s^n-1 for step-size control
            m_diff = state.m - m0

            # compute y^n-1 for step-size control
            dm = torch.linalg.cross(state.m, torch.linalg.cross(state.m, h))
            dm_diff = dm - dm0

            # compute dm_max as convergence indicator
            dm_max = dm.abs().max()
            if dm_max < dm_tol:
                break

            # next stepsize (alternate tau1 and tau2)
            if (steps % 2 == 0):
                tau = (m_diff*m_diff).sum() / (m_diff*dm_diff).sum()
            else:
                tau = (m_diff*dm_diff).sum() / (dm_diff*dm_diff).sum()
            tau = max(min(abs(tau), tau_max), tau_min) #* tau_sign

            logging.info_blue("[MinimizerBB] Step: %d, Tau: %.5g, dm_max: %.5g" % (steps, tau, dm_max))

            # increase step count
            steps += 1
            m0 = state.m.clone()
            h0 = h.clone()
            dm0 = dm.clone()

        return steps
