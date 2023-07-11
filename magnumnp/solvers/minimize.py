from magnumnp.common import logging, timedmethod, constants
import torch

__all__ = ["MinimizerBB"]

class MinimizerBB(object):
    def __init__(self, terms, tau_min = 1e-13, tau_max = 1e-5, dm_max = 1e-3, samples = 1):
        """
        This class implements the direct energy minimizing algorithm introduced in [Exl2014]_.

        .. note:: This feature is experimental.

        *Example*
          .. code:: python

            state = State(mesh)
            minimizer = MinimizerBB([ExchangeField()])
            minimizer.minimize(state)

        *Arguments*
          terms ([:class:`LLGTerm`])
            List of LLG contributions to be considered for energy minimization
          region (:class:`str`)
            region on which the energy is minimized
          tau_min (:class:`float`)
            minimum step size
          tau_max (:class:`float`)
            maximum step size
          dm_max (:class:`float`)
            stop criterion given as supremum norm of dm/dt
          sample (:class:`int`)
            number of subsequent steps the stop criterion has to be fulfilled
        """
        self._terms = terms
        self._tau_min = tau_min
        self._tau_max = tau_max
        self._dm_max = dm_max
        self._samples = samples
        self._tau = tau_min

    def E(self, state):
        return sum([term.E(state) for term in self._terms])

    def _h(self, state):
        return sum([term.h(state) for term in self._terms])

    def _dm(self, state, h):
        return torch.cross(state.m, torch.cross(state.m, h))

    def _midpoint(self, m, h, tau):
        """
        Solving the the semi-implicit midpoint scheme:
            m_i+1 = m_i + tau * (m_i + m_i+1)/2 x (m_i x Heff[m_i])

        see "Abert, 'Efficient Energyminimization in Finite-Difference Micromagnetics', 2014"
        see "Goldfarb, 'A Curvilinear Search Method for p-Harmonic Flows on Spheres', 2009"
        """
        mxh = torch.cross(m, h)
        mx, my, mz = m.unbind(-1)
        mxh_x, mxh_y, mxh_z = mxh.unbind(-1)

        N = 4 + tau*tau * (mxh_x*mxh_x + mxh_y*mxh_y + mxh_z*mxh_z)
        return torch.stack([(4*mx + 4*tau * (mxh_y*mz - mxh_z*my) + tau*tau*mx * (+ mxh_x*mxh_x - mxh_y*mxh_y - mxh_z*mxh_z) + 2*tau*tau*mxh_x * (mxh_y*my + mxh_z*mz)) / N,
                            (4*my + 4*tau * (mxh_z*mx - mxh_x*mz) + tau*tau*my * (- mxh_x*mxh_x + mxh_y*mxh_y - mxh_z*mxh_z) + 2*tau*tau*mxh_y * (mxh_z*mz + mxh_x*mx)) / N,
                            (4*mz + 4*tau * (mxh_x*my - mxh_y*mx) + tau*tau*mz * (- mxh_x*mxh_x - mxh_y*mxh_y + mxh_z*mxh_z) + 2*tau*tau*mxh_z * (mxh_x*mx + mxh_y*my)) / N], dim=-1)
        # mumax code
        #// m = 1 / (4 + τ²(m x H)²) [{4 - τ²(m x H)²} m - 4τ(m x m x H)]
        #// note: torque from LLNoPrecess has negative sign

    def minimize(self, state):
        tau = self._tau_min
        step = 0
        dm_max = 1e18
        last_dm_max = []
        E = sum([term.E(state) for term in self._terms])
        energy = [E]
        steps = []
        h = sum([term.h(state) for term in self._terms])
        logging.info_blue("Tau: %.5g, dm_max: %.5g, E: %g  " % (tau, dm_max, E))

        while len(last_dm_max) < self._samples or max(last_dm_max) > self._dm_max:
#        for i in range(1000):
            m_next = self._midpoint(state.m, h, tau)

            #m_next = state.m - tau*dm
            #m_next = m_next.normalize()

            # compute s^n-1 for step-size control
            m_diff = m_next - state.m

            m0 = torch.clone(state.m)
            dm = torch.cross(state.m, torch.cross(state.m, h))

            # update state
            state.m = torch.clone(m_next)
            E = sum([term.E(state) for term in self._terms])

#            if E > energy[-1]:
#                print("dE > 0: E= ", E.numpy())
#            #    tau = self._linesearch(state, dm_next)
#            #    #energy[-1] = E
#            #    print(f'E={E}')
#            #    print(f'tau={tau}')
#            #    continue
#
#                for dt in state.linspace(0, tau, 1000):
#                    #state.m = m0 - dt*dm
#                    state.m = self._midpoint(m0, h, dt)
#                    E = sum([term.E(state) for term in self._terms]).numpy()
#                    print(f'TEST E = {E}, dt/tau = {(dt/tau).numpy()}')
#                state.m = self._midpoint(m0, h, tau)
#                E = sum([term.E(state) for term in self._terms]).numpy()
#                exit(0)

            # compute y^n-1 for step-size control
            h = sum([term.h(state) for term in self._terms]) # TODO: calculate only once
            dm_next = torch.cross(state.m, torch.cross(state.m, h))
            dm_diff = dm_next - dm

            # compute dm_max as convergence indicator
            dm_max = dm_next.max()
            last_dm_max.append(dm_max)
            if len(last_dm_max) > self._samples: last_dm_max.pop(0)

            if m_diff.abs().max() == 0.:
                logging.info_blue("[MinimizerBB] dm == 0! Minimum reached.")
                break

            if dm_diff.abs().max() == 0.:
                logging.info_blue("[MinimizerBB] dm_diff == 0! Minimum reached.")
                break

            # next stepsize (alternate tau1 and tau2)
            try:
                if (step % 2 == 0):
                    tau = (m_diff*m_diff).sum() / (m_diff*dm_diff).sum()
                else:
                    tau = (m_diff*dm_diff).sum() / (dm_diff*dm_diff).sum()
            except ZeroDivisionError:
                tau = self._tau_max

            tau_sign = torch.sign(tau) #TODO: check tau_sign
            tau = max(min(abs(tau), self._tau_max), self._tau_min) #* tau_sign

            # increase step count
            step += 1
            steps.append(step)

            logging.info_blue("Tau: %.5g, dm_max: %.5g, E: %g  " % (tau, dm_max, E))

        return energy, steps


    def _linesearch(self, state, m0, h0, E0, dm0, tau):
        r = 0.5  # Reduction factor
        c = 0.5  # Sufficient decrease parameter
        m = -(constants.mu_0*state.material["Ms"]*state.cell_volumes*dm0*dm0).sum()
        t = -c*m

        for j in range(10):
            state.m = self._midpoint(m0, h0, tau)
            E = sum([term.E(state) for term in self._terms])
            if E0 - E >= tau * t:
                break
            tau *= r
            logging.info_blue("[MinimizerBB] Linesearch: %d, E=%g" % (j, E))

        return E

    def minimize2(self, state):
        tau = self._tau_min
        steps = 0
        dm_max = 1e18
        last_dm_max = []
        m0 = state.m.clone()
        h0 = sum([term.h(state) for term in self._terms])
        E0 = sum([term.E(state) for term in self._terms])
        dm0 = torch.cross(state.m, torch.cross(state.m, h0))

        #while len(last_dm_max) < self._samples or max(last_dm_max) > self._dm_max:
        while dm_max > self._dm_max:
            #E = self._linesearch(state, m0, h0, E0, dm0, tau)
            state.m = self._midpoint(m0, h0, tau)
            h = sum([term.h(state) for term in self._terms])
            E = sum([term.E(state) for term in self._terms])

            #if energy[-1] > energy[-2]:
            #    for dt in state._linspace(0, 1e-5, 1000):
            #        #state.m = m0 - dt*dm
            #        state.m = self._midpoint(m0, h, dt)
            #        E = sum([term.E(state) for term in self._terms]).numpy()
            #        print(f'{E} {dt}')
            #    exit(0)

            # compute s^n-1 for step-size control
            m_diff = state.m - m0

            # compute y^n-1 for step-size control
            dm = torch.cross(state.m, torch.cross(state.m, h))
            dm_diff = dm - dm0

            # compute dm_max as convergence indicator
            dm_max = dm.abs().max()

            # next stepsize (alternate tau1 and tau2)
            if (steps % 2 == 0):
                tau = (m_diff*m_diff).sum() / (m_diff*dm_diff).sum()
            else:
                tau = (m_diff*dm_diff).sum() / (dm_diff*dm_diff).sum()
            tau = max(min(abs(tau), self._tau_max), self._tau_min) #* tau_sign

            if torch.isnan(state.Tensor(tau)):
                logging.info_blue("[MinimizerBB] Minimum reached (isnan).")
                break

            logging.info_blue("[MinimizerBB] Step: %d, Tau: %.5g, dm_max: %.5g, E=%g, last: %g" % (steps, tau, dm_max, E, dm_max))

            # increase step count
            steps += 1
            m0 = state.m.clone()
            h0 = h.clone()
            E0 = E.clone()
            dm0 = dm.clone()

        return E, steps
