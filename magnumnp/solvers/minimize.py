from magnumnp.common import logging, timedmethod
import torch

__all__ = ["Minimizer_BB"]

class Minimizer_BB(object):
    def __init__(self, terms, tau_min = 1e-13, tau_max = 1e-5, dm_max = 1e4, samples = 10):
        """
        This class implements the direct energy minimizing algorithm introduced in [Exl2014]_.

        .. note:: This feature is experimental.

        *Example*
          .. code:: python

            state = State(mesh)
            minimizer = Minimizer_BB([ExchangeField()])
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

    def _dm(self, state):
        h = sum([term.h(state) for term in self._terms])
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
    
    def _linesearch(self, state, dm):
        alpha = 1.0 #initial step size
        tau = 0.5  # Reduction factor
        c = 0.1  # Sufficient decrease parameter   
        h = sum([term.h(state) for term in self._terms])
        dm = torch.cross(state.m, torch.cross(state.m, h))
        E = sum([term.E(state) for term in self._terms])
        m = (dm*dm).sum()
        t = -c*m
        j = 0
        maxIter = 1000
        
        while j < maxIter:
            h = sum([term.h(state) for term in self._terms])
            m_new = self._midpoint(state.m, h, alpha)
            state.m = m_new
            E_new = sum([term.E(state) for term in self._terms])
            sufficient_decrease = E - E_new - alpha * t
            
            if sufficient_decrease >= 0:
                break
            
            alpha *= tau
            j += 1
    
        return alpha

    def minimize(self, state):
        tau = self._tau_min
        step = 0
        dm_max = 1e18
        last_dm_max = []
        E = sum([term.E(state) for term in self._terms])
        energy = []
        #energy.append(E)
        steps = []
        h = sum([term.h(state) for term in self._terms])


        while len(last_dm_max) < self._samples or max(last_dm_max) > self._dm_max:
            m_next = self._midpoint(state.m, h, tau)
            
            #m_next = state.m - tau*dm
            #m_next = m_next.normalize()

            # compute s^n-1 for step-size control
            m_diff = m_next - state.m
            
            m0 =torch.clone(state.m)
            
            dm = torch.cross(state.m, torch.cross(state.m, h))
            
            # update state
            state.m = torch.clone(m_next)
            
            # compute y^n-1 for step-size control
            h = sum([term.h(state) for term in self._terms]) # TODO: calculate only once
            dm_next = torch.cross(state.m, torch.cross(state.m, h))
            dm_diff = dm_next - dm
                        
            E = sum([term.E(state) for term in self._terms])
            #if E > energy[-1]:
            #    tau = self._linesearch(state, dm_next)
            #    #energy[-1] = E
            #    print(f'E={E}')
            #    print(f'tau={tau}')
            #    continue
            energy.append(E)
            
                #for dt in state._linspace(0, tau, 1000):
                #    state.m = m0 - dt*dm
                #    #state.m = self._midpoint(m0, h, dt)
                #    E = sum([term.E(state) for term in self._terms]).numpy()
                #    print(f'E = {E}, dt = {-dt}')
                #exit(0)

            # compute dm_max as convergence indicator
            dm_max = dm_next.max()
            last_dm_max.append(dm_max)
            if len(last_dm_max) > self._samples: last_dm_max.pop(0)

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
            
            logging.info_blue("Tau: %.5g, dm_max: %.5g, E=%g" % (tau, dm_max, E))
            
        return state, energy, steps



#### magnum.fd
#   def minimize(self, max_dpns = 0.01, samples = 10, h_max = 1e-5, h_min = 1e-16):
#        # TODO make use of stephandlers for logging
#        h        = self.state.h
#        dpnslist = []
#        log      = ScreenLogMinimizer()
#
#        # Reset step
#        self.state.step = 0
#
#        while len(dpnslist) < samples or max(dpnslist) > max_dpns:
#            # Calculate next M and dM for minimization step
#            M_next = self.state.minimizer_M(h)
#            dM = self.state.minimizer_dM
#
#            # Get s^n-1 for step-size calculation
#            M_diff = VectorField(self.mesh)
#            M_diff.assign(M_next)
#            M_diff.add(self.state.M, -1.0)
#
#            # Set next M
#            self.state.y = M_next
#            self.state.finish_step() # normalize, TODO really need to do this every step?
#            self.state.flush_cache()
#
#            # Calculate deg per ns
#            # TODO M.absMax might be the wrong choice if different materials are in use
#            dp_timestep = (180.0 / math.pi) * math.atan2(M_diff.absMax(), self.state.M.absMax())
#            dpns = abs(1e-9 * dp_timestep / h)
#            dpnslist.append(dpns)
#            if len(dpnslist) > samples: dpnslist.pop(0)
#            self.state.deg_per_ns_minimizer = dpns
#
#            # Get y^n-1 for step-size calculation
#            dM_diff = VectorField(self.mesh)
#            dM_diff.assign(self.state.minimizer_dM)
#            dM_diff.add(dM, -1.0)
#
#            # Next stepsize (Alternate h1 and h2)
#            try:
#              if (self.state.step % 2 == 0):
#                h = M_diff.dotSum(M_diff) / M_diff.dotSum(dM_diff)
#              else:
#                h = M_diff.dotSum(dM_diff) / dM_diff.dotSum(dM_diff)
#            except ZeroDivisionError, ex:
#              h = h_max

