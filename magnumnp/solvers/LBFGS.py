from magnumnp.common import logging, timedmethod
import torch
import numpy as np

__all__ = ["Minimizer_LBFGS"]

class Minimizer_LBFGS(object):
    def __init__(self, terms, tolerance= 1e-8, max_iter=1000, memory_size = 5, dm_max = 1e4):
        self._terms = terms
        self._dm_max = dm_max
        self._memory_size = memory_size
        self._tolerance = tolerance
        self.max_iter = max_iter

    def _dm(self, state):
        h = sum([term.h(state) for term in self._terms])
        return torch.linalg.cross(state.m, torch.linalg.cross(state.m, h))
    
    def linesearch(self, state, p):
        alpha = 1.0 #initial step size
        tau = 0.5  # Reduction factor
        c = 0.1  # Sufficient decrease parameter   
        h = sum([term.h(state) for term in self._terms])
        dm = torch.linalg.cross(state.m, torch.linalg.cross(state.m, h))
        E = sum([term.E(state) for term in self._terms])
        m = (dm*p).sum()
        t = -c*m
        j = 0
        maxIter = 1000
        
        while j < maxIter:
            m_new = state.m + alpha * p
            state.m = m_new
            E_new = sum([term.E(state) for term in self._terms])
            sufficient_decrease = E - E_new - alpha * t
            
            if sufficient_decrease >= 0:
                break
            
            alpha *= tau
            j += 1
    
        return state

    def minimize(self, state, tolerance= 1e-8, max_iter=1000, memory_size = 5,):
        step = 0
        glob_step = 0
        glob_steps = []
        glob_steps.append(glob_step)
        E = sum([term.E(state) for term in self._terms]) #Energy
        energy = []
        energy.append(E)
        h = sum([term.h(state) for term in self._terms])
        s_vectors = [torch.zeros(state.m.size())]*memory_size
        y_vectors = [torch.zeros(state.m.size())]*memory_size
        dm = torch.linalg.cross(state.m, torch.linalg.cross(state.m, h)) #Gradient

        eps = 2.22e-16;
        eps2 = np.sqrt(eps)
        epsr = np.power(eps, 0.9)
        tolf2 = np.sqrt(tolerance)
        tolf3 = np.power(tolerance, 0.3333333333333333333333333)

        dm_max = dm.max()
        if dm_max < epsr*(1+torch.abs(E)):
            return state, energy, glob_steps
        
        alpha = np.zeros(memory_size, dtype=np.float64)
        q = torch.zeros(state.m.size(), dtype=torch.float64)
        s = torch.zeros(state.m.size(), dtype=torch.float64)
        y = torch.zeros(state.m.size(), dtype=torch.float64)
        E_old = 0.0
        dm_old = torch.zeros(dm.size(), dtype=torch.float64)
        m_old = state.m;
        H0k = 1

        while glob_step < max_iter:
            
            E_old = E
            m_old = state.m
            dm_old = dm
            
            q = dm.clone()
            k = min(memory_size, step)
            for i in range(k - 1, -1, -1):
                rho = 1.0 / (s_vectors[i]* y_vectors[i]).sum()
                alpha[i] = rho * (s_vectors[i]*q).sum()
                q -= alpha[i] * y_vectors[i]
                
            z = H0k * q
            for i in range(k):
                rho = 1.0 / (s_vectors[i]* y_vectors[i]).sum()
                beta = rho * (y_vectors[i]*z).sum()
                z += s_vectors[i] * (alpha[i] - beta)
                #z = -z for minimization?
                
            phiPrime0 = -(dm*z).sum()
            if phiPrime0 > 0:
                z = dm
                step = 0
 
            #linesearch and update of state.m
            state = self.linesearch(state, z)
            
            #update gradient
            h = sum([term.h(state) for term in self._terms])
            dm = torch.linalg.cross(state.m, torch.linalg.cross(state.m, h))

            dm_max = dm.max()
            if dm_max < epsr*(1+torch.abs(E)):
                return state, energy, glob_steps
            
            m_diff = state.m - m_old
            dm_diff = dm - dm_old
            
            ys = (dm_diff*m_diff).sum()
            #if ys <= eps2 * torch.norm(dm_diff) * torch.norm(m_diff):
            #    
            #else:
            if step < memory_size:
                s_vectors[step] = m_diff
                y_vectors[step] = dm_diff
            else:
                s_vectors[:-1] = s_vectors[1:]  # Cyclic shift of s_vectors
                s_vectors[memory_size - 1] = m_diff
                y_vectors[:-1] = y_vectors[1:]  # Cyclic shift of y_vectors
                y_vectors[memory_size - 1] = dm_diff
            H0k = ys / (dm_diff*dm_diff).sum()
            step += 1
            
            # increase step count
            glob_step += 1
            glob_steps.append(glob_step)
            
        return state, energy, glob_steps


