from magnumnp.common import timedmethod, constants
import torch

__all__ = ["ExchangeField", "ExchangeDMIField"]

class ExchangeField(object):
    def __init__(self, domain=None):
        self._domain = domain

    @timedmethod
    def h(self, state):
        h = state._zeros(state.mesh.n + (3,))

        A = state.material["A"]
        if self._domain != None:
            A = A * self._domain[:,:,:,None]
        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        for dim in range(3): # TODO: [] could be overloaded for constant Decorated Tensor in order to handle both cases
            if isinstance(A, torch.Tensor) and A.dim() == 4: # TODO: A could be a 1D tensor instead of a 4D tensor field
                A_avg = 2.*A[next]*A[current]/(A[next]+A[current])
                h[current] += A_avg * (state.m[next] - state.m[current]) / state.mesh.dx[dim]**2 # m_i+1 - m_i
                h[next]    += A_avg * (state.m[current] - state.m[next]) / state.mesh.dx[dim]**2 # m_i-1 - m_i
            else:
                h[current] += A * (state.m[next] - state.m[current]) / state.mesh.dx[dim]**2 # m_i+1 - m_i
                h[next]    += A * (state.m[current] - state.m[next]) / state.mesh.dx[dim]**2 # m_i-1 - m_i

            # rotate dimension
            current = current[-1:] + current[:-1]
            next = next[-1:] + next[:-1]

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))




class ExchangeDMIField(object):
    def __init__(self, domain=None):
        self._domain = domain

    @timedmethod
    def h(self, state):
        h = state._zeros(state.mesh.n + (3,))

        A = state.material["A"]
        Di = state.material["Di"]
        if self._domain != None:
            A = A * self._domain[:,:,:,None]
        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        # assue 1. dimension




# rhs = (A[current]*state.m[current] + A[next]*state.m[next])
#
# m_15 = torch.concat([+rhs[...,(0,)] *                       (A[current] + A[next]) / ((A[current] + A[next])**2 + (state.mesh.dx[0]/4.*(Di[current]-Di[next]))**2)
#                      -rhs[...,(2,)] * state.mesh.dx[dim]/4.*(Di[current]-Di[next]) / ((A[current] + A[next])**2 + (state.mesh.dx[0]/4.*(Di[current]-Di[next]))**2),
#
#                      rhs[...,(1,)] / (A[current] + A[next]),
#
#                      +rhs[...,(0,)] * state.mesh.dx[dim]/4.*(Di[current]-Di[next]) / ((A[current] + A[next])**2 + (state.mesh.dx[0]/4.*(Di[current]-Di[next]))**2)
#                      +rhs[...,(2,)] *                       (A[current] + A[next]) / ((A[current] + A[next])**2 + (state.mesh.dx[0]/4.*(Di[current]-Di[next]))**2)], dim=-1)
#
# # exchange field
# h[current] += 2.*A[current] * (m_15 - state.m[current]) / state.mesh.dx[dim]**2 # m_i+1 - m_i
# h[next]    += 2.*A[next]    * (m_15 - state.m[next])    / state.mesh.dx[dim]**2 # m_i-1 - m_i




        # calculate m_1.5:
        for dim in range(3): # TODO: omit singelton dimentions # TODO: material needs to be a TensorField (use expand)
            A0, A1 = A.pad(dim, -1), A.pad(dim, +1)
            Di0, Di1 = Di.pad(dim, -1), Di.pad(dim, +1)
            m0, m1 = state.m.pad(dim, -1), state.m.pad(dim, +1)

            rhs = A0*m0 + A1*m1

            m_15 = torch.concat([+rhs[...,(0,)] *                       (A0 + A1) / ((A0 + A1)**2 + (state.mesh.dx[0]/4.*(Di0-Di1))**2)
                                 -rhs[...,(2,)] * state.mesh.dx[dim]/4.*(Di0-Di1) / ((A0 + A1)**2 + (state.mesh.dx[0]/4.*(Di0-Di1))**2),

                                 rhs[...,(1,)] / (A0 + A1),

                                 +rhs[...,(0,)] * state.mesh.dx[dim]/4.*(Di0-Di1) / ((A0 + A1)**2 + (state.mesh.dx[0]/4.*(Di0-Di1))**2)
                                 +rhs[...,(2,)] *                       (A0 + A1) / ((A0 + A1)**2 + (state.mesh.dx[0]/4.*(Di0-Di1))**2)], dim=-1)

            # exchange field
            h += 2.*A * (m_15[current] - 2.*state.m + m_15[next]) / state.mesh.dx[dim]**2

#            # DMI field
#            #m:       0     1     2     3     4
#            #m15:  *     0     1     2     3      *   (*==0)???
#            dm = state._zeros(state.mesh.n + (3,))
#            dm[current] = m_15 - state.m[current]
#            dm[next] -=   m_15 - state.m[next]
#
#
#            #dm[0] = m_15[0] - m[0]
#            #      - 0
#            #dm[1] = m_15[1] - m[1]
#            #      - m_15[0] + m[1] =
#            #      = m_15[1] - m_15[0]
#            h += -state.material["Di"] * torch.stack((dm[...,2], 0.*dm[...,1], -dm[...,0]), dim=-1)

            # rotate dimension
            current = current[-1:] + current[:-1]
            next = next[-1:] + next[:-1]

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))
