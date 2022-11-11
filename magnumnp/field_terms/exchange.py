from magnumnp.common import timedmethod, constants
import torch

__all__ = ["ExchangeField", "ExchangeDMIField", "ExchangeDMIFieldOpenBC"]

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
        # calculate m_1.5:
        for dim in [0]: #range(3): # TODO: omit singelton dimentions # TODO: material needs to be a TensorField (use expand)
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

            # rotate dimension
            #current = current[-1:] + current[:-1]
            #next = next[-1:] + next[:-1]

        # DMI field
        dmdx = (m_15[next] - m_15[current]) / state.mesh.dx[dim]
        h += 2.*Di * torch.stack((dmdx[...,2], 0.*dmdx[...,1], -dmdx[...,0]), dim=-1)

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))


class ExchangeDMIFieldOpenBC(object):
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

        for dim in [0]: #range(3): # TODO: omit singelton dimentions # TODO: material needs to be a TensorField (use expand)
            # laplace m = (m_i+1 - m_i) + (m_i-1 - m_i) = m_i+1 - 2*m_i + m_i-1
            # dm_dx = (m_i+1 - m_i) - (m_i-1 - m_i) = m_i+1 - m_i-1

            # Exchange Field
            A_avg = 2.*A[next]*A[current]/(A[next]+A[current])
            dm_next    = state.m[next] - state.m[current]
            dm_current = state.m[current] - state.m[next]
            h[current] += A_avg * dm_next    / state.mesh.dx[dim]**2 # m_i+1 - m_i
            h[next]    += A_avg * dm_current / state.mesh.dx[dim]**2 # m_i-1 - m_i

            # DMI field
            dmdx = state._zeros(state.mesh.n + (3,))
            dmdx[current] += state.m[next]    / (2.*state.mesh.dx[dim]) #  m_i+1 - m_i
            dmdx[next]    -= state.m[current] / (2.*state.mesh.dx[dim]) # -m_i-1 + m_i
            h += Di * torch.stack((dmdx[...,2], 0.*dmdx[...,1], -dmdx[...,0]), dim=-1)
        
            # rotate dimension
            current = current[-1:] + current[:-1]
            next = next[-1:] + next[:-1]

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))
