from magnumnp.common import timedmethod, constants
import torch

__all__ = ["ExchangeField", "ExchangeDMI", "ExchangeDMIField2"]

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


class ExchangeDMI(object):
    def __init__(self, domain=None):
        self._domain = domain

    @timedmethod
    def h(self, state):
        #TODO: maybe pad m internally
        h = state._zeros(state.mesh.n + (3,))
        laplace = state._zeros(state.mesh.n + (3,))
        dm = state._zeros(state.mesh.n + (3,))

        # exchange contribution
        A = state.material["A"]
        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        # x - component
        dim = 0
        A_avg = 2.*A[next]*A[current]/(A[next]+A[current]) # assume A is tensor-field
        dm_next    = A_avg * (state.m[next] - state.m[current]) / state.mesh.dx[dim]**2 # m_i+1 - m_i
        dm_current = A_avg * (state.m[current] - state.m[next]) / state.mesh.dx[dim]**2 # m_i-1 - m_i
        laplace[...] = 0
        laplace[current] = dm_next
        laplace[next] += dm_current

        h += laplace

        dm[current] = dm_next
        dm[next] -= dm_current
        h[:,:,:,(0,)] += -2. * state.material["Di"] / constants.mu_0 / state.material["Ms"] * dm[:,:,:,(2,)]
        h[:,:,:,(2,)] -= -2. * state.material["Di"] / constants.mu_0 / state.material["Ms"] * dm[:,:,:,(0,)]

        current = current[-1:] + current[:-1] # rotate dimension
        next = next[-1:] + next[:-1]

        # y - component
        dim = 1
        A_avg = 2.*A[next]*A[current]/(A[next]+A[current]) # assume A is tensor-field
        dm_next    = A_avg * (state.m[next] - state.m[current]) / state.mesh.dx[dim]**2 # m_i+1 - m_i
        dm_current = A_avg * (state.m[current] - state.m[next]) / state.mesh.dx[dim]**2 # m_i-1 - m_i
        laplace[...] = 0
        laplace[current] = dm_next
        laplace[next] += dm_current

        h += laplace

        dm[current] = dm_next
        dm[next] -= dm_current
        h[:,:,:,(1,)] += -2. * state.material["Di"] / constants.mu_0 / state.material["Ms"] * dm[:,:,:,(2,)]
        h[:,:,:,(2,)] -= -2. * state.material["Di"] / constants.mu_0 / state.material["Ms"] * dm[:,:,:,(1,)]

        current = current[-1:] + current[:-1] # rotate dimension
        next = next[-1:] + next[:-1]

        # z - component
        dim = 2
        A_avg = 2.*A[next]*A[current]/(A[next]+A[current]) # assume A is tensor-field
        dm_next    = A_avg * (state.m[next] - state.m[current]) / state.mesh.dx[dim]**2 # m_i+1 - m_i
        dm_current = A_avg * (state.m[current] - state.m[next]) / state.mesh.dx[dim]**2 # m_i-1 - m_i
        laplace[...] = 0
        laplace[current] = dm_next
        laplace[next] += dm_current

        h += laplace

        current = current[-1:] + current[:-1] # rotate dimension
        next = next[-1:] + next[:-1]

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))






class ExchangeDMIField2(object):
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
        dim = 0
        rhs = (4.*A[current]*state.m[current] + 4.*A[next]*state.m[next])
        m_15 = torch.concat([+rhs[...,(0,)] * (4.*A[current] + 4.*A[next]) / ((4.*A[current] + 4.*A[next])**2 - (Di[current]-Di[next])**2)
                             -rhs[...,(2,)] * (Di[current]-Di[next])       / ((4.*A[current] + 4.*A[next])**2 - (Di[current]-Di[next])**2),

                             rhs[...,(1,)] / (4.*A[current] + 4.*A[next]),

                             -rhs[...,(0,)] * (Di[current]-Di[next])       / ((4.*A[current] + 4.*A[next])**2 - (Di[current]-Di[next])**2)
                             +rhs[...,(2,)] * (4.*A[current] + 4.*A[next]) / ((4.*A[current] + 4.*A[next])**2 - (Di[current]-Di[next])**2)], dim=-1)

        # exchange field
        h[current] += 2.*A[current] * (m_15 - state.m[current]) / state.mesh.dx[dim]**2 # m_i+1 - m_i
        h[next]    += 2.*A[next]    * (m_15 - state.m[next])    / state.mesh.dx[dim]**2 # m_i-1 - m_i

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))
