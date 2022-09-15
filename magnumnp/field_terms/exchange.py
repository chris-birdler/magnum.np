from magnumnp.common import timedmethod, constants
import torch

__all__ = ["ExchangeField", "ExchangeDMI"]

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
        h = state._zeros(state.mesh.n + (3,))

        # exchange contribution
        A = state.material["A"]
        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        for dim in range(3):
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

#        # DMI contribution
#        dim = [i for i in range(3) if state.mesh.n[i] > 1]
#        dx = [state.mesh.dx[i] for i in range(3) if state.mesh.n[i] > 1]
#
#        dmxdx = torch.gradient(state.m[:,:,:,0], spacing = dx[0], dim = 0)[0]
#        dmydy = torch.gradient(state.m[:,:,:,1], spacing = dx[1], dim = 1)[0]
#        dmzdx = torch.gradient(state.m[:,:,:,2], spacing = dx[0], dim = 0)[0]
#        dmzdy = torch.gradient(state.m[:,:,:,2], spacing = dx[1], dim = 1)[0]
#        h += -2. * state.material["Di"] / constants.mu_0 / state.material["Ms"] * torch.stack((dmzdx, dmzdy, -dmxdx-dmydy), dim=-1)
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))
