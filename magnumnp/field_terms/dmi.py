from magnumnp.common import timedmethod, constants
import torch
import numpy as np

__all__ = ["InterfaceDMIField", "BulkDMIField", "D2dDMIField"]

class DMIField(object):
    def __init__(self, dmi_vector, D = "D"):
        self._dmi_vector = dmi_vector
        self._D = D

    @timedmethod
    def h(self, state):
        h = state._zeros(state.mesh.n + (3,))
        D = state.material[self._D]

        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        for dim in range(3):
            v = state.Tensor(self._dmi_vector[dim]).expand(state.m[next].shape)
            if isinstance(D, torch.Tensor) and D.dim() == 4: # TODO: D could be a 1D tensor instead of a 4D tensor field
                D_avg = torch.where(D[next]*D[current] < 0,
                                    torch.sqrt(torch.sqrt(-D[next]*D[current])*torch.abs(D[next]+D[current]) / 2.),
                                    2.*D[next]*D[current]/(D[next]+D[current]))
                h[current] += D_avg * torch.linalg.cross(v, state.m[next]   ) / (2.*state.mesh.dx[dim])
                h[next]    -= D_avg * torch.linalg.cross(v, state.m[current]) / (2.*state.mesh.dx[dim])
            else:
                h[current] += D * torch.linalg.cross(v, state.m[next]   ) / (2.*state.mesh.dx[dim])
                h[next]    -= D * torch.linalg.cross(v, state.m[current]) / (2.*state.mesh.dx[dim])

            # rotate dimension
            current = current[-1:] + current[:-1]
            next = next[-1:] + next[:-1]

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * m * self.h(state))


class InterfaceDMIField(DMIField):
    def __init__(self, Di = "Di"):
        dmi_vector = [[ 0, 1, 0], # x
                      [-1, 0, 0], # y
                      [ 0, 0, 0]] # z
        super().__init__(dmi_vector, Di)


class BulkDMIField(object):
    @timedmethod
    def h(self, state):
        dx = state.mesh.dx

        if state.mesh.n[0] > 1:
            dmydx = torch.gradient(state.m[:,:,:,1], spacing = dx[0], dim = 0)[0]
            dmzdx = torch.gradient(state.m[:,:,:,2], spacing = dx[0], dim = 0)[0]
        else:
            dmydx = 0.
            dmzdx = 0.

        if state.mesh.n[1] > 1:
            dmxdy = torch.gradient(state.m[:,:,:,0], spacing = dx[1], dim = 1)[0]
            dmzdy = torch.gradient(state.m[:,:,:,2], spacing = dx[1], dim = 1)[0]
        else:
            dmxdy = 0.
            dmzdy = 0.

        if state.mesh.n[2] > 1:
            dmxdz = torch.gradient(state.m[:,:,:,0], spacing = dx[2], dim = 2)[0]
            dmydz = torch.gradient(state.m[:,:,:,1], spacing = dx[2], dim = 2)[0]
        else:
            dmxdz = 0.
            dmydz = 0.

        h = -2. * state.material["Db"] / constants.mu_0 / state.material["Ms"] * torch.stack((dmydz-dmzdy, dmzdx-dmxdz, dmxdy-dmydx), dim=-1)
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * m * self.h(state))


class D2dDMIField(object):
    @timedmethod
    def h(self, state):
        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        dx = [state.mesh.dx[i] for i in range(3) if state.mesh.n[i] > 1]

        dmxdy = torch.gradient(state.m[:,:,:,0], spacing = dx[1], dim = 1)[0]
        dmydx = torch.gradient(state.m[:,:,:,1], spacing = dx[0], dim = 0)[0]
        dmzdx = torch.gradient(state.m[:,:,:,2], spacing = dx[0], dim = 0)[0]
        dmzdy = torch.gradient(state.m[:,:,:,2], spacing = dx[1], dim = 1)[0]
        h = -2. * state.material["DD2d"] / constants.mu_0 / state.material["Ms"] * torch.stack((dmzdy, dmzdx, -dmxdy-dmydx), dim=-1)
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * m * self.h(state))
