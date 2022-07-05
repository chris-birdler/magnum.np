from magnumnp.common import timedmethod, constants
import torch

__all__ = ["InterfaceDMIField", "BulkDMIField"]

class InterfaceDMIField(object):
    def __init__(self, state):
        self._state = state
        self._mesh = state._mesh
        self._Js = constants.mu_0 * self._state._material["Ms"]
        self._Di = self._state._material["Di"]

    @timedmethod
    def h(self, t, m):
        dim = [i for i in range(3) if self._mesh.n[i] > 1]
        dx = [self._mesh.dx[i] for i in range(3) if self._mesh.n[i] > 1]

        dmxdx = torch.gradient(self._state.m[:,:,:,0], spacing = dx[0], dim = 0)[0]
        dmydy = torch.gradient(self._state.m[:,:,:,1], spacing = dx[1], dim = 1)[0]
        dmzdx = torch.gradient(self._state.m[:,:,:,2], spacing = dx[0], dim = 0)[0]
        dmzdy = torch.gradient(self._state.m[:,:,:,2], spacing = dx[1], dim = 1)[0]
        h = -2.*self._Di/self._Js * torch.stack((dmzdx, dmzdy, -dmxdx-dmydy), dim=-1)
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, t, m):
        return -0.5 * constants.mu_0 * self._mesh.cell_volume * torch.sum(self._Js * m * self.h(t, m))


class BulkDMIField(object):
    def __init__(self, state):
        self._state = state
        self._mesh = state._mesh
        self._Js = constants.mu_0 * self._state._material["Ms"]
        self._Db = self._state._material["Db"]

    @timedmethod
    def h(self, t, m):
        dx = self._mesh.dx

        if self._mesh.n[0] > 1:
            dmydx = torch.gradient(self._state.m[:,:,:,1], spacing = dx[0], dim = 0)[0]
            dmzdx = torch.gradient(self._state.m[:,:,:,2], spacing = dx[0], dim = 0)[0]
        else:
            dmydx = 0.
            dmzdx = 0.

        if self._mesh.n[1] > 1:
            dmxdy = torch.gradient(self._state.m[:,:,:,0], spacing = dx[1], dim = 1)[0]
            dmzdy = torch.gradient(self._state.m[:,:,:,2], spacing = dx[1], dim = 1)[0]
        else:
            dmxdy = 0.
            dmzdy = 0.

        if self._mesh.n[2] > 1:
            dmxdz = torch.gradient(self._state.m[:,:,:,0], spacing = dx[2], dim = 2)[0]
            dmydz = torch.gradient(self._state.m[:,:,:,1], spacing = dx[2], dim = 2)[0]
        else:
            dmxdz = 0.
            dmydz = 0.

        h = -2. * self._Db/self._Js * torch.stack((dmydz-dmzdy, dmzdx-dmxdz, dmxdy-dmydx), dim=-1)
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, t, m):
        return -0.5 * constants.mu_0 * self._mesh.cell_volume * torch.sum(self._Js * m * self.h(t, m))
