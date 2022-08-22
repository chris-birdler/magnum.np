from magnumnp.common import timedmethod, constants
import torch

__all__ = ["InterfaceDMIField", "BulkDMIField", "D2dDMIField"]

class InterfaceDMIField(object):
    @timedmethod
    def h(self, state):
        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        dx = [state.mesh.dx[i] for i in range(3) if state.mesh.n[i] > 1]

        dmxdx = torch.gradient(state.m[:,:,:,0], spacing = dx[0], dim = 0)[0]
        dmydy = torch.gradient(state.m[:,:,:,1], spacing = dx[1], dim = 1)[0]
        dmzdx = torch.gradient(state.m[:,:,:,2], spacing = dx[0], dim = 0)[0]
        dmzdy = torch.gradient(state.m[:,:,:,2], spacing = dx[1], dim = 1)[0]
        h = -2. * state.material["Di"] / constants.mu_0 / state.material["Ms"] * torch.stack((dmzdx, dmzdy, -dmxdx-dmydy), dim=-1)
        return torch.nan_to_num(h, posinf=0, neginf=0)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * m * self.h(state))


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
