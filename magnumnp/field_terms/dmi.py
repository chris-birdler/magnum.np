from magnumnp.common import timedmethod, constants
import torch

__all__ = ["InterfaceDMIField"]

class InterfaceDMIField(object):
    def __init__(self, state):
        self._state = state
        self._mesh = state._mesh
        self._Js = constants.mu_0 * self._state._material["Ms"]
        self._Di = self._state._material["Di"]
        self._Di_axis = state.zeros(state._mesh.n + (3,))
        self._Di_axis[:,:,:,:] = state.DoubleTensor([0, 0, 1])

    @timedmethod
    def h(self, t, m):
        dim = [i for i in range(3) if self._mesh.n[i] > 1]
        dx = [self._mesh.dx[i] for i in range(3) if self._mesh.n[i] > 1]

        dmxdx = torch.gradient(self._state.m[:,:,:,0], spacing = dx[0], dim = 0)[0]
        dmydy = torch.gradient(self._state.m[:,:,:,1], spacing = dx[1], dim = 1)[0]
        dmzdx = torch.gradient(self._state.m[:,:,:,2], spacing = dx[0], dim = 0)[0]
        dmzdy = torch.gradient(self._state.m[:,:,:,2], spacing = dx[1], dim = 1)[0]
        return -2.*self._Di/self._Js * torch.stack((dmzdx, dmzdy, -dmxdx-dmydy), dim=-1)
 

    def E(self, t, m):
        return -0.5 * constants.mu_0 * self._mesh.cell_volume * torch.sum(self._Js * m * self.h(t, m)) 


#needs correct boundary conditions
class InterfaceDMIFieldFULL(object):
    def __init__(self, state):
        self._state = state
        self._mesh = state._mesh
        self._Ms = self._state._material["Ms"]
        self._Di = self._state._material["Di"]
        self._Di_axis = state.zeros(state._mesh.n + (3,))
        self._Di_axis[:,:,:,:] = state.DoubleTensor(state._material["Di_axis"])

    @timedmethod
    def h(self, t, m):
        dim = [i for i in range(3) if self._mesh.n[i] > 1]
        dx = [self._mesh.dx[i] for i in range(3) if self._mesh.n[i] > 1]
        nablaedm = torch.stack(torch.gradient(torch.sum(self._Di_axis * self._state.m, dim = 3), spacing=dx, dim=dim), dim = -1)
        print(nablaedm)
        nablamed = torch.sum(torch.stack(torch.gradient(self._state.m, spacing = dx, dim = dim), dim = -1), dim = -1)
        print(nablamed)
        return - 2 * self._Di /(constants.mu_0 * self._Ms)*(nablaedm - nablamed)
 

    def E(self, t, m):
        raise -0.5 * constants.mu_0 * self._mesh.cell_volume * torch.sum(self._Js * m * self.h(t, m)) 
