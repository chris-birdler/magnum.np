from magnumnp.common import timedmethod
import torch

__all__ = ["DemagFieldPBC"]

class DemagFieldPBC(object):
    def __init__(self, state):
        self._state = state
        self._mesh = self._state._mesh
        self._Ms = self._state._material["Ms"]

    @timedmethod
    def h(self, t, m):
        m_fft = torch.fft.rfftn(self._Ms * m, dim = [i for i in range(3) if self._mesh.n[i] > 1]).squeeze(-1)
        dx, dy, dz = self._mesh.dx

        kx = (2. * np.pi * self._state._arange(self._mesh.n[0]) / self._mesh.n[0]).reshape(-1,1,1)
        ky = (2. * np.pi * self._state._arange(self._mesh.n[1]) / self._mesh.n[1]).reshape(1,-1,1)
        kz = (2. * np.pi * self._state._arange(self._mesh.n[2]) / self._mesh.n[2]).reshape(1,1,-1)

        div_fft = (1.-torch.exp(-1j*kx)) * m_fft[:,:,:,0] / dx \
                + (1.-torch.exp(-1j*ky)) * m_fft[:,:,:,1] / dy \
                + (1.-torch.exp(-1j*kz)) * m_fft[:,:,:,2] / dz

        u_fft = -div_fft / (4./dx**2*torch.sin(kx/2.)**2 + \
                            4./dy**2*torch.sin(ky/2.)**2 + \
                            4./dz**2*torch.sin(kz/2.)**2)
        u_fft[0,0,0] = 0
        
        h_fft = torch.empty_like(m_fft)
        h_fft[:,:,:,0] = (1.-torch.exp(1j*kx)) * u_fft / dx 
        h_fft[:,:,:,1] = (1.-torch.exp(1j*ky)) * u_fft / dy 
        h_fft[:,:,:,2] = (1.-torch.exp(1j*kz)) * u_fft / dz 

        h = torch.fft.irfftn(h_fft, dim = [i for i in range(3) if self._mesh.n[i] > 1])
        return h.real
