import numpy as np
import scipy.fftpack

__all__ = ["DemagFieldPBC_numpy", "DemagFieldPBC_scipy", "DemagFieldPBC_numpy_real"]

class DemagFieldPBC_numpy(object):
    def __init__(self, mesh, material):
        self._mesh = mesh
        self._Ms = material["Ms"]

    def h(self, t, m):
        m_fft = np.fft.fftn(self._Ms * m, axes = list(filter(lambda i: self._mesh.n[i] > 1, range(3))))
        dx, dy, dz = self._mesh.dx

        kx = (2. * np.pi * np.arange(self._mesh.n[0]) / self._mesh.n[0]).reshape(-1,1,1)
        ky = (2. * np.pi * np.arange(self._mesh.n[1]) / self._mesh.n[1]).reshape(1,-1,1)
        kz = (2. * np.pi * np.arange(self._mesh.n[2]) / self._mesh.n[2]).reshape(1,1,-1)

        div_fft = (1.-np.exp(-1j*kx)) * m_fft[:,:,:,0] / dx \
                  + (1.-np.exp(-1j*ky)) * m_fft[:,:,:,1] / dy \
                  + (1.-np.exp(-1j*kz)) * m_fft[:,:,:,2] / dz

        with np.errstate(divide='ignore', invalid='ignore'):
            u_fft = -div_fft / (4./dx**2*np.sin(kx/2.)**2 + \
                                4./dy**2*np.sin(ky/2.)**2 + \
                                4./dz**2*np.sin(kz/2.)**2)
            u_fft[0,0,0] = 0
        
        h_fft = np.empty_like(m_fft)
        h_fft[:,:,:,0] = (1.-np.exp(1j*kx)) * u_fft / dx 
        h_fft[:,:,:,1] = (1.-np.exp(1j*ky)) * u_fft / dy 
        h_fft[:,:,:,2] = (1.-np.exp(1j*kz)) * u_fft / dz 

        h = np.fft.ifftn(h_fft, axes = list(filter(lambda i: self._mesh.n[i] > 1, range(3))))
        return h.real

class DemagFieldPBC_scipy(object):
    def __init__(self, mesh, material):
        self._mesh = mesh
        self._Ms = material["Ms"]

    def h(self, t, m):
        m_fft = scipy.fftpack.fftn(self._Ms * m, axes = list(filter(lambda i: self._mesh.n[i] > 1, range(3))))
        dx, dy, dz = self._mesh.dx

        kx = (2. * np.pi * np.arange(self._mesh.n[0]) / self._mesh.n[0]).reshape(-1,1,1)
        ky = (2. * np.pi * np.arange(self._mesh.n[1]) / self._mesh.n[1]).reshape(1,-1,1)
        kz = (2. * np.pi * np.arange(self._mesh.n[2]) / self._mesh.n[2]).reshape(1,1,-1)

        div_fft = (1.-np.exp(-1j*kx)) * m_fft[:,:,:,0] / dx \
                  + (1.-np.exp(-1j*ky)) * m_fft[:,:,:,1] / dy \
                  + (1.-np.exp(-1j*kz)) * m_fft[:,:,:,2] / dz

        with np.errstate(divide='ignore', invalid='ignore'):
            u_fft = -div_fft / (4./dx**2*np.sin(kx/2.)**2 + \
                                4./dy**2*np.sin(ky/2.)**2 + \
                                4./dz**2*np.sin(kz/2.)**2)
            u_fft[0,0,0] = 0
        
        h_fft = np.empty_like(m_fft)
        h_fft[:,:,:,0] = (1.-np.exp(1j*kx)) * u_fft / dx 
        h_fft[:,:,:,1] = (1.-np.exp(1j*ky)) * u_fft / dy 
        h_fft[:,:,:,2] = (1.-np.exp(1j*kz)) * u_fft / dz 

        h = scipy.fftpack.fftn(h_fft, axes = list(filter(lambda i: self._mesh.n[i] > 1, range(3))))
        return h.real

class DemagFieldPBC_numpy_real(object):
    def __init__(self, mesh, material):
        self._mesh = mesh
        self._Ms = material["Ms"]

    def h(self, t, m):
        m_fft = np.fft.rfftn(self._Ms * m, axes = list(filter(lambda i: self._mesh.n[i] > 1, range(3)))[::-1])
        dx, dy, dz = self._mesh.dx

        kx = (2. * np.pi * np.arange(self._mesh.n[0] // 2+1) / self._mesh.n[0]).reshape(-1,1,1)
        ky = (2. * np.pi * np.arange(self._mesh.n[1]) / self._mesh.n[1]).reshape(1,-1,1)
        kz = (2. * np.pi * np.arange(self._mesh.n[2]) / self._mesh.n[2]).reshape(1,1,-1)

        div_fft = (1.-np.exp(-1j*kx)) * m_fft[:,:,:,0] / dx \
                  + (1.-np.exp(-1j*ky)) * m_fft[:,:,:,1] / dy \
                  + (1.-np.exp(-1j*kz)) * m_fft[:,:,:,2] / dz

        with np.errstate(divide='ignore', invalid='ignore'):
            u_fft = -div_fft / (4./dx**2*np.sin(kx/2.)**2 + \
                                4./dy**2*np.sin(ky/2.)**2 + \
                                4./dz**2*np.sin(kz/2.)**2)
            u_fft[0,0,0] = 0
        
        h_fft = np.empty_like(m_fft)
        h_fft[:,:,:,0] = (1.-np.exp(1j*kx)) * u_fft / dx 
        h_fft[:,:,:,1] = (1.-np.exp(1j*ky)) * u_fft / dy 
        h_fft[:,:,:,2] = (1.-np.exp(1j*kz)) * u_fft / dz 

        h = np.fft.irfftn(h_fft, axes = list(filter(lambda i: self._mesh.n[i] > 1, range(3)))[::-1])
        return h

