import numpy as np
from scipy import ndimage, constants

__all__ = ["SpinTorque"]

class SpinTorque(object):
    def __init__(self, mesh, material, j):
        self._stencil = np.zeros((3,3,3,3))
        self._stencil[0,:,1,1] = np.array([1,0,-1]) / (2 * mesh.dx[0])
        self._stencil[1,1,:,1] = np.array([1,0,-1]) / (2 * mesh.dx[1])
        self._stencil[2,1,1,:] = np.array([1,0,-1]) / (2 * mesh.dx[2])

        self._mesh = mesh
        self._j = j
        self._gamma = material["gamma"]
        self._b = material["b"]
        self._xi = material["xi"]

        # initialize scratch spaces
        self._jgradm = np.zeros(mesh.n + (3,))
        #self._h = np.zeros(mesh.n + (3,))

    def h(self, t, m):
        self._jgradm.fill(0.)
        for i in np.arange(3):
            for k in np.arange(3):
                if self._j[k] != 0:
                    self._jgradm[:, :, :, i] += self._j[k] * ndimage.convolve(m[:,:,:,i], self._stencil[k], mode = 'mirror')

        return self._b / self._gamma * (np.cross(m, self._jgradm) + self._xi * self._jgradm)

    def E(self, t, m):
        raise NotImplemented()
