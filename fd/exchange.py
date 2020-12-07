import numpy as np
from scipy import ndimage, constants

class ExchangeField(object):
    def __init__(self, mesh, material):
        self._mesh = mesh
        self._A = material["A"]
        self._Ms = material["Ms"]

        # initialize laplace kernel
        self._kernel = np.zeros((3,3,3))
        self._kernel[:,1,1] += np.array((1,-2,1)) / mesh.dx[0]**2
        self._kernel[1,:,1] += np.array((1,-2,1)) / mesh.dx[1]**2
        self._kernel[1,1,:] += np.array((1,-2,1)) / mesh.dx[2]**2

        # initialize scratch space
        self._h = np.zeros(mesh.n + (3,))

    def h(self, t, m):
        f = 2. * self._A / (constants.mu_0 * self._Ms)
        for i in range(3):
            self._h[:,:,:,i] = f * ndimage.convolve(m[:,:,:,i], self._kernel)

        return self._h

    def E(self, t, m):
        return - 0.5 * constants.mu_0 * self._mesh.cell_volume \
               * np.sum(self._Ms * m * self.h(t, m))
