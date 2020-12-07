import numpy as np
from scipy import constants

class ExternalField(object):
    def __init__(self, mesh, material, h):
        self._mesh = mesh
        self._Ms = material["Ms"]
        self._h = np.zeros(mesh.n + (3,))
        self._h[:,:,:,:] = np.array(h)

    def h(self, t, m):
        return self._h

    def E(self, t, m):
        return - constants.mu_0 * self._mesh.cell_volume \
               * np.sum(self._Ms * m * self.h(t, m))
