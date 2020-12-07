import numpy as np
from scipy import constants

class AnisotropyField(object):
    def __init__(self, mesh, material):
        self._mesh = mesh
        self._Ms = material["Ms"]
        self._K = material["K"]
        #self._K_axis = np.array(material["K_axis"])
        self._K_axis = np.tile(material["K_axis"], mesh.n).reshape(mesh.n + (3,))

    def h(self, t, m):
        h = np.zeros(m.shape)
        return 2. * self._K * self._K_axis / (constants.mu_0 * self._Ms) * np.sum(self._K_axis * m, axis=3)
        #return 2. * self._K * self._K_axis / (constants.mu_0 * self._Ms) * \
        #       np.sum(m * self._K_axis, axis=3).repeat(3).reshape(m.shape)


    def E(self, t, m):
        return - 0.5 * constants.mu_0 * self._mesh.cell_volume \
               * np.sum(self._Ms * m * self.h(t, m))
