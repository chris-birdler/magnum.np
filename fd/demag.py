import os
import numpy as np
from math import asinh, atan, sqrt, pi
from scipy import constants

# newell f
def newell_f(p):
    x, y, z = abs(p[0]), abs(p[1]), abs(p[2])

    result = 1.0 / 6.0 * (2*x**2 - y**2 - z**2) * sqrt(x**2 + y**2 + z**2)

    if x**2 + z**2 > 0:
        result += y / 2.0 * (z**2 - x**2) * asinh(y / (sqrt(x**2 + z**2)))
  
    if x**2 + y**2 > 0:
        result += z / 2.0 * (y**2 - x**2) * asinh(z / (sqrt(x**2 + y**2)))
  
    if x * (x**2 + y**2 + z**2) > 0:
        result -= x*y*z * atan(y*z / (x * sqrt(x**2 + y**2 + z**2)))

    return result

# newell g
def newell_g(p):
    x, y, z = p[0], p[1], abs(p[2])
    
    result = - x*y * sqrt(x**2 + y**2 + z**2) / 3.0

    if x**2 + y**2 > 0:
        result += x*y*z * asinh(z / (sqrt(x**2 + y**2)))

    if y**2 + z**2 > 0:
        result += y / 6.0 * (3.0 * z**2 - y**2) * asinh(x / (sqrt(y**2 + z**2)))

    if x**2 + z**2 > 0:
        result += x / 6.0 * (3.0 * z**2 - x**2) * asinh(y / (sqrt(x**2 + z**2)))

    if z * (x**2 + y**2 + z**2) != 0:
        result -= z**3 / 6.0 * atan(x*y / (z * sqrt(x**2 + y**2 + z**2)))

    if y * (x**2 + y**2 + z**2) != 0:
        result -= z * y**2 / 2.0 * atan(x*z / (y * sqrt(x**2 + y**2 + z**2)))

    if x * (x**2 + y**2 + z**2) != 0:
        result -= z * x**2 / 2.0 * atan(y*z / (x * sqrt(x**2 + y**2 + z**2)))

    return result


class DemagField(object):
    def __init__(self, mesh, material):
        self._mesh = mesh
        self._Ms = material["Ms"]
        self._init_N()

    def _init_N_component(self, c, permute, func):
        it = np.nditer(self._N[:,:,:,c], flags=['multi_index'], op_flags=['writeonly'])
        while not it.finished:
            value = 0.0
            for i in np.rollaxis(np.indices((2,)*6), 0, 7).reshape(64, 6):
                idx = map(lambda k: (it.multi_index[k] + self._mesh.n[k]) % \
                                    (2*self._mesh.n[k]) - self._mesh.n[k], range(3))
                value += (-1)**sum(i) * func(map(lambda j: (idx[j] + i[j] - i[j+3]) * self._mesh.dx[j], permute))
            it[0] = - value / (4 * pi * np.prod(self._mesh.dx))
            it.iternext()

    def _init_N(self):
        if os.path.isfile("cache/%s.npy" % self._mesh):
            self._N = np.load("cache/%s.npy" % self._mesh)
        else:
            self._N = np.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [6])
            for i, t in enumerate(((newell_f,0,1,2),
                                   (newell_g,0,1,2),
                                   (newell_g,0,2,1),
                                   (newell_f,1,2,0),
                                   (newell_g,1,2,0),
                                   (newell_f,2,0,1))):
                self._init_N_component(i, t[1:], t[0])

            np.save("cache/%s" % self._mesh, self._N)

        self._N_fft = np.fft.rfftn(self._N, axes = filter(lambda i: self._mesh.n[i] > 1, range(3)))

        # init scratch spaces
        self._m_pad = np.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [3])
        self._h_fft = np.zeros(self._N_fft.shape[:3] + (3,), dtype=self._N_fft.dtype)


    def h(self, t, m):
        self._m_pad[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2],:] = self._Ms * m
        m_pad_fft = np.fft.rfftn(self._m_pad, axes = filter(lambda i: self._mesh.n[i] > 1, range(3)))

        self._h_fft[:,:,:,0] = (self._N_fft[:,:,:,(0, 1, 2)]*m_pad_fft).sum(axis = 3)
        self._h_fft[:,:,:,1] = (self._N_fft[:,:,:,(1, 3, 4)]*m_pad_fft).sum(axis = 3)
        self._h_fft[:,:,:,2] = (self._N_fft[:,:,:,(2, 4, 5)]*m_pad_fft).sum(axis = 3)

        h_pad = np.fft.irfftn(self._h_fft, axes = filter(lambda i: self._mesh.n[i] > 1, range(3)))
        return h_pad[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2],:]

    def E(self, t, m):
        return - 0.5 * constants.mu_0 * self._mesh.cell_volume \
               * np.sum(self._Ms * m * self.h(t, m))
