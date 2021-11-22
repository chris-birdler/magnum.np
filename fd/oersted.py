import os
import numpy as np
from math import atan, sqrt, pi, log
from scipy import constants

def krueger_g(p):
    x, y, z = p
    R = sqrt(x**2 + y**2 + z**2)

    res = (3.*x**2 + 3.*y**2 - 2.*z**2)*z*R/24. 
    res += pi*z/4.*abs(x*y*z)

    if x != 0. or y != 0:
        res += (x**4 - 6.*x**2*y**2 + y**4)/24. * log(z+R)
    if y != 0:
        res += x*y/6. * (y**2 - 3.*z**2) * atan(x*z/(y*R))
    if x != 0:
        res += x*y/6. * (x**2 - 3.*z**2) * atan(y*z/(x*R))
    if y != 0. or z != 0:
        res += z/6. * x * (z**2 - 3.*y**2) * log(x+R) 
    if x != 0. or z != 0:
        res += z/6. * y * (z**2 - 3.*x**2) * log(y+R)

    return res

class OerstedField(object):
    def __init__(self, mesh):
        self._mesh = mesh
        self._init_K()

    def _init_K_component(self, c, permute, func):
        it = np.nditer(self._K[:,:,:,c], flags=['multi_index'], op_flags=['writeonly'])
        while not it.finished:
            value = 0.0
            for i in np.rollaxis(np.indices((3,)*3), 0, 4).reshape(27, -1) - 1:
                idx = [(it.multi_index[k] + self._mesh.n[k]) % (2*self._mesh.n[k]) - self._mesh.n[k] for k in range(3)]
                value += np.prod(2.-3*abs(i)) * func([(idx[j] + i[j]) * self._mesh.dx[j] for j in permute])
            it[0] = value / (4 * pi * np.prod(self._mesh.dx))
            it.iternext()

    def _init_K(self):
        if os.path.isfile("cache/K%s.npy" % self._mesh):
            self._K = np.load("cache/K%s.npy" % self._mesh)
        else:
            self._K = np.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [3])
            for i, t in enumerate(((krueger_g,0,1,2),
                                   (krueger_g,1,2,0),
                                   (krueger_g,2,0,1))):
                self._init_K_component(i, t[1:], t[0])

            np.save("cache/K%s" % self._mesh, self._K)

        self._K_fft = np.fft.rfftn(self._K, axes = [i for i in range(3) if self._mesh.n[i] > 1])

        # init scratch spaces
        self._j_pad = np.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [3])
        self._h_fft = np.zeros(self._K_fft.shape[:3] + (3,), dtype=self._K_fft.dtype)


    def h(self, t, m):
        self._j_pad[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2],:] = m
        j_pad_fft = np.fft.rfftn(self._j_pad, axes = [i for i in range(3) if self._mesh.n[i] > 1])

        self._h_fft[:,:,:,0] =                                           - self._K_fft[:,:,:,0]*j_pad_fft[:,:,:,1] + self._K_fft[:,:,:,2]*j_pad_fft[:,:,:,2]
        self._h_fft[:,:,:,1] = + self._K_fft[:,:,:,0]*j_pad_fft[:,:,:,0]                                           - self._K_fft[:,:,:,1]*j_pad_fft[:,:,:,2]
        self._h_fft[:,:,:,2] = - self._K_fft[:,:,:,2]*j_pad_fft[:,:,:,0] + self._K_fft[:,:,:,1]*j_pad_fft[:,:,:,1]

        h_pad = np.fft.irfftn(self._h_fft, axes = [i for i in range(3) if self._mesh.n[i] > 1])
        return h_pad[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2],:]
