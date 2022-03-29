import os
from scipy import constants
import numpy as np
import torch
import torch.fft
from torch import asinh, atan, sqrt, log, abs
from torch.cuda import IntTensor, DoubleTensor
import os
CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
cuda = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")
from time import time

import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

def krueger_g(points):
    x = points[:,:,:,0]
    y = points[:,:,:,1]
    z = points[:,:,:,2]

    R = sqrt(x**2 + y**2 + z**2)

    res = (3.*x**2 + 3.*y**2 - 2.*z**2)*z*R/24. 
    res += np.pi*z/4.*abs(x*y*z)

    #x**2 + y**2 > 0
    temp = torch.zeros_like(x, device = cuda)
    mask = (x**2 + y**2).gt(0)
    temp[mask] = ((x**4 - 6.*x**2*y**2 + y**4)/24. * log(z+R))[mask]
    res += temp

    #y**2 > 0
    temp = torch.zeros_like(x, device = cuda)
    mask = (y**2).gt(0)
    temp[mask] = (x*y/6. * (y**2 - 3.*z**2) * atan(x*z/(y*R)))[mask]
    res += temp

    #x**2 > 0
    temp = torch.zeros_like(x, device = cuda)
    mask = (x**2).gt(0)
    temp[mask] = (x*y/6. * (x**2 - 3.*z**2) * atan(y*z/(x*R)))[mask]
    res += temp

    #y**2 + z**2 > 0
    temp = torch.zeros_like(x, device = cuda)
    mask = (y**2+z**2).gt(0)
    temp[mask] = (z/6. * x * (z**2 - 3.*y**2) * log(x+R))[mask]
    res += temp

    #x**2 + z**2 > 0
    temp = torch.zeros_like(x, device = cuda)
    mask = (x**2+z**2).gt(0)
    temp[mask] = (z/6. * y * (z**2 - 3.*x**2) * log(y+R))[mask]
    res += temp

    return res


def dipole_g(points):
    x = points[:,:,:,0]
    y = points[:,:,:,1]
    z = points[:,:,:,2]

    R = sqrt(x**2 + y**2 + z**2)
    res = -z/R**3
    res[0,0,0] = 0.
    return res


class OerstedField(object):
    def __init__(self, mesh, p = 15):
        self._mesh = mesh
        self._p = p
        self._init_K()

    def _init_K_component(self, c, perm, func_near, func_far):
        # dipole far-field
        ij = [torch.fft.fftshift(torch.arange(n, dtype=torch.float64, device=cuda)) - n//2 for n in self._K.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        r = torch.stack([ij[ind]*self._mesh.dx[ind] for ind in perm], dim=-1)
        self._K[:,:,:,c] = func_far(r) * np.prod(self._mesh.dx) / (4.*np.pi)

        # newell near-field
        K_near = torch.zeros([1 if i==1 else 2*i for i in np.minimum(self._mesh.n, self._p)], dtype=torch.float64, device=cuda)
        ij = [torch.fft.fftshift(torch.arange(n, dtype=torch.float64, device=cuda)) - n//2 for n in K_near.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        for k in np.rollaxis(np.indices((3,)*3), 0, 4).reshape(27, -1) - 1:
            r = torch.stack([(ij[ind] + k[ind])*self._mesh.dx[ind] for ind in perm], dim=-1)
            K_near[:,:,:] += np.prod(2.-3*np.abs(k)) * func_near(r) / (4.*np.pi*np.prod(self._mesh.dx))

        n_near = np.minimum(self._mesh.n, self._p)
        self._K[-n_near[0]:,-n_near[1]:,-n_near[2]:,c] = K_near[-n_near[0]:,-n_near[1]:,-n_near[2]:]
        self._K[:n_near[0],:n_near[1],:n_near[2],c] = K_near[:n_near[0],:n_near[1],:n_near[2]]

    def _init_K(self):
        if os.path.isfile("cache/K%s.pt" % self._mesh):
            self._K = torch.load("cache/K%s.pt" % self._mesh, map_location=cuda)
            logging.info("[DEMAG]: Use cached oersted kernel")
        else:
            self._K = torch.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [3], dtype=torch.float64, device=cuda)

            time_kernel = time()
            for i, t in enumerate(((krueger_g,dipole_g,0,1,2),
                                   (krueger_g,dipole_g,1,2,0),
                                   (krueger_g,dipole_g,2,0,1))):
                self._init_K_component(i, t[2:], t[0], t[1])

            logging.info(f"[OERSTED]: Time calculation of oersted kernel = {time() - time_kernel} s")
            if not os.path.isdir("cache"):
                os.makedirs("cache")
            torch.save(self._K, "cache/K%s.pt" % self._mesh)
            logging.info("[OERSTED]: Saved oersted kernel")

        self._K_fft = torch.fft.rfftn(self._K, dim = [i for i in range(3) if self._mesh.n[i] > 1])

        # init scratch spaces
        self._j_pad = torch.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [3], dtype=torch.float64, device=cuda)
        self._h_fft = torch.zeros(list(self._K_fft.shape[:3]) + [3], dtype=self._K_fft.dtype, device=cuda)

    def h(self, t, j):
        self._j_pad[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2],:] = j
        j_pad_fft = torch.fft.rfftn(self._j_pad, dim = [i for i in range(3) if self._mesh.n[i] > 1])

        self._h_fft[:,:,:,0] =                                           - self._K_fft[:,:,:,0]*j_pad_fft[:,:,:,1] + self._K_fft[:,:,:,2]*j_pad_fft[:,:,:,2]
        self._h_fft[:,:,:,1] = + self._K_fft[:,:,:,0]*j_pad_fft[:,:,:,0]                                           - self._K_fft[:,:,:,1]*j_pad_fft[:,:,:,2]
        self._h_fft[:,:,:,2] = - self._K_fft[:,:,:,2]*j_pad_fft[:,:,:,0] + self._K_fft[:,:,:,1]*j_pad_fft[:,:,:,1]

        h_pad = torch.fft.irfftn(self._h_fft, dim = [i for i in range(3) if self._mesh.n[i] > 1])

        return h_pad[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2],:]

#    def E(self, t, m):
#        return - 0.5 * constants.mu_0 * self._mesh.cell_volume \
#               * torch.sum(self._Ms * m * self.h(t, m))

    def __str__(self):
        return "oersted"
