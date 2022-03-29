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

def dipole_g(points):
    x = points[:,:,:,0]
    y = points[:,:,:,1]
    z = points[:,:,:,2]

    R = sqrt(x**2 + y**2 + z**2)
    res = -z/R**3
    res[0,0,0] = 0.
    return res


class OerstedFieldDipole(object):
    def __init__(self, mesh):
        self._mesh = mesh
        self._init_K()

    def _init_K_component(self, c, perm, func):
        ij = [torch.fft.fftshift(torch.arange(n, dtype=torch.float64, device=cuda)) - n//2 for n in self._K.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        r = torch.stack([ij[ind]*self._mesh.dx[ind] for ind in perm], dim=-1)
        self._K[:,:,:,c] += func(r) * np.prod(self._mesh.dx) / (4.*np.pi)

    def _init_K(self):
        if os.path.isfile("cache/K%s.pt" % self._mesh):
            self._K = torch.load("cache/K%s.pt" % self._mesh, map_location=cuda)
            logging.info("[DEMAG]: Use cached oersted kernel")
        else:
            self._K = torch.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [3], dtype=torch.float64, device=cuda)

            time_kernel = time()
            for i, t in enumerate(((dipole_g,0,1,2),
                                   (dipole_g,1,2,0),
                                   (dipole_g,2,0,1))):
                self._init_K_component(i, t[1:], t[0])

            logging.info(f"[OERSTED]: Time calculation of demag kernel = {time() - time_kernel} s")
            if not os.path.isdir("cache"):
                os.makedirs("cache")
            torch.save(self._K, "cache/K%s.pt" % self._mesh)
            logging.info("[OERSTED]: Saved demag kernel")

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
