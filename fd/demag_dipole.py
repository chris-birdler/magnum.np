from scipy import constants
import numpy as np
import torch
import torch.fft
from torch import asinh, atan, sqrt, log, abs
import os
CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
cuda = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")
from time import time

import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

def dipole_f(points):
    x = points[:,:,:,0]
    y = points[:,:,:,1]
    z = points[:,:,:,2]

    result = (2.*x**2 - y**2 - z**2) * pow(x**2 + y**2 + z**2, -5./2.)
    result[0,0,0] = 0.
    return result

def dipole_g(points):
    x = points[:,:,:,0]
    y = points[:,:,:,1]
    z = points[:,:,:,2]

    result = 3.*x*y * pow(x**2 + y**2 + z**2, -5./2.)
    result[0,0,0] = 0.
    return result


class DemagFieldDipole(object):
    def __init__(self, mesh, material):
        self._mesh = mesh
        self._Ms = material["Ms"]
        self._init_N()

    def _init_N_component(self, c, perm, func):
        ij = [torch.fft.fftshift(torch.arange(n, dtype=torch.float64, device=cuda)) - n//2 for n in self._N.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        r = torch.stack([(ij[ind])*self._mesh.dx[ind] for ind in perm], dim=-1)
        self._N[:,:,:,c] = func(r) * np.prod(self._mesh.dx) / (4.*np.pi)

    def _init_N(self):
        if os.path.isfile("cache/Ndipole_%s.pt" % self._mesh):
            self._N = torch.load("cache/Ndipole_%s.pt" % self._mesh, map_location=cuda)
            logging.info("[DEMAG]: Use cached demag kernel")
        else:
            self._N = torch.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [6], dtype=torch.float64, device=cuda)

            time_kernel = time()
            for i, t in enumerate(((dipole_f,0,1,2),
                                   (dipole_g,0,1,2),
                                   (dipole_g,0,2,1),
                                   (dipole_f,1,2,0),
                                   (dipole_g,1,2,0),
                                   (dipole_f,2,0,1))):
                self._init_N_component(i, t[1:], t[0])

            logging.info(f"[DEMAG]: Time calculation of demag kernel = {time() - time_kernel} s")
            if not os.path.isdir("cache"):
                os.makedirs("cache")
            torch.save(self._N, "cache/Ndipole_%s.pt" % self._mesh)
            logging.info("[DEMAG]: Saved demag kernel")

        self._N_fft = torch.fft.rfftn(self._N, dim = [i for i in range(3) if self._mesh.n[i] > 1])

        # init scratch spaces
        self._m_pad = torch.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [3], dtype=torch.float64, device=cuda)
        self._h_fft = torch.zeros(list(self._N_fft.shape[:3]) + [3], dtype=self._N_fft.dtype, device=cuda)

    def h(self, t, m):
        self._m_pad[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2],:] = self._Ms * m
        m_pad_fft = torch.fft.rfftn(self._m_pad, dim = [i for i in range(3) if self._mesh.n[i] > 1])

        self._h_fft[:,:,:,0] = (self._N_fft[:,:,:,(0, 1, 2)]*m_pad_fft).sum(axis = 3)
        self._h_fft[:,:,:,1] = (self._N_fft[:,:,:,(1, 3, 4)]*m_pad_fft).sum(axis = 3)
        self._h_fft[:,:,:,2] = (self._N_fft[:,:,:,(2, 4, 5)]*m_pad_fft).sum(axis = 3)

        h_pad = torch.fft.irfftn(self._h_fft, dim = [i for i in range(3) if self._mesh.n[i] > 1])

        return h_pad[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2],:]

    def E(self, t, m):
        return - 0.5 * constants.mu_0 * self._mesh.cell_volume * torch.sum(self._Ms * m * self.h(t, m))

    def __str__(self):
        return "demag"
