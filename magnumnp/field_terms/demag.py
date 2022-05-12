from magnumnp.common import logging, timedmethod, constants
import numpy as np
import torch
import torch.fft
from torch import asinh, atan, sqrt, log, abs
import os
CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
cuda = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")
from time import time

__all__ = ["DemagField"]

def newell_f(points):
    x = abs(points[:,:,:,0])
    y = abs(points[:,:,:,1])
    z = abs(points[:,:,:,2])

    result = 1.0 / 6.0 * (2*x**2 - y**2 - z**2) * sqrt(x**2 + y**2 + z**2)

    # x**2 + z**2 > 0:
    temp = torch.zeros_like(x, dtype=torch.float64, device = cuda)
    mask = (x**2 + z**2).gt(0)
    temp[mask] = (y / 2.0 * (z**2 - x**2) * asinh(y / sqrt(x**2 + z**2)))[mask]
    result += temp

    # x**2 + y**2 > 0:
    temp = torch.zeros_like(x, dtype=torch.float64, device = cuda)
    mask = (x**2 + y**2).gt(0)
    temp[mask] = (z / 2.0 * (y**2 - x**2) * asinh(z / sqrt(x**2 + y**2)))[mask]
    result += temp

    # x * (x**2 + y**2 + z**2) > 0:
    temp = torch.zeros_like(x, dtype=torch.float64, device = cuda)
    mask = (x * (x**2 + y**2 + z**2)).gt(0)
    temp[mask] = (x * y * z * atan(y*z / (x * sqrt(x**2 + y**2 + z**2))))[mask]
    result -= temp

    return result

def newell_g(points):
    x = points[:,:,:,0]
    y = points[:,:,:,1]
    z = abs(points[:,:,:,2])

    result = - x*y * sqrt(x**2 + y**2 + z**2) / 3.0

    # x**2 + y**2 > 0
    temp = torch.zeros_like(x, dtype=torch.float64, device = cuda)
    mask = (x**2 + y**2).gt(0)
    temp[mask] = (x*y*z * asinh(z / sqrt(x**2 + y**2)))[mask]
    result += temp

    # y**2 + z**2 > 0
    temp = torch.zeros_like(x, dtype=torch.float64, device = cuda)
    mask = (y**2 + z**2).gt(0)
    temp[mask] = (y / 6.0 * (3.0 * z**2 - y**2) * asinh(x / sqrt(y**2 + z**2)))[mask]
    result += temp

    # x**2 + z**2 > 0
    temp = torch.zeros_like(x, dtype=torch.float64, device = cuda)
    mask = (x**2 + z**2).gt(0)
    temp[mask] = (x / 6.0 * (3.0 * z**2 - x**2) * asinh(y / sqrt(x**2 + z**2)))[mask]
    result += temp

    # z * (x**2 + y**2 + z**2) != 0
    temp = torch.zeros_like(x, dtype=torch.float64, device = cuda)
    mask = (z * (x**2 + y**2 + z**2)).ne(0)
    temp[mask] = ( z**3 / 6.0 * atan(x*y / (z * sqrt(x**2 + y**2 + z**2))))[mask]
    result -= temp

    # y * (x**2 + y**2 + z**2) != 0
    temp = torch.zeros_like(x, dtype=torch.float64, device = cuda)
    mask = (y * (x**2 + y**2 + z**2)).ne(0)
    temp[mask] = (z * y**2 / 2.0 * atan(x*z / (y * sqrt(x**2 + y**2 + z**2))))[mask]
    result -= temp

    # x * (x**2 + y**2 + z**2) != 0
    temp = torch.zeros_like(x, dtype=torch.float64, device = cuda)
    mask = (x * (x**2 + y**2 + z**2)).ne(0)
    temp[mask] = (z * x**2 / 2.0 * atan(y*z / (x * sqrt(x**2 + y**2 + z**2))))[mask]
    result -= temp

    return result


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


class DemagField(object):
    def __init__(self, mesh, material, p = 20):
        self._mesh = mesh
        self._Ms = material["Ms"]
        self._p = p
        self._init_N()
        torch.cuda.empty_cache()

    def _init_N_component(self, N, c, perm, func_near, func_far):
        # dipole far-field
        ij = [torch.fft.fftshift(torch.arange(n, dtype=torch.float64, device=cuda)) - n//2 for n in N.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        r = torch.stack([ij[ind]*self._mesh.dx[ind] for ind in perm], dim=-1)
        N[:,:,:,c] = func_far(r) * np.prod(self._mesh.dx) / (4.*np.pi)

        # newell near-field
        n_near = np.minimum(self._mesh.n, self._p)
        N_near = torch.zeros([1 if i==1 else 2*i for i in n_near], dtype=torch.float64, device=cuda)
        ij = [torch.fft.fftshift(torch.arange(n, dtype=torch.float64, device=cuda)) - n//2 for n in N_near.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        for kl in np.rollaxis(np.indices((2,)*6), 0, 7).reshape(64, 6):
            k, l = kl[:3], kl[3:]
            r = torch.stack([(ij[ind] + k[ind] - l[ind])*self._mesh.dx[ind] for ind in perm], dim=-1)
            N_near[:,:,:] -= (-1)**np.sum(kl) * func_near(r) / (4.*np.pi*np.prod(self._mesh.dx))

        N[ :n_near[0], :n_near[1], :n_near[2],c] = N_near[ :n_near[0], :n_near[1], :n_near[2]]
        N[ :n_near[0], :n_near[1],-n_near[2]:,c] = N_near[ :n_near[0], :n_near[1],-n_near[2]:]
        N[ :n_near[0],-n_near[1]:, :n_near[2],c] = N_near[ :n_near[0],-n_near[1]:, :n_near[2]]
        N[ :n_near[0],-n_near[1]:,-n_near[2]:,c] = N_near[ :n_near[0],-n_near[1]:,-n_near[2]:]
        N[-n_near[0]:, :n_near[1], :n_near[2],c] = N_near[-n_near[0]:, :n_near[1], :n_near[2]]
        N[-n_near[0]:, :n_near[1],-n_near[2]:,c] = N_near[-n_near[0]:, :n_near[1],-n_near[2]:]
        N[-n_near[0]:,-n_near[1]:, :n_near[2],c] = N_near[-n_near[0]:,-n_near[1]:, :n_near[2]]
        N[-n_near[0]:,-n_near[1]:,-n_near[2]:,c] = N_near[-n_near[0]:,-n_near[1]:,-n_near[2]:]

    def _init_N(self):
        if os.path.isfile("cache/Ndipole_%s.pt" % self._mesh):
            N = torch.load("cache/Ndipole_%s.pt" % self._mesh, map_location=cuda)
            logging.info("[DEMAG]: Use cached demag kernel")
        else:
            N = torch.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [6], dtype=torch.float64, device=cuda)

            time_kernel = time()
            for i, t in enumerate(((newell_f,dipole_f,0,1,2),
                                   (newell_g,dipole_g,0,1,2),
                                   (newell_g,dipole_g,0,2,1),
                                   (newell_f,dipole_f,1,2,0),
                                   (newell_g,dipole_g,1,2,0),
                                   (newell_f,dipole_f,2,0,1))):
                self._init_N_component(N, i, t[2:], t[0], t[1])

            logging.info(f"[DEMAG]: Time calculation of demag kernel = {time() - time_kernel} s")
            if not os.path.isdir("cache"):
                os.makedirs("cache")
            torch.save(N, "cache/Ndipole_%s.pt" % self._mesh)

        self._N_fft = torch.fft.rfftn(N, dim = [i for i in range(3) if self._mesh.n[i] > 1])

        # init scratch spaces
        self._m_pad = torch.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [3], dtype=torch.float64, device=cuda)
        self._h_fft = torch.zeros(list(self._N_fft.shape[:3]) + [3], dtype=self._N_fft.dtype, device=cuda)

    @timedmethod
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
