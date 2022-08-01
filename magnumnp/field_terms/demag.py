from magnumnp.common import logging, timedmethod, constants, Timer
import numpy as np
import torch
import torch.fft
from torch import asinh, atan, sqrt, log, abs
from time import time
import os

__all__ = ["DemagField", "newell_f", "newell_g", "dipole_f", "dipole_g"]

def newell_f(points):
    x = abs(points[:,:,:,0])
    y = abs(points[:,:,:,1])
    z = abs(points[:,:,:,2])

    result = 1.0 / 6.0 * (2*x**2 - y**2 - z**2) * sqrt(x**2 + y**2 + z**2)

    mask = (x**2 + z**2).gt(0)
    result[mask] += (y / 2.0 * (z**2 - x**2) * asinh(y / sqrt(x**2 + z**2)))[mask]

    mask = (x**2 + y**2).gt(0)
    result[mask] += (z / 2.0 * (y**2 - x**2) * asinh(z / sqrt(x**2 + y**2)))[mask]

    mask = (x * (x**2 + y**2 + z**2)).gt(0)
    result[mask] -= (x * y * z * atan(y*z / (x * sqrt(x**2 + y**2 + z**2))))[mask]

    return result

def newell_g(points):
    x = points[:,:,:,0]
    y = points[:,:,:,1]
    z = abs(points[:,:,:,2])

    result = - x*y * sqrt(x**2 + y**2 + z**2) / 3.0

    mask = (x**2 + y**2).gt(0) # x**2 + y**2 > 0
    result[mask] += (x*y*z * asinh(z / sqrt(x**2 + y**2)))[mask]

    mask = (y**2 + z**2).gt(0)
    result[mask] += (y / 6.0 * (3.0 * z**2 - y**2) * asinh(x / sqrt(y**2 + z**2)))[mask]

    mask = (x**2 + z**2).gt(0)
    result[mask] += (x / 6.0 * (3.0 * z**2 - x**2) * asinh(y / sqrt(x**2 + z**2)))[mask]

    mask = (z * (x**2 + y**2 + z**2)).ne(0)
    result[mask] -= ( z**3 / 6.0 * atan(x*y / (z * sqrt(x**2 + y**2 + z**2))))[mask]

    mask = (y * (x**2 + y**2 + z**2)).ne(0)
    result[mask] -= (z * y**2 / 2.0 * atan(x*z / (y * sqrt(x**2 + y**2 + z**2))))[mask]

    mask = (x * (x**2 + y**2 + z**2)).ne(0)
    result[mask] -= (z * x**2 / 2.0 * atan(y*z / (x * sqrt(x**2 + y**2 + z**2))))[mask]

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
    @timedmethod
    def __init__(self, state, p = 20):
        self._state = state
        self._mesh = state._mesh
        self._Ms = state._material["Ms"]
        self._p = p
        self._init_N()
        torch.cuda.empty_cache()

    def _init_N_component(self, N, c, perm, func_near, func_far):
        # dipole far-field
        ij = [torch.fft.fftshift(self._state.arange(n)) - n//2 for n in N.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        r = torch.stack([ij[ind]*self._mesh.dx[ind] for ind in perm], dim=-1)
        N[:,:,:,c] = func_far(r) * np.prod(self._mesh.dx) / (4.*np.pi)

        # newell near-field
        n_near = np.minimum(self._mesh.n, self._p)
        N_near = self._state.zeros([1 if i==1 else 2*i for i in n_near])
        ij = [torch.fft.fftshift(self._state.arange(n)) - n//2 for n in N_near.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        for kl in np.rollaxis(np.indices((2,)*6), 0, 7).reshape(64, 6):
            k, l = kl[:3], kl[3:]
            r = torch.stack([(ij[ind] + k[ind] - l[ind])*self._mesh.dx[ind] for ind in perm], dim=-1)
            N_near[:,:,:] -= (-1)**np.sum(kl) * func_near(r) / (4.*np.pi*np.prod(self._mesh.dx))

        N[:n_near[0]   ,:n_near[1]   ,:n_near[2]   ,c] = N_near[:n_near[0]   ,:n_near[1]   ,:n_near[2]   ]
        N[:n_near[0]   ,:n_near[1]   ,-n_near[2]+1:,c] = N_near[:n_near[0]   ,:n_near[1]   ,-n_near[2]+1:]
        N[:n_near[0]   ,-n_near[1]+1:,:n_near[2]   ,c] = N_near[:n_near[0]   ,-n_near[1]+1:,:n_near[2]   ]
        N[:n_near[0]   ,-n_near[1]+1:,-n_near[2]+1:,c] = N_near[:n_near[0]   ,-n_near[1]+1:,-n_near[2]+1:]
        N[-n_near[0]+1:,:n_near[1]   ,:n_near[2]   ,c] = N_near[-n_near[0]+1:,:n_near[1]   ,:n_near[2]   ]
        N[-n_near[0]+1:,:n_near[1]   ,-n_near[2]+1:,c] = N_near[-n_near[0]+1:,:n_near[1]   ,-n_near[2]+1:]
        N[-n_near[0]+1:,-n_near[1]+1:,:n_near[2]   ,c] = N_near[-n_near[0]+1:,-n_near[1]+1:,:n_near[2]   ]
        N[-n_near[0]+1:,-n_near[1]+1:,-n_near[2]+1:,c] = N_near[-n_near[0]+1:,-n_near[1]+1:,-n_near[2]+1:]


    def _init_N(self):
        N = self._state.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [6])

        time_kernel = time()
        for i, t in enumerate(((newell_f,dipole_f,0,1,2),
                               (newell_g,dipole_g,0,1,2),
                               (newell_g,dipole_g,0,2,1),
                               (newell_f,dipole_f,1,2,0),
                               (newell_g,dipole_g,1,2,0),
                               (newell_f,dipole_f,2,0,1))):
            self._init_N_component(N, i, t[2:], t[0], t[1])

        logging.info(f"[DEMAG]: Time calculation of demag kernel = {time() - time_kernel} s")

        N = torch.fft.rfftn(N, dim = [i for i in range(3) if self._mesh.n[i] > 1])
        Nxx = N[:,:,:,0].real.clone()
        Nxy = N[:,:,:,1].real.clone()
        Nxz = N[:,:,:,2].real.clone()
        Nyy = N[:,:,:,3].real.clone()
        Nyz = N[:,:,:,4].real.clone()
        Nzz = N[:,:,:,5].real.clone()

        self._N = [[Nxx, Nxy, Nxz],
                   [Nxy, Nyy, Nyz],
                   [Nxz, Nyz, Nzz]]

    @timedmethod
    def h(self, t, m):
        hx = self._state.zeros(list(self._N[0][0].shape[:3]), dtype=torch.complex128)
        hy = self._state.zeros(list(self._N[0][0].shape[:3]), dtype=torch.complex128)
        hz = self._state.zeros(list(self._N[0][0].shape[:3]), dtype=torch.complex128)
        for ax in range(3):
            with Timer("fft"):
                m_pad_fft1D = torch.fft.rfftn(self._Ms[:,:,:,0] * m[:,:,:,ax], dim = [i for i in range(3) if self._mesh.n[i] > 1], s = [2*self._mesh.n[i] for i in range(3) if self._mesh.n[i] > 1])

            with Timer("multiply"):
                hx[:,:,:] += self._N[0][ax] * m_pad_fft1D[:,:,:]
                hy[:,:,:] += self._N[1][ax] * m_pad_fft1D[:,:,:]
                hz[:,:,:] += self._N[2][ax] * m_pad_fft1D[:,:,:]


        with Timer("ifft"):
            hx = torch.fft.irfftn(hx, dim = [i for i in range(3) if self._mesh.n[i] > 1])
            hy = torch.fft.irfftn(hy, dim = [i for i in range(3) if self._mesh.n[i] > 1])
            hz = torch.fft.irfftn(hz, dim = [i for i in range(3) if self._mesh.n[i] > 1])

        with Timer("stack"):
            h_fft = torch.stack([hx[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2]], \
                                 hy[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2]], \
                                 hz[:self._mesh.n[0],:self._mesh.n[1],:self._mesh.n[2]]], dim=3)
        return h_fft

    def E(self, t, m):
        return - 0.5 * constants.mu_0 * self._mesh.cell_volume * torch.sum(self._Ms * m * self.h(t, m))

    def __str__(self):
        return "demag"
