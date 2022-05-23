from magnumnp.common import logging, timedmethod, constants
import numpy as np
import torch
import torch.fft
from torch import asinh, atan, sqrt, log, abs
import os
from time import time

__all__ = ["OerstedField"]

def krueger_g(points):
    x = points[:,:,:,0]
    y = points[:,:,:,1]
    z = points[:,:,:,2]

    R = sqrt(x**2 + y**2 + z**2)

    res = (3.*x**2 + 3.*y**2 - 2.*z**2)*z*R/24.
    res += np.pi*z/4.*abs(x*y*z)

    mask = (x**2 + y**2).gt(0)
    res[mask] += ((x**4 - 6.*x**2*y**2 + y**4)/24. * log(z+R))[mask]

    mask = (y**2).gt(0)
    res[mask] += (x*y/6. * (y**2 - 3.*z**2) * atan(x*z/(y*R)))[mask]

    mask = (x**2).gt(0)
    res[mask] += (x*y/6. * (x**2 - 3.*z**2) * atan(y*z/(x*R)))[mask]

    mask = (y**2+z**2).gt(0)
    res[mask] += (z/6. * x * (z**2 - 3.*y**2) * log(x+R))[mask]

    mask = (x**2+z**2).gt(0)
    res[mask] += (z/6. * y * (z**2 - 3.*x**2) * log(y+R))[mask]

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
    def __init__(self, state, p = 15):
        self._state = state
        self._mesh = state._mesh
        self._p = p
        self._init_K()
        torch.cuda.empty_cache()

    def _init_K_component(self, K, c, perm, func_near, func_far):
        # dipole far-field
        ij = [torch.fft.fftshift(torch.arange(n, dtype=torch.float64, device=cuda)) - n//2 for n in K.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        r = torch.stack([ij[ind]*self._mesh.dx[ind] for ind in perm], dim=-1)
        K[:,:,:,c] = func_far(r) * np.prod(self._mesh.dx) / (4.*np.pi)

        # newell near-field
        n_near = np.minimum(self._mesh.n, self._p)
        K_near = self._state.zeros([1 if i==1 else 2*i for i in n_near])
        ij = [torch.fft.fftshift(self._state.arange(n)) - n//2 for n in K_near.shape[:3]]
        ij = torch.meshgrid(*ij,indexing='ij')

        for k in np.rollaxis(np.indices((3,)*3), 0, 4).reshape(27, -1) - 1:
            r = torch.stack([(ij[ind] + k[ind])*self._mesh.dx[ind] for ind in perm], dim=-1)
            K_near[:,:,:] += np.prod(2.-3*np.abs(k)) * func_near(r) / (4.*np.pi*np.prod(self._mesh.dx))

        K[ :n_near[0], :n_near[1], :n_near[2],c] = K_near[ :n_near[0], :n_near[1], :n_near[2]]
        K[ :n_near[0], :n_near[1],-n_near[2]:,c] = K_near[ :n_near[0], :n_near[1],-n_near[2]:]
        K[ :n_near[0],-n_near[1]:, :n_near[2],c] = K_near[ :n_near[0],-n_near[1]:, :n_near[2]]
        K[ :n_near[0],-n_near[1]:,-n_near[2]:,c] = K_near[ :n_near[0],-n_near[1]:,-n_near[2]:]
        K[-n_near[0]:, :n_near[1], :n_near[2],c] = K_near[-n_near[0]:, :n_near[1], :n_near[2]]
        K[-n_near[0]:, :n_near[1],-n_near[2]:,c] = K_near[-n_near[0]:, :n_near[1],-n_near[2]:]
        K[-n_near[0]:,-n_near[1]:, :n_near[2],c] = K_near[-n_near[0]:,-n_near[1]:, :n_near[2]]
        K[-n_near[0]:,-n_near[1]:,-n_near[2]:,c] = K_near[-n_near[0]:,-n_near[1]:,-n_near[2]:]

    def _init_K(self):
        if os.path.isfile("cache/K%s.pt" % self._mesh):
            K = torch.load("cache/K%s.pt" % self._mesh, map_location=self._state._device)
            logging.info("[DEMAG]: Use cached oersted kernel")
        else:
            K = self._state.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [3], dtype=torch.float64, device=cuda)

            time_kernel = time()
            for i, t in enumerate(((krueger_g,dipole_g,0,1,2),
                                   (krueger_g,dipole_g,1,2,0),
                                   (krueger_g,dipole_g,2,0,1))):
                self._init_K_component(K, i, t[2:], t[0], t[1])

            logging.info(f"[OERSTED]: Time calculation of oersted kernel = {time() - time_kernel} s")
            if not os.path.isdir("cache"):
                os.makedirs("cache")
            torch.save(K, "cache/K%s.pt" % self._mesh)

        self._K_fft = torch.fft.rfftn(K, dim = [i for i in range(3) if self._mesh.n[i] > 1])

        # init scratch spaces
        self._j_pad = self._state.zeros([1 if i==1 else 2*i for i in self._mesh.n] + [3])
        self._h_fft = self._state.zeros(list(self._K_fft.shape[:3]) + [3], dtype=self._K_fft.dtype)

    @timedmethod
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
