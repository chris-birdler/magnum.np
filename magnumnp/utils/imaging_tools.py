#
# This file is part of the magnum.np distribution
# (https://gitlab.com/magnum.np/magnum.np).
# Copyright (c) 2023 magnum.np team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

from magnumnp.common import constants
from magnumnp.common.mesh import Mesh
from magnumnp.common.state import State
from magnumnp.field_terms import DemagField
import torch
from math import sin, cos, pi, sqrt
import numpy as np

__all__ = ["LTEM", "MFM"]

class LTEM(object): # TODO: document and improve interface
    def __init__(self, state, method = "beleggia", comp = 2, voltage = 300e3, kcx = 0.1, kcy = 0.1, theta = 0.):
        self._state = state
        self._mesh = state.mesh
        self._Ms = self._state.material["Ms"]
        self._m = (self._state.m(self._state.t)*self._state.material["Ms"])
        self._h = self._mesh.n[comp] * self._mesh.dx[comp]
        self._dim = [0, 1]
        self._volt = voltage
        self._kcx, self._kcy = kcx, kcy
        self._comp = comp
        self._theta = np.deg2rad(theta)
        self._method = method
        if self._state._dtype == torch.float32:
            self._dtypecompl = torch.complex64
        else:
            self._dtypecompl = torch.complex128

    def _lambda_el(self):
        return 2.*pi*constants.hbar/sqrt(2 * constants.me * self._volt * constants.e * (1 + constants.e * self._volt / (2 * constants.me * constants.c**2)))

    def _k(self):
        dim = [ni for ni in self._mesh.n]
        dim[self._comp] = 1
        ij = [torch.fft.fftshift(torch.fft.fftfreq(dim[ind], self._mesh.dx[ind])) for ind in [0,1,2]]
        dk = [abs(ij[0][1]-ij[0][0]), abs(ij[1][1]-ij[1][0])]
        ij = torch.meshgrid(*ij, indexing = 'ij')
        k = torch.stack(ij, dim = -1).squeeze(dim = 2).to(device = self._state._device, dtype = self._dtypecompl)
        return k, dk

    def _Phim(self):
        m_int = torch.sum(self._m, dim = self._comp)*self._mesh.dx[self._comp]
        Mmn = torch.fft.fftn(m_int, dim = self._dim)

        k, dk = self._k()
        kx = k.real[:,:,0]
        ky = k.real[:,:,1]

        denom = (kx**2. +  ky**2.)/(kx**2. + ky**2. + dk[0]**2 * self._kcx**2. + self._kcy**2 * dk[1]**2. )**2.
        denom = denom.unsqueeze(-1)
        prefactor = 1j*constants.mu_0*constants.e / (2.*pi*constants.hbar)
        Fphim = denom * prefactor * torch.linalg.cross(Mmn, k)[:,:,2].unsqueeze(-1)
        Phim = torch.fft.ifftn(Fphim, dim = self._dim).real
        return Phim

    def Defocus(self, df, cs):
        k, dk = self._k()
        kx = k.real[:,:,0]
        ky = k.real[:,:,1]
        ks = kx**2. + ky**2.
        if self._method == "mansuripur":
            phim = self._Phim_Mansuripur()
        else:
            phim = self._Phim()

        cts = - df + 0.5 * self._lambda_el() ** 2. * cs * ks
        exp = torch.exp(pi * cts * 1j * ks * self._lambda_el()).unsqueeze(dim=-1)
        def_wf_cts = torch.fft.ifft2(torch.fft.fftshift(torch.fft.fft2(torch.exp(phim * 1j), dim = self._dim)) * exp, dim = self._dim)
        return (torch.conj(def_wf_cts) * def_wf_cts).real

    def MagneticPhaseShift(self):
        return torch.rad2deg(self._Phim())

    def InductionMap(self):
        if self._method == "mansuripur":
            phim = self._Phim_Mansuripur()
        else:
            phim = self._Phim()
        dphidx = torch.gradient(phim, dim = 0, spacing = self._mesh.dx[0])[0]
        dphidy = torch.gradient(phim, dim = 1, spacing = self._mesh.dx[1])[0]
        return constants.hbar/(constants.e*self._h)*torch.stack([-dphidy, dphidx], dim = -1)

    def _Phim_Mansuripur(self):
        m_int = torch.sum(self._m, dim = self._comp)*self._mesh.dx[self._comp]
        Mmn = torch.fft.fftn(m_int, dim = self._dim)
        p = self._state.Tensor([0., sin(self._theta), cos(self._theta)], dtype = self._dtypecompl).unsqueeze(dim=0).unsqueeze(dim=0)
        ez = self._state.Tensor([0, 0, 1], dtype = self._dtypecompl).unsqueeze(dim=0).unsqueeze(dim=0)
        k, dk = self._k()
        ks = torch.linalg.norm(k, dim = -1, keepdim = True)# + (dk[0]**2 * self._kcx**2. + dk[1]**2 * self._kcy**2.)
        ks[0,0, :] = 1e-15
        kbool = (ks == 0)
        ks[kbool] = 1e-15

        pxpxM = torch.linalg.cross(p, torch.linalg.cross(p, Mmn))
        khat = k/ks
        Gp = 1./(torch.sum(p * khat, dim = -1, keepdim = True)**2. + p[:,:,2]**2.) * torch.sinc(self._h * ks * torch.sum(p * khat, dim = -1, keepdim = True)/ p[:,:,2])
        func = self._h / ks * Gp * (torch.sum(torch.linalg.cross(khat, ez) * pxpxM, dim = -1, keepdim = True))  * 1j
        phim = 2. * constants.e / constants.hbar * torch.fft.ifftn(func, dim = self._dim).real
        return phim


class MFM(object): # TODO: document and improve interface
    def __init__(self, state, height = 100e-9, q = 1, k =1):
        self._state = state
        self._m = self._state.m(self._state.t)
        self._t = self._state.t
        self._h = height
        self._q = q
        self._k = k

        n = self._state.mesh.n
        dx = self._state.mesh.dx

        nz = int(n[2] + self._h/dx[2])

        n_new = [n[0], n[1], nz]

        mesh_new = Mesh(n_new, dx)
        self._state0 = State(mesh_new)
        self._state0.material["Ms"] = self._state.material["Ms"]
        self._state0.m = self._state0.Constant([0, 0, 0])
        self._state0.m[:,:,:n[2],:] = self._m

    def _Hdemag(self):
        demag = DemagField()
        return demag.h(self._state0)

    def PhaseShift(self, mm_tip = None, dm_tip = None):
        dHszdz = torch.gradient(self._Hdemag(), dim = [2], spacing = self._state.mesh.dx[2])[0][:,:,:,2]
        d2Hsdz2 = torch.gradient(torch.stack(torch.gradient(self._Hdemag(), dim = [2], spacing = self._state.mesh.dx[2]), dim = -1), dim = 2, spacing = self._state.mesh.dx[2])[0]

        if mm_tip is not None:
            phi = self._q*constants.mu_0/self._k*mm_tip*dHszdz
        elif dm_tip is not None:
            phi = self._q*constants.mu_0/self._k*(dm_tip*d2Hsdz2).sum(dim = -1, keepdim = True).squeeze(-1)
        else:
            raise TypeError("Provide Monopole Moment mm_tip, or Dipolar Moment Vector dm_tip")

        return phi[:,:,-1,:].unsqueeze(2)
