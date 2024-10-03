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

from magnumnp.common import constants, complex_dtype
from magnumnp.common.mesh import Mesh
from magnumnp.common.state import State
from magnumnp.field_terms import DemagField
import torch
from math import sin, cos, pi, sqrt
import numpy as np

__all__ = ["LTEM", "MFM", "to_discretisedfield"]

class LTEM(object): # TODO: document and improve interface
    def __init__(self, state, method = "beleggia", comp = 2, voltage = 300e3, kcx = 0.1, kcy = 0.1, theta = 0.):
        self._state = state
        self._mesh = state.mesh
        self._Ms = self._state.material["Ms"]
        self._M = self._state.m * self._state.material["Ms"]
        self._h = self._mesh.n[comp] * self._mesh.dx[comp]
        self._dim = [0, 1]
        self._volt = voltage
        self._kcx, self._kcy = kcx, kcy
        self._comp = comp
        self._theta = np.deg2rad(theta)
        self._method = method

    def _lambda_el(self):
        return 2.*pi*constants.hbar/sqrt(2 * constants.me * self._volt * constants.e * (1 + constants.e * self._volt / (2 * constants.me * constants.c**2)))

    def _k(self):
        dim = [ni for ni in self._mesh.n]
        dim[self._comp] = 1
        ij = [torch.fft.fftshift(torch.fft.fftfreq(dim[ind], self._mesh.dx[ind])) for ind in [0,1,2]]
        dk = [abs(ij[0][1]-ij[0][0]), abs(ij[1][1]-ij[1][0])]
        ij = torch.meshgrid(*ij, indexing = 'ij')
        k = torch.stack(ij, dim = -1).squeeze(dim = 2).to(dtype = complex_dtype[self._state.dtype])
        return k, dk

    def _Phim(self):
        m_int = torch.sum(self._M, dim = self._comp)*self._mesh.dx[self._comp]
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
        m_int = torch.sum(self._M, dim = self._comp)*self._mesh.dx[self._comp]
        Mmn = torch.fft.fftn(m_int, dim = self._dim)
        p = torch.tensor([0., sin(self._theta), cos(self._theta)], dtype = complex_dtype[state.dtype]).unsqueeze(dim=0).unsqueeze(dim=0)
        ez = torch.tensor([0.,0.,1.], dtype = complex_dtype[state.dtype]).unsqueeze(dim=0).unsqueeze(dim=0)
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


class MFM(object):
    def __init__(self, height = 100e-9, Q = 1, k = 1, mm_tip = 0., dm_tip = None):
        r""" Calculation of the phase shift of an MFM tip.

        The contrast in MFM images originates from the magnetic interaction between
        the magnetic tip of an oscillating cantilever and the samples stray field.
        As MFM is based on the stray field outside of a sample, an 'airbox' method
        is used. The saturation magnetisation has to be set to zero outside of the sample
        in which we wish to perform these MFM measurements.

        The magnetic cantilever is driven to oscillate near its resonant frequency
        when there is no stray field. In the presence of a stray magnetic field the
        phase shift of MFM cantilever is given by

        .. math::
            \Delta \phi = \frac{Q\mu_0}{k} \left( q \frac{\partial
                          {\bf H}_{sz}}{\partial z} + {\bf M}_t \cdot
                          \frac{\partial^2{\bf H}_{s}}{\partial z^2} \right),

        NOTE: this code is based on the MFM implementation of Ubermag
              (see also: https://ubermag.github.io/documentation/notebooks/mag2exp/Magnetic_Force_Microscopy.html)

        *Example*
          .. code:: python

            state = State(mesh)
            state.material = {"Ms": 8e5}
            state.material["Ms"][~magnetic] = 0. # zero Ms outside of sample

            x, y, z = meste.SpatialCoordinate()
            state.m = torch.stack([y, -x, 0*z], dim=-1)
            state.m.normalize()

            demag = DemagField()
            mfm = MFM(height=10e-9, mm_tip = 10e-9, dm_tip=20e-9)

            mfm_logger = FieldLogger("data/phi.pvd", [mfm.PhaseShift])
            mfm_logger << mfm.extended_state
        """
        self._prefactor = Q*constants.mu_0/k
        self._mm_tip = mm_tip
        self._dm_tip = dm_tip
        self._height = height
        self._xdemag = DemagField()

    def extend_demag(self, state):
        n = state.mesh.n
        if not hasattr(self, "extended_state"):
            dx = state.mesh.dx

            nz = int(n[2] + self._height/dx[2])
            n_new = [n[0], n[1], nz]

            mesh_new = Mesh(n_new, dx)
            self.extended_state = State(mesh_new)
            self.extended_state.material["Ms"] = 0.
            self.extended_state.m = self.extended_state.Constant([0.,0.,00.])

        # copy data
        self.extended_state.material["Ms"][:,:,:n[2]] = state.material["Ms"]
        self.extended_state.m[:,:,:n[2],:] = state.m
        return self._xdemag.h(self.extended_state)

    def PhaseShift(self, state, mm_tip = None, dm_tip = None):
        h_demag = self.extend_demag(state)
        spacing = state.mesh.dx[2]
        mm_tip = mm_tip or self._mm_tip
        dm_tip = dm_tip or self._dm_tip

        # monopole
        dHsdz = torch.gradient(h_demag, dim = [2], spacing = spacing)[0]
        phi = mm_tip * dHsdz[:,:,:,2]

        # dipole
        if dm_tip != None:
            d2Hsdz2 = torch.gradient(dHsdz, dim = [2], spacing = spacing)[0]
            phi += (dm_tip * d2Hsdz2).sum(dim=-1)

        return phi.unsqueeze(-1)


def to_discretisedfield(data):
    '''
    Utility function to convert vector fields to Ubermag's Field object.

    This allows using all sorts of post-processing features provided by Ubermag.
    The code has been contributed by @AbsoluteVacuum (see #33 in gitlab)
    '''
    from discretisedfield import Field
    from xarray import DataArray
    xarr = DataArray(
        data.detach().cpu().numpy(),
        dims=["x", "y", "z", "comp"],
        coords=dict(
            x=torch.arange(nx).numpy()*dx,
            y=torch.arange(ny).numpy()*dy,
            z=torch.arange(nz).numpy()*dz,
            comp=["x", "y", "z"],
        ),
        name="mag",
        attrs=dict(cell=mesh.dx, p1=mesh.origin, p2=[a*b+c for a,b,c in zip(mesh.n, mesh.dx, mesh.origin)])
    )
    return Field.from_xarray(xarr)
