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

from magnumnp.common import logging, constants, write_vti, complex_dtype, timedmethod
import os
import torch
import numpy as np
from scipy.sparse.linalg import LinearOperator, aslinearoperator, eigs
from scipy.interpolate import interp2d
from scipy.linalg import eig
from xml.etree import cElementTree
from xml.dom import minidom

__all__ = ["EigenSolver", "EigenResult"]

class EigenSolver(object):
    def __init__(self, state, linear_terms, constant_terms, domain=Ellipsis):
        """
        This class implements solution of the linearized LLG in the frequncy domain [dAquino2009]_.
        The corresponding eigenvalue problem is solved using Scipy (CPU only).

        *Example*
            .. code:: python

                llg = LLGSolver([demag, exchange, external])
                logger = Logger("data", ['t', 'm'])
                while state.t < 1e-9-eps:
                    llg.step(state, 1e-11)
                    logger << state

        *Arguments*
            terms ([:class:`LLGTerm`])
                List of LLG contributions to be considered for energy minimization
            state ([:class:`State`])
                the State object containing the equilibrium magnetization state.m
            linear_terms (list)
                list of linear field terms
            constant_terms (list)
                list of constant field terms
            domain ([:class:`torch.Tensor`])
                integrate without precession term (default: False)
        """
        self._domain = domain
        self._linear_terms = linear_terms
        self._state = state
        self._m0 = state.m
        self._h0 = torch.sum(sum([term.h(state) for term in self._linear_terms + constant_terms])*self._m0, dim=-1, keepdim=True)

        ez = state.Constant([1e-15,0.,1.])
        self._e1 = torch.linalg.cross(ez, self._m0)
        self._e1 = self._e1 / torch.linalg.norm(self._e1, axis=3, keepdim=True)
        self._e0 = -torch.linalg.cross(self._e1, self._m0)
        self._e0 = self._e0 / torch.linalg.norm(self._e0, axis=3, keepdim=True)
        self._vv = state.Constant([0.,0.], dtype=torch.complex128)

        self._it = 0

    def _C(self, m):
        self._state.m = m
        return sum([term.h(self._state) for term in self._linear_terms])

    def _D0(self, vv):
        self._it += 1
        if self._it % 500 == 0:
            logging.info_blue("[Eigensolver] it= %d" % self._it)

        vv = torch.from_numpy(vv).to(dtype=self._vv.dtype, device=self._state.device)
        vv = vv.reshape(self._vv[self._domain].shape)
        self._vv[...] = 0.
        self._vv[self._domain] = vv

        # apply R
        vvv = self._vv[:,:,:,(0,)]*self._e0 + self._vv[:,:,:,(1,)]*self._e1

        # calculate A0 v = (C+I H0)*v
        hr = self._C(vvv.real)
        hi = self._C(vvv.imag)
        h = hr + 1j*hi
        h -= self._h0*vvv

        # apply B0
        rrr = -constants.gamma * torch.linalg.cross(self._m0.to(dtype = h.dtype), h)

        # apply R^T (reuse vv tensor)
        self._vv[:,:,:,0] = (rrr*self._e0).sum(dim=-1)
        self._vv[:,:,:,1] = (rrr*self._e1).sum(dim=-1)

        return self._vv[self._domain].reshape(-1).detach().cpu().numpy()

    def B0(self, vv):
        return torch.stack([-1j*vv[...,1,:], 1j*vv[...,0,:]], dim=-2)

    @timedmethod
    def solve(self, k=10, tol=1e-6):
        N = np.prod(self._m0[self._domain].shape[:-1])
        D0 = LinearOperator((2*N,2*N), self._D0, dtype=np.complex128)

        evals, evecs2D = eigs(D0, k = 2*k, which = 'SM', tol = tol, v0 = np.ones(2*N))
        #evals, evecs2D = eigs(D0, k = 2*k, sigma = 0, which = 'LM', tol = tol, v0 = np.ones(2*N))

        evalvecs_sorted = sorted(zip(evals.imag,evecs2D.T), key=lambda x: np.abs(x[0]))
        evals = np.array([x[0] for x in evalvecs_sorted if x[0] > 1000.])
        evecs2D = np.array([x[1] for x in evalvecs_sorted if x[0] > 1000.]).transpose()
        evecs2D = torch.from_numpy(evecs2D).reshape(-1,2,evecs2D.shape[-1])

        omega = torch.tensor(evals)

        res = torch.zeros(self._m0.shape[:3] + (2,evecs2D.shape[-1]), dtype=torch.complex128)
        res[self._domain] = evecs2D.reshape(res[self._domain].shape)
        evecs2D = res

        # NOTE: our mode normalization (phi_i,phi_j)_l2 = \delta_ij seems to differ from the published version!
        #       thus, re-normalize w_i (phi_i,phi_j)_B0 = \delta_ij
        # NOTE: Using volume-averaged scalar product: <f,g> = mean_i(f_i * g_i)
        norm = omega * (evecs2D.conj() * self.B0(evecs2D)).sum(dim=3).mean(dim=(0,1,2)).real
        evecs2D /= torch.sqrt(norm).reshape(1, 1, 1, 1, -1)

        return EigenResult(omega, evecs2D, self._state, m0 = self._m0, e0 = self._e0, e1 = self._e1, D0 = D0)


class EigenResult(object):
    def __init__(self, omega, evecs2D, state, **kwargs):
        self._omega = omega
        self._evecs2D = evecs2D
        self._state = state
        self.__dict__.update(kwargs)

    def store(self, filename):
        torch.save({"m0":self.m0, "omega":self._omega, "evecs2D":self._evecs2D}, filename)
        logging.info_green("[Eigensolver] Stored %d eigenvalues to '%s'" % (len(self._omega), filename))

    @staticmethod
    def load(state, filename):
        stored = torch.load(filename, map_location=state.device)
        m0, omega, evecs2D = stored['m0'], stored['omega'], stored['evecs2D']
        state.m = m0

        ez = state.Constant([1e-15,0.,1.])
        e1 = torch.linalg.cross(ez, m0)
        e1 = e1 / torch.linalg.norm(e1, axis=3, keepdim=True)
        e0 = -torch.linalg.cross(e1, m0)
        e0 = e0 / torch.linalg.norm(e0, axis=3, keepdim=True)

        logging.info_green("[Eigensolver] Loaded %d eigenvalues from '%s'" % (len(omega), filename))
        return EigenResult(omega, evecs2D, state, m0 = m0, e0 = e0, e1 = e1)

    @property
    def omega(self):
        return self._omega

    @property
    def freq(self):
        return self._omega/2./torch.pi

    def evecs(self, N = slice(None)):
        vvv = self._evecs2D[:,:,:,(0,),:]*self.e0[:,:,:,:,None] + self._evecs2D[:,:,:,(1,),:]*self.e1[:,:,:,:,None]
        return vvv[...,N]

    def save_evecs3D(self, filename, which = "abs", N = slice(None)):
        if which == "abs":
            op = torch.abs
        elif which == "real":
            op = torch.real
        elif which == "imag":
            op = torch.imag
        else:
            op = which

        freq = self.freq
        evecs = [op(v.squeeze(-1)) for v in self.evecs(N = N).split(1,dim=-1)]
        xmlroot = cElementTree.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
        cElementTree.SubElement(xmlroot, "Collection")
        for i, vvv in enumerate(evecs):
            filename_vti = "%s_%04d.vti" % (os.path.splitext(filename)[0], i)
            write_vti(vvv, filename_vti, self._state)
            cElementTree.SubElement(xmlroot[0], "DataSet", timestep=str(self.freq[i].cpu().numpy()), file=os.path.basename(filename_vti))

        with open(filename, 'w') as fd:
            fd.write(minidom.parseString(" ".join(cElementTree.tostring(xmlroot).decode().replace("\n","").split()).replace("> <", "><")).toprettyxml(indent="  "))

    @property
    def domega(self):
        ### domega_k = alpha * omega_k^2 * ||phi_k||^2   # TODO: add reference!
        return self._omega**2 * (self._state.material["alpha"][...,None] * (self._evecs2D.conj()*self._evecs2D).real).sum(axis=3).mean(axis=(0,1,2))

    def _B0(self, vv):
        return torch.stack([-1j * vv[..., 1, :], 1j * vv[..., 0, :]], dim=-2)


    def projection(self, omega, delta_m):
        """Compute volume-averaged PSD using simplified Lorentzian formula.

        Parameters
        ----------
        delta_m : torch.Tensor
            Real-valued deviation field used to obtain modal amplitudes.
        omega : array_like
            Angular frequency axis (rad/s) at which to evaluate the spectrum.

        Returns
        -------
        torch.Tensor
            Volume-averaged PSD evaluated on ``omega``.
        """
        w = torch.tensor(omega)
        w_k = self.omega.unsqueeze(-1)
        dw_k = self.domega.unsqueeze(-1)

        m2d = self._state.Constant([0., 0.])
        m2d[:,:,:,0] = (delta_m * self.e0).sum(axis=-1)
        m2d[:,:,:,1] = (delta_m * self.e1).sum(axis=-1)
        a_k = (self._omega * (self._evecs2D.conj() * self._B0(m2d.unsqueeze(-1))).sum(axis=3).mean(axis=(0,1,2))).unsqueeze(-1)

        phi2 = (self._state.material["Ms"]**2 * (self._evecs2D.conj()*self._evecs2D).real.sum(axis=3)).mean(axis=(0,1,2)).unsqueeze(-1)

        lorentz_pos = 1.0 / ((w - w_k)**2 + dw_k**2)
        lorentz_neg = 1.0 / ((w + w_k)**2 + dw_k**2)
        lorentz = torch.abs(a_k)**2 * w_k**2 * phi2 * (lorentz_pos + lorentz_neg)
        return lorentz.sum(axis=0) / 2 # P(ω) = ⟨Ms² |δm̂|²⟩/2  [A²/m²]


    def spectrum(self, omega, h_excite):
        """Compute volume-averaged absorbed power using Eq. (39) of d'Aquino & Hertel (JAP 133, 033902 (2023)).

        Parameters
        ----------
        omega : array_like
            Angular frequencies (rad/s) at which to evaluate the absorbed power.
        h_excite : torch.Tensor
            Real-valued excitation field dh^ac(x) (same spatial shape as state.m).

        Returns
        -------
        torch.Tensor
            Volume-averaged absorbed power P_abs(omega)
        """
        # calculate h_k = phi_k^H * R^T * P_m0 * h_excite
        h2d = self._state.Constant([0.,0.])
        h2d[:,:,:,0] = (h_excite*self.e0).sum(axis=-1)
        h2d[:,:,:,1] = (h_excite*self.e1).sum(axis=-1)
        h_k = (self._evecs2D.conj() * h2d[...,None]).sum(axis=3).mean(axis=(0,1,2)).unsqueeze(-1)

        # calculate coefficiencs a_k = (omega_k/(omega_k-omega+i*domega_k) phi_k^H * R^T * P_m0 * h_excite)
        w = torch.tensor(omega)
        w_k = self.omega.unsqueeze(-1)
        w_k_prime = (self.omega + 1j * self.domega).unsqueeze(-1)
        a_k_pos = constants.gamma * h_k * w_k_prime / (w_k_prime - w)
        a_k_neg = constants.gamma * h_k * w_k_prime / (w_k_prime + w)

        # phi2 = Ms² ||phi_k||^2: sum over 2 components, Ms²-weighted mean over space
        phi2 = (self._state.material["Ms"]**2 * (self._evecs2D.conj()*self._evecs2D).real.sum(axis=3)).mean(axis=(0,1,2)).unsqueeze(-1)
        p = 0.5 * phi2 * ((a_k_pos.conj()*a_k_pos).real + (a_k_neg.conj()*a_k_neg).real)
        return p.sum(axis=0) # P(ω) = ⟨Ms² |δm̂|²⟩/2  [A²/m²]


    def absorption(self, omega, h_excite):
        """Compute volume-averaged absorbed power density using Eq. (40) of d'Aquino & Hertel (JAP 133, 033902 (2023)).

        Parameters
        ----------
        omega : array_like
            Angular frequencies (rad/s) at which to evaluate the absorbed power.
        h_excite : torch.Tensor
            Real-valued excitation field dh^ac(x) (same spatial shape as state.m).

        Returns
        -------
        torch.Tensor
            Volume-averaged absorbed power density P_abs(omega) [W/m³]
        """
        # calculate h_k = phi_k^H * R^T * P_m0 * h_excite
        h2d = self._state.Constant([0., 0.])
        h2d[:,:,:,0] = (h_excite*self.e0).sum(axis=-1)
        h2d[:,:,:,1] = (h_excite*self.e1).sum(axis=-1)
        h_k = (self._evecs2D.conj() * h2d[...,None]).sum(axis=3).mean(axis=(0,1,2)).unsqueeze(-1)
        ###h_k2 = (h_k.conj() * h_k).real

        # Ms-weighted projection: g_k = mean_x{ Ms(x) * phi_k^H(x) * h2d(x) }
        Ms = self._state.material["Ms"]
        print("Ms:", Ms.shape)
        print("evecs:", self._evecs2D.conj().shape)
        print("h2d:", h2d[...,None].shape)
        g_k = (Ms[...,None] * self._evecs2D.conj() * h2d[...,None]).sum(axis=3).mean(axis=(0,1,2)).unsqueeze(-1)

        # calculate Pabs = mu_0/2 * Re{i*omega * sum_k g_k^* a_k}
        w = torch.tensor(omega)
        w_k = self.omega.unsqueeze(-1)
        w_k_prime = (self.omega + 1j * self.domega).unsqueeze(-1)
        a_k = constants.gamma * h_k * w_k / (w_k_prime - w)

        Pabs = 0.5 * constants.mu_0 * (1j * w * g_k.conj() * a_k)

        # correct for eigenvector normalization using total-mesh mean (line 130)
        f = (self._state.material["Ms"] > 0).float().mean()
        return Pabs.sum(axis=0).squeeze(0).real / f


    def dispersion(self, points, dx, num_omega=1000):
        vvv = self.evecs()
        m = points(vvv)
        N = m.shape[0]
        kk = 2.*torch.pi * torch.fft.fftshift(torch.fft.fftfreq(N, dx))

        window = torch.hann_window(N, device=m.device)[:,None]
        mfft = torch.abs(torch.fft.fftshift(torch.fft.fft(m * window, dim=0), dim=0))
        mfft = mfft.T

        w = self._omega
        ww = torch.linspace(w.abs().min(), w.abs().max(), num_omega, device=w.device)

        # inline 1D linear interpolation along w
        idx = torch.clamp(torch.searchsorted(w, ww), 1, w.numel() - 1)
        t = (ww - w[idx-1]) / (w[idx] - w[idx-1])
        mm = mfft[idx-1] + (mfft[idx] - mfft[idx-1]) * t.unsqueeze(1)

        return kk, ww, mm
