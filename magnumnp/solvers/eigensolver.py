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
        #evals, evecs2D = eigs(D0, k = 2*k, sigma = 0, which = 'LM', tol = tol)

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
        norm = omega * (evecs2D.conj() * self.B0(evecs2D)).sum(dim=(0,1,2,3)).real
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
        return self._omega**2 * (self._state.material["alpha"][...,None] * (self._evecs2D.conj()*self._evecs2D).real).sum(axis=(0,1,2,3))

    def spectrum(self, omega, h_excite):
        """Compute absorbed power using Eq. (39) of d'Aquino & Hertel (JAP 133, 033902 (2023)).

        Parameters
        ----------
        omega : array_like
            Angular frequencies (rad/s) at which to evaluate the absorbed power.
        h_excite : torch.Tensor
            Real-valued excitation field dh^ac(x) (same spatial shape as state.m).

        Returns
        -------
        torch.Tensor
            Absorbed power P_abs(omega)
        """
        # calculate h_k = phi_k^H * R^T * P_m0 * h_excite
        h2d = self._state.Constant([0.,0.])
        h2d[:,:,:,0] = (h_excite*self.e0).sum(axis=-1)
        h2d[:,:,:,1] = (h_excite*self.e1).sum(axis=-1)
        h_k = (self._evecs2D.conj() * h2d[...,None]).sum(axis=(0,1,2,3)).unsqueeze(-1)

        # calculate coefficiencs a_k = (omega_k/(omega_k-omega+i*domega_k) phi_k^H * R^T * P_m0 * h_excite)
        w = torch.tensor(omega)
        w_k = self.omega.unsqueeze(-1)
        w_k_prime = (self.omega + 1j * self.domega).unsqueeze(-1)
        a_k = constants.gamma * h_k * w_k_prime / (w_k_prime - w)

        phi2 = ((self._evecs2D.conj()*self._evecs2D).real).sum(axis=(0,1,2,3)).unsqueeze(-1)
        p = 0.5 * self._state.mesh.cell_volumes * phi2 * (a_k.conj()*a_k).real
        return 2.*p.sum(axis=0) # consider factor of 2 since only positive eigenfrequencies are stored

    def _B0(self, vv):
        return torch.stack([-1j * vv[..., 1, :], 1j * vv[..., 0, :]], dim=-2)

    def project(self, field):
        """Project a spatial vector field onto the eigenmode basis.

        This is useful to express magnetization deviations or excitation fields in the modal
        coordinates used by the eigensolver.

        Parameters
        ----------
        field : torch.Tensor
            Real-valued vector field with the same spatial shape as ``state.m``.

        Returns
        -------
        torch.Tensor
            Complex modal coefficients c_k satisfying δm ≈ Σ_k c_k φ_k.
        """
        proj2d = self._state.Constant([0., 0.])
        proj2d[:,:,:,0] = (field * self.e0).sum(axis=-1)
        proj2d[:,:,:,1] = (field * self.e1).sum(axis=-1)

        proj2d = proj2d.unsqueeze(-1)
        proj2d = self._B0(proj2d)
        proj2d = proj2d.expand_as(self._evecs2D)

        coeffs = (self._evecs2D.conj() * proj2d).sum(axis=(0,1,2,3))
        return self._omega * coeffs

    def modal_projection_psd(self, delta_m, times, volume_scale=1.0, coeffs=None):
        """Compute the PSD of a modal reconstruction that matches a time-domain ring-down.

        Parameters
        ----------
        delta_m : torch.Tensor
            Real-valued deviation field. Ignored when ``coeffs`` is provided.
        times : torch.Tensor
            1D tensor with uniform time samples.
        volume_scale : float, optional
            Scale factor (typically ``num_cells * cell_volume``) to match the PSD normalization.
        coeffs : torch.Tensor, optional
            Precomputed modal coefficients ``a_k`` to avoid re-projection.

        Returns
        -------
        tuple(np.ndarray, np.ndarray)
            Frequency axis (Hz) and PSD averaged over space with components along the last axis.
        """
        coeffs = self.project(delta_m)

        phi = self.evecs().to(dtype=torch.complex128)
        omega = self.omega.to(dtype=torch.complex128)
        damping = self.domega.to(dtype=torch.complex128)

        tt = times.to(dtype=torch.complex128, device=phi.device)
        time_phase = torch.exp(((-damping) - 1j * omega).unsqueeze(-1) * tt)
        amplitudes = coeffs.to(dtype=torch.complex128, device=phi.device)[:, None] * time_phase

        modal_delta = 2.0 * torch.tensordot(phi, amplitudes, dims=([4], [0])).real
        modal_delta = modal_delta.permute(4, 0, 1, 2, 3).contiguous()
        modal_fft = torch.fft.rfft(modal_delta, dim=0)
        modal_power = (modal_fft.abs()**2).mean(dim=(1,2,3)) * volume_scale

        num_steps = times.shape[0]
        dt = float((times[1] - times[0]).detach().cpu().item())
        freq = np.fft.rfftfreq(num_steps, d=dt)
        return freq, modal_power.detach().cpu().numpy()

    def simple_modal_projection(self, delta_m, freq, volume_scale=1.0, eps=1e-30, coeffs=None):
        """Build a Lorentzian sum directly from modal amplitudes ``a_k``.

        Parameters
        ----------
        delta_m : torch.Tensor
            Real-valued deviation field used to obtain modal amplitudes. Ignored when ``coeffs`` is provided.
        freq : array_like
            Frequency axis in Hz at which to evaluate the spectrum.
        volume_scale : float, optional
            Additional scaling to match PSD normalization (e.g. ``cell_volume``).
        eps : float, optional
            Small positive number to avoid division by zero for undamped modes.
        coeffs : torch.Tensor, optional
            Precomputed modal coefficients ``a_k``.

        Returns
        -------
        np.ndarray
            Scalar PSD evaluated on ``freq``.
        """
        freq = np.asarray(freq, dtype=float)
        omega_axis = 2.0 * np.pi * freq

        if coeffs is None:
            if delta_m is None:
                raise ValueError("Either delta_m or coeffs must be provided")
            coeffs_np = self.project(delta_m).detach().cpu().numpy()
        else:
            coeffs_np = coeffs.detach().cpu().numpy()
        omega = self.omega.detach().cpu().numpy()
        domega = self.domega.detach().cpu().numpy()

        mode_weights = np.abs(coeffs_np)**2
        widths = domega[:, None] + eps
        lorentz = mode_weights[:, None] * widths / ((omega_axis - omega[:, None])**2 + widths**2)
        return (lorentz.sum(axis=0)) * volume_scale

    def absorption(self, omega, h_excite):
        """Compute absorbed power using Eq. (40) of d'Aquino & Hertel (JAP 133, 033902 (2023)).

        Parameters
        ----------
        omega : array_like
            Angular frequencies (rad/s) at which to evaluate the absorbed power.
        h_excite : torch.Tensor
            Real-valued excitation field dh^ac(x) (same spatial shape as state.m).

        Returns
        -------
        torch.Tensor
            Absorbed power P_abs(omega)
        """
        # calculate h_k = phi_k^H * R^T * P_m0 * h_excite
        h2d = self._state.Constant([0., 0.])
        h2d[:,:,:,0] = (h_excite*self.e0).sum(axis=-1)
        h2d[:,:,:,1] = (h_excite*self.e1).sum(axis=-1)
        h_k = (self._evecs2D.conj() * h2d[...,None]).sum(axis=(0,1,2,3)).unsqueeze(-1)
        h_k2 = (h_k.conj() * h_k).real

        # calculate Pabs = 1/(2N) * sum(i*omega*|h_k|^2 / ((w_k - w) + i alpha omega_k |phi_k|^2)
        w = torch.tensor(omega)
        w_k = self.omega.unsqueeze(-1)
        w_k_prime = (self.omega + 1j * self.domega).unsqueeze(-1)

        Pabs_complex = 0.5 / constants.mu_0 * self._state.mesh.cell_volumes * (1j * w * h_k2 * w_k / (w_k_prime - w))
        return (Pabs_complex.sum(axis=0).squeeze(0)).real


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
