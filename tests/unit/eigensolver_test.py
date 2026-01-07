import pytest
import torch
import numpy as np
from math import sqrt
from magnumnp import *

def test_singlespin_hext():
    alpha = 0.008
    hext = 1./constants.mu_0
    n  = (10, 10, 10)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms":1., "alpha":alpha}
    state.m = state.Constant([0.,1./sqrt(2.),1./sqrt(2.)])
    external = ExternalField([0.,hext/sqrt(2.),hext/sqrt(2.)])
    eigen = EigenSolver(state, [], [external])
    res = eigen.solve(k=20)

    w0 = constants.gamma*hext
    dw0 = alpha * w0
    torch.testing.assert_close(res.omega.abs(), torch.full_like(res.omega, w0), atol=1e-10, rtol=1e-10)
    torch.testing.assert_close(res.domega.abs(), torch.full_like(res.omega, dw0), atol=1e-10, rtol=1e-10)

def test_singlespin_exchange():
    hext = 1./constants.mu_0
    alpha = 0.008
    n  = (10, 10, 10)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": 1.,
                      "A": 1.3e-11,
                      "alpha": alpha}
    state.m = state.Constant([0.,1./sqrt(2.),1./sqrt(2.)])
    external = ExternalField([0.,hext/sqrt(2.),hext/sqrt(2.)])
    exchange = ExchangeField()
    eigen = EigenSolver(state, [exchange], [external])
    res = eigen.solve(k=1)

    w0 = torch.tensor(constants.gamma * hext)
    dw0 = alpha * w0
    torch.testing.assert_close(res.omega[0], w0, atol=1e-10, rtol=1e-10)
    torch.testing.assert_close(res.domega[0], dw0, atol=1e-10, rtol=1e-10)

    # check spectrum
    h_amp = 1e-3
    prefactor0 = h_amp**2 / (2.0 * w0**2)
    omega = torch.linspace(0.9 * w0, 1.1 * w0, 5)

    h_excite = state.Constant([h_amp, 0., 0.])
    spectrum_numeric = res.spectrum(omega, h_excite)

    w_prime = torch.complex(w0, dw0)
    w = torch.complex(omega, torch.zeros_like(omega))
    spectrum_expected = prefactor0 * (torch.abs(w_prime)**2 / torch.abs(w_prime - w)**2)

    print("spectrum_numeric:", spectrum_numeric)
    print("spectrum_expected:", spectrum_expected)
    print("diff:", spectrum_numeric-spectrum_expected)
    torch.testing.assert_close(spectrum_numeric, spectrum_expected, atol=0, rtol=1e-10)

def test_singlespin_aniso():
    alpha = 0.008
    Ms = 1./constants.mu_0
    n  = (10, 10, 10)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {
            "Ms": Ms,
            "Ku": 0.5/constants.mu_0,
            "Ku_axis": [1,0,0],
            "alpha": alpha
            }
    aniso = UniaxialAnisotropyField()
    state.m = state.Constant([1.,0.,0.])
    eigen = EigenSolver(state, [aniso], [])
    res = eigen.solve(k=20)

    w0 = constants.gamma*Ms
    dw0 = alpha * w0
    torch.testing.assert_close(res.omega.abs(), torch.full_like(res.omega, w0), atol=1e-10, rtol=1e-10)
    torch.testing.assert_close(res.domega.abs(), torch.full_like(res.omega, dw0), atol=1e-10, rtol=1e-10)

def test_saturated_thinfilm():
    lex = 5.71e-9
    Js = 1.
    n  = (40, 40, 1)
    dx = (2.5e-9, 2.5e-9, 3e-9)

    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {
            "A": state.Constant(lex**2*Js**2/(2.*constants.mu_0)),
            "Ms": state.Constant(1./constants.mu_0),
            }
    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField(state.Constant([0.,0.,1.2/constants.mu_0]))

    state.m = state.Constant([0.,0.,1.])
    eigen = EigenSolver(state, [demag, exchange], [external])
    res = eigen.solve(k=20, tol=1e-6)
    #print("evals[GHz]:", res.omega.numpy()/2./torch.pi*1e-9)
    #res.save_evecs3D("data/evecs.vti")
    torch.testing.assert_close(res.omega.abs()[:5]/2./torch.pi*1e-9, torch.tensor([8.23468553,10.29218845,10.36022486,12.28226803,13.60508024]), atol=1e-3, rtol=1e-2)


def test_absorption():
    hext = 1./constants.mu_0
    alpha = 0.008
    n  = (10, 10, 10)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": 1., "A": 1.3e-11, "alpha": alpha}
    state.m = state.Constant([0.,1./sqrt(2.),1./sqrt(2.)])

    external = ExternalField([0.,hext/sqrt(2.),hext/sqrt(2.)])
    exchange = ExchangeField()
    eigen = EigenSolver(state, [exchange], [external])
    res = eigen.solve(k=1)

    w0 = torch.tensor(constants.gamma * hext)
    dw0 = alpha * w0

    torch.testing.assert_close(res.omega[0], w0, atol=1e-10, rtol=1e-10)
    torch.testing.assert_close(res.domega[0], dw0, atol=1e-10, rtol=1e-10)

    h_ac = 1e-3
    omega = torch.linspace(0.9 * w0, 1.1 * w0, 200)
    h_excite = state.Constant([h_ac, 0.0, 0.0])
    absorption = res.absorption(omega.numpy(), h_excite)

    # |h_k|^2 must only depend on the physical excitation amplitude (Eq. 40 in d'Aquino & Hertel)
    h_k2 = h_ac**2 / (2.0 * w0)
    absorption_analytic = (0.5 * state.mesh.volume * (1j * omega * h_k2 * w0) / ((w0 + 1j * dw0) - omega)).real

    torch.testing.assert_close(absorption, absorption_analytic, atol=0.0, rtol=1e-10)
    assert torch.max(torch.abs(absorption - absorption_analytic)) / absorption_analytic.max() < 1e-10

    # Plotting code (kept for reference but disabled during tests):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(omega / 2.0 / np.pi * 1e-9, absorption, label="numerical", linewidth=2)
    ax.plot(omega / 2.0 / np.pi * 1e-9, absorption_analytic, "--", label="analytic")
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("Absorbed power [J/s]")
    ax.set_title("FMR absorption of a tiny macrospin")
    ax.legend()
    ax.grid(True, linestyle=":", linewidth=0.5)
    fig.tight_layout()
    fig.savefig("result_absorption.png", dpi=150)


#def test_vortex():
#    lex = 5.71e-9
#    Js = 1.
#    size = (200e-9, 200e-9, 3e-9)
#    n  = (40, 40, 1)
#    dx = tuple(np.array(size) / np.array(n))
#
#    mesh = Mesh(n, dx)
#    material = {
##            "A":lex**2*Js**2/(2.*constants.mu_0),
#            "A":lex**2/constants.mu_0,
#            "Ms":Js/constants.mu_0,
#            "gamma":2.211e5,
#            "alpha":0.
#            }
#    print("A:", material["A"])
#    demag    = DemagField(mesh, material)
#    exchange = ExchangeField(mesh, material)
#
#    m0 = np.zeros(n + (3,))
#    m0[:20,:,:,1]   = -1.0
#    m0[20:,:,:,1]   = 1.0
#    m0[20,20,:,1]   = 0.
#    m0[20,20,:,2]   = 1.
#
#    try:
#        m = m0.copy()
#        m[:,:,:,:] = np.load("data/m0_vortex.npy")
#    except:
#        minimizer = Minimizer([demag, exchange])
#        m = minimizer.minimize(m0, 1e-8, 1e-5)
#        write_vtr(m, "data/m0_vortex")
#        np.save("data/m0_vortex.npy", m)
#
#    #eigen = EigenSolver([demag, exchange], [], material, m)
#    eigen = EigenSolver2([demag, exchange], [], material, m)
#    #evals, evecs = eigen.solve(k=20, method="eig")
#    evals, evecs = eigen.solve(k=20, method="eig_finite")
#    print("evals[GHz]:", evals.imag/2./np.pi*1e-9)
#    print("evals[GHz]:", evals.real/2./np.pi*1e-9)
#    e = np.sort(np.abs(evals.imag))
#    print("evals.min[GHz]:", e[:10]/2./np.pi*1e-9)
#    print("evals.min[GHz]:", (e[e>1000.])[:10]/2./np.pi*1e-9)
#    #eigen.save_evecs3D(evecs, "data/evecs_vortex_iterative")
#    #eigen.save_evecs3D(evecs, "data/evecs_vortex_finite")
#    #eigen.save_evecs3D(evecs, "data/evecs_vortex_analytic")
#    return evals, evecs
