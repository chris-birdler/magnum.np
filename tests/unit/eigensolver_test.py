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
    omega0 = constants.gamma*hext
    domega0 = alpha * omega0

    torch.testing.assert_close(res.omega.abs(), torch.full_like(res.omega, omega0), atol=1e-10, rtol=1e-10)
    torch.testing.assert_close(res.domega.abs(), torch.full_like(res.omega, domega0), atol=1e-10, rtol=1e-10)

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
    omega0 = torch.tensor(constants.gamma * hext)
    domega0 = alpha * omega0

    torch.testing.assert_close(res.omega[0], omega0, atol=1e-10, rtol=1e-10)
    torch.testing.assert_close(res.domega[0], domega0, atol=1e-10, rtol=1e-10)

    # check spectrum
    h_amp = 1e-3
    prefactor0 = h_amp**2 / (2.0 * omega0**2)
    omega = torch.linspace(0.9 * omega0, 1.1 * omega0, 5)

    h_excite = state.Constant([h_amp, 0., 0.])
    spectrum_numeric = res.spectrum(omega, h_excite)

    w_prime = torch.complex(omega0, domega0)
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
    omega0 = constants.gamma*Ms
    torch.testing.assert_close(res.omega.abs(), torch.full_like(res.omega, omega0), atol=1e-10, rtol=1e-10)
    torch.testing.assert_close(res.domega.abs(), torch.full_like(res.omega, alpha * omega0), atol=1e-10, rtol=1e-10)

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
    n = (4, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    origin = (-n[0] * dx[0] / 2.0, -n[1] * dx[1] / 2.0, -n[2] * dx[2] / 2.0)
    mesh = Mesh(n, dx, origin=origin)
    state = State(mesh)

    Ms = 800e3
    A = 13e-12
    alpha = 0.01
    state.material = {"Ms": Ms, "A": A, "alpha": alpha}
    state.m = state.Constant([0.0, 0.0, 1.0])

    H_bias = 1e5
    bias = ExternalField([0.0, 0.0, H_bias])
    exchange = ExchangeField()

    magnetic = state.material["Ms"].squeeze(-1) > 0.0
    num_cells = int(magnetic.sum().item())
    volume = num_cells * state.mesh.cell_volumes

    solver = EigenSolver(state, [exchange], [bias])
    result = solver.solve(k=1, tol=1e-10)

    omega_num = result.omega[0]
    domega_num = result.domega[0]
    omega_analytic = torch.tensor(constants.gamma * H_bias, dtype=omega_num.dtype)
    domega_analytic = alpha * omega_analytic

    torch.testing.assert_close(omega_num, omega_analytic, atol=0.0, rtol=1e-12)
    torch.testing.assert_close(domega_num, domega_analytic, atol=0.0, rtol=1e-12)

    freq = np.linspace(1.0e9, 6.0e9, 200)
    omega = 2.0 * np.pi * freq
    h_ac = 1e3
    h_excite = state.Constant([h_ac, 0.0, 0.0])

    abs_numeric = result.absorption(omega, h_excite, magnetic).detach().cpu().numpy()
    h_k_sq = h_ac**2 * num_cells / (2.0 * omega_analytic.item())
    abs_analytic = 0.5 * volume * (1j * omega * h_k_sq * omega_analytic.item()) / (
        (omega_analytic.item() + 1j * domega_analytic.item()) - omega
    )
    abs_analytic = abs_analytic.real

    np.testing.assert_allclose(abs_numeric, abs_analytic, rtol=1e-10, atol=0.0)
    rel_err = np.max(np.abs(abs_numeric - abs_analytic)) / abs_analytic.max()
    assert rel_err < 1e-10

    # Plotting code (kept for reference but disabled during tests):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(freq * 1e-9, abs_numeric, label="numerical", linewidth=2)
    ax.plot(freq * 1e-9, abs_analytic, "--", label="analytic")
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
