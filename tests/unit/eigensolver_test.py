import pytest
import torch
from math import sqrt
from magnumnp import *

def test_singlespin_hext():
    hext = 1./constants.mu_0
    n  = (10, 10, 10)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": 1.,}
    state.m = state.Constant([0.,1./sqrt(2.),1./sqrt(2.)])

    external = ExternalField([0.,hext/sqrt(2.),hext/sqrt(2.)])

    eigen = EigenSolver(state, [], [external])
    res = eigen.solve(k=20)
    torch.testing.assert_close(res.omega.abs(), torch.full_like(res.omega, constants.gamma*hext), atol=1e-10, rtol=1e-10)

def test_singlespin_exchange():
    hext = 1./constants.mu_0
    n  = (10, 10, 10)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": 1., "A": 1.3e-11}
    state.m = state.Constant([0.,1./sqrt(2.),1./sqrt(2.)])

    external = ExternalField([0.,hext/sqrt(2.),hext/sqrt(2.)])
    exchange = ExchangeField()

    eigen = EigenSolver(state, [exchange], [external])
    res = eigen.solve(k=20)

    torch.testing.assert_close(res.omega[0].abs(), state.Tensor(constants.gamma*hext), atol=0, rtol=1e-6)

# this test only succeeds if ran as a single test (seems to depend on the random intial value of eigs)
##def test_singlespin_aniso():
##    Ms = 1./constants.mu_0
##    n  = (10, 10, 10)
##    dx = (1e-9, 1e-9, 1e-9)
##    mesh = Mesh(n, dx)
##    state = State(mesh)
##    state.material = {
##            "Ms": Ms,
##            "Ku": 0.5/constants.mu_0,
##            "Ku_axis": (1,0,0),
##            }
##    aniso = UniaxialAnisotropyField()
##
##    state.m = state.Constant([1.,0.,0.])
##
##    eigen = EigenSolver(state, [aniso], [])
##    res = eigen.solve(k=20)
##    print("omega:", res.omega)
##    print("omega:", res.omega[0].abs().numpy(), state.Tensor(constants.gamma*Ms).numpy())
##    torch.testing.assert_close(res.omega.abs(), torch.full_like(res.omega, constants.gamma*Ms), atol=1e-10, rtol=1e-10)

def test_saturated_thinfilm():
    lex = 5.71e-9
    Js = 1.
    n  = (40, 40, 1)
    dx = (2.5e-9, 2.5e-9, 3e-9)

    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {
            "A":lex**2*Js**2/(2.*constants.mu_0),
            "Ms":1./constants.mu_0,
            }
    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField([0.,0.,1.2/constants.mu_0])

    state.m = state.Constant([0.,0.,1.])
    eigen = EigenSolver(state, [demag, exchange], [external])
    res = eigen.solve(k=20, tol=1e-6)
    #print("evals[GHz]:", res.omega.numpy()/2./torch.pi*1e-9)
    #res.save_evecs3D("data/evecs.vti")
    torch.testing.assert_close(res.omega.abs()[:5]/2./torch.pi*1e-9, state.Tensor([8.23468553,10.29218845,10.36022486,12.28226803,13.60508024]), atol=1e-3, rtol=1e-2)

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
