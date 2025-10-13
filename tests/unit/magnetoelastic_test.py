import pytest
import torch
from magnumnp import *
import numpy as np

def test_magnetoelastic_field_regression():
    n = (200, 10, 10)
    dx = (4e-9, 2e-9, 2e-9)
    mesh = Mesh(n, dx)

    state = State(mesh)

    rho = 7000.
    C = C_cubic(2e11, 1.5e11, 1e11)
    l_m = 5e-5

    c = np.sqrt(C[0,0]/rho)
    k = 4e9/c

    state.material = {
        "Ms" : state.Constant(1./constants.mu_0),
        "C": state.Constant(C),
        "A" : state.Constant(1.),
        "lambda_100" : state.Constant(l_m),
        "lambda_111" : state.Constant(l_m)
        }

    x, y, z = mesh.SpatialCoordinate()
    state.ud = state.Constant((0,0,0))
    state.ud[:,:,:,0] = torch.sin(k*x)

    state.m = state.Constant((1,0,0))
    state.m[:,:,:,1] = 0.1*torch.sin(k*x)
    state.m[:,:,:,0] = 0.1*torch.cos(k*x)
    normalize(state.m)

    # Field regression test
    # ---------------------

    elastic = MagnetoElasticField()

    h = elastic.h(state).abs().sum()
    assert h.cpu().numpy() == pytest.approx(4.883955064603125e+17, rel=1e-5)

def test_magnetoelastic_field_equivalence():
    # Test the equivalence of the magnetoelastic field terms for isotropic materials

    n = (1,1,1)
    dx = (2e-9, 2e-9, 2e-9)

    mesh = Mesh(n,dx)
    state = State(mesh)

    eps0 = 1e-10 / 1e-6 # 10 nm on a distance of one micron
    state.eps = state.Constant((0.,0.,0.,0.,0.,0.))
    state.eps[...,0] = eps0
    state.eps[...,1] = 0.7*eps0
    state.eps[...,3] = 0.1*eps0
    state.eps[...,5] = 0.3*eps0

    # material is assumed to be isotropic
    nu = 0.3
    E = 2e11
    C = C_isotropic(E, nu)

    l100 = 20e-6
    l111 = l100

    state.material = {"Ms" : state.Constant(1./constants.mu_0),
                      "C" : state.Constant(C),
                      "lambda_100" : state.Constant(l100),
                      "lambda_111" : state.Constant(l111)}

    # init magnetization
    phi = 35.*np.pi / 180.
    theta = -60.*np.pi / 180.
    mx = np.cos(phi)*np.sin(theta)
    my = np.sin(phi)*np.sin(theta)
    mz = np.cos(theta)
    state.m = state.Constant((mx, my, mz))

    # get fields
    hFull = MagnetoElasticField(mechanical_strain=state.eps)
    hLinear = LinearMagnetoElasticField(mechanical_strain=state.eps)

    # compute torques
    T1 = torch.linalg.cross(hFull.h(state), state.m)
    T2 = torch.linalg.cross(hLinear.h(state), state.m)

    # compare torques
    torch.testing.assert_close(T1, T2, rtol=1e-4, atol=0.)

def test_set_constant():
    # Set Strain / Set Displacement / Compare fields
    # ----------------------------------------------

    n = (35, 3, 1)
    dx = (3e-9, 5e-9, 2e-9)

    Lx = n[0]*dx[0]
    Ly = n[1]*dx[1]
    kx = 2.*np.pi / Lx

    mesh = Mesh(n, dx, pbc=(1,0,1))
    state = State(mesh)

    x,y,z = mesh.SpatialCoordinate()

    u0 = 1e-11
    ud = torch.zeros(n+(3,))
    ud[...,0] = u0*torch.sin(kx*x)
    ud[...,2] = u0*y/Ly

    eps = torch.zeros(n+(6,))
    eps[...,0] = kx*u0*torch.cos(kx*x)
    eps[...,3] = u0/Ly

    state.m = state.Constant((1,0,0))
    state.m[:,:,:,1] = 0.1*torch.sin(kx*x)
    state.m[:,:,:,0] = 0.4*torch.cos(kx*x)
    normalize(state.m)

    C = C_cubic(2e11, 1.5e11, 1e11)
    l_m = 5e-5

    state.material = {
        "Ms" : state.Constant(1./constants.mu_0),
        "C": state.Constant(C),
        "A" : state.Constant(1.),
        "lambda_100" : state.Constant(l_m),
        "lambda_111" : state.Constant(0.4*l_m)
    }

    # set displacement in the constructor
    term_1_set_ud = LinearMagnetoElasticField(ud=ud)
    term_2_set_ud = MagnetoElasticField(ud=ud)
    h1_set_ud = torch.sum(torch.abs(term_1_set_ud.h(state)))
    h2_set_ud = torch.sum(torch.abs(term_2_set_ud.h(state)))

    # set strain in the constructor
    # ... testing also the strain computation by comparison to this case
    term_1_set_strain = LinearMagnetoElasticField(mechanical_strain=eps)
    term_2_set_strain = MagnetoElasticField(mechanical_strain=eps)
    h1_set_strain = torch.sum(torch.abs(term_1_set_strain.h(state)))
    h2_set_strain = torch.sum(torch.abs(term_2_set_strain.h(state)))

    # displacement as state member
    state.ud = ud
    term_1 = LinearMagnetoElasticField()
    term_2 = MagnetoElasticField()
    h1 = torch.sum(torch.abs(term_1.h(state)))
    h2 = torch.sum(torch.abs(term_2.h(state)))

    torch.testing.assert_close(h1, h1_set_ud, rtol=1e-6, atol=0)
    torch.testing.assert_close(h2, h2_set_ud, rtol=1e-6, atol=0)

    torch.testing.assert_close(h1, h1_set_strain, rtol=1e-2, atol=0)
    torch.testing.assert_close(h2, h2_set_strain, rtol=1e-2, atol=0)

def test_set_function():
    # Set Strain and Displacement as Function of State
    # ------------------------------------------------

    eps0 = 1e-5

    n = (4,1,1)
    dx = (3e-9, 3e-9, 3e-9)
    mesh = Mesh(n, dx, pbc=(0,1,1))
    state = State(mesh)

    # isotropic material
    # (inconsequential since linear field and more general version are NOT compared directly)
    C = C_isotropic(2e11, 0.3)
    C11 = C[0,0]
    C12 = C[0,1]
    l_m = 5e-5

    state.material = {
        "Ms" : state.Constant(1./constants.mu_0),
        "C": state.Constant(C),
        "A" : state.Constant(1.),
        "lambda_100" : state.Constant(l_m),
        "lambda_111" : state.Constant(l_m)
    }

    # magnetic strain only in x-direction
    state.m = state.Constant((1.,0.,0.))

    # functions passed as argument during the constructor call
    def f_ud(state):
        n = state.mesh.n
        x = state.mesh.SpatialCoordinate()[0]

        v = torch.zeros(n+(3,))
        v[...,0] = eps0*x
        return state.t* v

    def f_strain(state):
        n = state.mesh.n
        v = torch.zeros(n+(6,))
        v[...,0] = eps0

        return state.t * v

    # init field terms
    term_1_ud = LinearMagnetoElasticField(ud=f_ud)
    term_2_ud = MagnetoElasticField(ud=f_ud)

    term_1_eps = LinearMagnetoElasticField(mechanical_strain=f_strain)
    term_2_eps = MagnetoElasticField(mechanical_strain=f_strain)

    # test correct function referencing
    for i in range(3):
        state.t = i

        t1 = torch.sum(term_1_ud.h(state))
        t2 = torch.sum(term_2_ud.h(state))

        t3 = torch.sum(term_1_eps.h(state))
        t4 = torch.sum(term_2_eps.h(state))

        B1 = -3.*l_m*(C11-C12)/2.
        f1 = -2.*B1*i*eps0*n[0]
        f1 = torch.tensor(f1)

        sig_m = (C11-C12)*l_m
        sig_el = (C11*i*eps0 - sig_m)*n[0]
        f2 = 3.*l_m*sig_el
        f2 = torch.tensor(f2)

        torch.testing.assert_close(t1, f1, rtol=1e-6, atol=0)
        torch.testing.assert_close(t2, f2, rtol=1e-6, atol=0)
        torch.testing.assert_close(t3, f1, rtol=1e-6, atol=0)
        torch.testing.assert_close(t4, f2, rtol=1e-6, atol=0)