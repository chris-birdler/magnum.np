import pytest
import pathlib
import torch
from math import pi, cos, sin
from magnumnp import *

def test_rkky():
    J_rkky = -0.001
    n  = (10, 10, 2)
    dx = (1e-9, 1e-9, 1e-9)
    L  = (n[0]*dx[0], n[1]*dx[1], n[2]*dx[2])
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {"Ms": state.Constant(1./constants.mu_0),
                      "A": state.Constant(1e-11)}
    state.m = state.Constant([1,0,0])

    domain1 = state.Constant(False, dtype=torch.bool)
    domain1[:,:,0] = True

    domain2 = state.Constant(False, dtype=torch.bool)
    domain2[:,:,1] = True

    exchange1 = ExchangeField(domain1)
    exchange2 = ExchangeField(domain2)
    rkky = RKKYField(J_rkky, "z", 0, 1)

    for phi in torch.linspace(0,2*pi,36):
        state.m[domain2] = torch.tensor([cos(phi), sin(phi), 0])
        h = rkky.h(state)

        h1 = torch.linalg.cross(state.avg(state.m[domain1]), state.avg(h[domain1]))
        h2 = torch.linalg.cross(state.avg(state.m[domain2]), state.avg(h[domain2]))

        h1_analytic = torch.linalg.cross(state.avg(state.m[domain1]), J_rkky/dx[2]*state.avg(state.m[domain2]))
        h2_analytic = torch.linalg.cross(state.avg(state.m[domain2]), J_rkky/dx[2]*state.avg(state.m[domain1]))

        E = rkky.E(state).detach().cpu().numpy()
        E_analytic = -J_rkky*L[0]*L[1]*cos(phi)

        torch.testing.assert_close(h1, h1_analytic, atol=1e-10, rtol=1e-10)
        torch.testing.assert_close(h2, h2_analytic, atol=1e-10, rtol=1e-10)
        assert E == pytest.approx(E_analytic, rel=1e-10, abs=1e-25)

def test_rkky_biquadratic():
    J_rkky = 0.001
    n  = (10, 10, 2)
    dx = (1e-9, 1e-9, 1e-9)
    L  = (n[0]*dx[0], n[1]*dx[1], n[2]*dx[2])
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {"Ms":1./constants.mu_0, "A":1e-11}
    state.m = state.Constant([1,0,0])

    domain1 = state.Constant(False, dtype=torch.bool)
    domain1[:,:,0] = True

    domain2 = state.Constant(False, dtype=torch.bool)
    domain2[:,:,1] = True

    exchange1 = ExchangeField(domain1)
    exchange2 = ExchangeField(domain2)
    rkky = BiquadraticRKKYField(J_rkky, "z", 0, 1)

    for phi in torch.linspace(0,2*pi,36):
        state.m[domain2] = torch.tensor([cos(phi), sin(phi), 0])
        h = rkky.h(state)

        h1 = state.avg(h[domain1])
        h2 = state.avg(h[domain2])

        m1 = state.avg(state.m[domain1])
        m2 = state.avg(state.m[domain2])

        h1_analytic = 2.*J_rkky/dx[2]*(m1*m2).sum()*state.avg(state.m[domain2])
        h2_analytic = 2.*J_rkky/dx[2]*(m1*m2).sum()*state.avg(state.m[domain1])

        E = rkky.E(state).detach().cpu().numpy()
        E_analytic = -J_rkky*L[0]*L[1]*cos(phi)**2

        torch.testing.assert_close(h1, h1_analytic, atol=1e-10, rtol=1e-10)
        torch.testing.assert_close(h2, h2_analytic, atol=1e-10, rtol=1e-10)
        assert E == pytest.approx(E_analytic, rel=1e-10, abs=1e-25)

def test_rkky_nonequi():
    J_rkky = -0.001
    n  = (10, 10, 2)
    dx = (1e-9, 1e-9, [5e-9, 5e-9])
    L  = (n[0]*dx[0], n[1]*dx[1], n[2]*dx[2])
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {"Ms": state.Constant(1./constants.mu_0),
                      "A": state.Constant(1e-11)}
    state.m = state.Constant([1,0,0])

    domain1 = state.Constant(False, dtype=torch.bool)
    domain1[:,:,0] = True

    domain2 = state.Constant(False, dtype=torch.bool)
    domain2[:,:,1] = True

    exchange1 = ExchangeField(domain1)
    exchange2 = ExchangeField(domain2)
    rkky = RKKYField(J_rkky, "z", 0, 1)

    for phi in torch.linspace(0,2*pi,36):
        state.m[domain2] = torch.tensor([cos(phi), sin(phi), 0])
        h = rkky.h(state)

        h1 = torch.linalg.cross(state.m[0,0,0,:], h[0,0,0,:])
        h2 = torch.linalg.cross(state.m[0,0,1,:], h[0,0,1,:])

        h1_analytic = torch.linalg.cross(state.m[0,0,0,:], J_rkky/dx[2][0]*state.m[0,0,1,:])
        h2_analytic = torch.linalg.cross(state.m[0,0,1,:], J_rkky/dx[2][1]*state.m[0,0,0,:])

        E = rkky.E(state).detach().cpu().numpy()
        E_analytic = -J_rkky*L[0]*L[1]*cos(phi)

        torch.testing.assert_close(h1, h1_analytic, atol=1e-10, rtol=1e-10)
        torch.testing.assert_close(h2, h2_analytic, atol=1e-10, rtol=1e-10)
        assert E == pytest.approx(E_analytic, rel=1e-10, abs=1e-25)

def test_rkky_biquadratic_nonequi():
    J_rkky = 0.001
    n  = (10, 10, 2)
    dx = (1e-9, 1e-9, [5e-9, 5e-9])
    L  = (n[0]*dx[0], n[1]*dx[1], n[2]*dx[2])
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {"Ms":1./constants.mu_0, "A":1e-11}
    state.m = state.Constant([1,0,0])

    domain1 = state.Constant(False, dtype=torch.bool)
    domain1[:,:,0] = True

    domain2 = state.Constant(False, dtype=torch.bool)
    domain2[:,:,1] = True

    exchange1 = ExchangeField(domain1)
    exchange2 = ExchangeField(domain2)
    rkky = BiquadraticRKKYField(J_rkky, "z", 0, 1)

    for phi in torch.linspace(0,2*pi,36):
        state.m[domain2] = torch.tensor([cos(phi), sin(phi), 0])
        h = rkky.h(state)

        h1 = h[0,0,0,:]
        h2 = h[0,0,1,:]

        m1 = state.m[0,0,0,:]
        m2 = state.m[0,0,1,:]

        h1_analytic = 2.*J_rkky/dx[2][0]*(m1*m2).sum()*m2
        h2_analytic = 2.*J_rkky/dx[2][1]*(m1*m2).sum()*m1

        E = rkky.E(state).detach().cpu().numpy()
        E_analytic = -J_rkky*L[0]*L[1]*cos(phi)**2

        torch.testing.assert_close(h1, h1_analytic, atol=1e-10, rtol=1e-10)
        torch.testing.assert_close(h2, h2_analytic, atol=1e-10, rtol=1e-10)
        assert E == pytest.approx(E_analytic, rel=1e-10, abs=1e-25)
