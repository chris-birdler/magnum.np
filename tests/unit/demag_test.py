import pytest
import pathlib
import torch
from magnumnp import *


@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_energy_cube(dtype):
    n  = (20, 20, 20)
    dx = (1e-9, 1e-9, 1e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh, dtype=dtype)
    state.material = {"Ms": Ms}
    demag = DemagField()

    state.m = state.Constant([1,0,0])

    h = demag.h(state)
    torch.testing.assert_close(constants.mu_0*h[10,10,10], state.Tensor([-1./3., 0.0, 0.0]), atol=1e-3, rtol=1e-3)

    E = demag.E(state)
    assert E.cpu() == pytest.approx(1./6.*mesh.volume*constants.mu_0*Ms**2, abs=0, rel=1e-3)

def test_Ms_domain():
    n  = (4, 4, 4)
    dx = (1e-9, 1e-9, 1e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant([Ms])}
    state.m = state.Constant([1,0,0])
    demag = DemagField()
    h = demag.h(state)
    E = demag.E(state)
    assert E.cpu() == pytest.approx(1./6.*mesh.volume*constants.mu_0*Ms**2, abs=0, rel=1e-3)

def test_PBC():
    n  = (4, 4, 4)
    dx = (1e-9, 1e-9, 1e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant([Ms])}
    state.m = state.Constant([1,0,0])
    demag = DemagFieldPBC()

    h = demag.h(state)
    torch.testing.assert_close(h, state.Constant([0.0, 0.0, 0.0]), atol=1e-15, rtol=0)
    E = demag.E(state)
    assert E.cpu() == pytest.approx(0., abs=1e-15, rel=0.)

def test_nonequi_vs_equi():
    n  = (4, 4, 4)
    dx = (1., 2., 5.)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant([1.])}

    x, y, z = state.SpatialCoordinate()
    state.m = torch.stack([x*y, y*z, z*x], dim=-1)

    demag = DemagField()
    h1 = demag.h(state)

    demag2 = DemagFieldNonEquidistant()
    h2 = demag2.h(state)

    torch.testing.assert_close(h1, h2, atol=1e-10, rtol=1e-10)

def test_nonequidistant():
    n  = (20, 4, 6)
    dx1 = (5., 3., 1.)
    mesh1 = Mesh(n, dx1)
    state1 = State(mesh1)
    state1.material = {'Ms': state1.Constant([1.])}

    state1.m = state1.Constant([0.,0.,0.])
    state1.m[0,0,0:2,0] = 1.

    demag1 = DemagField()
    h1 = demag1.h(state1)

    n  = (20, 4, 4)
    dx2 = torch.ones(n[2]) * 2.
    dx2[2:] = 1.
    dx = (5., 3., dx2)
    mesh2 = Mesh(n, dx)
    state2 = State(mesh2)
    state2.material = {'Ms': state2.Constant([1.])}

    state2.m = state2.Constant([0.,0.,0.])
    state2.m[0,0,0,0] = 1.

    demag2 = DemagFieldNonEquidistant()
    h2 = demag2.h(state2)

    torch.testing.assert_close(h1[:,:,-2:,:], h2[:,:,-2:,:], atol=1e-10, rtol=1e-10)
