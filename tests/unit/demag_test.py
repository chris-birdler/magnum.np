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
    n  = (4,4,4)
    dx = (1e-9, 1e-9, 1e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant([Ms])}
    state.m = state.Constant([1,0,0])
    demag = DemagField()
    h = demag.h(state)
    E = demag.E(state)
    print("E1:", E)
    print("E2:", 1./6.*mesh.volume*constants.mu_0*Ms**2)
    assert E.cpu() == pytest.approx(1./6.*mesh.volume*constants.mu_0*Ms**2, abs=0, rel=1e-3)

def test_PBC():
    n  = (4,4,4)
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
