import pytest
import pathlib
import torch
from magnumnp import *

def test_different_tensors():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh_volume = n[0] * n[1] * n[2] * dx[0] * dx[1] * dx[2]
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Constant([Ms])}
    state.m = state.Constant([1.,0.,0.])

    ## test vector [3] -> deprecated
    #external = ExternalField([1.,0.,0.])
    #torch.testing.assert_close(avg(external.h(state)), torch.tensor([1.,0.,0.]))
    #assert external.E(state).cpu() / (-mesh_volume*constants.mu_0*Ms) == pytest.approx(1., abs=1e-6, rel=1e-6)

    # test vector [nx,ny,nz,3]
    external = ExternalField(state.Constant([1.,0.,0.]))
    torch.testing.assert_close(avg(external.h(state)), torch.tensor([1.,0.,0.]))
    assert external.E(state).cpu() / (-mesh_volume*constants.mu_0*Ms) == pytest.approx(1., abs=1e-6, rel=1e-6)

    ## test lambda [3] -> deprecated
    #h_func = lambda t: torch.tensor([t,0.,0.])
    #external = ExternalField(h_func)
    #state.t = 0.
    #torch.testing.assert_close(avg(external.h(state)), torch.tensor([0.,0.,0.]), atol=1e-15, rtol=1e-15)
    #state.t = 1.
    #torch.testing.assert_close(avg(external.h(state)), torch.tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)

    # test lambda [nx,ny,nz,3]
    h_func = lambda t: state.Constant([t,0.,0.])
    external = ExternalField(h_func)
    state.t = 0.
    torch.testing.assert_close(avg(external.h(state)), torch.tensor([0.,0.,0.]), atol=1e-15, rtol=1e-15)
    state.t = 1.
    torch.testing.assert_close(avg(external.h(state)), torch.tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)

def test_setter():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Constant([Ms])}
    state.m = state.Constant([1.,0.,0.])

    external = ExternalField(state.Constant([1.,0.,0.]))
    torch.testing.assert_close(avg(external.h(state)), torch.tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)

    external.h = state.Constant([0.,1.,0.])
    torch.testing.assert_close(avg(external.h(state)), torch.tensor([0.,1.,0.]), atol=1e-15, rtol=1e-15)
