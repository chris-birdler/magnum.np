import pytest
import pathlib
import torch
from magnumnp import *

def test_different_tensors():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": Ms}
    state.m = state.Constant([1,0,0])

    # test vector [3]
    external = ExternalField([1,0,0])
    torch.testing.assert_close(external.h(state), state.tensor([1,0,0]))
    assert external.E(state) / (-mesh.volume*constants.mu_0*Ms) == pytest.approx(1.)

    # test vector [nx,ny,nz,3]
    external = ExternalField(state.Constant([1,0,0]))
    torch.testing.assert_close(external.h(state).avg(), state.tensor([1,0,0]))
    assert external.E(state) / (-mesh.volume*constants.mu_0*Ms) == pytest.approx(1.)

    # test lambda [3]
    h_func = lambda t: state.Tensor([t,0,0])
    external = ExternalField(h_func)
    state.t = 0.
    torch.testing.assert_close(external.h(state).avg(), state.tensor([0,0,0]))
    state.t = 1.
    torch.testing.assert_close(external.h(state), state.tensor([1,0,0]))

    # test lambda [nx,ny,nz,3]
    h_func = lambda t: state.Constant([t,0,0])
    external = ExternalField(h_func)
    state.t = 0.
    torch.testing.assert_close(external.h(state).avg(), state.tensor([0,0,0]))
    state.t = 1.
    torch.testing.assert_close(external.h(state).avg(), state.tensor([1,0,0]))

def test_setter():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": Ms}
    state.m = state.Constant([1,0,0])

    external = ExternalField([1,0,0])
    torch.testing.assert_close(external.h(state), state.tensor([1,0,0]))

    external.h = [0,1,0]
    torch.testing.assert_close(external.h(state), state.tensor([0,1,0]))
