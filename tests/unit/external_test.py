import pytest
import pathlib
import torch
from magnumnp import *

def test_different_tensors():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    material = {"Ms": Ms}
    state = State(mesh, material)
    state.m = state.Constant([1,0,0])

    # test vector [3]
    external = ExternalField(state, [1,0,0])
    torch.testing.assert_close(external.h(0,state.m), state.tensor([1,0,0]))
    assert external.E(0, state.m) / (-mesh.volume*constants.mu_0*Ms) == pytest.approx(1.)

    # test vector [nx,ny,nz,3]
    external = ExternalField(state, state.Constant([1,0,0]))
    torch.testing.assert_close(external.h(0,state.m).avg(), state.tensor([1,0,0]))
    assert external.E(0, state.m) / (-mesh.volume*constants.mu_0*Ms) == pytest.approx(1.)

    # test lambda [3]
    h_func = lambda t: state.Tensor([t,0,0])
    external = ExternalField(state, h_func)
    torch.testing.assert_close(external.h(0,state.m).avg(), state.tensor([0,0,0]))
    torch.testing.assert_close(external.h(1.,state.m), state.tensor([1,0,0]))

    # test lambda [nx,ny,nz,3]
    h_func = lambda t: state.Constant([t,0,0])
    external = ExternalField(state, h_func)
    torch.testing.assert_close(external.h(0,state.m).avg(), state.tensor([0,0,0]))
    torch.testing.assert_close(external.h(1.,state.m).avg(), state.tensor([1,0,0]))

def test_setter():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    material = {"Ms": Ms}
    state = State(mesh, material)
    state.m = state.Constant([1,0,0])

    external = ExternalField(state, [1,0,0])
    torch.testing.assert_close(external.h(0,state.m), state.tensor([1,0,0]))

    external.h = [0,1,0]
    torch.testing.assert_close(external.h(0,state.m), state.tensor([0,1,0]))
