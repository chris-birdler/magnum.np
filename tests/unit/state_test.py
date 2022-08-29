import pytest
import pathlib
import torch
from magnumnp import *

def test_decorated_function():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": Ms}
    state.m = state.Constant([1,1,1])
    torch.testing.assert_close(state.m.avg(), state.Tensor([1,1,1]))

    state.m.normalize()
    avg = state.m.avg()
    assert avg[0]**2+avg[1]**2+avg[2]**2 == pytest.approx(1.)

    x, y, z = state.SpatialCoordinates()
    domain1 = x < 4e-9
    state.m[domain1] = state.Tensor([0,0,1])
    torch.testing.assert_close(state.m[domain1].avg(), state.Tensor([0,0,1]))
    torch.testing.assert_close(state.m.avg(), state.Tensor([0.28867513, 0.28867513, 0.78867513]))

def test_Constant():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": Ms}
    state.m = state.Constant([0,0,1])
    torch.testing.assert_close(state.m.avg(), state.Tensor([0,0,1]))
