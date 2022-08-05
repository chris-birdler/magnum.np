import pytest
import pathlib
import torch
from magnumnp import *

def test_decorated_function():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    material = {"Ms": Ms}
    state = State(mesh, material)
    state.m = state.Constant([1,1,1])
    torch.testing.assert_close(state.m.avg(), state.tensor([1,1,1]))

    state.m.normalize()
    avg = state.m.avg()
    assert avg[0]**2+avg[1]**2+avg[2]**2 == pytest.approx(1.)
    
def test_Constant():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    material = {"Ms": Ms}
    state = State(mesh, material)
    state.m = state.Constant([0,0,1])
    
