import pytest
import pathlib
import torch
from magnumnp import *

def test_zhangli():
    n  = (2,1,2)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {
        "Ms": 8e5,
        "A": 1.3e-11,
        "alpha": 0.1,
        "xi": 0.05,
        "b": 72.17e-12
        }
    
    state.m = state.Constant([0,1,0])
    state.m[:1,:,:,1] = -1.
    
    state.j = state.Tensor((1e12, 0, 0))
    torque = SpinTorqueZhangLi()
    h1 = torque.h(state)

    state.j = state.Constant([1e12,0,0])
    h2 = torque.h(state)

    assert torch.allclose(h1, h2)

