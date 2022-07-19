import pytest
import pathlib
import torch
from magnumnp import *

def test_zhangli():
    n  = (2,2,2)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    material = {
            "Ms": 8e5,
            "A": 1.3e-11,
            "alpha": 0.1,
            "xi": 0.05,
            "b": 72.17e-12
            }
    
    state = State(mesh, material)
    state.m = state.zeros(n + (3,))
    state.m[:1,:,:,1] = -1.
    state.m[1:,:,:,1] = 1.
    
    state.j = state.tensor((1e12, 0, 0))
    torque = SpinTorqueZhangLi(state)
    h1 = torque.h(0, state.m)

    state.j = state.zeros(n + (3,))
    state.j[:,:,:,0] = 1e12
    torque = SpinTorqueZhangLi(state)
    h2 = torque.h(0, state.m)

    assert torch.allclose(h1, h2)

