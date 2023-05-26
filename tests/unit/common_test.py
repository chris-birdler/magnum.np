import pytest
import torch
from magnumnp import *

def test_timeinterpolator():
    n  = (5,5,5)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.m = state.Constant((0,0,1))

    interpolator = TimeInterpolator(state, {
        0.00e-9: [0.0, 0.0, 0.0],
        1.00e-9: [0.0, 0.0, 1.0],
        2.00e-9: [0.0, 0.0, 3.0]})

    state.t = 0.7e-9
    torch.testing.assert_close(interpolator(-1.0e-9), state.Tensor([0.0,0.0,-1.0]))
    torch.testing.assert_close(interpolator(0.0e-9), state.Tensor([0.0,0.0,0.0]))
    torch.testing.assert_close(interpolator(0.5e-9), state.Tensor([0.0,0.0,0.5]))
    torch.testing.assert_close(interpolator(1.1e-9), state.Tensor([0.0,0.0,1.2]))
    torch.testing.assert_close(interpolator(1.5e-9), state.Tensor([0.0,0.0,2.0]))
    torch.testing.assert_close(interpolator(2.5e-9), state.Tensor([0.0,0.0,4.0]))
    torch.testing.assert_close(interpolator(state.t), state.Tensor([0.0,0.0,0.7]))

    external = ExternalField(interpolator)
    torch.testing.assert_close(external.h(state).avg(), state.Tensor([0.0,0.0,0.7]), atol=1e-6, rtol=1e-6)

def test_timeinterpolator_field():
    n  = (1, 1, 10)
    dx = (2e-9, 2e-9, 2e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.m = state.Constant([0,0,1])

    zeeman = ExternalField(TimeInterpolator(state, {0e-0: [0, 0, -1],
                                                    1e-9: [0, 0,  1]}))

    zeeman.h(state)
