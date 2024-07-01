import pytest
import torch
from magnumnp import *

def test_timeinterpolator():
    n  = (5,5,5)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.m = state.Constant([0,0,1])

    interpolator = TimeInterpolator(state, {
        0.00e-9: [0.0, 0.0, 0.0],
        1.00e-9: [0.0, 0.0, 1.0],
        2.00e-9: [0.0, 0.0, 3.0]})

    state.t = -1.0e-9
    torch.testing.assert_close(interpolator(state), state.Constant([0.0,0.0,-1.0]))
    state.t = 0.0e-9
    torch.testing.assert_close(interpolator(state), state.Constant([0.0,0.0,0.0]))
    state.t = 0.5e-9
    torch.testing.assert_close(interpolator(state), state.Constant([0.0,0.0,0.5]))
    state.t = 1.1e-9
    torch.testing.assert_close(interpolator(state), state.Constant([0.0,0.0,1.2]))
    state.t = 1.5e-9
    torch.testing.assert_close(interpolator(state), state.Constant([0.0,0.0,2.0]))
    state.t = 2.5e-9
    torch.testing.assert_close(interpolator(state), state.Constant([0.0,0.0,4.0]))
    state.t = 0.7e-9
    torch.testing.assert_close(interpolator(state), state.Constant([0.0,0.0,0.7]))

    external = ExternalField(interpolator)
    torch.testing.assert_close(state.avg(external.h(state)), torch.tensor([0.0,0.0,0.7]), atol=1e-6, rtol=1e-6)

def test_expression():
    n  = (1, 1, 10)
    dx = (2e-9, 2e-9, 2e-9)
    mesh = Mesh(n, dx)
    x, y, z = mesh.SpatialCoordinate()

    Ms = Expression(x)
    assert Ms.shape == torch.Size([1,1,10,1])

    Ku_axis = Expression((x,y,z))
    assert Ku_axis.shape == torch.Size([1,1,10,3])

def test_randomize():
    n  = (10,10,10)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.m = state.RandM()

    state.m2 = state.Constant([0.,0.,1.])
    randomize(state.m2)
