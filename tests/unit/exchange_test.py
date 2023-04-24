import pytest
import torch
from magnumnp import *
from helpers import *

def test_call():
    n  = (2, 3, 4)
    dx = (1, 2, 5)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"A": 1.3e-11,
                      "Ms": 800e3}
    state.m = state.Constant([1,0,0])
    exchange = ExchangeField()
    exchange.h(state)
    exchange.E(state)

def test_PBC(simple_state):
    exchange = ExchangeFieldPBC()
    simple_state.m = simple_state.Constant([0,0,1])
    simple_state.m[:50,:,:,2] = -1
    h = exchange.h(simple_state)
    torch.testing.assert_close(h[:50,...], -h[50:,...], atol=1e-15, rtol=1e-15)

def test_nonequi_vs_equi():
    n  = (10, 10, 10)
    dx1 = (1e-9, 2e-9, 5e-9)
    mesh1 = Mesh(n, dx1)
    state1 = State(mesh1)
    state1.material = {"A": 1., #1.3e-11,
                      "Ms": 1.} #800e3}

    dx2 = (torch.ones(n[0]) * 1e-9, 2e-9, 5e-9)
    mesh2 = Mesh(n, dx2)
    state2 = State(mesh2)
    state2.material = {"A": 1., #1.3e-11,
                      "Ms": 1.} #800e3}
    x, y, z = state1.SpatialCoordinate()
    state1.m = torch.stack([x*y, y*z, z*x], dim=-1)
    state2.m = torch.stack([x*y, y*z, z*x], dim=-1)
    
    exchange = ExchangeField()
    h1 = exchange.h(state1)
    h2 = exchange.h(state2)
    
    torch.testing.assert_close(h1/h1.max(), h2/h1.max(), atol=1e-15, rtol=1e-15)

def test_nonequidistant():
    n  = (9, 2, 3)
    dx0 = torch.ones(n[0]) * 1e-9
    dx0[4:] = 2e-9
   
    dx = (dx0, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"A": 1.3e-11,
                      "Ms": 800e3}
    state.m = state.Constant([1,0,0])
    state.m[3:6,:,:,2] = -1
    
    exchange = ExchangeField()
    h = exchange.h(state)
    
    write_vtr({"h": h, "m": state.m}, "data/h.vtr", state)
    torch.testing.assert_close(h[2,...], 4*h[6,...], atol=1e-15, rtol=1e-15)
