import pytest
import torch
from magnumnp import *
from helpers import *

def test_call():
    n  = (2, 3, 4)
    dx = (1, 2, 5)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"A": state.Constant(1.3e-11),
                      "Ms": state.Constant(800e3)}
    state.m = state.Constant([1,0,0])
    exchange = ExchangeField()
    exchange.h(state)
    exchange.E(state)

def test_PBC():
    n  = (100, 25, 1)
    dx = (5e-9, 5e-9, 3e-9)
    mesh = Mesh(n, dx, pbc="xyz")
    state = State(mesh)
    state.material = {"A": state.Constant(1.3e-11),
                      "Ms": state.Constant(800e3)}
    exchange = ExchangeField()
    state.m = state.Constant([0,0,1])
    state.m[:50,:,:,2] = -1
    h = exchange.h(state)
    write_vti({"m":state.m, "h":h}, "data/h_exchangePBC.vti")
    torch.testing.assert_close(h[:50,...], -h[50:,...], atol=1e-15, rtol=1e-15)

def test_nonequi_vs_equi():
    n  = (10, 10, 10)
    dx1 = (1e-9, 2e-9, 5e-9)
    mesh1 = Mesh(n, dx1)
    state1 = State(mesh1)
    state1.material = {"A": state1.Constant(1.), #1.3e-11,
                      "Ms": state1.Constant(1.)} #800e3}

    dx2 = (torch.ones(n[0]) * 1e-9, 2e-9, 5e-9)
    mesh2 = Mesh(n, dx2)
    state2 = State(mesh2)
    state2.material = {"A": state2.Constant(1.), #1.3e-11,
                      "Ms": state2.Constant(1.)} #800e3}
    x, y, z = state1.SpatialCoordinate()
    state1.m = Expression([x*y, y*z, z*x])
    state2.m = Expression([x*y, y*z, z*x])

    exchange = ExchangeField()
    h1 = exchange.h(state1).cpu()
    h2 = exchange.h(state2).cpu()
    torch.testing.assert_close(h1/h1.max(), h2/h1.max(), atol=1e-10, rtol=0)

def test_nonequidistant():
    n  = (9, 2, 3)
    dx0 = torch.ones(n[0]) * 1e-9
    dx0[4:] = 2e-9

    dx = (dx0, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"A": state.Constant(1.3e-11),
                      "Ms": state.Constant(800e3)}
    state.m = state.Constant([1,0,0])
    state.m[3:6,:,:,2] = -1

    exchange = ExchangeField()
    h = exchange.h(state)

    write_vtr({"h": h, "m": state.m}, "data/h.vtr", state)
    torch.testing.assert_close(h[2,...], 4*h[6,...], atol=1e-15, rtol=1e-15)
