import pytest
import pathlib
import torch
from magnumnp import *

def test_call():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.m = state.Constant([1,0,0])
    state.material = {"Ms": Ms,
                      "Di": 1.,
                      "Db": 1.,
                      "DD2d": 1.}

    dmi = InterfaceDMIField()
    dmi.h(state)
    dmi = BulkDMIField()
    dmi.h(state)
    dmi = D2dDMIField()
    dmi.h(state)

def test_exchange_only():
    n  = (20, 20, 1)
    dx = (0.5e-9, 0.5e-9, 1e-9)
    mesh = Mesh(n, dx, origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., 0.0))
    state = State(mesh)
    state.material["Ms"] = state.Constant([1./constants.mu_0])
    state.material["A"] = state.Constant([1e-12])
    state.material["Di"] = state.Constant([0.])

    x, y, z = state.SpatialCoordinates()
    domain = x**2 + y**2 < 5e-9**2

    state.m = torch.stack([-y,x,0*z], dim=-1)
    state.m[~domain] = 0

    exchange1 = ExchangeField()
    exchange2 = ExchangeDMIField()

    h1 = exchange1.h(state)
    h2 = exchange2.h(state)

    torch.testing.assert_close(h1, h2)
