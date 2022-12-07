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
    assert torch.allclose(h[:50,...], -h[50:,...])
