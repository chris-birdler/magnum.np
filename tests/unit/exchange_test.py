import pytest
import torch
from math import pi, cos, sin
from magnumnp import *

def test_call():
    n  = (1, 1, 1)
    dx = (1, 2, 5)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"A": 1.3e-11,
                      "Ms": 800e3}
    state.m = state.Constant([1,0,0])
    exchange = ExchangeField()
    exchange.h(state)
    exchange.E(state)
