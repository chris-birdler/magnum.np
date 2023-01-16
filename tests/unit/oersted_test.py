import pytest
import torch
from magnumnp import *

def test_call():
    n  = (1, 2, 2)
    dx = (1, 2, 5)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.j = state.Constant([0,0,1])

    oersted = OerstedField()
    h = oersted.h(state)
