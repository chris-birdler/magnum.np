import pytest
import torch
from magnumnp import *

@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_call(dtype):
    n  = (1, 2, 2)
    dx = (1, 2, 5)
    mesh = Mesh(n, dx)
    state = State(mesh, dtype=dtype)
    state.j = state.Constant([0,0,1])

    oersted = OerstedField()
    h = oersted.h(state)
