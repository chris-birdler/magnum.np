import pytest
import pathlib
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

def test_regression():
    n  = (100, 100, 10)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.j = state.Constant([1,0,0])

    oersted = OerstedField()
    h_oersted = oersted.h(state)

    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "ref" / "h_oersted_regression.vti"
    ## Uncomment to updated reference data
    #write_vti({"h_oersted":h_oersted}, filename)
    mesh, ref = read_vti(filename)

    torch.testing.assert_close(h_oersted, ref["h_oersted"], atol=1e-10, rtol=1e-6)
