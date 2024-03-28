import pytest
import pathlib
from magnumnp import *
import os

#added comment
def test_read_vti():
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "ref" / "fields.vti"

    mesh, fields = read_vti(filename)
    assert len(fields) == 3
    assert fields['m'].shape == (100, 25, 1, 3)
    assert fields['h'].shape == (100, 25, 1, 3)
    assert fields['domain1'].shape == (100, 25, 1)
    torch.testing.assert_close(avg(fields['m']), torch.tensor([0.98,0.02,0.]))
    assert fields['m'].dtype == torch.float64

    torch.set_default_dtype(torch.float32)
    mesh, fields = read_vti(filename)
    assert fields['m'].dtype == torch.float32

def test_write_vti(tmp_path):
    os.chdir(tmp_path)

    n  = (2, 2, 2)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.m = state.Constant([1,0,0])

    write_vti({"m": state.m}, "m.vti", state)
    write_vti({"m": state.m}, "data/m.vti", state)

    for f in [tmp_path / "m.vti", tmp_path / "data" / "m.vti"]:
        mesh, fields = read_vti(f)
        assert len(fields) == 1
        assert fields['m'].shape == (2, 2, 2, 3)

