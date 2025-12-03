import pytest
import pathlib
from magnumnp import *
import os

#added comment
def test_read_vti():
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "ref" / "fields.vti"

    mesh, fields = read_vti(filename)
    state = State(mesh)
    assert len(fields) == 3
    assert fields['m'].shape == (100, 25, 1, 3)
    assert fields['h'].shape == (100, 25, 1, 3)
    assert fields['domain1'].shape == (100, 25, 1)
    torch.testing.assert_close(state.avg(fields['m']), torch.tensor([0.98,0.02,0.]))
    assert fields['m'].dtype == torch.float64

    torch.set_default_dtype(torch.float32)
    mesh, fields = read_vti(filename)
    assert fields['m'].dtype == torch.float32
    torch.set_default_dtype(torch.float64)

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

def test_read_image():
    """
    Reads in a greyscale image ("init_mag_greyscale.png", 8bit). Calls the read_image function,
    which should be tested. The values between 0 and 256 are changed to only -1 and 1 values. The
    array is taken as initial mz magnetization.
    """
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "ref" / "stripe_domains.png"
    n  = (256, 256, 1)
    dx = (16e-9, 16e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    field = read_image(mesh, filename, fix_aspect_ratio = True)

    state.m = state.Constant([0,0,0])
    state.m[:,:,:,2] = Expression((field > 100)*2. - 1.)
    torch.testing.assert_close(state.avg(state.m), torch.tensor([0.,0.,0.42355347]), rtol=1e-7, atol=1e-7)

def test_read_neper():
    """
    reads Neper .tesr file create using:

    .. code::
        neper -T -n 1000 -domain "cube(500,250,40)" -group "id<500?1:2" -o  test -format tess,tesr -tesrsize 500:125:10
    """
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "ref" / "test.tesr"

    mesh, domains, groups = read_neper(filename, scale = 1e-9)

    assert mesh.n == (500, 125, 10)
    assert mesh.dx == (1e-9, 2e-9, 4e-9)
    assert domains.max().cpu().numpy() == 1000
    assert groups.max().cpu().numpy() == 2

    state = State(mesh, scale = 1e9)
    state.write_vtk(domains, "domains.vti")
    state.material['Ms'] = state.Constant(8e5)
    state.material['Ms'][groups==2] = 0.
    state.write_vtk(state.material, "material.vti")


