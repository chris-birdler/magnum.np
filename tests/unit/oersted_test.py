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
    h_oersted = oersted.h(state).cpu()

    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "ref" / "h_oersted_regression.vti"
    ## Uncomment to updated reference data
    #write_vti({"h_oersted":h_oersted}, filename)
    mesh, ref = read_vti(filename)

    torch.testing.assert_close(h_oersted, ref["h_oersted"], atol=1e-10, rtol=1e-6)

def test_wire():
    N = 101
    n  = (N,N,101)
    dx = (1e-9, 1e-9, 5e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2.,-n[1]*dx[1]/2.,-n[2]*dx[2]/2.))
    state = State(mesh)
    state.j = state.Constant([0,0,0])
    state.j[N//2,N//2,:,2] = 1

    oersted = OerstedField()
    h1 = oersted.h(state)
    h1 = h1[:,n[1]//2,n[2]//2,1]

    x,y,z = state.SpatialCoordinate()
    h2 = dx[0]**2/(2.*torch.pi*x[:,n[1]//2,n[2]//2])
    h2[n[0]//2] = 0

    #import matplotlib.pyplot as plt
    #fig, ax = plt.subplots()
    #ax.plot(x[:,n[1]//2,n[2]//2], h1, '-')
    #ax.plot(x[:,n[1]//2,n[2]//2], h2, '-')
    #ax.set_ylim([-2e-10, 2e-10])
    #ax.grid()
    #fig.savefig("data/results.png")

    torch.testing.assert_close(h1[:n[0]//2-5]/h1.max(), h2[:n[0]//2-5]/h1.max(), atol=1e-3, rtol=1e-3)
