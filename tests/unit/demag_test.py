import pytest
import pathlib
import torch
from magnumnp import *
from magnumnp.field_terms.demag import newell, f, g
import numpy as np


@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_energy_cube(dtype):
    torch.set_default_dtype(dtype)
    n  = (20, 20, 20)
    dx = (1e-9, 1e-9, 1e-9)
    mesh_volume = n[0] * n[1] * n[2] * dx[0] * dx[1] * dx[2]
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Constant(Ms)}
    demag = DemagField()

    state.m = state.Constant([1,0,0])

    h = demag.h(state)
    torch.testing.assert_close(constants.mu_0*h[10,10,10], torch.tensor([-1./3., 0.0, 0.0]), atol=1e-3, rtol=1e-3)

    E = demag.E(state)
    assert E.cpu() == pytest.approx(1./6.*mesh_volume*constants.mu_0*Ms**2, abs=0, rel=1e-3)

def test_single():
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    demag = DemagField()

    state.m = state.Constant([1,0,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([-1./3.,0.,0.]), atol=1e-10, rtol=1e-10)

def test_truePBC():
    n  = (4, 4, 4)
    dx = (1e-9, 1e-9, 1e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant(Ms)}
    state.m = state.Constant([1,0,0])
    demag = DemagFieldPBC()

    h = demag.h(state)
    torch.testing.assert_close(h, state.Constant([0.0, 0.0, 0.0]), atol=1e-15, rtol=0)
    E = demag.E(state)
    assert E.cpu() == pytest.approx(0., abs=1e-15, rel=0.)

def test_pseudoPBC():
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)

    # cube
    mesh = Mesh(n, dx, pbc = (0,0,0))
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    demag = DemagField()
    state.m = state.Constant([1,0,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([-1./3.,0.,0.]), atol=1e-10, rtol=1e-10)
    state.m = state.Constant([0,1,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,-1./3.,0.]), atol=1e-10, rtol=1e-10)
    state.m = state.Constant([0,0,1])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,-1./3.]), atol=1e-10, rtol=1e-10)

    # long cylinder
    mesh = Mesh(n, dx, pbc = (10,0,0))
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    demag = DemagField()
    state.m = state.Constant([1,0,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,0.]), atol=1e-2, rtol=1e-2)
    state.m = state.Constant([0,1,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,-1./2.,0.]), atol=1e-2, rtol=1e-2)
    state.m = state.Constant([0,0,1])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,-1./2.]), atol=1e-2, rtol=1e-2)

    # thin film
    mesh = Mesh(n, dx, pbc = (10,10,0))
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    demag = DemagField()
    state.m = state.Constant([1,0,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,0.]), atol=5e-2, rtol=5e-2)
    state.m = state.Constant([0,1,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,0.]), atol=5e-2, rtol=5e-2)
    state.m = state.Constant([0,0,1])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,-1.]), atol=5e-2, rtol=5e-2)

@pytest.mark.parametrize("nx", [5,6])
def test_pseudoPBC_single(nx):
    n  = (nx, 4, 3)
    dx = (1e-9, 2e-9, 5e-9)

    mesh = Mesh(n, dx, pbc = (1,0,0))
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    demag = DemagField()
    state.m = state.Constant([0,0,0])
    state.m[0,0,0,0] = 1.

    h = demag.h(state)

    h2xa = newell(f, torch.tensor(2*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2xb = newell(f, torch.tensor((2+n[0])*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2xc = newell(f, torch.tensor((2-n[0])*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2x = h2xa + h2xb + h2xc
    torch.testing.assert_close(h[2,2,0,0]/h.max(), h2x/h.max(), atol=1e-10, rtol=1e-10)

    h2ya = newell(g, torch.tensor(2*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2yb = newell(g, torch.tensor((2+n[0])*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2yc = newell(g, torch.tensor((2-n[0])*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2y = h2ya + h2yb + h2yc

    torch.testing.assert_close(h[2,2,0,1]/h.max(), h2y/h.max(), atol=1e-10, rtol=1e-10)

def test_pseudoPBC_vs_repeat():
    n  = (108, 1, 50)
    dx = (1e-9, 16e-9, 1e-9)

    # pbc
    mesh = Mesh(n, dx, pbc = (1,0,0))
    state = State(mesh)
    state.material = {'Ms': state.Constant(2.1/constants.mu_0)}

    demag = DemagField()
    state.m = state.Constant([0.0, 1.0, 0.0])
    state.m[:20,:,:,:] = 0
    state.m[-20:,:,:,:] = 0
    h = demag.h(state)

    # nopbc (repeat)
    mesh2 = Mesh((3*n[0], n[1], n[2]), dx, origin=(-n[0]*dx[0],0.,0.))
    state2 = State(mesh2)
    state2.material = {'Ms': state2.Constant(2.1/constants.mu_0)}
    demag2 = DemagField()
    state2.m = torch.vstack([state.m,state.m,state.m])
    write_vti(state.m, "data/m.vti", state)
    write_vti(state2.m, "data/m2.vti", state2)
    h2 = demag2.h(state2)

    torch.testing.assert_close(h/h.abs().max(), h2[n[0]:2*n[0],:,:,:]/h.abs().max(), atol=1e-3, rtol=1e-3)

def test_nonequi_vs_equi():
    n  = (4, 4, 4)
    dx = (1., 2., 5.)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}

    x, y, z = mesh.SpatialCoordinate()
    state.m = Expression([x*y, y*z, z*x])

    demag = DemagField()
    h1 = demag.h(state)

    demag2 = DemagFieldNonEquidistant()
    h2 = demag2.h(state)

    torch.testing.assert_close(h1, h2, atol=1e-10, rtol=1e-10)

def test_nonequi_vs_equi_far_layers():
    # layer pairs beyond the Newell radius (p = 20 cell diagonals): the
    # far-field on-axis coupling must be kept. Previous versions zeroed the
    # local [0,0,0] entry of the dipole formula positionally, which dropped
    # the dominant interlayer coupling for distant layers; the field of a
    # single magnetized layer then vanished in layers further away than the
    # Newell radius.
    n  = (4, 4, 50)
    dx = (1., 1., 1.)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    state.m = state.Constant([0., 0., 0.])
    state.m[:,:,0,2] = 1. # magnetize only the bottom layer

    h1 = DemagField().h(state)
    h2 = DemagFieldNonEquidistant().h(state)

    # the far layers carry a small but essential dipole field
    assert h1[:,:,-1,:].abs().max() > 0.
    torch.testing.assert_close(h1, h2, atol=1e-10, rtol=1e-6)

def test_nonequidistant():
    n  = (20, 4, 6)
    dx1 = (5., 3., 1.)
    mesh1 = Mesh(n, dx1)
    state1 = State(mesh1)
    state1.material = {'Ms': state1.Constant(1.)}

    state1.m = state1.Constant([0.,0.,0.])
    state1.m[0,0,0:2,0] = 1.

    demag1 = DemagField()
    h1 = demag1.h(state1)

    n  = (20, 4, 4)
    dx2 = np.ones(n[2]) * 2.
    dx2[2:] = 1.
    dx = (5., 3., dx2)
    mesh2 = Mesh(n, dx)
    state2 = State(mesh2)
    state2.material = {'Ms': state2.Constant(1.)}

    state2.m = state2.Constant([0.,0.,0.])
    state2.m[0,0,0,0] = 1.

    demag2 = DemagFieldNonEquidistant()
    h2 = demag2.h(state2)

    torch.testing.assert_close(h1[:,:,-2:,:], h2[:,:,-2:,:], atol=1e-10, rtol=1e-10)

def test_regression():
    n  = (101, 101, 101)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Constant(1.)}
    state.m = state.Constant([0,0,0])
    state.m[50,50,50] = 1.

    demag = DemagField()
    h_demag = demag.h(state)

    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "ref" / "h_demag_regression.vti"
    ## Uncomment to updated reference data
    #write_vti({"h_demag":h_demag}, filename)
    mesh, ref = read_vti(filename)

    torch.testing.assert_close(h_demag, ref["h_demag"], atol=1e-10, rtol=1e-6)

def test_precision():
    ''' if Newell formula is evaluated with single precision large fluctuations are expected '''
    n = (100,1,1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    dtype = torch.get_default_dtype()
    torch.set_default_dtype(torch.float32)
    state = State(mesh)
    state.material = {"Ms": state.Constant(1.)}
    state.m = state.Constant([0,0,0])
    state.m[0,0,0,2] = 1.

    demag = DemagField()
    h_demag = demag.h(state).abs().log()[:,0,0,2]

    ref = torch.arange(n[0])
    ref = torch.log(1/4./torch.pi/ref**3)

    torch.set_default_dtype(dtype)
    torch.testing.assert_close(h_demag[20:], ref[20:], atol=1e-2, rtol=1e-2)
