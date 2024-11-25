import pytest
import torch
from magnumnp import *

def test_voronoi():
    n  = (100, 100, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2.,-n[1]*dx[1]/2.,-n[2]*dx[2]/2.))
    state = State(mesh)
    voi = Voronoi(mesh, 10)
    domains = voi.relax()
    voi.add_intergrain_phase(2)
    assert voi.domains.unique().size() == torch.Size([11])

def test_voronoi2D():
    n  = (100, 100, 100)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2.,-n[1]*dx[1]/2.,-n[2]*dx[2]/2.))
    state = State(mesh)

    num_grains = 10
    L = torch.tensor(mesh.dx)*torch.tensor(mesh.n)
    offset = torch.tensor(mesh.origin)
    seed_points = L * torch.rand((num_grains, 3)) + offset
    seed_points[:,2] = 0

    voi = Voronoi(mesh, seed_points = seed_points)
    domains = voi.relax()
    assert voi.domains.unique().size() == torch.Size([10])

def test_MFM():
    n = (100, 20, 10)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.))
    state = State(mesh)
    state.material = {"Ms": 8e5}
    x, y, z = mesh.SpatialCoordinate()

    state.m = torch.stack([y, -x, 0*z], dim=-1)
    normalize(state.m)

    demag = DemagField()
    mfm = MFM(height=10e-9, mm_tip = 10e-9, dm_tip=20e-9)

    phi = mfm.PhaseShift(state)
    assert phi.shape == torch.Size([100,20,12,1])
    assert mfm.extended_state.mesh.n == (100,20,12)

def test_LTEM():
    n  = (100, 100, 3)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2.,-n[1]*dx[1]/2.,-n[2]*dx[2]/2.))
    state = State(mesh)
    state.material["Ms"] = 1.
    x, y, z = mesh.SpatialCoordinate()
    state.m = Expression([-y, x, 1e-3+0*z])
    normalize(state.m)

    ltem = LTEM(state)
    ltem.InductionMap()
    ltem.Defocus(1., 1.)
    ltem.MagneticPhaseShift()

def test_nsk():
    n  = (100, 100, 3)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2.,-n[1]*dx[1]/2.,-n[2]*dx[2]/2.))
    state = State(mesh)
    state.material["Ms"] = 1.
    x, y, z = mesh.SpatialCoordinate()
    state.m = Expression([-y, x, 1e-3+0*z])
    normalize(state.m)

    x = nsk(state)
    torch.testing.assert_close(x, torch.tensor(7.9577e-10), atol=0, rtol=1e-4)
