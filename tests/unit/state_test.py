import pytest
import torch
import numpy as np
from magnumnp import *
from helpers import *

def test_domain():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Constant(Ms)}
    state.m = state.Constant([1,1,1])
    torch.testing.assert_close(state.avg(state.m), torch.tensor([1.,1.,1.]))

    normalize(state.m)
    m_avg = state.avg(state.m)
    torch.testing.assert_close(m_avg[0]**2+m_avg[1]**2+m_avg[2]**2, torch.tensor(1.))

    x, y, z = mesh.SpatialCoordinate()
    domain1 = x < 4e-9
    state.m[domain1] = torch.tensor([0.,0.,1.])
    torch.testing.assert_close(state.avg(state.m[domain1]), torch.tensor([0.,0.,1.]))
    torch.testing.assert_close(state.avg(state.m), torch.tensor([0.28867513, 0.28867513, 0.78867513]))


def test_normal(simple_state):
    x = torch.normal(0.,1., size = simple_state.m.shape)
    assert isinstance(x, torch.Tensor)
    torch.testing.assert_close(x.shape, simple_state.m.shape)


def test_spatial_coordinate():
    n = (100, 20, 10)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.))
    state = State(mesh)
    x, y, z = mesh.SpatialCoordinate()
    torch.testing.assert_close(x.min(), torch.tensor(-n[0]*dx[0]/2.+dx[0]/2.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(x.max(), torch.tensor( n[0]*dx[0]/2.-dx[0]/2.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(y.min(), torch.tensor(-n[1]*dx[1]/2.+dx[1]/2.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(y.max(), torch.tensor( n[1]*dx[1]/2.-dx[1]/2.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(z.min(), torch.tensor(-n[2]*dx[2]/2.+dx[2]/2.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(z.max(), torch.tensor( n[2]*dx[2]/2.-dx[2]/2.), atol=1e-15, rtol=1e-15)


def test_spatial_coordinate_nonequi():
    n = (1, 1, 4)
    dx = (5e-9, 5e-9, (np.arange(n[2])+1.)*1e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., 0) )
    state = State(mesh)
    x, y, z = mesh.SpatialCoordinate()

    torch.testing.assert_close(z.reshape(-1), torch.tensor([0.5e-9, 2.0e-9, 4.5e-9, 8.0e-9]), atol=1e-15, rtol=1e-15)


def test_average():
    n = (2, 3, 5)
    dx = (5e-9, 5e-9, 5e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.) )
    state = State(mesh)
    state.m = state.Constant([1.,0.,0.])
    state.material["A"] = state.Constant(1.)
    x, y, z = mesh.SpatialCoordinate()

    torch.testing.assert_close(state.avg(state.m), torch.tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(state.avg(state.material["A"]), torch.tensor([1.]), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(state.avg(z), torch.tensor(0.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(state.avg(state.m[:,:,:2,:]), torch.tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)


def test_average_nonequi():
    n = (1, 1, 4)
    dx = (5e-9, 5e-9, (np.arange(n[2])+1.)*1e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., 0))
    state = State(mesh)
    state.m = state.Constant([1.,0.,0.])
    state.material["A"] = state.Constant(1.)
    x, y, z = mesh.SpatialCoordinate()

    torch.testing.assert_close(state.avg(state.m, state.mesh.cell_volumes), torch.tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(state.avg(state.material["A"], state.mesh.cell_volumes), torch.tensor([1.]), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(state.avg(z, state.mesh.cell_volumes), torch.tensor(5e-9), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(state.avg(state.m[:,:,:2,:], state.mesh.cell_volumes[:,:,:2]), torch.tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)

def test_time():
    n = (2, 3, 5)
    dx = (5e-9, 5e-9, 5e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.) )
    state = State(mesh)
    assert state.t.dtype == torch.float64

    state.t = 1
    assert state.t.dtype == torch.float64

    state.t = 2e-9
    assert state.t.dtype == torch.float64
