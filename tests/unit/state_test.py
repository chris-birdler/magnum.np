import pytest
import torch
from magnumnp import *
from helpers import *

def test_domain():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Constant([Ms])}
    state.m = state.Constant([1,1,1])
    torch.testing.assert_close(avg(state.m), state.Tensor([1,1,1]))

    normalize(state.m)
    m_avg = avg(state.m)
    torch.testing.assert_close(m_avg[0]**2+m_avg[1]**2+m_avg[2]**2, state.Tensor(1.))

    x, y, z = state.SpatialCoordinate()
    domain1 = x < 4e-9
    state.m[domain1] = state.Tensor([0,0,1])
    torch.testing.assert_close(avg(state.m[domain1]), state.Tensor([0,0,1]))
    torch.testing.assert_close(avg(state.m), state.Tensor([0.28867513, 0.28867513, 0.78867513]))


def test_normal(simple_state):
    x = simple_state._normal(0.,1., size = simple_state.m.shape)
    assert isinstance(x, torch.Tensor)
    torch.testing.assert_close(x.shape, simple_state.m.shape)


def test_spatial_coordinate():
    n = (100, 20, 10)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.))
    state = State(mesh)
    x, y, z = state.SpatialCoordinate()
    torch.testing.assert_close(x.min(), state.Tensor(-n[0]*dx[0]/2.+dx[0]/2.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(x.max(), state.Tensor( n[0]*dx[0]/2.-dx[0]/2.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(y.min(), state.Tensor(-n[1]*dx[1]/2.+dx[1]/2.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(y.max(), state.Tensor( n[1]*dx[1]/2.-dx[1]/2.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(z.min(), state.Tensor(-n[2]*dx[2]/2.+dx[2]/2.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(z.max(), state.Tensor( n[2]*dx[2]/2.-dx[2]/2.), atol=1e-15, rtol=1e-15)


def test_spatial_coordinate_nonequi():
    n = (1, 1, 4)
    dx = (5e-9, 5e-9, (torch.arange(n[2])+1.)*1e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., 0) )
    state = State(mesh)
    x, y, z = state.SpatialCoordinate()

    torch.testing.assert_close(z.reshape(-1), state.Tensor([0.5e-9, 2.0e-9, 4.5e-9, 8.0e-9]), atol=1e-15, rtol=1e-15)


def test_average():
    n = (2, 3, 5)
    dx = (5e-9, 5e-9, 5e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.) )
    state = State(mesh)
    state.m = state.Constant([1.,0.,0.])
    state.material["A"] = state.Constant([1.])
    x, y, z = state.SpatialCoordinate()

    torch.testing.assert_close(avg(state.m), state.Tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(avg(state.material["A"]), state.Tensor([1.]), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(avg(z), state.Tensor(0.), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(avg(state.m[:,:,:2,:]), state.Tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)


def test_average_nonequi():
    n = (1, 1, 4)
    dx = (5e-9, 5e-9, (torch.arange(n[2])+1.)*1e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., 0))
    state = State(mesh)
    state.m = state.Constant([1.,0.,0.])
    state.material["A"] = state.Constant([1.])
    x, y, z = state.SpatialCoordinate()

    torch.testing.assert_close(avg(state.m, state.cell_volumes), state.Tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(avg(state.material["A"], state.cell_volumes), state.Tensor([1.]), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(avg(z, state.cell_volumes), state.Tensor(5e-9), atol=1e-15, rtol=1e-15)
    torch.testing.assert_close(avg(state.m[:,:,:2,:], state.cell_volumes[:,:,:2]), state.Tensor([1.,0.,0.]), atol=1e-15, rtol=1e-15)
