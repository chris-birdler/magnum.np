import pytest
import torch
from magnumnp import *
from helpers import *

def test_float():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Constant(1./constants.mu_0),
                      "A": state.Constant(1e-11)}
    Ms = state.material["Ms"]

    state.m = state.Constant([1,0,0])
    demag = DemagField()
    exchange = ExchangeField()
    assert Ms.cpu() == pytest.approx(1./constants.mu_0)
    assert isinstance(state.t, float)
    assert isinstance(Ms, torch.Tensor)
    assert avg(Ms).cpu() == pytest.approx(1./constants.mu_0)
    demag.h(state)
    exchange.h(state)

def test_setter(simple_state):
    ## set float -> deprecated
    #simple_state.material = {"Ms": simple_state.Constant(1./constants.mu_0), "A": 1e-11}
    #Ms = simple_state.material["Ms"]
    #assert Ms.shape == (100, 25, 1, 1)
    #assert avg(Ms).cpu() == pytest.approx(1./constants.mu_0)

    ## set tensor -> deprecated
    #simple_state.material = {"Ms": simple_torch.tensor([1./constants.mu_0]), "A": simple_torch.tensor([1e-11])}
    #Ms = simple_state.material["Ms"]
    #assert Ms.shape == (100, 25, 1, 1)
    #assert avg(Ms).cpu() == pytest.approx(1./constants.mu_0)

    ## set tensor (no brackets) -> deprecated
    #simple_state.material = {"Ms": simple_torch.tensor(1./constants.mu_0), "A": simple_torch.tensor(1e-11)}
    #Ms = simple_state.material["Ms"]
    #assert Ms.shape == (100, 25, 1, 1)
    #assert avg(Ms).cpu() == pytest.approx(1./constants.mu_0)

    # set tensor-field
    simple_state.material = {"Ms": simple_state.Constant(1./constants.mu_0),
                             "A": simple_state.Constant(1e-11)}
    Ms = simple_state.material["Ms"]
    assert Ms.shape == (100, 25, 1, 1)
    assert avg(Ms).cpu() == pytest.approx(1./constants.mu_0)

    # set tensor-field (no brackets)
    simple_state.material = {"Ms": simple_state.Constant(1./constants.mu_0),
                             "A": simple_state.Constant(1e-11)}
    Ms = simple_state.material["Ms"]
    assert Ms.shape == (100, 25, 1, 1)
    assert avg(Ms).cpu() == pytest.approx(1./constants.mu_0)

    ## set lambda (float) -> deprecated
    #simple_state.material["Ms"] = lambda t: 1./constants.mu_0 * (t+1.)
    #simple_state.t = 0.
    #Ms = simple_state.material["Ms"]
    #assert Ms.shape == (100, 25, 1, 1)
    #assert avg(Ms).cpu() == pytest.approx(1./constants.mu_0)

    #simple_state.t = 1.
    #Ms = simple_state.material["Ms"]
    #assert Ms.shape == (100, 25, 1, 1)
    #assert avg(Ms).cpu() == pytest.approx(2./constants.mu_0)

    # set lambda (tensorfield)
    simple_state.material["Ms"] = lambda t: simple_state.Constant(1./constants.mu_0 * (t+1.))
    simple_state.t = 0.
    Ms = simple_state.material["Ms"]
    assert Ms.shape == (100, 25, 1, 1)
    assert avg(Ms).cpu() == pytest.approx(1./constants.mu_0)

    simple_state.t = 1.
    Ms = simple_state.material["Ms"]
    assert Ms.shape == (100, 25, 1, 1)
    assert avg(Ms).cpu() == pytest.approx(2./constants.mu_0)


def test_material_as_dict():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"A": state.Constant(1.), "B": state.Constant(2.)}

    assert len(state.material.items()) == 2


def test_domain():
    n  = (8,1,1)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    x,y,z = mesh.SpatialCoordinate()
    domain1 = x < 2e-9
    domain2 = x > 6e-9
    state.material["Ms"] = state.Constant(1.)
    state.material["Ms"][domain1] = 0.
    assert avg(state.material["Ms"]).cpu() == pytest.approx(0.75)
    assert avg(state.material["Ms"][~domain1]).cpu() == pytest.approx(1.0)

    state.material["Ku_axis"] = state.Constant([0.,0.,1.])
    state.material["Ku_axis"][domain1] = torch.tensor([1.,0.,0.])

    torch.testing.assert_close(avg(state.material["Ku_axis"]), torch.tensor([0.25,0.00,0.75]))
    torch.testing.assert_close(avg(state.material["Ku_axis"][~domain1]), torch.tensor([0.,0.,1.]))


def test_spatial():
    n  = (8,1,1)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    x,y,z = mesh.SpatialCoordinate()
    state.material["Ms"] = 5.3e5 * x.unsqueeze(-1)
    assert avg(state.material["Ms"]).cpu() == pytest.approx(0.00212)

def test_set():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    x,y,z = mesh.SpatialCoordinate()
    domain1 = x < 2e-9
    domain2 = x > 6e-9
    state.material.set({"Ms": 1.}, domain1)
    assert avg(state.material["Ms"]).cpu() == pytest.approx(0.25)

    state.material.set({"Ms": 2.})
    assert avg(state.material["Ms"]).cpu() == pytest.approx(2.)

    state.material.set({"Ms": 1.}, domain2)
    assert avg(state.material["Ms"]).cpu() == pytest.approx(1.75)
