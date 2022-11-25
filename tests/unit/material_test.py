import pytest
import torch
from magnumnp import *

def test_float():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": 1./constants.mu_0, "A":1e-11}
    Ms = state.material["Ms"]
    state.m = state.Constant([1,0,0])
    demag = DemagField()
    exchange = ExchangeField()
    assert Ms.cpu() == pytest.approx(1./constants.mu_0)
    assert isinstance(state.t, torch.Tensor)
    assert isinstance(Ms, torch.Tensor)
    assert Ms.avg().cpu() == pytest.approx(1./constants.mu_0)
    demag.h(state)
    exchange.h(state)

def test_tensor():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Tensor([1./constants.mu_0]), "A":state.Tensor([1e-11])}
    Ms = state.material["Ms"]
    assert Ms.cpu() == pytest.approx(1./constants.mu_0)
    assert isinstance(state.t, torch.Tensor)
    assert isinstance(Ms, torch.Tensor)
    assert Ms.avg().cpu() == pytest.approx(1./constants.mu_0)

    state.m = state.Constant([1,0,0])
    demag = DemagField()
    exchange = ExchangeField()
    demag.h(state)
    exchange.h(state)

def test_tensorfield():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Constant([1./constants.mu_0]), "A":state.Constant([1e-11])} # TODO: state.Constant(1./constants.mu_0) without [] brackets does not work!
    Ms = state.material["Ms"]
    assert isinstance(state.t, torch.Tensor)
    assert isinstance(Ms, torch.Tensor)
    assert Ms.avg().cpu() == pytest.approx(1./constants.mu_0)

    state.m = state.Constant([1,0,0])
    demag = DemagField()
    exchange = ExchangeField()
    demag.h(state)
    exchange.h(state)

def test_lambda():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"A":state.Constant([1e-11])}
    state.material["Ms"] = lambda t: 1./constants.mu_0 * (t+1.)
    Ms = state.material["Ms"]
    assert isinstance(state.t, torch.Tensor)
    assert isinstance(Ms, torch.Tensor)
    assert Ms.avg().cpu() == pytest.approx(1./constants.mu_0)

    state.m = state.Constant([1,0,0])
    demag = DemagField()
    h1 = demag.h(state).avg()
    state.t = 1.
    h2 = demag.h(state).avg()
    torch.testing.assert_close(h2, 2*h1)

def test_material_as_dict():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"A": 1, "B": 2}

    assert len(state.material.items()) == 2

def test_slice_on_constant():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Tensor([1./constants.mu_0])}
    Ms = state.material["Ms"]
    x = Ms[1:,:,:]
    assert x.dim() == Ms.dim()
