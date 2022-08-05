import pytest
import pathlib
import torch
from magnumnp import *

def test_simple():
    n  = (10, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    material = {
            "Ms": 8e5,
            "A": 1.3e-11,
            "alpha": 0.02
            }
    state = State(mesh, material)

    demag    = DemagField(state)
    exchange = ExchangeField(state)
    external = ExternalField(state, [-24.6e-3/constants.mu_0,
                                     +4.3e-3/constants.mu_0,
                                     0.0])

    state.m = state.Constant([1,0,0])
    state.m[5:,:,:,0] = -1.0

    llg = LLGSolver(state, [demag, exchange, external])
    llg.step(1e-11)
    assert state.t == pytest.approx(1e-11)


def test_material_tensors():
    n  = (10, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    material = {}
    state = State(mesh, material)

    material["alpha"] = state.Constant([0.02])
    material["Ms"] = state.Constant([8e5])
    material["A"] = state.Constant([1.3e-11])
    material["K"] = state.Constant([1e5])
    material["K_axis"] = state.Constant([0,1,0])

    demag    = DemagField(state)
    exchange = ExchangeField(state)
    aniso    = AnisotropyField(state)
    external = ExternalField(state, [-24.6e-3/constants.mu_0,
                                     +4.3e-3/constants.mu_0,
                                     0.0])

    state.m = state.Constant([1,0,0])
    state.m[5:,:,:,0] = -1.0

    llg = LLGSolver(state, [demag, exchange, aniso, external])
    llg.step(1e-11)
    assert state.t == pytest.approx(1e-11)
