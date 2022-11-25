import pytest
import pathlib
import torch
from magnumnp import *

def test_simple():
    n  = (10, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {
        "Ms": 8e5,
        "A": 1.3e-11,
        "alpha": 0.02
        }

    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField([-24.6e-3/constants.mu_0,
                              +4.3e-3/constants.mu_0,
                              0.0])

    state.m = state.Constant([1,0,0])
    state.m[5:,:,:,0] = -1.0

    llg = LLGSolver([demag, exchange, external])
    llg.step(state, 1e-11)
    assert state.t.cpu() == pytest.approx(1e-11)


def test_material_tensors():
    n  = (10, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)

    state.material["alpha"] = state.Constant([0.02])
    state.material["Ms"] = state.Constant([8e5])
    state.material["A"] = state.Constant([1.3e-11])
    state.material["Ku"] = state.Constant([1e5])
    state.material["Ku_axis"] = state.Constant([0,1,0])

    demag    = DemagField()
    exchange = ExchangeField()
    aniso    = UniaxialAnisotropyField()
    external = ExternalField([-24.6e-3/constants.mu_0,
                              +4.3e-3/constants.mu_0,
                              0.0])

    state.m = state.Constant([1,0,0])
    state.m[5:,:,:,0] = -1.0

    llg = LLGSolver([demag, exchange, aniso, external])
    llg.step(state, 1e-11)
    assert state.t.cpu() == pytest.approx(1e-11)
