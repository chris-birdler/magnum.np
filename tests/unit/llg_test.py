import pytest
import torch
from magnumnp import *
from helpers import *

def test_simple(simple_state):
    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField([-24.6e-3/constants.mu_0,
                              +4.3e-3/constants.mu_0,
                              0.0])

    simple_state.m = simple_state.Constant([1,0,0])
    simple_state.m[5:,:,:,0] = -1.0

    llg = LLGSolver([demag, exchange, external])
    llg.step(simple_state, 1e-11)
    assert simple_state.t.cpu() == pytest.approx(1e-11, abs=0, rel=1e-6)


# TODO: move to fieldterm test
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
    assert state.t.cpu() == pytest.approx(1e-11, abs=0, rel=1e-6)
