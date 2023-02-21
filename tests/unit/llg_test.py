import pytest
import torch
import numpy as np
from magnumnp import *
from helpers import *

@pytest.mark.parametrize("solver", [RKF45, ScipyODE, ScipyOdeint, TorchDiffEq])
def test_step(simple_state, solver):
    if torch.cuda.is_available():
        if solver == ScipyODE or solver == ScipyOdeint:
           pytest.skip()

    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField([-24.6e-3/constants.mu_0,
                              +4.3e-3/constants.mu_0,
                              0.0])

    simple_state.m = simple_state.Constant([1,0,0])
    simple_state.m[5:,:,:,0] = -1.0

    llg = LLGSolver([demag, exchange, external], solver = solver)
    llg.step(simple_state, 1e-11)
    assert simple_state.t.cpu() == pytest.approx(1e-11, abs=0, rel=1e-6)


@pytest.mark.parametrize("solver", [RKF45, ScipyODE, ScipyOdeint, TorchDiffEq])
def test_precession(solver):
    if torch.cuda.is_available():
        if solver == ScipyODE or solver == ScipyOdeint:
           pytest.skip()

    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {
        "Ms": 8e5,
        "alpha": 0.00
        }
    state.m = state.Constant([0.1,0,1])
    state.m.normalize()

    hz = 0.1/constants.mu_0
    external = ExternalField([0.0, 0.0, hz])

    f0 = constants.gamma * hz  / (2. * np.pi) # lamour frequency

    llg = LLGSolver([external], solver = solver, rtol = 1e-5, atol = 1e-5)
    llg.step(state, 1/f0)

    assert state.m.average()[1].cpu() == pytest.approx(0., abs=1e-4, rel=0)

@pytest.mark.parametrize("solver", [RKF45, ScipyODE, ScipyOdeint, TorchDiffEq])
def test_relax(solver):
    if torch.cuda.is_available():
        if solver == ScipyODE or solver == ScipyOdeint:
           pytest.skip()

    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {
        "Ms": 8e5,
        "Ku": 1e5,
        "Ku_axis": state.Tensor([0,0,1]),
        "alpha": 0.01
        }
    state.m = state.Constant([1,0,0.1])
    state.m.normalize()

    aniso = UniaxialAnisotropyField()

    llg = LLGSolver([aniso], solver = solver)
    llg.relax(state, rtol = 1e-8)

    torch.testing.assert_close(state.m.avg(), state.Tensor([0,0,1]), atol=1e-3, rtol=1e-3)

def test_stochastic():
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {
        "Ms": 8e5,
        "Ku": 1e5,
        "Ku_axis": state.Tensor([0,0,1]),
        "alpha": 0.01
        }
    state.m = state.Constant([0,0,1])
    state.m.normalize()
    state.T = 300

    aniso = UniaxialAnisotropyField()

    llg = LLGSolver([aniso])
    logger = Logger("data", ['t', 'm'])
    while state.t < 1e-9:
        llg.step(state, 1e-11)
        logger << state

