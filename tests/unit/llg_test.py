import pytest
import torch
from magnumnp import *
from helpers import *

@pytest.mark.parametrize("solver", [RKF45, ScipyODE, ScipyOdeint, TorchDiffEq])
def test_step(simple_state, solver):
    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField(simple_state.Constant([-24.6e-3/constants.mu_0,
                                                    +4.3e-3/constants.mu_0,
                                                     0.0]))

    simple_state.m = simple_state.Constant([1,0,0])
    simple_state.m[5:,:,:,0] = -1.0

    llg = LLGSolver([demag, exchange, external], solver = solver)
    llg.step(simple_state, 1e-11)
    assert simple_state.t.cpu() == pytest.approx(1e-11, abs=0, rel=1e-6)


@pytest.mark.parametrize("solver", [RKF45, ScipyODE, ScipyOdeint, TorchDiffEq])
def test_precession(solver):
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {
        "Ms": state.Constant([8e5]),
        "alpha": state.Constant([0.00])
        }
    state.m = state.Constant([0.1,0,1])
    normalize(state.m)

    hz = 0.1/constants.mu_0
    external = ExternalField(state.Constant([0.0, 0.0, hz]))

    f0 = constants.gamma * hz  / (2. * torch.pi) # lamour frequency

    llg = LLGSolver([external], solver = solver, rtol = 1e-5, atol = 1e-5)
    llg.step(state, 1/f0)

    assert avg(state.m)[1].cpu() == pytest.approx(0., abs=1e-4, rel=0)

@pytest.mark.parametrize("solver", [RKF45, ScipyODE, ScipyOdeint, TorchDiffEq])
def test_relax(solver):
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {
        "Ms": state.Constant([8e5]),
        "Ku": state.Constant([1e5]),
        "Ku_axis": state.Constant([0,0,1]),
        "alpha": state.Constant([0.01])
        }
    state.m = state.Constant([1,0,0.1])
    normalize(state.m)

    aniso = UniaxialAnisotropyField()

    llg = LLGSolver([aniso], solver = solver)
    llg.relax(state, rtol = 1e-8)

    torch.testing.assert_close(avg(state.m), state.Tensor([0,0,1]), atol=1e-2, rtol=1e-2)

def test_stochastic():
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {
        "Ms": state.Constant([8e5]),
        "Ku": state.Constant([1e5]),
        "Ku_axis": state.Constant([0,0,1]),
        "alpha": state.Constant([0.01])
        }
    state.m = state.Constant([0,0,1])
    normalize(state.m)
    state.T = 300

    aniso = UniaxialAnisotropyField()

    llg = LLGSolver([aniso])
    logger = Logger("data", ['t', 'm'])
    while state.t < 1e-9:
        llg.step(state, 1e-11)
        logger << state

