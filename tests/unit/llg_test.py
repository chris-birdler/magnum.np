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
    llg.step(simple_state, 5e-11)
    assert simple_state.t.cpu() == pytest.approx(5e-11, abs=0, rel=1e-6)
    torch.testing.assert_close(simple_state.avg(simple_state.m), torch.tensor([-0.8912,0.1284,0.0115]), atol=1e-3, rtol=1e-3)

def test_step_adjoint(simple_state):
    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField(simple_state.Constant([-24.6e-3/constants.mu_0,
                                                    +4.3e-3/constants.mu_0,
                                                     0.0]))

    simple_state.m = simple_state.Constant([1,0,0])
    simple_state.m[5:,:,:,0] = -1.0

    llg = LLGSolver([demag, exchange, external], solver = TorchDiffEqAdjoint, adjoint_parameters = [])
    llg.step(simple_state, 5e-11)
    assert simple_state.t.cpu() == pytest.approx(5e-11, abs=0, rel=1e-6)
    torch.testing.assert_close(simple_state.avg(simple_state.m), torch.tensor([-0.8912,0.1284,0.0115]), atol=1e-3, rtol=1e-3)


@pytest.mark.parametrize("solver", [ScipyOdeint, TorchDiffEq])
def test_solve(simple_state, solver):
    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField(simple_state.Constant([-24.6e-3/constants.mu_0,
                                                    +4.3e-3/constants.mu_0,
                                                     0.0]))

    simple_state.m = simple_state.Constant([1,0,0])
    simple_state.m[5:,:,:,0] = -1.0

    llg = LLGSolver([demag, exchange, external], solver = solver)
    tt = torch.linspace(0., 5e-11, steps=5)
    res = llg.solve(simple_state, tt)
    assert simple_state.t.cpu() == pytest.approx(5e-11, abs=0, rel=1e-6)
    torch.testing.assert_close(simple_state.avg(simple_state.m), torch.tensor([-0.8912,0.1284,0.0115]), atol=1e-3, rtol=1e-3)

def test_solve_adjoint(simple_state):
    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField(simple_state.Constant([-24.6e-3/constants.mu_0,
                                                    +4.3e-3/constants.mu_0,
                                                     0.0]))

    simple_state.m = simple_state.Constant([1,0,0])
    simple_state.m[5:,:,:,0] = -1.0

    llg = LLGSolver([demag, exchange, external], solver = TorchDiffEqAdjoint, adjoint_parameters = [])
    tt = torch.linspace(0., 5e-11, steps=5)
    res = llg.solve(simple_state, tt)
    assert simple_state.t.cpu() == pytest.approx(5e-11, abs=0, rel=1e-6)
    torch.testing.assert_close(simple_state.avg(simple_state.m), torch.tensor([-0.8912,0.1284,0.0115]), atol=1e-3, rtol=1e-3)


@pytest.mark.parametrize("solver", [RKF45, ScipyODE, ScipyOdeint, TorchDiffEq])
def test_precession(solver):
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {
        "Ms": state.Constant(8e5),
        "alpha": state.Constant(0.00)
        }
    state.m = state.Constant([0.1,0,1])
    normalize(state.m)

    hz = 0.1/constants.mu_0
    external = ExternalField(state.Constant([0.0, 0.0, hz]))

    f0 = constants.gamma * hz  / (2. * torch.pi) # lamour frequency

    llg = LLGSolver([external], solver = solver, rtol = 1e-5, atol = 1e-5)
    llg.step(state, 1/f0)

    assert state.avg(state.m)[1].cpu() == pytest.approx(0., abs=1e-4, rel=0)

@pytest.mark.parametrize("solver", [RKF45, ScipyODE, ScipyOdeint, TorchDiffEq])
def test_relax(solver):
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {
        "Ms": state.Constant(8e5),
        "Ku": state.Constant(1e5),
        "Ku_axis": state.Constant([0,0,1]),
        "alpha": state.Constant(0.01)
        }
    state.m = state.Constant([1,0,0.1])
    normalize(state.m)

    aniso = UniaxialAnisotropyField()

    llg = LLGSolver([aniso], solver = solver)
    llg.relax(state)

    torch.testing.assert_close(state.avg(state.m), torch.tensor([0.,0.,1.]), atol=1e-2, rtol=1e-2)

def test_stochastic():
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)

    state = State(mesh)
    state.material = {
        "Ms": state.Constant(8e5),
        "Ku": state.Constant(1e5),
        "Ku_axis": state.Constant([0,0,1]),
        "alpha": state.Constant(0.01)
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



@pytest.mark.parametrize("dtype, tol", [(torch.float64, 1e-12), (torch.float32, 1e-5)])
@pytest.mark.parametrize("n_fields", [1, 3, 5])
@pytest.mark.parametrize("alpha_field", [True, False])
@pytest.mark.parametrize("no_precession", [False, True])
def test_llg_rhs_sum_against_reference(dtype, tol, n_fields, alpha_field, no_precession):
    """The fused right-hand side (field sum and cross products in one pass, audit item A1) gives the result of
    accumulate_h and _llg_rhs within round-off (PLAN 11, S1)."""
    from magnumnp.solvers.llg import _llg_rhs, _llg_rhs_sum
    torch.set_default_dtype(dtype)
    g = torch.Generator().manual_seed(2)
    n = (6, 5, 4)
    r = torch.randn(n + (3,), generator=g, dtype=torch.float64)
    m = (r / r.norm(dim=-1, keepdim=True)).to(dtype)
    hs = [(1e5 * torch.randn(n + (3,), generator=g, dtype=torch.float64)).to(dtype) for _ in range(n_fields)]
    alpha = (0.05 + 0.1 * torch.rand(n + (1,), generator=g, dtype=torch.float64)).to(dtype) if alpha_field else 0.1
    ref = _llg_rhs(m, accumulate_h(iter(hs)), alpha, no_precession)
    new = _llg_rhs_sum(m, hs, alpha, no_precession)
    assert float((new - ref).abs().max() / ref.abs().max()) < tol
