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
    default = torch.get_default_dtype()
    torch.set_default_dtype(dtype)
    try:
        _check_rhs_sum(dtype, tol, n_fields, alpha_field, no_precession)
    finally:
        torch.set_default_dtype(default)


def _check_rhs_sum(dtype, tol, n_fields, alpha_field, no_precession):
    from magnumnp.solvers.llg import _llg_rhs, _llg_rhs_sum
    g = torch.Generator().manual_seed(2)
    dev = torch.empty(0).device                     # the default device of magnum.np (a GPU with CUDA_DEVICE set)
    n = (6, 5, 4)
    r = torch.randn(n + (3,), generator=g, dtype=torch.float64, device="cpu")
    m = (r / r.norm(dim=-1, keepdim=True)).to(dtype=dtype, device=dev)
    hs = [(1e5 * torch.randn(n + (3,), generator=g, dtype=torch.float64, device="cpu")).to(dtype=dtype, device=dev)
          for _ in range(n_fields)]
    alpha = (0.05 + 0.1 * torch.rand(n + (1,), generator=g, dtype=torch.float64, device="cpu")).to(dtype=dtype, device=dev) \
        if alpha_field else 0.1
    ref = _llg_rhs(m, accumulate_h(iter(hs)), alpha, no_precession)
    new = _llg_rhs_sum(m, hs, alpha, no_precession)
    assert new.shape == m.shape
    assert float((new - ref).abs().max() / ref.abs().max()) < tol
    # alpha as a material value of a State (a tensor subclass, [n, 1]) and a drive field that is expanded
    state = State(Mesh(n, (1e-9,) * 3))
    state.material = {"alpha": 0.1}
    hs[-1] = hs[-1][0, 0, 0].expand(n + (3,))
    ref = _llg_rhs(m, accumulate_h(iter(hs)), state.material["alpha"], no_precession)
    new = _llg_rhs_sum(m, hs, state.material["alpha"], no_precession)
    assert new.shape == m.shape
    assert float((new - ref).abs().max() / ref.abs().max()) < tol


@pytest.mark.parametrize("dtype, tol", [(torch.float64, 1e-12), (torch.float32, 1e-5)])
def test_dm_compiles_the_rhs_with_a_cow_alpha(dtype, tol):
    """PLAN 15 F3: a homogeneous alpha is a CoWTensor in the State; dm gives _llg_rhs_sum a plain view of it, so that
    torch.compile traces the whole function (with the subclass it traced only the field sum). The result equals the
    eager function within tol of the largest value."""
    import torch._dynamo as dyn
    from magnumnp.common.cow import CoWTensor
    from magnumnp.solvers.llg import _llg_rhs_sum
    default = torch.get_default_dtype()
    torch.set_default_dtype(dtype)
    try:
        n = (6, 5, 4)
        state = State(Mesh(n, (1e-9,) * 3))
        state.material = {"alpha": 0.1, "Ms": 8e5}
        assert type(state.material["alpha"]) is CoWTensor          # the case of production (run_loops)
        g = torch.Generator().manual_seed(3)
        r = torch.randn(n + (3,), generator=g, dtype=torch.float64).to(dtype)
        m = r / r.norm(dim=-1, keepdim=True)
        hs = [torch.randn(n + (3,), generator=g, dtype=torch.float64).to(dtype) * 1e5 for _ in range(3)]

        class Fixed:                                               # a field term that returns a fixed field
            def __init__(self, h):
                self._h = h

            def h(self, state):
                return self._h
        llg = LLGSolver([Fixed(h) for h in hs])
        dm = llg.dm(0.0, m, state)
        eager = _llg_rhs_sum._torchdynamo_orig_callable
        ref = eager(m, hs, state.material["alpha"].as_subclass(torch.Tensor), False)
        assert float((dm - ref).abs().max()) <= tol * float(ref.abs().max())
        alpha = state.material["alpha"]
        for a, traced in ((alpha, False), (alpha.as_subclass(torch.Tensor), True)):
            dyn.reset()
            e = dyn.explain(eager)(m, hs, a, False)
            assert (e.graph_break_count == 0 and e.op_count >= 8) == traced
    finally:
        torch.set_default_dtype(default)
