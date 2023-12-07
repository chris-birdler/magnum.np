import pytest
import magnumnp

__all__ = ["simple_mesh", "simple_state"]

@pytest.fixture
def simple_mesh():
    n  = (100, 25, 1)
    dx = (5e-9, 5e-9, 3e-9)
    return magnumnp.Mesh(n, dx)

@pytest.fixture
def simple_state(simple_mesh):
    state = magnumnp.State(simple_mesh)

    state.material = {
        "Ms": state.Constant(8e5),
        "A": state.Constant(1.3e-11),
        "Ku": state.Constant(1e3),
        "Ku_axis": state.Constant([0, 0, 1]),
        "Di": state.Constant(1e-3),
        "Db": state.Constant(1e-3),
        "DD2d": state.Constant(1e-3),
        "eta_damp": state.Constant(-0.1),
        "eta_field": state.Constant(0.3),
        "p": state.Constant([0, -1, 0]),
        "d": state.Constant(3e-9),
        "je": state.Constant(6.9e10),
        "xi": state.Constant(0.05),
        "b": state.Constant(72.17e-12),
        "alpha": state.Constant(0.02)
        }

    state.m = state.Constant([1,0,0])
    return state
