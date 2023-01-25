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
        "Ms": 8e5,
        "A": 1.3e-11,
        "Ku": 1e3,
        "Ku_axis": [0, 0, 1],
        "Di": 1e-3,
        "Db": 1e-3,
        "DD2d": 1e-3,
        "eta_damp": -0.1,
        "eta_field": 0.3,
        "p": [0, -1, 0],
        "d": 3e-9,
        "je": 6.9e10,
        "xi": 0.05,
        "b": 72.17e-12,
        "alpha": 0.02
        }

    state.m = state.Constant([1,0,0])
    return state
