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
        "alpha": 0.02
        }

    state.m = state.Constant([1,0,0])
    return state
