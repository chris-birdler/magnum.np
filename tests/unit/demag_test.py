import pytest
import pathlib
from magnumnp import *

def test_energy_cube():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1.
    mesh = Mesh(n, dx)
    material = {"Ms": Ms}
    state = State(mesh, material)
    demag    = DemagField(state)

    state.m = state.zeros(n + (3,))
    state.m[:,:,:,0] = 1.0

    h = demag.h(0,state.m)
    E = demag.E(0,state.m)
    assert E == pytest.approx(1./6.*mesh.volume*constants.mu_0*Ms**2)
