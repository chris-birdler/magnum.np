import pytest
import pathlib
import torch
from magnumnp import *

def test_energy_cube():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": Ms}
    demag = DemagField()

    state.m = state.Constant([1,0,0])

    h = demag.h(state)
    E = demag.E(state)
    assert E == pytest.approx(1./6.*mesh.volume*constants.mu_0*Ms**2)

def test_Ms_domain():
    n  = (4,4,4)
    dx = (1e-9, 1e-9, 1e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant([Ms])}
    state.m = state.Constant([1,0,0])
    demag = DemagField()

    h = demag.h(state)
    E = demag.E(state)
    assert E == pytest.approx(1./6.*mesh.volume*constants.mu_0*Ms**2)
