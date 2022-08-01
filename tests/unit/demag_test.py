import pytest
import pathlib
import torch
from magnumnp import *

def test_energy_cube():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    material = {"Ms": Ms}
    state = State(mesh, material)
    demag = DemagField(state)

    state.m = state.zeros(n + (3,))
    state.m[:,:,:,0] = 1.0

    h = demag.h(0,state.m)
    E = demag.E(0,state.m)
    assert E == pytest.approx(1./6.*mesh.volume*constants.mu_0*Ms**2)

def test_Ms_domain():
    n  = (4,4,4)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    material = {}
    state = State(mesh, material)

    state._material['Ms'] = state.zeros(n + (1,))
    state._material['Ms'][:,:,:,0] = 1.

    state.m = state.zeros(n + (3,))
    state.m[:,:,:,0] = 1.0

    demag = DemagField(state)

    h = demag.h(0,state.m)
    E = demag.E(0,state.m)

def test_symmetry():
    n  = (30, 1, 1)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    material = {"Ms": 3e8}
    state = State(mesh, material)

    demag = DemagField(state)
    N = state.zeros([1 if i==1 else 2*i for i in mesh.n] + [6])

    for i, t in enumerate(((newell_f,dipole_f,0,1,2),
                           (newell_g,dipole_g,0,1,2),
                           (newell_g,dipole_g,0,2,1),
                           (newell_f,dipole_f,1,2,0),
                           (newell_g,dipole_g,1,2,0),
                           (newell_f,dipole_f,2,0,1))):
        demag._init_N_component(N, i, t[2:], t[0], t[1])

    N_fft = torch.fft.rfftn(N, dim = [i for i in range(3) if mesh.n[i] > 1])
    assert torch.allclose(N_fft.imag, state.tensor([0.]))
    
