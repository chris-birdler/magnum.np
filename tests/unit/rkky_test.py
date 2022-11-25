import pytest
import pathlib
import torch
from math import pi, cos, sin
from magnumnp import *

def test_call():
    n  = (10, 10, 2)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    
    state = State(mesh)
    state.material = {"Ms":1./constants.mu_0, "A":1e-11}
    state.m = state.Constant([1,0,0])
    
    domain1 = state.Constant(False, dtype=torch.bool)
    domain1[:,:,0] = True
    
    domain2 = state.Constant(False, dtype=torch.bool)
    domain2[:,:,1] = True
    
    exchange1 = ExchangeField(domain1)
    exchange2 = ExchangeField(domain2)
    rkky = RKKYField(-0.002, "z", 0, 1)
    
    for phi in torch.linspace(0, 2*pi, 100):
        state.m[domain2] = state.Tensor([cos(phi), sin(phi), 0])
        E_rkky = rkky.E(state).detach().cpu().numpy()
        E_ex1 = exchange1.E(state).detach().cpu().numpy()
        E_ex2 = exchange1.E(state).detach().cpu().numpy()
        h = rkky.h(state)
        print("phi:", phi.cpu().numpy(), "E_rkky:", E_rkky, "E_ex1:", E_ex1, "E_ex2:", E_ex2, "E_tot:", E_rkky+E_ex1+E_ex2, "m:", state.m[:,:,:,0].mean().cpu().numpy(), state.m[:,:,:,1].mean().cpu().numpy(), state.m[:,:,:,2].mean().cpu().numpy())
    #assert E == pytest.approx(1./6.*mesh.volume*constants.mu_0*Ms**2)

