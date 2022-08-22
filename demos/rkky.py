from magnumnp import *
import torch
import numpy as np

Timer.enable()
A = 10e-12
J_rkky = -1e-3
Ms = 1e6

# initialize mesh
n  = (10, 10, 2)
dx = (1e-9, 1e-9, 1e-9)
mesh = Mesh(n, dx)
state = State(mesh)
state.material = {"Ms": Ms, "A": A}
state.m = state.Constant((1,0,0))

domain1 = state.Constant(False, dtype=torch.bool)
domain1[:,:,0] = True

domain2 = state.Constant(False, dtype=torch.bool)
domain2[:,:,1] = True

exchange1 = ExchangeField(domain1)
exchange2 = ExchangeField(domain2)
rkky = RKKYField(J_rkky, "z", 0, 1)

for phi in torch.linspace(0, 2*np.pi, 100):
    state.m[domain2] = state._tensor([np.cos(phi), np.sin(phi), 0])
    E_rkky = rkky.E(state).detach().cpu().numpy()
    E_ex1 = exchange1.E(state).detach().cpu().numpy()
    E_ex2 = exchange1.E(state).detach().cpu().numpy()
    h = rkky.h(state)
    print("phi:", phi.numpy(), "E_rkky:", E_rkky, "E_ex1:", E_ex1, "E_ex2:", E_ex2, "E_tot:", E_rkky+E_ex1+E_ex2, "m:", state.m[:,:,:,0].mean().numpy(), state.m[:,:,:,1].mean().numpy(), state.m[:,:,:,2].mean().numpy())
