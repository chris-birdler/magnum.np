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
material = {"Ms":Ms, "A":A}

state = State(mesh, material)
state.m = state.zeros(n + (3,))
state.m[:,:,:,0] = 1.0

domain1 = state.zeros(n, dtype=torch.bool)
domain1[:,:,0] = True

domain2 = state.zeros(n, dtype=torch.bool)
domain2[:,:,1] = True

exchange1 = ExchangeField(state, domain1)
exchange2 = ExchangeField(state, domain2)
rkky = RKKYField(state, J_rkky, "z", 0, 1)

for phi in torch.linspace(0, 2*np.pi, 100):
    state.m[domain2] = torch.DoubleTensor([np.cos(phi), np.sin(phi), 0])
    E_rkky = rkky.E(0.,state.m).detach().cpu().numpy()
    E_ex1 = exchange1.E(0.,state.m).detach().cpu().numpy()
    E_ex2 = exchange1.E(0.,state.m).detach().cpu().numpy()
    h = rkky.h(0., state.m)
    print("phi:", phi.numpy(), "E_rkky:", E_rkky, "E_ex1:", E_ex1, "E_ex2:", E_ex2, "E_tot:", E_rkky+E_ex1+E_ex2, "m:", state.m[:,:,:,0].mean().numpy(), state.m[:,:,:,1].mean().numpy(), state.m[:,:,:,2].mean().numpy())
