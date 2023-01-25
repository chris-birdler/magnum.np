from magnumnp import *
import torch
from torch import sin, cos

torch.set_default_dtype(torch.float64)

# initialize mesh
n  = (101, 51, 1)
dx = (5e-9, 5e-9, 5e-9)
mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2.,-n[1]*dx[1]/2.,-n[2]*dx[2]/2.))
state = State(mesh)
state.material = {"Ms": 1.}

state.m = state.Constant([0,1,0], requires_grad = True)
demag = DemagField()

h = demag.h(state)
J = h[n[0]//2, 0, n[2]//2, 1]
#J = -0.5*h[n[0]//2, 0, n[2]//2, 0]**2

J.backward()
write_vti(state.m.grad, "data/m_grad.vti")
m_opt = state.m.grad[...,1] > 0.
write_vti(m_opt, "data/m_opt.vti")
