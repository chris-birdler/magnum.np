from magnumnp import *
import torch
from torch import sin, cos

# initialize mesh
n  = (50, 50, 1)
dx = (5e-9, 5e-9, 5e-9)
mesh = Mesh(n, dx)
h = 300e-3/constants.mu_0
tf = 2e-9
dt = 1e-12

# initialize material
material = {
           "Ms": 8e5,
           "A": 1.3e-11,
           "alpha": 0.08,
           "K": 1e5,
           "K_axis": [0,0,1]
           }

state = State(mesh, material)

m0 = state.Constant([0,0,0])
m0[:, :, :, 2] = 1
m0[:2, :, :, 0] = 0.1
m0[2:, :, :, 0] = -0.1
m0.normalize()
write_vti(m0, "data/m0.vti", state)

m_t = state.Constant([0,0,0])
m_t[:, :, :, 0] = 1
m_t[:, :, :, 2] = 1
m_t.normalize()

phi = state.Tensor([0.], requires_grad = True)
theta = state.Tensor([1.5708], requires_grad = True)
h_ext = lambda t: torch.stack([sin(theta)*cos(phi),
                               sin(theta)*sin(phi),
                               cos(theta)], dim=-1)

#h_ext = state.Tensor([1,0,0], requires_grad = True)

demag    = DemagField(state)
exchange = ExchangeField(state)
aniso    = AnisotropyField(state)
external = ExternalField(state, h_ext)

optimizer = torch.optim.Adam([phi, theta], lr=0.05)
#optimizer = torch.optim.Adam([h_ext], lr=0.05)
my_loss = torch.nn.L1Loss()

for epoch in range(100):
    print("epoch:", epoch)
    optimizer.zero_grad()
    with Timer("forward"):
        llg = LLGSolver(state, [demag, exchange, external, aniso])
        state.m = m0
        state.t = 0
        llg.step(1e-11)#tf)
        L = my_loss(state.m, m_t)
        write_vti(state.m, 'data/m1_%04d.vti' % epoch)

    with Timer("backward"):
        L.backward()
    with Timer("optimizer"):
        optimizer.step()
    print("phi:", phi)
    print("theta:", theta)
