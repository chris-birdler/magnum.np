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

# initialize state
state = State(mesh)
state.material = {
    "Ms": 8e5,
    "A": 1.3e-11,
    "alpha": 0.08,
    "Ku": 1e5,
    #"Ku_axis": [0,0,1]
    "Ku_axis": state.Tensor([0,0,1])
    }

m0 = state.Constant([0,0,0])
m0[:, :, :, 2] = 1
m0[:2, :, :, 0] = 0.1
m0[2:, :, :, 0] = -0.1
m0.normalize()
write_vti(m0, "data/m0.vti", state)

m_t = state.Constant([1,0,1])
m_t.normalize()

phi = state.Tensor([0.], requires_grad = True)
theta = state.Tensor([1.5708], requires_grad = True)
h_ext = lambda t: torch.stack([sin(theta)*cos(phi),
                               sin(theta)*sin(phi),
                               cos(theta)], dim=-1)

demag    = DemagField()
exchange = ExchangeField()
aniso    = UniaxialAnisotropyField()
external = ExternalField(h_ext)

optimizer = torch.optim.Adam([phi, theta], lr=0.05)
my_loss = torch.nn.L1Loss()
llg = LLGSolver([demag, exchange, external, aniso])

for epoch in range(100):
    print("epoch:", epoch)
    optimizer.zero_grad()
    with Timer("forward"):
        state.m = m0
        state.t = 0
        llg.step(state, 1e-11)#tf)
        L = my_loss(state.m, m_t)
        write_vti(state.m, 'data/m1_%04d.vti' % epoch)

    with Timer("backward"):
        L.backward()
    with Timer("optimizer"):
        optimizer.step()
    print("phi:", phi)
    print("theta:", theta)
