from magnumnp import *
import torch
from torch import sin, cos

Timer.enable()

# initialize mesh
n  = (10, 10, 10)
dx = (5e-9, 5e-9, 5e-9)
mesh = Mesh(n, dx)
h = 300e-3/constants.mu_0
tf = 1e-9
dt = 1e-12

# initialize state
state = State(mesh)
state.material = {
    "Ms": 8e5,
    "A": 1.3e-11,
    "alpha": 1,
    "Ku": 1e5,
    "Ku_axis": state.Tensor([0,0,1])
    }

m0 = state.Constant([0,0,0])
m0[:, :, :, 2] = 1
m0[:n[0]//2, :, :, 0] = 0.1
m0[n[0]//2:, :, :, 0] = -0.1
m0.normalize()

m_t = state.Constant([1,0,1])
m_t.normalize()

state.phi = state.Tensor([1.5708], requires_grad = True)
state.theta = state.Tensor([1.5708], requires_grad = True)

h_ext = lambda t: torch.stack([h*sin(state.theta)*cos(state.phi),
                               h*sin(state.theta)*sin(state.phi),
                               h*cos(state.theta)], dim=-1)

demag    = DemagField()
exchange = ExchangeField()
aniso    = UniaxialAnisotropyField()

# relax system
with torch.no_grad():
    state.m = m0
    llg = LLGSolver([demag, exchange, aniso])
    llg.relax(state)
    m0 = state.m
    write_vti(m0, "data/m0.vti")


# start optimization
external = ExternalField(h_ext)
llg = LLGSolver([demag, exchange, aniso, external], solver = TorchDiffEqAdjoint, adjoint_parameters = (state.phi, state.theta))

phi_grad = ("phi_grad", lambda state: state.phi.grad)
theta_grad = ("theta_grad", lambda state: state.theta.grad)
epoch_logger = ScalarLogger("data/epoch.dat", ["phi", "theta", "loss", phi_grad, theta_grad])

optimizer = torch.optim.Adam([state.phi, state.theta], lr=0.05)
my_loss = torch.nn.L1Loss()

for epoch in range(100):
    print("epoch: ", epoch)

    optimizer.zero_grad()
    state.m = m0
    state.t = 0.

    logger = ScalarLogger("data/m_%03d.dat" % epoch, ['t', 'm'])
    with Timer("forward"):
        while state.t < tf:
            llg.step(state, dt)
            logger << state

    with Timer("backward"):
        state.loss = my_loss(state.m, m_t)
        state.loss.backward()

    optimizer.step()

    epoch_logger << state

Timer.print_report()
