from magnumnp import *
import torch

Timer.enable()

# initialize mesh
eps = 1e-15
n  = (1, 1, 1)
dx = (1.0e-9, 1.0e-9, 1e-9)
material = {}

mesh = Mesh(n, dx)
state = State(mesh, material)

# initialoze polarization, p, and charge current amplitude
# thickness of thin film on which the SOT acts
p = state.tensor((0, -1, 0))
je = 6.9e10 
d = n[2] * dx[2]
Keff = 1200e3*constants.mu_0*0.4/2./constants.mu_0

# initialize material
state._material = {
        "Ms": 1200e3,
        "A": 15e-12,
        "K": Keff,
        "K_axis": [0, 0, 1], 
        "gamma": 2.211e5,
        "alpha": 0.048,
        "eta_damp": -0.1,# both eta with opposite sign as magnum.af, same as magnum.pi
        "eta_field": 0.3,
        "p": p,
        "d": d,
        "je": je
        }

state.m = state.zeros(n + (3,))
state.m[:,:,:,2]   = 1.

# initialize field terms
exchange = ExchangeField(state)
aniso = AnisotropyField(state)
torque = SpinOrbitTorque(state)

#relax magnetization state before applying SOT
llg = LLGSolver(state, [aniso, exchange])
while state.t < 2e-9-eps:
    llg.step(1e-9)

# initialize LLG, and solve for 10 ns
Timer.enable()
llg = LLGSolver(state, [aniso, exchange, torque])
state.t = 0.0
with open('data/log.dat', 'w') as f:
    i = 0
    while state.t < 2e-9:
        llg.step(1e-12)
        f.write("%g %g %g %g\n" % ((state.t,) + tuple(torch.mean(state.m, axis=(0,1,2)))))
        i += 1

Timer.print_report()
