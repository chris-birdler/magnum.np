from magnumnp import *
import torch

Timer.enable()

# initialize mesh
eps = 1e-15
n  = (1, 1, 1)
dx = (1.0e-9, 1.0e-9, 1e-9)
mesh = Mesh(n, dx)

# initialoze polarization, p, and charge current amplitude
# thickness of thin film on which the SOT acts
p = torch.DoubleTensor((0, -1, 0)).to(cuda)
je = 6.9e10 
d = n[2] * dx[2]

# Calculate effective anisotorpy constant to account for shape anisotropy
Keff = 1200e3*constants.mu_0*0.4/2./constants.mu_0

# initialize material
material = {
        "Ms": 1200e3,
        "A": 15e-12,
        "K": Keff,
        "K_axis": [0, 0, 1], 
        "gamma": 2.211e5,
        "alpha": 0.004,
        "eta_damp": -0.1,# both eta with opposite sign as magnum.af, same as magnum.pi
        "eta_field": 0.3,
        "p": p,
        "d": d,
        "je": je
        }


#initialize state class
state = State(mesh, material)


# initialize magnetization parallel to [0, 0, 1] 
state.m = state.zeros(n + (3,))
state.m[:,:,:,2]   = 1.



# initialize field terms
exchange = ExchangeField(state)
aniso = AnisotropyField(state)
torque = SpinOrbitTorque(state)



#relax magnetization state before applying SOT
llg = LLGSolver(state, [aniso, exchange])
with open('data/log_relax.dat', 'w') as f:
    i = 0
    while state.t < 2e-9-eps:
        llg.step(1e-9)
        f.write("%g %g %g %g\n" % ((state.t,) + tuple(torch.mean(state.m, axis=(0,1,2)))))
        i += 1

# print timing for relaxation
Timer.print_report()



# initialize LLG, and solve for 10 ns
llg = LLGSolver(state, [aniso, exchange, torque])
state.t = 0.0
with open('data/log.dat', 'w') as f:
    i = 0
    while state.t < 10e-9:
        llg.step(1e-12)
        f.write("%g %g %g %g\n" % ((state.t,) + tuple(torch.mean(state.m, axis=(0,1,2)))))
        i += 1

# print second time report for SOT
Timer.print_report()
