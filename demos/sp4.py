from magnumnp import *
import torch

Timer.enable()

# initialize mesh
eps = 1e-15
n  = (100, 25, 1)
dx = (5e-9, 5e-9, 3e-9)
mesh = Mesh(n, dx)

# initialize material
material = {
        "Ms": 8e5,
        "A": 1.3e-11,
        "gamma": 2.211e5,
        "alpha": 0.02
        }

state = State(mesh, material)

# initialize field terms
demag    = DemagField(state)
exchange = ExchangeField(state)
external = ExternalField(state, [-24.6e-3/constants.mu_0,
                                 +4.3e-3/constants.mu_0,
                                 0.0])

# initialize magnetization that relaxes into s-state
state.m = state.zeros(n + (3,))
state.m[1:-1,:,:,0]   = 1.0
state.m[(-1,0),:,:,1] = 1.0

# initialize sstate
minimizer = Minimizer(state, [demag, exchange])
#m = minimizer.minimize(m0, 1e-2, 1e-4)
minimizer.minimize_event()
write_vtr(state.m, "data/sp4_m0")

# perform integration with external field
llg = LLGSolver(state, [demag, exchange, external])
with open('data/sp4.dat', 'w') as f:
    while state.t < 1e-9-eps:
        llg.step(1e-11)
        f.write("%g %g %g %g\n" % ((state.t,) + tuple(torch.mean(state.m, axis=(0,1,2)))))

Timer.print_report()
