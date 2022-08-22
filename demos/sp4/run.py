from magnumnp import *
import torch

set_log_level(10)
Timer.enable()

# initialize mesh
eps = 1e-15
n  = (100, 25, 1)
dx = (5e-9, 5e-9, 3e-9)
mesh = Mesh(n, dx)
state = State(mesh)

state.material = {
    "Ms": 8e5,
    "A": 1.3e-11,
    "alpha": 1.00
    }


# initialize field terms
demag    = DemagField()
exchange = ExchangeField()
external = ExternalField(state, [-24.6e-3/constants.mu_0,
                                 +4.3e-3/constants.mu_0,
                                 0.0])

# initialize magnetization that relaxes into s-state
state.m = state.Constant([0,0,0])
state.m[1:-1,:,:,0]   = 1.0
state.m[(-1,0),:,:,1] = 1.0

# relax without external field
#llg = LLGSolver(state, [demag, exchange])
#llg.relax()
#write_vti(state.m, "data/m0.vti", state)

# perform integration with external field
state.material["alpha"] = 0.02
llg = LLGSolver([demag, exchange, external])
with open('data/m.dat', 'w') as f:
    while state.t < 1e-9-eps:
        llg.step(state, 1e-11)
        f.write("%g %g %g %g\n" % ((state.t,) + tuple(torch.mean(state.m, axis=(0,1,2)))))
        f.flush()

Timer.print_report()
