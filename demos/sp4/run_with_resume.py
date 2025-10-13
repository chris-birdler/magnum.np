from magnumnp import *
import torch

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
    "alpha": 0.02
    }

# initialize field terms
demag    = DemagField()
exchange = ExchangeField()
external = ExternalField([-24.6e-3/constants.mu_0,
                          +4.3e-3/constants.mu_0,
                          0.0])

# initialize magnetization that relaxes into s-state
state.m = state.Constant([0,0,0])
state.m[1:-1,:,:,0]   = 1.0
state.m[(-1,0),:,:,1] = 1.0

logger = Logger("data", ['t', 'm'], ["m"], fields_every=10)
if not logger.is_resumable():
    # relax without external field
    minimizer = MinimizerBB([demag, exchange])
    minimizer.minimize(state)
    state.write_vtk(state.m, "data/m0")

# perform integration with external field
logger.resume(state)
llg = LLGSolver([demag, exchange, external])
while state.t < 1e-9-eps:
    logger << state
    llg.step(state, 1e-11)

Timer.print_report()
