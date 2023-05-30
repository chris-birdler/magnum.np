from magnumnp import *
import numpy as np

set_log_level(10)
Timer.enable()

n  = (64, 64, 64)
dx = (1e-9, 1e-9, 1e-9)
mesh = Mesh(n, dx)

eta = 30
Js = 1.

state = State(mesh)
state.material = {
    "Ms": Js/constants.mu_0,
    "alpha": 1.0
    }
state.m = state.Constant([0,0,1])
state.m.normalize()
state.T = Js * state.cell_volumes.max() * np.prod(mesh.n) / (constants.kb * eta)
print("T:", state.T.numpy())
#print("cell_volume:", state.cell_volumes)

aniso = UniaxialAnisotropyField()
external = ExternalField([0,0,0])

llg = LLGSolver([external], solver=Heun)
logger = Logger("data", ['t', external.h, 'm'])

# relax
for i in range(10000):
    llg.step(state, 1e-13)

for h in np.linspace(0,0.1, num=11):
    external.h = [h / constants.mu_0, 0, 0]
    for i in range(10000):
        llg.step(state, 1e-13)
    logger << state
        
Timer.print_report()
