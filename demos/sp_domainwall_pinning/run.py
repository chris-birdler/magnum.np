# %% [markdown]
# # Domainwall Pinning Standard Problem

# %% [markdown]
# ## Run Simulation

# %%
from magnumnp import *
import numpy as np
import pathlib

Timer.enable()
try:
    this_dir = pathlib.Path(__file__).resolve().parent
except:
    this_dir = pathlib.Path().resolve()

Hextmax=1.8/constants.mu_0
Hextmin=1.4/constants.mu_0
tfinal = 20e-9

n  = (80, 1, 1)
dx = (1e-9, 1e-9, 1e-9)
origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.,)
mesh = Mesh(n, dx, origin)
state = State(mesh)
state.material = {"alpha": state.Constant(1.)}
x, y, z = mesh.SpatialCoordinate()

soft = (x < 0)
hard = (x >= 0)

Ms = state.Constant(0.0)
Ms[soft] = 0.25/constants.mu_0
Ms[hard] = 1./constants.mu_0

Ku = state.Constant(0.0)
Ku[soft] = 1e5
Ku[hard] = 1e6

A = state.Constant(0.0)
A[soft] = 0.25e-11
A[hard] = 1.e-11

state.material['Ms'] = Ms
state.material['A'] = A
state.material['Ku'] = Ku
state.material['Ku_axis'] = state.Constant([0,1,0])

state.m = state.Constant([np.sin(0.3), np.cos(0.3), 0.0])
state.m[hard] = torch.tensor([0.0, -1.0, 0.0])
normalize(state.m)

exchange = ExchangeField()
aniso    = UniaxialAnisotropyField()
external = ExternalField(state.Constant([0, 0, 0]))

minimizer = MinimizerBB([exchange, aniso, external])
logger = ScalarLogger(this_dir / "data" / "m.dat", ['t', external.h, 'm'])
for h in torch.linspace(Hextmin, Hextmax, steps=100):
    external.h = state.Constant([0, h, 0])
    minimizer.minimize(state)
    logger << state

Timer.print_report()
