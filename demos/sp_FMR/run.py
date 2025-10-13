# %% [markdown]
# # FMR Standard Problem

# %% [markdown]
# ## Run Simulation

# %%
from magnumnp import *
import torch
import pathlib
from tqdm import tqdm

set_log_level(25) # show info_green, but hide info_blue
Timer.enable()
try:
    this_dir = pathlib.Path(__file__).resolve().parent
except:
    this_dir = pathlib.Path().resolve()

# initialize state
dt = 5e-12
n = (24, 24, 2)
l = (120e-9, 120e-9, 10e-9)
dx = (l[0]/n[0], l[1]/n[1], l[2]/n[2])
origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.)
mesh = Mesh(n, dx, origin=origin)

state = State(mesh)
state.material = {
    "alpha": 0.008,
    "Ms": 800e3,
    "A": 13e-12
    }
state.m = state.Constant([0, 0, 1])

# initialize field terms
demag    = DemagField()
exchange = ExchangeField()
bias     = ExternalField(80e3*normalize(state.Constant([1., 0.715, 0])))

#relax state
minimizer = MinimizerBB([demag, exchange, bias])
minimizer.minimize(state)

#integrate
bias = ExternalField(80e3*normalize(state.Constant([1., 0.7, 0])))

llg = LLGSolver([demag, exchange, bias])
logger = Logger(this_dir / "data", ['t', 'm'])

for i in tqdm(torch.arange(0, 10e-9, dt)):
    llg.step(state, dt)
    logger << state

Timer.print_report()
