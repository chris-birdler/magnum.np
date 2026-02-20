# %% [markdown]
# # Slonczewski STT Demo 2
#
# Example 2 taken from the Mumax3 Website (https://mumax.github.io/examples.html)
# Flipping the magnetization of a elliptical nanodisk with Slonczewski STT.
# Implemented by Jed Cheng (jed.chengXmag.ed.kyushu-u.ac.jp, replace X with @)

# %% [markdown]
# ## Run Simulation

# %%
from magnumnp import *
import torch
import numpy as np
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
n = (64, 32, 1)
dx = (2.5e-9, 2.5e-9, 5e-9)
L = (n[0]*dx[0], n[1]*dx[1], n[2]*dx[2])
origin = (-L[0]/2., -L[1]/2., -L[2]/2.)
mesh = Mesh(n, dx, origin = origin)
state = State(mesh)

state.material = {
    "Ms": 8e5,
    "A": 1.3e-11,
    "alpha": 0.01,
    "P": 0.5669,
    "Lambda": 1,
    "epsilon_prime": 0,
    "mp": [np.cos(np.radians(20)), np.sin(np.radians(20)), 0],
    "d": L[2],
    "J": -8e11,
    }

x, y, z = mesh.SpatialCoordinate()
disk = (x/80e-9)**2 + (y/40e-9)**2 < 1
write_vti(disk, "data/domain.vti")

state.m = state.Constant([0.,0.,0.])
state.m[disk] = torch.tensor([1.,0.,0.])

# initialize field terms
demag    = DemagField()
exchange = ExchangeField()
torque   = SpinTorqueSlonczewski()

# perform integration with spin torque
m_disk = lambda state: state.m[disk] # only consider disk region
llg = LLGSolver([demag, exchange, torque])
logger = Logger(this_dir / "data", ['t', 'm', m_disk], ["m"], fields_every = 10)

for i in tqdm(torch.arange(0, 1e-9, dt)):
    llg.step(state, dt)
    logger << state

Timer.print_report()
