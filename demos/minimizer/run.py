# %% [markdown]
# # Minimizer Standard Problem
#
# This micromagnetic problem tests the calculation of the demagnetization curve for hard magnetic 
# grains. The grains have slightly different orientation of the anisotropy axes. The critical 
# fields at which the grains revers are the Stoner-Wohlfarth switching fields when the 
# demagnetizing fields are not taken into account.

# | grain | K(MJ/m3) | kx          | ky          | kz         | Js(T) | A(pJ/m) | angle    | Bsw(mT)|
# |-------|----------|-------------|-------------|------------|-------|---------|----------|--------|
# | 1     | 2.37     | -0.04202257 | -0.12252924 | 0.99157485 | 1.61  | 9.2     | 0.129900 | 2646   |
# | 2     | 2.37     | 0.04902541  | 0.0919583   | 0.99455527 | 1.61  | 9.2     | 0.104400 | 2753   |
# | 3     | 2.37     | -0.07325518 | 0.0437713   | 0.99635222 | 1.61  | 9.2     | 0.085440 | 2845   |
# | 4     | 2.37     | 0.02283161  | -0.08166655 | 0.99639816 | 1.61  | 9.2     | 0.084900 | 2848   |


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
N = 10
n = 7*N, N, N 
dx = 7*20e-9 / n[0], 20e-9 / n[1], 20e-9 / n[2]

mesh = Mesh(n, dx)
state = State(mesh)

cube1 = state.Constant(False, dtype=torch.bool)
cube1[0*N:1*N, :, :] = True
cube2 = state.Constant(False, dtype=torch.bool)
cube2[2*N:3*N, :, :] = True
cube3 = state.Constant(False, dtype=torch.bool)
cube3[4*N:5*N, :, :] = True
cube4 = state.Constant(False, dtype=torch.bool)
cube4[6*N:7*N, :, :] = True
magnetic = cube1 | cube2 | cube3 | cube4

state.material["Ms"] = state.Constant([0.])
state.material["Ms"][magnetic] = 1.61/constants.mu_0

state.material["A"] = state.Constant([0.])
state.material["A"][magnetic] = 9.2e-12

state.material["Ku"] = state.Constant([0.])
state.material["Ku"][magnetic] = 2.37e6

state.material["Ku_axis"] = state.Constant([0, 0, 0])
state.material["Ku_axis"][cube1] = torch.tensor([-0.04202257, -0.12252924, 0.99157485])
state.material["Ku_axis"][cube2] = torch.tensor([0.04902541, 0.0919583, 0.99455527])
state.material["Ku_axis"][cube3] = torch.tensor([-0.07325518, 0.0437713, 0.99635222])
state.material["Ku_axis"][cube4] = torch.tensor([0.02283161, -0.08166655, 0.99639816])

state.write_vtk(state.material, "data/m0")

# initialize field terms
demag    = DemagField()
exchange = ExchangeField()
aniso    = UniaxialAnisotropyField()
external = ExternalField([0.0, 0.0, 0.0])

# initialize magnetization
state.m = state.Constant([0,0,1])
state.m[~magnetic] = 0.

#minimizer = MinimizerBB([demag, exchange, aniso, external])
minimizer = MinimizerBB([exchange, aniso, external])
m_magnetic = lambda state: state.m[magnetic].sum() / magnetic.sum()
logger = Logger(this_dir / "data", [external.h, 'm', m_magnetic])

for i in tqdm(range(4000)):
    external.h = [0.0, 0.0, -i * 1e-3 / constants.mu_0]
    steps = minimizer.minimize(state)
    logger << state

Timer.print_report()
