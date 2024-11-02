from magnumnp import *
# %% [markdown]
# # DMI Standard Problem

# %% [markdown]
# ## Run Simulation

# %%
from magnumnp import *
import torch
import numpy as np
import pathlib

Timer.enable()
try:
    this_dir = pathlib.Path(__file__).resolve().parent
except:
    this_dir = pathlib.Path().resolve()

n  = (100, 1, 1)
dx = (1e-9, 1e-9, 1e-9)
origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.,)

mesh = Mesh(n, dx, origin)
state = State(mesh)
state.material.set({
    "alpha": 1.,
    "A": 13e-12,
    "Ms": 800e3,
    "Di": -3e-3,
    "Ku": 0.4e6,
    })
state.material['Ku_axis'] = state.Constant([0,0,1])
x, y, z = mesh.SpatialCoordinate()

state.write_vtk(state.material, "data/material")

state.m = state.Constant([-0.1, 0.0, 0.9])
normalize(state.m)

exchange = ExchangeField()
aniso    = UniaxialAnisotropyField()
dmi      = InterfaceDMIField()

minimizer = MinimizerBB([exchange, aniso, dmi])
minimizer.minimize(state)
np.savetxt(this_dir / "data" / "m0_magnumnp.dat", torch.concat((x[:,0,0,None], state.m[:,0,0,:]), axis=1).cpu().numpy())

Timer.print_report()
