# %% [markdown]
# # MuMag Standard Problem #4

# %% [markdown]
# ## Run Simulation

# %%
from magnumnp import *
import torch
import pathlib

set_log_level(25) # show info_green, but hide info_blue
Timer.enable()
try:
    this_dir = pathlib.Path(__file__).resolve().parent
except:
    this_dir = pathlib.Path().resolve()

# initialize state
n  = (101, 51, 1)
dx = (5e-9, 5e-9, 5e-9)
mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2.,-n[1]*dx[1]/2.,-n[2]*dx[2]/2.))
state = State(mesh)
Ms = state.Constant([1.], requires_grad = True)
state.material = {"Ms": Ms}
state.m = state.Constant([0,1,0])
demag = DemagField()

h = demag.h(state)
J = h[n[0]//2, 0, n[2]//2, 1]

J.backward()
state.write_vtk(Ms.grad, this_dir / "data" / "m_grad.vti")
m_opt = Ms.grad > 0.
state.write_vtk(m_opt, this_dir / "data" / "m_opt.vti")
torch.save([torch.tensor(n), torch.tensor(dx), Ms.grad], "data/grad.pt")

Timer.print_report()
