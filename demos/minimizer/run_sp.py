import torch
from magnumnp import *

#torch.set_default_dtype(torch.float32)
Timer.enable(log_mem = True)

N = 10

# initialize state
eps = 1e-15
nx, ny, nz = 7*N, N, N 
dx, dy, dz = 7*20e-9 / nx, 20e-9 / ny, 20e-9 / nz

mesh = Mesh((nx,ny,nz), (dx,dy,dz))
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
state.material["Ku_axis"][cube1] = state.Tensor([-0.04202257, -0.12252924, 0.99157485])
state.material["Ku_axis"][cube2] = state.Tensor([0.04902541, 0.0919583, 0.99455527])
state.material["Ku_axis"][cube3] = state.Tensor([-0.07325518, 0.0437713, 0.99635222])
state.material["Ku_axis"][cube4] = state.Tensor([0.02283161, -0.08166655, 0.99639816])

write_vti(state.material, "data/material.vti", state)

# initialize field terms
demag    = DemagField()
exchange = ExchangeField()
aniso    = UniaxialAnisotropyField()
external = ExternalField([0.0, 0.0, 0.0])

# initialize magnetization
state.m = state.Constant([0,0,1])
state.m[~magnetic] = 0.

#minimizer = MinimizerBB([demag, exchange, aniso, external])
minimizer = MinimizerBB([exchange, aniso, external], samples=10)
m_magnetic = lambda state: state.m[magnetic].sum() / magnetic.sum()
logger = Logger("data", [external.h, 'm', m_magnetic], ["m"])

for i in range(4000):
    external.h = [0.0, 0.0, -i * 1e-3 / constants.mu_0]
    E, steps = minimizer.minimize2(state)
    print("i:", i, "E:", E, "steps:", steps)
    logger << state

Timer.print_report()
