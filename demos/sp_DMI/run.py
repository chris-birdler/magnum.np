from magnumnp import *
import torch
import numpy as np
import os
import math

Timer.enable()

n  = (102, 1, 1)
dx = (1e-9, 1e-9, 1e-9)
origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.,)

mesh = Mesh(n, dx, origin)
state = State(mesh)
state.material = {"alpha": 1.}
x, y, z = state.SpatialCoordinates()

A = state.Constant([13e-12]) 
Ms = state.Constant([0.86e6]) 
Ku = state.Constant([0.4e6]) 
Di = state.Constant([3e-3]) 
A[0,...] = 0
A[-1,...] = 0
Ms[0,...] = 0
Ms[-1,...] = 0
Di[0,...] = 0
Di[-1,...] = 0
Ku[0,...] = 0
Ku[-1,...] = 0

state.material['Ms'] = Ms
state.material['A'] = A
state.material['Ku'] = Ku
state.material['Ku_axis'] = [0,0,1]
state.material['Di'] = -Di

write_vti(state.material, "data/material.vti", state)

state.m = state.Constant([-0.1, 0.0, 0.9])
state.m.normalize()
state.m[0,...] = 0
state.m[-1,...] = 0

exchange = ExchangeDMIField2()
aniso    = UniaxialAnisotropyField()

llg = LLGSolver([exchange, aniso])
logger = ScalarLogger("data/m_magnumnp.dat", ['t', 'm'])
while state.t <= 2e-9:
    logger << state
    llg.step(state, 1e-11)

write_vti(state.m, "data/m0_magnumnp.vti", state)
np.savetxt("data/m0_magnumnp.dat", torch.concat((x[:,0,0,None], state.m[:,0,0,:]), axis=1).cpu().numpy())
Timer.print_report()
