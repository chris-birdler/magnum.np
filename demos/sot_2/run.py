from magnumnp import *
import torch
import argparse
import itertools
import os
import numpy as np
from infineon_params import InfineonParameters
from math import sqrt

torch.set_default_dtype(torch.float64)
#torch.set_float32_matmul_precision('high')


parser = argparse.ArgumentParser(description='Run simple simulation.')
parser.add_argument('-r', '--rank', default=0, type=int)
args = parser.parse_args()

Hkarr = np.linspace(-10, 10, 51)*1e-3
EFarr = np.linspace(-0.05, 0.05, 51)
EDarr = np.linspace(-0.05, 0.05, 51)

params = list(itertools.product(Hkarr, EFarr, EDarr))
var = params[args.rank]

Hk = var[0] 
eta_field = var[1] 
eta_damp = var[2] 
je = 2.4e11
Ms = 0.75/constants.mu_0 
Keff = Ms*Hk/2.

thickness = 1 

outdir = "data/Hk_%g_mT/"%(Hk*1e3)

if not os.path.exists(outdir):
    os.makedirs(outdir)

Timer.enable()

# initialize mesh
eps = 1e-15
n  = (1, 1, 1)
dx = (1e-9, 1e-9, thickness*1e-9)
origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.) #shift origin of the mesh

#initialize mesh and state
mesh = Mesh(n, dx, origin)
state = State(mesh, dtype = torch.float64, device="cpu" )
x, y, z = state.SpatialCoordinate()

# thickness of thin film on which the SOT acts
p = state.Constant([0, -1, 0])#make it to function
d = n[2] * dx[2]#define thickness of magnetic layer for SOT


state.material = {
        "Ms": Ms,
        "Ku": Keff,
        "Ku_axis": state.Tensor([0,0,1]),
        "alpha": 1.0,
        "eta_damp": eta_field,# both eta with opposite sign as magnum.af, same as magnum.pi
        "eta_field": eta_damp,
        "p": p,
        "d": d,
        "je": je,
        }
torch.manual_seed(30032023)
state.m = state.Constant([0, 0, 0]) 
mrand = torch.empty_like(state.m).normal_(mean = 0., std = 1.0)
state.m += mrand
state.m.normalize()
minit = state.m.clone()

# initialize field terms
aniso = UniaxialAnisotropyField()
zeeman = ExternalField(state.Tensor([0, 0, 0]))
logger = ScalarLogger("%s/Eta_Field_%.2f_Eta_Damp_%.2f.dat"%(outdir, eta_field, eta_damp), ['t', 'm', zeeman.h])

def simstep(m0, Hinput, jinput, logger= logger):
    state.m = m0
    state.material["je"] = jinput
    torque = SpinOrbitTorque()
    zeeman = ExternalField(state.Tensor([Hinput/sqrt(2), 0, Hinput/sqrt(2)]))
    relaxtime = 100e-9    
    llg = LLGSolver([aniso, torque, zeeman], solver = ScipyODE)
    while state.t <= relaxtime:
        llg.step(state, 20e-9)
    logger << state
    state.t=0

for hi in np.arange(-20, 20.25, 0.25)*1e-3/constants.mu_0:
    for ji in [-1*je,je]:
        simstep(minit, hi, ji)

Timer.print_report()