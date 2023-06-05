from magnumnp import *
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Run simple simulation.')
parser.add_argument('-o', '--output', default="Output", type=str)
args = parser.parse_args()

folder = args.output
if not os.path.exists(folder):
    os.makedirs(folder)

A1 = 3.5
A2 = 4.5
A3 = 1.5

nx, ny, nz = 10, 10, 10
nzsp = 2


B1 = 0.0
B2 = 0.0
B3 = 0.0

J_funcFeFe=lambda d: A1/(d-B1)**2.
J_funcCoFe=lambda d: A2/(d-B2)**2.
J_funcCoCo=lambda d: A3/(d-B3)**2.

Fe_ratio = 0.8
a = 2.35e-10
dx = (a, a, a)
phi = 1.0
# initialize mesh
n  = (nx, ny, nz)
mesh = Mesh(n, dx)
state = State(mesh, device = "cpu")

J_FeFe = 7.e-21
J_CoCo = 6.e-21
J_CoFe = 4.5e-21

# 0 = Ru, 1 = X (e.g. Fe), 2 = Y (e.g. Co) ...
mat = state.Constant([1], dtype = int)
while torch.sum(mat) > nx*ny*nz * (Fe_ratio):
    mat[np.random.randint(0,n[0]),
        np.random.randint(0,n[1]),
        np.random.randint(0,n[2])] = 0


J = state.zeros((3, 4))
J[2,:] = state.Tensor([J_FeFe] + [- float(J_funcFeFe(i*a*1e9))*1e-3 * a**2 for i in range(1,4)]) # Fe Fe

Ms = 1500e3
state.material = {
        "RuxJij": J,
        "Ms": Ms,
        "RuxDistribution": mat,
        "alpha": 1.0
}


# initialize mag with 90 degree angle
x,y,z = state.SpatialCoordinate()
rux = (mat.squeeze(-1) == 0)
magnetic = (mat.squeeze(-1) == 1)

state.m = state.Constant([0,0,0])
state.m += torch.rand((nx, ny, nz, 3))
state.m[~magnetic] = 0.
state.m.normalize()

h = 2. / constants.mu_0

external = ExternalField([h / constants.mu_0, 0., 0.])
rux = AtomisticRuXExchangeField(Jij = J)
h = rux.h(state)
write_vti(h, "Hrux.vti")

# minimize energy
logger = Logger(folder, ['m', rux.E, external.E], ['m', rux.h])
llg = LLGSolver([rux, external])
for i in range(1000):
    llg.step(state, 1e-15)
    logger << state


