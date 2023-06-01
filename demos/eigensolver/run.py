from magnumnp import *
from torch import sqrt

n = (20, 20, 1)
dx = (10e-9, 10e-9, 1e-9)
origin = (-1./2.*n[0]*dx[0], -1./2.*n[1]*dx[1], -1./2.*n[2]*dx[2])

mesh = Mesh(n, dx, origin)
state = State(mesh)

x,y,z = state.SpatialCoordinate()
r = x**2 + y**2

state.m = state.Tensor(torch.stack([y/r, -x/r, 1e-9/r], dim = -1))
state.m.normalize()

state.material = {"Ms": 800e3,
                  "A": 13e-12,
                  "Ku": -5e3,
                  "Ku_axis": (0,0,1)}

exchange = ExchangeField()
demag = DemagField()
external = ExternalField([10e3, 0, 0])
eigs = EigenSolver(state, [demag, exchange], [external])
res = eigs.solve(k = 5)
#omega, evecs2D, state, m0, D0 = eigs.solve(k = 5)

evecs = res.evecs()
print("f[GHz]:",  res.omega/2./torch.pi*1e-9)
res.save_evecs3D("data/evecs.vti")
#print(res[0])
#print(res[1].shape)
