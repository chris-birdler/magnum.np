from magnumnp import *
import torch
import numpy as np

Timer.enable()

# initialize mesh
eps = 1e-15
n  = (50, 50, 50)
dc = np.sqrt(2)*2e-8
dx = (dc, dc, dc)
mesh = Mesh(n, dx)

# initialize state
state = State(mesh)
state.material = {
        "Ms": 1.07 / constants.mu_0,
        "A": 1.0e-11,
        "alpha": 0.10,
        "Kc1": 6e3,
        "Kc2": 5e2,
        "Kc_alpha": 0.,
        "Kc_beta ": 0.,
        "Kc_gamma": 0.
        }

x, y, z = state.SpatialCoordinates()
r = 16 * dc
gap = 50*dc - 2. * r
l = 2.*r + gap

magnetic1  = ((x  )**2 + (y  )**2 + (z  )**2 <= r**2)
magnetic2  = ((x-l)**2 + (y  )**2 + (z  )**2 <= r**2)
magnetic3  = ((x  )**2 + (y-l)**2 + (z  )**2 <= r**2)
magnetic4  = ((x  )**2 + (y  )**2 + (z-l)**2 <= r**2)
magnetic5  = ((x-l)**2 + (y  )**2 + (z-l)**2 <= r**2)
magnetic6  = ((x-l)**2 + (y-l)**2 + (z  )**2 <= r**2)
magnetic7  = ((x  )**2 + (y-l)**2 + (z-l)**2 <= r**2)
magnetic8  = ((x-l)**2 + (y-l)**2 + (z-l)**2 <= r**2)
magnetic9  = ((x-l/2.)**2 + (y     )**2 + (z-l/2.)**2 <= r**2)
magnetic10 = ((x-l/2.)**2 + (y-l/2.)**2 + (z     )**2 <= r**2)
magnetic11 = ((x     )**2 + (y-l/2.)**2 + (z-l/2.)**2 <= r**2)
magnetic12 = ((x-l/2.)**2 + (y-l   )**2 + (z-l/2.)**2 <= r**2)
magnetic13 = ((x-l/2.)**2 + (y-l/2.)**2 + (z-l   )**2 <= r**2)
magnetic14 = ((x-l   )**2 + (y-l/2.)**2 + (z-l/2.)**2 <= r**2)

magnetic = magnetic1 | magnetic2 | magnetic3 | magnetic4 | magnetic5 | magnetic6 | magnetic7 | magnetic8 | magnetic9 | magnetic10 | magnetic11 | magnetic12 | magnetic13 | magnetic14
write_vti(magnetic, "data/domain.vti")

Ms = state.Constant([0.0])
Ms[magnetic] = 1.07 / constants.mu_0

Kc1 = state.Constant([0.0])
Kc1[magnetic] = 6e3

Kc2 = state.Constant([0.0])
Kc2[magnetic] = 5e2

A = state.Constant([0.0])
A[magnetic] = 1e-11

state.material["Ms"] = Ms
state.material["A"] = A
state.material["Kc1"] = Kc1
state.material["Kc2"] = Kc2

state.material["Kc_alpha"] =  state.Constant([0.])
state.material["Kc_beta"] =  state.Constant([0.])
state.material["Kc_gamma"] =  state.Constant([0.])

# initialize field terms
exchange = ExchangeField()
cubicaniso = CubicAnisotropyField()
demagPBC = DemagFieldPBC()

# initialize magnetization that relaxes into s-state
state.m = state.Constant([1.0, 0.5, 1.0])
state.m[~magnetic] = 0.0
state.m.normalize()

llg = LLGSolver([exchange, cubicaniso, demagPBC])
with open('data/m.dat', 'w') as f:
    i = 0
    while state.t < 50e-9-eps:
        llg.step(state, 1e-10)
        f.write("%g %g %g %g %g %g %g\n" % ((state.t,) + tuple(state.m.avg()) + tuple(state.m[magnetic].avg())))
        f.flush()
        i += 1

Timer.print_report()
