from magnumnp import *
import torch
import itertools
from math import sin, cos, pi
import numpy as np

Timer.enable()

# initialize state
eps = 1e-15
n  = (100, 40, 1)
dx = (5e-9, 5e-9, 4e-9)
mesh = Mesh(n, dx, origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., 0.0))

state = State(mesh)
state.material = {
    "Ms": 800e3,
    "A": 1.3e-11,
    "alpha": 0.2
    }

x, y, z = state.SpatialCoordinate()
a = 470e-9/2.
b = 170e-9/2.
magnetic = ((x/a)**2. + (y/b)**2. <= 1.)

# initialize field terms
demag = DemagField()
exchange = ExchangeField()

# initialize magnetization that relaxes into SAF Vortex
m0 = state.Constant([1, 0, 0])
m0[~magnetic] = 0.
m1 = state.Constant([0, 1, 0])
m1[~magnetic] = 0.
m2 = state.Constant([-1, 0, 0])
m2[~magnetic] = 0.

images = [m0, m1, m2]
string = StringSolver([demag, exchange], num_images = 21)

with open("data/E.dat", "w") as efile:
    for i in range(100):
        images = string.step(state, images)

        E = string.E(state, images)
        for k, ek, in enumerate(E):
            efile.write("%g %g %g\n" % (i, k, ek))
            write_vti(images[k], "data/m_%.4d.vti" % k)
        efile.write("\n")
        efile.flush()

print("min(E):", min(E), "max(E):", max(E), "dE:", max(E)-min(E))
efile.close()
Timer.print_report()
