from fd import *
import numpy as np
from scipy import constants

# initialize mesh
n  = (10, 10, 10)
dx = (1e-9, 1e-9, 1e-9)
mesh = Mesh(n, dx)

# initialize material
material = {
        "Ms": 1.,
        "A": 1.,
        "gamma": 1.,
        }

# initialize field terms
#demag    = DemagField(mesh, material)
#exchange = ExchangeField(mesh, material)
external = ExternalField(mesh, material, [1.,0.,0.])

# initialize magnetization that relaxes into s-state
m0 = np.zeros(n + (3,))
m0[:,:,:,0]   = 1.0

# initialize sstate
#eigen = EigenSolver([demag, exchange, external], material, m0)
eigen = EigenSolver([], [external], material, m0)
vv = np.zeros(n + (2,))
vv[:,:,:,0] = 1.
vv = vv.reshape(-1)
print("D0:", eigen._D0(vv))

eigen.solve(k=5)
